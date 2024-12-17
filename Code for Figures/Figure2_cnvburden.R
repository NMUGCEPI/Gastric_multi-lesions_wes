##########################################################################################

library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(optparse)
library(ggpubr)
library(patchwork)

###########################################################################################

info <- data.frame(fread(sample_file))
seg <- data.frame(fread(seg_file))
class_order <- data.frame(fread(class_order_file))

info_public <- data.frame(fread(sample_public_file))
seg_public <- data.frame(fread(seg_public_file))

t_logR <- 0.2

###########################################################################################
col <- c(
  brewer.pal(9,"YlGnBu")[6],
  rgb(234,106,79,alpha=255,maxColorValue=255),
  rgb(203,24,30,alpha=255,maxColorValue=255),
  rgb(255,0,0,alpha=255,maxColorValue=255)
  )

names(col) <- c("IM" , "IGC" , "DGC" , "GC")
col <- col[1:4]

col_im <- brewer.pal(9,"YlGnBu")[6:8]
names(col_im) <- c("IM(IGC)" , "IM(DGC)" , "IM(IGC_DGC)")

###########################################################################################

dat <- seg
dat$seg.mean <- dat$Median_logR
dat$loc.start <- dat$Start
dat$loc.end <- dat$End
dat$Tumor <- sapply( strsplit( dat$Sample , "_" ) , "[" , 1)
dat$Normal <- sapply( strsplit( dat$Sample , "_" ) , "[" , 2)

###########################################################################################

dat <- merge( dat , info[,c("Tumor" , "Class" , "Type" , "TCGA_Class")] , by = "Tumor" )

###########################################################################################

dat_public <- seg_public
dat_public <- merge( dat_public , info_public[,c("Tumor" , "Class" , "From" , "Molecular.subtype")] , by.x = "Sample" , by.y = "Tumor" )
dat_public$Normal <- dat_public$Sample

###########################################################################################

use_col <- c("Tumor" , "Chromosome" , "Start" , "End" , "Length.snp." , "Median_logR" , "Class" , "Normal" , "TCGA_Class")
dat_im <- subset(dat , Class=="IM")[,use_col]
dat_gc <- subset(dat , Class!="IM")[,use_col]
dat_nmu <- rbind(dat_im , dat_gc)
dat_nmu$From <- "NJMU"
dat_nmu <- dat_nmu[,c("Tumor" , "Chromosome" , "Start" , "End" , "Length.snp." , "Median_logR" , "Class" , "From" , "TCGA_Class" , "Normal")]
colnames(dat_nmu) <- colnames(dat_public)

###########################################################################################

dat_combine <- rbind(dat_nmu , dat_public)
dat_combine$seg.mean <- dat_combine$Segment_Mean
dat_combine$loc.start <- dat_combine$Start
dat_combine$loc.end <- dat_combine$End
dat_combine$Tumor <- dat_combine$Sample
dat_combine <- subset( dat_combine , Chromosome %in% c(1:22) )

###########################################################################################
dat <- dat_combine
dat$TCN <- 2*2^dat$seg.mean
dat$Length <- as.numeric(abs(dat$loc.start - dat$loc.end))

loss_region <- dat[which(dat$seg.mean <= -t_logR),]
gain_region <- dat[which(dat$seg.mean >= t_logR),]

result_loss <- c()
result_gain <- c()
for(Tumor in unique(dat$Tumor)){
  print(Tumor)
  class_s <- unique(dat[which(dat$Tumor==Tumor),'Class'])
  normal <- unique(dat[which(dat$Tumor==Tumor),'Normal'])
  Molecular.subtype <- unique(dat[which(dat$Tumor==Tumor),'Molecular.subtype'])
  
  loss_rate <- sum(loss_region[which(loss_region$Tumor==Tumor),'Length'])/sum(dat[which(dat$Tumor==Tumor),'Length'])
  result_loss <- rbind(result_loss,data.frame(Tumor=Tumor,Normal=normal,Rate=loss_rate,Class=class_s , Molecular.subtype = Molecular.subtype))

  gain_rate <- sum(gain_region[which(gain_region$Tumor==Tumor),'Length'])/sum(dat[which(dat$Tumor==Tumor),'Length'])
  result_gain <- rbind(result_gain,data.frame(Tumor=Tumor,Normal=normal,Rate=gain_rate,Class=class_s , Molecular.subtype = Molecular.subtype))
}

min_burden <- -0.0001
max_burden <- 0.7

result_gain$Class <- factor(result_gain$Class,levels= unique(class_order$Class), ordered=TRUE)
result_loss$Class <- factor(result_loss$Class,levels= unique(class_order$Class), ordered=TRUE)

###########################################################################################
dat <- result_loss
dat_loss<- c() 
for(Normal in unique(dat$Normal)){
  print(Normal)
  tmp <- dat[which(dat$Normal==Normal),]

  tmp <- tmp %>%
  group_by( Normal , Class , Molecular.subtype) %>%
  summarize( CNV_Burden_Loss = median(Rate) )

  tmp <- na.omit(data.frame(tmp))

  dat_loss <- rbind(dat_loss,tmp)
}

## gain
dat <- result_gain
dat_gain <- c() 
for(Normal in unique(dat$Normal)){
  print(Normal)
  tmp <- dat[which(dat$Normal==Normal),]

  tmp <- tmp %>%
  group_by( Normal , Class , Molecular.subtype) %>%
  summarize( CNV_Burden_Gain = median(Rate) )

  tmp <- na.omit(data.frame(tmp))
  dat_gain <- rbind(dat_gain,tmp)
}

###########################################################################################
trans <- function(num){
    up <- floor(log10(num))
    down <- round(num*10^(-up),2)
    text <- paste("p == ",down," %*% 10","^",up)
    return(text)
}

plotFunction <- function(dat_plot = dat_plot , y_lab = y_lab , my_comparisons_1 = my_comparisons_1 , title = title){

  sample_num <- dat_plot %>%
  group_by( Class ) %>%
  summarize( nums = length(unique(Normal)) )

  dat_plot <- merge( dat_plot , sample_num , by = "Class" )
  dat_plot$Class_num <- paste0( dat_plot$Class , "\n(" , dat_plot$nums , ")" )
  dat_plot$Class_num <- factor( dat_plot$Class_num , levels = unique(dat_plot$Class_num)[2:1] )


  a <- dat_plot[dat_plot$Class==unique(dat_plot$Class)[1],"Burden_use"]
  b <- dat_plot[dat_plot$Class==unique(dat_plot$Class)[2],"Burden_use"]
  p <- wilcox.test( a , b )$p.value

  if( p < 0.01 ){
      p_text <- trans(p)
  }else{
      p_text <- paste0( "p == " , round(as.numeric(p) , 3) ) 
  }
  dat_plot$p_text <- ""
  dat_plot$p_text[1] <- p_text

  col_tmp <- c(
    rgb(red=247,green=184,blue=71,alpha=255,max=255) ,
    rgb(red=2,green=100,blue=190,alpha=255,max=255) 
  )

  if(1!=1){
    plot <- ggplot(data=dat_plot,mapping = aes(x=Class_num,y=Burden_use))+
    geom_boxplot(lwd=1.5,aes(color=Class) , outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.2),aes(color=Class)) +
    scale_color_manual(values=col_tmp) +
    xlab(NULL) +
    ylab(y_lab)+
    labs(title=title) +
    theme_bw() +
    geom_text(aes(label=p_text , y = 0.7 , x = 1.5),parse = TRUE,size=4 , color = "black") +
    ylim(min_burden,max_burden) +
    theme(panel.background = element_blank(),#设置背影为白色#清除网格线
          legend.position ='none', # 隐藏图例
          legend.title = element_blank() ,
          panel.grid.major=element_line(colour=NA),
          legend.text = element_text(size = 8,color="black",face='bold'),
          axis.text.x = element_text(size = 10,color="black",face='bold'),
          axis.text.y = element_text(size = 10,color="black",face='bold'),
          axis.title.x = element_text(size = 10,color="black",face='bold'),
          axis.title.y = element_text(size = 10,color="black",face='bold'),
          axis.line = element_line(size = 0.5)) 
  }

  ebtop<-function(x){
    return(quantile(x)[4])
  }
  ebbottom<-function(x){
    return(quantile(x)[2])
  }

  y_max <- 0.93
  plot <- ggplot( data = dat_plot , aes( x = Class_num , y = Burden_use , fill = Class ))+
  stat_summary(geom = "bar",fun = "median",
             position = position_dodge(0.9))+
  stat_summary(geom = "errorbar",
             fun.min = ebbottom,
             fun.max = ebtop,
             position = position_dodge(0.9),
             width=0.4)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  geom_text(aes(label=p_text , y = 0.6 , x = 1.5) , size=4.5 , color = "black" , face = "bold" , parse = TRUE) +
  theme_bw()+
  xlab(NULL) +
  ylab(y_lab)+
  ylim(-0.02 , 0.65) +
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=col_tmp) +
  theme(
      legend.position = 'none',
      legend.title = element_blank() ,
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(size = 12,color="black",face='bold'),
      legend.text = element_text(size = 12,color="black",face='bold'),
      axis.text.y = element_text(size = 12,color="black",face='bold'),
      axis.title.x = element_text(size = 12,color="black",face='bold'),
      axis.title.y = element_text(size = 12,color="black",face='bold'),
      axis.text.x = element_text(size = 12,color="black",face='bold') ,
      axis.ticks.length = unit(0.2, "cm") ,
      strip.text.x = element_text(size = 15, colour = "black",face='bold') ,
      axis.line = element_line(size = 0.5)) 
  return(plot)
}

my_comparisons_1 <- list( c(1, 2) )
title <- "All"

y_lab <- "Fraction of CNV altered with\n(gain)"
dat_plot$Burden_use <- dat_plot$CNV_Burden_Gain

dat_plot <- dat_plot %>%
group_by( Normal , Class ) %>%
summarize( Burden_use = median(Burden_use) )
dat_plot$Class <- factor(dat_plot$Class,levels=c("IGC" , "DGC") , ordered=T)

p1 <- plotFunction(dat_plot = dat_plot , y_lab = y_lab , my_comparisons_1 = my_comparisons_1 , title = "All")

y_lab <- "Fraction of CNV altered with\n(loss)"
dat_plot$Burden_use <- dat_plot$CNV_Burden_Loss

dat_plot <- dat_plot %>%
group_by( Normal , Class ) %>%
summarize( Burden_use = median(Burden_use) )
dat_plot$Class <- factor(dat_plot$Class,levels=c("IGC" , "DGC") , ordered=T)
p2 <- plotFunction(dat_plot = dat_plot , y_lab = y_lab , my_comparisons_1 = my_comparisons_1 , title = "All")

images_name <- paste(images_path,"/Fig2_cnvburden.All.pdf",sep="")
result_p <- p1 + p2
ggsave(file=images_name,plot=result_p,width=7.3/1.6,height=5.18/1.6)

###########################################################################################
my_comparisons_1 <- list( c(1, 2) )

for(molecular in c("GS" , "CIN" , "MSI")){

  title <- molecular
  y_lab <- "Fraction of CNV altered with\n(gain)"
  dat_plot <- subset( dat_gain , Class != "IM" & Molecular.subtype == molecular )
  dat_plot$Burden_use <- dat_plot$CNV_Burden_Gain
  dat_plot$Class <- factor(dat_plot$Class,levels=c("IGC" , "DGC") , ordered=T)
  p1 <- plotFunction(dat_plot = dat_plot , y_lab = y_lab , my_comparisons_1 = my_comparisons_1 , title = title)

  y_lab <- "Fraction of CNV altered with\n(loss)"
  dat_plot <- subset( dat_loss , Class != "IM" & Molecular.subtype == molecular )
  dat_plot$Burden_use <- dat_plot$CNV_Burden_Loss
  dat_plot$Class <- factor(dat_plot$Class,levels=c("IGC" , "DGC") , ordered=T)
  p2 <- plotFunction(dat_plot = dat_plot , y_lab = y_lab , my_comparisons_1 = my_comparisons_1 , title = title)

  images_name <- paste(images_path,"/Fig2_cnvburden.",molecular,".pdf",sep="")
  result_p <- p1 + p2
  ggsave(file=images_name,plot=result_p,width=7.3/1.6,height=5.18/1.6)
}
