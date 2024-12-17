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

dat <- merge( dat , info[,c("Tumor" , "Class" , "Type" , "TCGA_Class")] , by = "Tumor" )
colnames(dat)[ncol(dat)] <- "Molecular.subtype"

###########################################################################################

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
  type <- unique(dat[which(dat$Tumor==Tumor),'Type'])
  Molecular.subtype <- unique(dat[which(dat$Tumor==Tumor),'Molecular.subtype'])

  
  loss_rate <- sum(loss_region[which(loss_region$Tumor==Tumor),'Length'])/sum(dat[which(dat$Tumor==Tumor),'Length'])
  result_loss <- rbind(result_loss,data.frame(Tumor=Tumor,Normal=normal,Rate=loss_rate,Class=class_s,Type= type,Molecular.subtype=Molecular.subtype))

  gain_rate <- sum(gain_region[which(gain_region$Tumor==Tumor),'Length'])/sum(dat[which(dat$Tumor==Tumor),'Length'])
  result_gain <- rbind(result_gain,data.frame(Tumor=Tumor,Normal=normal,Rate=gain_rate,Class=class_s,Type= type,Molecular.subtype=Molecular.subtype))
}

min_burden <- -0.0001
max_burden <- 0.6

result_gain$Class <- factor(result_gain$Class,levels= unique(class_order$Class), ordered=TRUE)
result_loss$Class <- factor(result_loss$Class,levels= unique(class_order$Class), ordered=TRUE)

###########################################################################################
dat <- result_loss
dat_loss<- c() 
for(Normal in unique(dat$Normal)){
  print(Normal)
  tmp <- dat[which(dat$Normal==Normal),]

  tmp <- tmp %>%
  group_by( Normal , Class , Type , Molecular.subtype ) %>%
  summarize( CNV_Burden_Loss = median(Rate) )

  tmp <- na.omit(data.frame(tmp))

  dat_loss <- rbind(dat_loss,tmp)
}

dat <- result_gain
dat_gain <- c() 
for(Normal in unique(dat$Normal)){
  print(Normal)
  tmp <- dat[which(dat$Normal==Normal),]

  tmp <- tmp %>%
  group_by( Normal , Class , Type , Molecular.subtype ) %>%
  summarize( CNV_Burden_Gain = median(Rate) )

  tmp <- na.omit(data.frame(tmp))
  dat_gain <- rbind(dat_gain,tmp)
}

###########################################################################################
trans <- function(num){
    up <- floor(log10(num))
    down <- round(num*10^(-up),2)
    text <- paste("P == ",down," %*% 10","^",up)
    return(text)
}

plotFunction <- function(dat_plot = dat_plot , y_lab = y_lab , my_comparisons_1 = my_comparisons_1 , title = title){

  dat_tmp <- c()

  for(type in unique(dat_plot$Type) ){

    dat_plot_tmp <- subset( dat_plot , Type == type )

    a <- dat_plot_tmp[dat_plot_tmp$Class==unique(dat_plot_tmp$Class)[1],"Burden_use"]
    b <- dat_plot_tmp[dat_plot_tmp$Class==unique(dat_plot_tmp$Class)[2],"Burden_use"]
    p <- wilcox.test( a , b )$p.value

    if( p < 0.01 ){
        p_text <- trans(p)
    }else{
        p_text <- paste0( "P == " , round(as.numeric(p) , 3) ) 
    }
    dat_plot_tmp$p_text <- ""
    dat_plot_tmp$p_text[1] <- p_text
    dat_tmp <- rbind( dat_plot_tmp , dat_tmp )
  }
  col_tmp <- c(
    rgb(red=247,green=184,blue=71,alpha=255,max=255) ,
    rgb(red=2,green=100,blue=190,alpha=255,max=255) ,
    rgb(red=2,green=100,blue=190,alpha=255,max=255) 
  )

  plot <- ggplot(data=dat_tmp , mapping = aes(x=Class,y=Burden_use))+
  geom_line( aes( group = Normal ) , size = 0.3 , color = "gray" ) +
  geom_boxplot(lwd=1.5,aes(color=Class) , outlier.shape = NA , size = 0.2) +
  geom_jitter(position=position_jitter(0.2),aes(color=Class) , size = 1) +
  scale_color_manual(values=col_tmp) +
  facet_grid(.~Type , scales="free_x" ) +
  xlab(NULL) +
  ylab(y_lab)+
  theme_bw() +
  geom_text(aes(label=p_text , y = 0.6 , x = 1.5),parse = TRUE,size=3 , color = "black") +
  ylim(min_burden,max_burden) +
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
      strip.text.x = element_text(size = 14, colour = "black",face='bold') ,
      axis.line = element_line(size = 0.5)) 
  return(plot)
}

my_comparisons_1 <- list( c(1, 2) )
molecular <- "All"
dat_plot_tmp_gain <- dat_gain
dat_plot_tmp_loss <- dat_loss
title <- molecular

y_lab <- "Fraction of CNV altered\n(gain)"
dat_plot <- dat_plot_tmp_gain
dat_plot$Burden_use <- dat_plot$CNV_Burden_Gain
dat_plot <- dat_plot %>%
group_by( Normal , Class , Type ) %>%
summarize( Burden_use = median(Burden_use) )
dat_plot$Class <- factor(dat_plot$Class,levels=c( "IM" , "IGC" , "DGC") , ordered=T)
sample_num <- dat_plot %>%
group_by( Type ) %>%
summarize( nums = length(unique(Normal)) )
dat_plot <- merge( dat_plot , sample_num , by = "Type" )
dat_plot$Type <- factor( dat_plot$Type , levels = c("IM + IGC" , "IM + DGC") , order = T )
p1 <- plotFunction(dat_plot = dat_plot , y_lab = y_lab , my_comparisons_1 = my_comparisons_1 , title = title)

y_lab <- "Fraction of CNV altered\n(loss)"
dat_plot <- dat_plot_tmp_loss
dat_plot$Burden_use <- dat_plot$CNV_Burden_Loss
dat_plot <- dat_plot %>%
group_by( Normal , Class , Type  ) %>%
summarize( Burden_use = median(Burden_use) )
dat_plot$Class <- factor(dat_plot$Class,levels=c( "IM" , "IGC" , "DGC") , ordered=T)
sample_num <- dat_plot %>%
group_by( Type ) %>%
summarize( nums = length(unique(Normal)) )
dat_plot <- merge( dat_plot , sample_num , by = "Type" )
p2 <- plotFunction(dat_plot = dat_plot , y_lab = y_lab , my_comparisons_1 = my_comparisons_1 , title = title)

## 合并
images_name <- paste(images_path,"/Fig3B.pdf",sep="")
result_p <- p1 + p2
ggsave(file=images_name,plot=result_p,width=6.1/0.8,height=4.8/0.8)

