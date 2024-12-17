##########################################################################################

library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(optparse)

###########################################################################################

col <- c(
  brewer.pal(9,"YlGnBu")[6],
  rgb(234,106,79,alpha=255,maxColorValue=255),
  rgb(203,24,30,alpha=255,maxColorValue=255),
  rgb(255,0,0,alpha=255,maxColorValue=255)
  )

names(col) <- c("IM" , "IGC" , "DGC" , "GC")
col_im_gc <- col[c(1,4)]
col_igc_dgc <- col[c(2,3)]

###########################################################################################

dat_mutRateGene <- data.frame(fread( mut_rate_gene_file ))
smg <- data.frame(fread(smg_file))$Gene_Symbol

###########################################################################################

dat_mutRateGene <- subset(dat_mutRateGene , Hugo_Symbol %in% smg)

dat_plot <- rbind( dat_mutRateGene )
dat_plot$Class <- factor( dat_plot$Class , levels = c("IM" , "GC" , "IGC" , "DGC") , order = T )
dat_plot$value_text <- paste0( round(dat_plot$MutRate , 2) * 100 , "%") 

###########################################################################################

dat_plot <- subset( dat_plot , From == from )

im_num <- unique(subset(dat_plot , Class=="IM")$SampleNum)
gc_num <- unique(subset(dat_plot , Class=="GC")$SampleNum)
igc_num <- unique(subset(dat_plot , Class=="IGC")$SampleNum)
dgc_num <- unique(subset(dat_plot , Class=="DGC")$SampleNum)

###########################################################################################
dat_plot$p.value = ""
dat_plot$p_text = ""

result_tmp <- c()

for(geneN in unique(dat_plot$Hugo_Symbol)){

    print(geneN)

    tmp_1 <- subset( dat_plot , Hugo_Symbol == geneN & Class %in% c("GC") )
    tmp_2 <- subset( dat_plot , Hugo_Symbol == geneN & Class %in% c("IM") )

    if(nrow(tmp_1)==0){
        tmp_1 <- tmp_2
        tmp_1$Class <- "GC"
        tmp_1$SampleNum <- gc_num
        tmp_1$MutNum <- 0
        tmp_1$MutRate <- 0
        tmp_1$value_text <- "0%"
    }

    if(nrow(tmp_2)==0){
        tmp_2 <- tmp_1
        tmp_2$Class <- "IM"
        tmp_2$SampleNum <- im_num
        tmp_2$MutNum <- 0
        tmp_2$MutRate <- 0
        tmp_2$value_text <- "0%"
    }

    tmp <- rbind( tmp_1 , tmp_2 )

    tmp_fisher <- matrix(c(tmp$MutNum , tmp$SampleNum - tmp$MutNum) , ncol = 2)
    p <- fisher.test(tmp_fisher)$p.value

    tmp$p.value <- p

    result_tmp <- rbind( result_tmp , tmp )

}

for(geneN in unique(dat_plot$Hugo_Symbol)){

    print(geneN)

    tmp_1 <- subset( dat_plot , Hugo_Symbol == geneN & Class %in% c("IGC") )
    tmp_2 <- subset( dat_plot , Hugo_Symbol == geneN & Class %in% c("DGC") )

    if(nrow(tmp_1)==0){
        tmp_1 <- tmp_2
        tmp_1$Class <- "IGC"
        tmp_1$SampleNum <- igc_num
        tmp_1$MutNum <- 0
        tmp_1$MutRate <- 0
        tmp_1$value_text <- "0%"
    }

    if(nrow(tmp_2)==0){
        tmp_2 <- tmp_1
        tmp_2$Class <- "DGC"
        tmp_2$SampleNum <- dgc_num
        tmp_2$MutNum <- 0
        tmp_2$MutRate <- 00
        tmp_2$value_text <- "0%"
    }

    tmp <- rbind( tmp_1 , tmp_2 )
    tmp_fisher <- matrix(c(tmp$MutNum , tmp$SampleNum - tmp$MutNum) , ncol = 2)
    p <- fisher.test(tmp_fisher)$p.value

    tmp$p.value <- p

    result_tmp <- rbind( result_tmp , tmp )
}

result_tmp_gc <- subset( result_tmp , Class == "GC" )
result_tmp_im <- subset( result_tmp , Class == "IM" )

result_tmp_order <- merge( result_tmp_gc[,c("Hugo_Symbol" , "MutRate")] , result_tmp_im[,c("Hugo_Symbol" , "MutRate")] , by = "Hugo_Symbol" )

gene_order <- c("TP53" , "ARID1A" , "CDH1" , "APC" , "SMAD4" , "DNAH3" , "MUC6" , "PIK3CA" , "CTNNB1" , 
    "RHOA" , "ERBB2" , "CFTR" , "KRAS" , "MAP2K7" , "ARID2" , "RNF43" , "TGFBR2" ,
    "BMP6" , "FBXW7" , "CDKN2A" , "MTRR" , "MICAL2")

result_tmp$Hugo_Symbol <- factor( result_tmp$Hugo_Symbol , 
    levels = gene_order , order = T)


###########################################################################################

result_use <- result_tmp

result_use$p <- result_use$p.value
result_use$p_text=ifelse(result_use$p>=0.05,"","*")
result_use$p_text=ifelse(result_use$p<0.05 & result_use$p>0.01,"*",result_use$p_text)
result_use$p_text=ifelse(result_use$p<0.01 & result_use$p>0.001,"**",result_use$p_text)
result_use$p_text=ifelse(result_use$p<0.001 ,"***",result_use$p_text)

result_use$p_pos <- ""
result_use[result_use$Class %in% c("IGC" , "DGC"),"p_pos"] <- 0.8
result_use[result_use$Class %in% c("IM" , "GC"),"p_pos"] <- -0.8
result_use$p_pos <- as.numeric(result_use$p_pos)

result_use$label_pos <- ""
result_use[result_use$Class %in% c("IGC" , "DGC"),"label_pos"] <- 0.5
result_use[result_use$Class %in% c("IM" , "GC"),"label_pos"] <- -0.5

result_use$Type <- ""
result_use[result_use$Class %in% c("IGC" , "DGC"),"Type"] <- "IGC vs DGC"
result_use[result_use$Class %in% c("IM" , "GC"),"Type"] <- "IM vs GC"

result_use$MutRate <- ifelse( result_use$Type == "IGC vs DGC" , result_use$MutRate , -result_use$MutRate )
result_use$percent_pos <- ifelse( result_use$Type == "IGC vs DGC" , result_use$MutRate + 0.03 , result_use$MutRate-0.06 )
result_use$percent_pos <- as.numeric(result_use$percent_pos)

p <- ggplot(result_use,mapping = aes(x = Hugo_Symbol , y = MutRate , fill = Class)) +
geom_bar(stat = 'identity', position = 'dodge' , width = 0.8 , color = 'black') + 
facet_grid(vars(Type) , scales = "free")+
theme_bw() +
scale_y_continuous(
    breaks = seq(-0.8,0.8,0.2) , 
    label = c( "80%" , "60%" , "40%" , "20%" , 
        0 ,
        "20%" , "40%" , "60%" , "80%"
        )
    ) +
geom_text(aes(label= value_text , x = Hugo_Symbol, y = percent_pos ), position=position_dodge(0.9),size=4, vjust=0 ,color="black", face='bold' , family="Helvetica")+
geom_text(aes(label=p_text , x = Hugo_Symbol , y = p_pos ),size=7,family="Helvetica")+
#geom_text(aes(label=Type , x = 18 , y = label_pos ),size=5,family="Helvetica")+
xlab("") +
ylab('Percentage of samples with mutations') +
scale_fill_manual(values = col) +
theme(
    title =element_text(size=4, face='bold'),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.position = c(0.90,0.8) ,
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1 , size = 14  , color = 'black' , family="Helvetica" ) ,
    axis.title.y =  element_text(size = 16 , color = 'black' ) ,
    axis.text.y =  element_text(size = 14 , color = 'black' , family="Helvetica") ,
    axis.ticks.length = unit(0.2, "cm") ,
    panel.grid=element_blank() 
)

width <- 14
height <- 5.3
images_name <- paste0(images_path , "/Fig4B.pdf")
ggsave( images_name , p , width = width , height = height )
