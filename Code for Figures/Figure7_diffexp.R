##########################################################################################
library(data.table)
library(optparse)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggsci)
library(ggrastr)

##########################################################################################

data <- data.frame(fread(input_file))
data2 <- data
colnames(data)[3] <- "-log10_P_Value"

##########################################################################################

q_t <- 0.05
foldchange_t <- 1
logFC_cutoff <- log2(foldchange_t)
log10_P_Value_cutoff <- -log10(q_t)

data$gene <- ifelse( data$gene %in% c("GKN1" , "GKN2" , "DSC2" , "DSG2") , data$gene , "" )

col <- c('#376B6D','#9F353A' , '#F7C242','#BDC0BA' )
names(col) <- c('DOWN', 'UP', 'TARGET', "GREY" )

alpha <- 1

options(ggrepel.max.overlaps = 200)
plot1 <- ggplot(data = data,aes(x = log_FC,y = `-log10_P_Value`))+
        geom_point(data = subset(data,abs(log_FC)<logFC_cutoff),
                aes(size = abs(log_FC)),col = col["GREY"],alpha = alpha)+
        geom_point(data = subset(data,abs(`-log10_P_Value`)<log10_P_Value_cutoff & abs(log_FC)>logFC_cutoff),
                aes(size = abs(log_FC)),col = col["GREY"],alpha = alpha)+
        geom_point(data = subset(data,abs(`-log10_P_Value`)>log10_P_Value_cutoff & log_FC>logFC_cutoff),
                aes(size = abs(log_FC)),col = col["UP"],alpha = alpha)+
        geom_point(data = subset(data,abs(`-log10_P_Value`)>log10_P_Value_cutoff & log_FC< -logFC_cutoff),
                aes(size = abs(log_FC)),col = col["DOWN"],alpha = alpha)+
        geom_point(data = subset(data,gene!=""),
                aes(size = abs(log_FC)),col = col["TARGET"],alpha = alpha)+
        theme_bw()+
        labs(x='log fold change\n(Scissor+ cells/All other cells)',y='-log10(adjusted p-value)')+
        geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff),lty = 3,col = 'black',lwd = 0.4)+
        geom_hline(yintercept = log10_P_Value_cutoff,lty = 3,col = 'black',lwd = 0.4) +
        geom_text_repel(data = subset(data,abs(`-log10_P_Value`)>log10_P_Value_cutoff & abs(log_FC)>logFC_cutoff), 
            fontface="bold", nudge_x = 0.65, nudge_y = 1, 
            aes(label = gene),size = 5.5,col = 'black' , face = 'bold') +
        theme(panel.background = element_blank(),
                legend.position ='none',
                legend.title = element_blank() ,
                panel.grid.major=element_line(colour=NA),
                legend.text = element_text(size = 12,color="black",face='bold'),
                axis.text.x = element_text(size = 12,color="black",face='bold'),
                axis.text.y = element_text(size = 10,color="black",face='bold'),
                axis.title.x = element_text(size = 14,color="black",face='bold'),
                axis.title.y = element_text(size = 14,color="black",face='bold'),
                strip.text.x = element_text(size = 12,color="black",face='bold'),
                axis.ticks.length = unit(0.2, "cm") ,
                axis.line = element_line(size = 0.5)
                )

image_name <- paste0( out_path , "/Fig7G.pdf" )     
ggsave( image_name , plot1 , width = 5.5 , height = 6 )
