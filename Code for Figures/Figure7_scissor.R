##########################################################################################

library(dplyr)
library(data.table)
library(optparse)
library(ggplot2)
library(Scissor)

###########################################################################################

dat_im <- data.frame(fread(input_im_file))

cell_order <- dat_im %>%
group_by( Cell_Type ) %>%
summarize( Cells = sum(Cells) )

tmp <- load(single_cell_file)

###########################################################################################

cell_order <- unique(data.frame(cell_order)[order(cell_order$Cells , decreasing=T),"Cell_Type"])

tmp_data <- dat_im
tmp_data$Scissor_Type <- as.character(tmp_data$Scissor_Type)
tmp_data$Cell_Type <- factor( tmp_data$Cell_Type , levels = cell_order , order = T )
col <- c('grey','indianred1','royalblue')
names(col) <- c("Background" , "Scissor+" , "Scissor-" )

p1 <- ggplot(tmp_data,aes(x=Cell_Type,y=Cells,fill=factor(Scissor_Type))) +
    geom_bar(stat="identity") +
    ylab("Cell counts") +
    xlab(NULL) +
    theme_bw() +
    theme(panel.background = element_blank(),
        legend.title = element_blank() ,
        panel.grid.major=element_line(colour=NA),
        legend.text = element_text(size = 10,color="black",face='bold'),
        legend.position = c(0.7,0.7) ,
        axis.text.x = element_text(size = 12,color="black",face='bold' , angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10,color="black",face='bold'),
        axis.title.x = element_text(size = 10,color="black",face='bold'),
        axis.title.y = element_text(size = 13,color="black",face='bold'),
        strip.text.x = element_text(size = 15,color="black",face='bold'),
        axis.ticks.length = unit(0.2, "cm") ,
        axis.line = element_line(size = 0.5))  +
    scale_fill_manual(values=c(col))

out_name <- paste0( out_path , "/Fig7E.pdf"  )
ggsave(file=out_name,plot=p1,width=4,height=4)


###########################################################################################
p1 <- DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1)) +
    theme(panel.background = element_blank(),
        legend.title = element_blank() ,
        panel.grid.major=element_line(colour=NA),
        legend.text = element_text(size = 12,color="black",face='bold'),
        legend.position = 'right' ,
        axis.text.x = element_text(size = 10,color="black",face='bold'),
        axis.text.y = element_text(size = 10,color="black",face='bold'),
        axis.title.x = element_text(size = 10,color="black",face='bold'),
        axis.title.y = element_text(size = 13,color="black",face='bold'),
        strip.text.x = element_text(size = 15,color="black",face='bold'),
        axis.ticks.length = unit(0.2, "cm") ,
        axis.line = element_line(size = 0.5))

p2 <- DimPlot(sc_dataset, reduction = 'umap', group.by = 'celltype' , pt.size = 1.2, order = c(2,1)) +
    theme(panel.background = element_blank(),
        legend.title = element_blank() ,
        panel.grid.major=element_line(colour=NA),
        legend.text = element_text(size = 12,color="black",face='bold'),
        legend.position = 'right' ,
        axis.text.x = element_text(size = 10,color="black",face='bold'),
        axis.text.y = element_text(size = 10,color="black",face='bold'),
        axis.title.x = element_text(size = 10,color="black",face='bold'),
        axis.title.y = element_text(size = 13,color="black",face='bold'),
        strip.text.x = element_text(size = 15,color="black",face='bold'),
        axis.ticks.length = unit(0.2, "cm") ,
        axis.line = element_line(size = 0.5))

plot <- p1 + p2
out_name <- paste0( out_path , "/Fig7D.pdf"  )
ggsave(file=out_name,plot=plot,width=10,height=4.5)
