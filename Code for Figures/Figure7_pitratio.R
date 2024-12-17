##########################################################################################

library(Seurat)
library(data.table)
library(optparse)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

##########################################################################################

info_singlecell <- data.frame(fread(singleCell_sample_file))
info <- data.frame(fread(sample_list_file))

sc_dataset <- load(single_cell_file, verbose = F)
sc_dataset_all <- epiall_nor_PCA_50_RE0.5
Idents(sc_dataset_all) <- "sample"   

sc_dataset_scissor <- load(single_cell_scissor_file, verbose = F)
sc_dataset_scissor <- sc_dataset

##########################################################################################

muc6_pit <- names(which(sc_dataset_scissor$scissor=="1" & sc_dataset_scissor$celltype=="Pit"))
other_pit <- names(which(sc_dataset_scissor$scissor!="1" & sc_dataset_scissor$celltype=="Pit"))

##########################################################################################

result_tmp <- c()
sc_dataset <- sc_dataset_all 

sc_dataset <- subset(sc_dataset_all , idents=c(class))
Idents(sc_dataset) <- "celltype"   
sc_dataset <- subset(sc_dataset , idents=c("Pit"))

info_singlecell_use <- info_singlecell

expression_matrix <- as.matrix(GetAssayData(sc_dataset, assay = "RNA"))

index1 <- grep( info_singlecell_use$singlecell_ID[1] , colnames(expression_matrix) )
index2 <- grep( info_singlecell_use$singlecell_ID[2] , colnames(expression_matrix) )
index3 <- grep( info_singlecell_use$singlecell_ID[3] , colnames(expression_matrix) )
index4 <- grep( info_singlecell_use$singlecell_ID[4] , colnames(expression_matrix) )

expression_matrix <- expression_matrix[,c(index1,index2,index3,index4)]

he_t1 <- as.numeric(quantile(expression_matrix["GKN1",])["50%"])
he_t2 <- as.numeric(quantile(expression_matrix["GKN2",])["50%"])

tmp_expression <- expression_matrix[c("GKN1" , "GKN2"),muc6_pit]
high_ratio <- length(which(tmp_expression["GKN1",] > he_t1 & tmp_expression["GKN2",] > he_t2))/ncol(tmp_expression)
tmp <- data.frame( Sample = "MUC6\nMut" , High_Ratio = high_ratio )
result_tmp <- rbind(result_tmp , tmp)

tmp_expression <- expression_matrix[c("GKN1" , "GKN2"),colnames(expression_matrix)[!(colnames(expression_matrix) %in% muc6_pit)]]
high_ratio <- length(which(tmp_expression["GKN1",] > he_t1 & tmp_expression["GKN2",] > he_t2))/ncol(tmp_expression)
tmp <- data.frame( Sample = "Other" , High_Ratio = high_ratio )
result_tmp <- rbind(result_tmp , tmp)

result_tmp_lowratio <- data.frame( Sample = result_tmp$Sample  , High_Ratio = 1 - result_tmp$High_Ratio)
result_tmp_lowratio$Type <- "Low"
result_tmp$Type <- "High"

result_final <- rbind( result_tmp , result_tmp_lowratio )
colnames(result_final)[2] <- "Ratio"

p <- fisher.test(round(matrix(c(result_final$Ratio),ncol=2) * 100))$p.value
trans <- function(num){
    up <- floor(log10(num))
    down <- round(num*10^(-up),2)
    text <- paste("p == ",down," %*% 10","^",up)
    return(text)
}
if( p < 0.01 ){
    p_text <- trans(p)
}else{
    p_text <- paste0( "p == " , round(as.numeric(p) , 3) ) 
}

result_final$value_text <- paste0( round(result_final$Ratio , 2) * 100 , "%") 
result_final$value_text <- paste0( round(result_final$Ratio , 2) * 100 , "%") 

col <- c(
    rgb(red=179,green=34,blue=35,alpha=255,max=255), 
    rgb(red=2,green=100,blue=190,alpha=255,max=255) 
    )

names(col) <- c("High" , "Low" )

p1 <- ggplot(result_final,aes(x=Sample,y=Ratio,fill=factor(Type))) +
	geom_bar(stat="identity") +
	ylab(paste0( "Proportion of cells with high expression \nGKN1 and GKN2 in pit cells")) +
	geom_text(aes(label=value_text) , position=position_stack(vjust = 0.5) , size=4 , color="white")+
	geom_text(aes(label=p_text , y = 1.05 , x = 1.5),parse = TRUE,size=5 , color = "black") +
	xlab(NULL) +
	theme_bw() +
  	theme(panel.background = element_blank(),
        legend.position ='right',
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
        )  +
	scale_fill_manual(values=c(col))

print(result_final)
out_name <- paste0( out_path , "/Fig7H.pdf"  )
ggsave(file=out_name,plot=p1,width=3.5,height=4.5)
