##########################################################################################

library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(optparse)
library(ggpubr)
library(ggsci)

###########################################################################################

col <- c(
  brewer.pal(9,"YlGnBu")[6],
  rgb(234,106,79,alpha=255,maxColorValue=255),
  rgb(203,24,30,alpha=255,maxColorValue=255),
  rgb(255,0,0,alpha=255,maxColorValue=255)
  )

names(col) <- c("IM" , "IGC" , "DGC" , "GC")
col <- col[1:3]
class_order <- c("IGC" , "DGC" , "IM")
molecular_order <- c("GS" , "CIN" , "MSI")

Variant_Types <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation")
base_col <- c("Gender" , "Age_divide" , "Tobacco" , "Alcohol" , "HP")

###########################################################################################

dat_input <- data.frame(fread( input_file ))

###########################################################################################
dat_plot2 <- dat_input %>%
group_by( Patient , Alcohol , Class , Type , MS_Type , TCGA_Class ) %>%
summarize( BurdenAll = median(BurdenAll) ,BurdenExon = median(BurdenExon) )

dat_plot2$Class <- factor( dat_plot2$Class , levels = class_order , order = T )

###########################################################################################

gs_sample <- subset(dat_plot2 , TCGA_Class== "GS")$Patient
cin_sample <- subset(dat_plot2 , TCGA_Class== "CIN")$Patient
msi_sample <- subset(dat_plot2 , TCGA_Class== "MSI")$Patient

###########################################################################################
trans <- function(num){
    up <- floor(log10(num))
    down <- round(num*10^(-up),2)
    text <- paste("P == ",down," %*% 10","^",up)
    return(text)
}

plotBurden <- function(dat_plot_tmp = dat_plot_tmp , out_name = out_name , width = width , height = height ){

    dat_tmp <- c()
    for(class in unique(dat_plot_tmp$Class) ){

        dat_plot_tmp_use <- subset( dat_plot_tmp , Class == class )

        a <- dat_plot_tmp_use[dat_plot_tmp_use$useCol_nums==unique(dat_plot_tmp_use$useCol_nums)[1],"MutBurden_use"]
        b <- dat_plot_tmp_use[dat_plot_tmp_use$useCol_nums==unique(dat_plot_tmp_use$useCol_nums)[2],"MutBurden_use"]

        if(is.na(a)){
            a <- 0
        }

        if(is.na(b)){
            b <- 0
        }
        
        
        p <- wilcox.test( a , b )$p.value

        if( p < 0.01 ){
            p_text <- trans(p)
        }else if( p > 0.05 ){
            p_text <- paste0( "P == " , round(as.numeric(p) , 2) ) 
        }else{
            p_text <- paste0( "P == " , round(as.numeric(p) , 3) ) 
        }

        dat_plot_tmp_use$p_text <- ""
        dat_plot_tmp_use$p_text[1] <- p_text

        dat_tmp <- rbind( dat_plot_tmp_use , dat_tmp )
    }

    y_max <- max(dat_tmp$MutBurden_use) + 0.1
    y_lab <- "Mutation rate per MB"
    
    col_tmp <- c(
        rgb(red=179,green=60,blue=59,alpha=255,max=255) ,
        rgb(red=14,green=90,blue=170,alpha=255,max=255)
    )

    plot <- ggplot( dat_tmp , aes( x = useCol , y = MutBurden_use , color = useCol ) ) +
        geom_boxplot(size = 1.2 , outlier.shape = NA ) + ## 去除散点，加粗线
        facet_grid( .~ Class , scales = "free_x" ) +
        scale_color_manual(values=col_tmp) +
        scale_fill_manual(values =col_tmp) +
        geom_text(aes(label=p_text , y = y_max , x = 1.5),parse = TRUE,size=5 , color = "black", face='bold') +
        xlab(NULL) +
        ylab(y_lab)+
        theme_bw() +
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
            axis.title.y = element_text(size = 14,color="black",face='bold'),
            axis.text.x = element_text(size = 14,color="black",face='bold') ,
            axis.ticks.length = unit(0.2, "cm") ,
            strip.text.x = element_text(size = 17, colour = "black",face='bold') ,
            axis.line = element_line(size = 0.5)) 
    ggsave(file=out_name,plot=plot,width=width,height=height)

}

###########################################################################################
baseuse <- "Alcohol"

dat_plot_tmp <- data.frame(dat_plot2)
dat_plot_tmp$useCol <- dat_plot2[[baseuse]]
dat_plot_tmp <- subset( dat_plot_tmp , !is.na(useCol) )
dat_plot_tmp <- subset( dat_plot_tmp , Type != "IM + IGC + DGC")
dat_plot_tmp <- subset( dat_plot_tmp , Class == "IM" )
dat_plot_tmp$Class <- ifelse( dat_plot_tmp$Patient %in% gs_sample , "GS" , dat_plot_tmp$Class )
dat_plot_tmp$Class <- ifelse( dat_plot_tmp$Patient %in% cin_sample , "CIN" , dat_plot_tmp$Class )
dat_plot_tmp$Class <- ifelse( dat_plot_tmp$Patient %in% msi_sample , "MSI/POLE" , dat_plot_tmp$Class )

dat_plot_tmp$id <- paste0( dat_plot_tmp$Class , "_" , dat_plot_tmp$useCol )
sample_num <- dat_plot_tmp %>%
group_by( id ) %>%
summarize( nums = length(unique( Patient )) )

dat_plot_tmp <- merge( dat_plot_tmp , sample_num , by = "id" )
dat_plot_tmp$useCol_nums <- paste0( dat_plot_tmp$useCol , "\n(" , dat_plot_tmp$nums , ")" )

dat_plot_tmp$useCol <- ifelse( dat_plot_tmp$useCol == "Drink" , "Drinking" , dat_plot_tmp$useCol )

mutType <- "cds"
dat_plot_tmp$MutBurden_use <- dat_plot_tmp$BurdenExon
out_name <- paste0( images_path , "/Fig1.imburden.pdf" ) 

width <- 4.8/1
height <- 4.47/1
plotBurden(dat_plot_tmp = dat_plot_tmp , out_name = out_name , width = width , height = height)
