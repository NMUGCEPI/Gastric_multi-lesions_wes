###################################################################################

library(data.table)
library(optparse)
library(dplyr)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(ggpubr)

##############################################################################

col <- c(
    brewer.pal(9,"YlGnBu")[6],
    rgb(red=247,green=184,blue=71,alpha=255,max=255) ,
    rgb(red=2,green=100,blue=190,alpha=255,max=255) ,
    rgb(255,0,0,alpha=255,maxColorValue=255)
    )

col <- c(
    brewer.pal(9,"YlGnBu")[6],
    rgb(red=179,green=60,blue=59,alpha=255,max=255) ,
    rgb(red=14,green=90,blue=170,alpha=255,max=255) ,
    rgb(255,0,0,alpha=255,maxColorValue=255)
    )

names(col) <- c("IM" , "IGC" , "DGC" , "GC")

##############################################################################

dat <- data.frame(fread(input_file))
dat_all <- dat
dat_all$From <- "All"
dat_all <- rbind(dat_all , dat)
dat_all <- subset( dat_all , Molecular.subtype != "unknown" )

dat <- subset( dat_all , From == "All" )

dat <- subset( dat , Molecular.subtype %in% c("CIN" , "POLE" , "GS" , "MSI") )
dat <- subset( dat , Class!="IGC + DGC" )
dat$Molecular.subtype <- ifelse( dat$Molecular.subtype == "POLE" , "MSI" , dat$Molecular.subtype )

trans <- function(num){
    up <- floor(log10(num))
    down <- round(num*10^(-up),2)
    text <- paste("P == ",down," %*% 10","^",up)
    return(text)
}

dat_tmp <- c()
for( type in unique(unique(dat$Molecular.subtype)) ){

    dat_plot_tmp_use <- subset( dat , Molecular.subtype == type )

    sample_class_num <- dat_plot_tmp_use %>%
    group_by(Class) %>%
    summarize( nums_class = length(unique(Tumor)) )

    dat_plot_tmp_use <- merge( dat_plot_tmp_use , sample_class_num , by = "Class")

    a <- dat_plot_tmp_use[dat_plot_tmp_use$Class==unique(dat_plot_tmp_use$Class)[1],"BurdenExon"]
    b <- dat_plot_tmp_use[dat_plot_tmp_use$Class==unique(dat_plot_tmp_use$Class)[2],"BurdenExon"]
    
    p <- wilcox.test( a , b )$p.value

    if( p < 0.01 ){
        p_text <- trans(p)
    }else{
        p_text <- paste0( "P == " , round(as.numeric(p) , 2) ) 
    }
    
    dat_plot_tmp_use$p_text <- ""
    dat_plot_tmp_use$p_text[1] <- p_text
    dat_plot_tmp_use$Class_use <- paste0( dat_plot_tmp_use$Class , "\n" , "(" , dat_plot_tmp_use$nums_class , ")" )
    dat_tmp <- rbind(dat_tmp , dat_plot_tmp_use)
}

y_max <- max(dat_tmp$BurdenExon) + 280
dat_tmp$Class_use <- factor( dat_tmp$Class_use , levels = unique(dat_tmp$Class_use)[order(unique(dat_tmp$Class_use) , decreasing=T)] , order = T )
dat_tmp$Class <- factor( dat_tmp$Class , levels = unique(dat_tmp$Class)[order(unique(dat_tmp$Class) , decreasing=T)] , order = T )
y_lab <- "Mutation rate per MB"

dat_tmp$Molecular.subtype <- ifelse( dat_tmp$Molecular.subtype == "MSI" , "MSI/POLE" , dat_tmp$Molecular.subtype )
dat_tmp$Molecular.subtype <- factor( dat_tmp$Molecular.subtype , levels = c("CIN" , "GS" , "MSI/POLE") , order = T )

plot <- ggplot(dat_tmp, aes(x = Class , y = BurdenExon, color = Class , fill = Class)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.2,position=position_dodge(0.9),fill="white",color="black")+ 
        scale_y_log10() +
        facet_grid( .~ Molecular.subtype , scales = "free_x" ) +
        scale_color_manual(values=col[c("IGC" , "DGC")]) +
        scale_fill_manual(values = col[c("IGC" , "DGC")]) +
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
    
out_name <- paste0(out_path , "/Fig1.gcburden.pdf")
ggsave(file=out_name,plot=plot,width=4.8/1,height=4.47/1)
