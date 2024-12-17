##########################################################################################

library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(optparse)
library(ggpubr)

###########################################################################################

col <- c(
    brewer.pal(9,"YlGnBu")[6],
    rgb(red=247,green=184,blue=71,alpha=255,max=255) ,
    rgb(red=2,green=100,blue=190,alpha=255,max=255) ,
    rgb(255,0,0,alpha=255,maxColorValue=255)
    )

names(col) <- c("IM" , "IGC" , "DGC" , "GC")
col_im <- brewer.pal(9,"YlGnBu")[6:8]
names(col_im) <- c("IM(IGC)" , "IM(DGC)" , "IM(IGC_DGC)")

Variant_Types <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation")

###########################################################################################

dat_sample <- data.frame(fread( info_file ))
dat_maf <- data.frame(fread( maf_file ))

###########################################################################################

dat_maf <- subset( dat_maf , Variant_Classification %in% Variant_Types )
exon_mut <- dat_maf %>%
group_by( Tumor_Sample_Barcode ) %>%
summarize( exonMutNum = length(Start_position) )

dat_plot <- merge( dat_sample , exon_mut , by.x = "Tumor" , by.y = "Tumor_Sample_Barcode" , all.x = TRUE )
dat_plot$BurdenAll <- dat_plot$mutNum/(dat_plot$coverage_All/1024/1024)
dat_plot$BurdenExon <- dat_plot$exonMutNum/(dat_plot$coverage_CDS/1024/1024)

###########################################################################################
dat_plot2 <- dat_plot %>%
group_by( Patient , Class , Type ) %>%
summarize( BurdenAll = median(BurdenAll) ,BurdenExon = median(BurdenExon) )

###########################################################################################

trans <- function(num){
    up <- floor(log10(num))
    down <- round(num*10^(-up),2)
    text <- paste("P == ",down," %*% 10","^",up)
    return(text)
}

###########################################################################################
dat_plot_tmp <- data.frame(dat_plot2)
dat_plot_tmp$Class_use <- ifelse(dat_plot_tmp$Class == "IM" , "IM" , "GC")
dat_plot_tmp$Class_use <- factor( dat_plot_tmp$Class_use , levels = c("IM" , "GC") , order = T )
dat_plot_tmp$Type <- factor( dat_plot_tmp$Type , levels = c("IM + IGC" , "IM + DGC") , order = T )
my_comparisons_1 <- list( c(1, 2) )

dat_plot_tmp_use <- dat_plot_tmp
type_num <- dat_plot_tmp_use %>%
group_by( Type ) %>%
summarize( type_nums = length(unique(Patient)) )
dat_plot_tmp_use <- merge( dat_plot_tmp_use , type_num , by = "Type" )

dat_tmp <- c()

y_max <- max(dat_plot_tmp_use$BurdenExon) + 3

for( type in unique(dat_plot_tmp_use$Type) ){

    dat <- subset( dat_plot_tmp_use , Type == type )
    dat <- dat[order(dat$Patient),]

    a <- dat[dat$Class==unique(dat$Class)[1],"BurdenExon"]
    b <- dat[dat$Class==unique(dat$Class)[2],"BurdenExon"]
    
    p <- wilcox.test( a , b , paired = T )$p.value

    if( p < 0.01 ){
        p_text <- trans(p)
    }else{
        p_text <- paste0( "P == " , round(as.numeric(p) , 3) ) 
    }
    dat$p_text <- ""
    dat$p_text[1] <- p_text
    dat_tmp <- rbind(dat_tmp , dat)
}
dat_tmp$Class <- factor( dat_tmp$Class , levels = c("IM" , "IGC" , "DGC") , order = T )
dat_tmp$Type <- factor( dat_tmp$Type , levels = c("IM + IGC" , "IM + DGC") , order = T )

col_tmp <- c(
    rgb(red=247,green=184,blue=71,alpha=255,max=255) ,
    rgb(red=2,green=100,blue=190,alpha=255,max=255) ,
    rgb(red=2,green=100,blue=190,alpha=255,max=255) 
)

names(col_tmp) <- c("IM" , "IGC" , "DGC")

plot <- ggplot( dat_tmp , aes( x = Class , y = BurdenExon , color = Class , fill = Class ) ) +
    geom_line( aes( group = Patient ) , size = 0.4 , color = "gray" ) +  ## 配对样本加线
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.2,position=position_dodge(0.9),fill="white",color="black")+ #绘制箱线图
    scale_color_manual(values=col_tmp) +
    scale_fill_manual(values=col_tmp) +
    facet_grid(.~Type,space='free_x',scales='free_x') +
    xlab(NULL) +
    ylab("Mutation rate per MB")+
    theme_bw() +
    geom_text(aes(label=p_text , y = y_max , x = 1.5),parse = TRUE,size=5 , color = "black" , face='bold') +
    theme(
        legend.position = 'none',
        legend.title = element_blank() ,
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 12,color="black",face='bold'),
        legend.text = element_text(size = 12,color="black",face='bold'),
        axis.text.y = element_text(size = 15,color="black",face='bold'),
        axis.title.x = element_text(size = 15,color="black",face='bold'),
        axis.title.y = element_text(size = 12,color="black",face='bold'),
        axis.text.x = element_text(size = 15,color="black",face='bold') ,
        axis.ticks.length = unit(0.2, "cm") ,
        strip.text.x = element_text(size = 15, colour = "black",face='bold') ,
        axis.line = element_line(size = 0.5)) 
out_name <- paste0( images_path , "Fig3A.pdf" )  
ggsave(file=out_name,plot=plot,width=4.7/1.2,height=4.8/1.2)

###########################################################################################

dat_sample <- data.frame(fread( info_file ))

burden_type <- "All"
dat_plot2 <- dat_sample %>%
    group_by( Patient , Class , Type , TCGA_Class ) %>%
    summarize( BurdenExon = median(BurdenExon) )

msi_sample <- unique(subset( dat_sample , TCGA_Class == "MSI" )$Patient)
gs_sample <- unique(subset( dat_sample , TCGA_Class == "GS" )$Patient)
cin_sample <- unique(subset( dat_sample , TCGA_Class == "CIN" )$Patient)

trans <- function(num){
    up <- floor(log10(num))
    down <- round(num*10^(-up),2)
    text <- paste("p == ",down," %*% 10","^",up)
    return(text)
}

plotBurden <- function(dat_plot_tmp = dat_plot_tmp , type = type ){

    dat_plot_tmp$title <- type
    dat_plot_tmp$TCGA_Class <- ifelse( dat_plot_tmp$Patient %in% msi_sample , "MSI/POLE" , dat_plot_tmp$TCGA_Class)
    dat_plot_tmp$TCGA_Class <- ifelse( dat_plot_tmp$Patient %in% gs_sample , "GS" , dat_plot_tmp$TCGA_Class)
    dat_plot_tmp$TCGA_Class <- ifelse( dat_plot_tmp$Patient %in% cin_sample , "CIN" , dat_plot_tmp$TCGA_Class)
    dat_plot_tmp$TCGA_Class <- factor( dat_plot_tmp$TCGA_Class , levels = c("MSI/POLE" , "CIN" , "GS") , order = T )
    sample_num <- dat_plot_tmp %>%
    group_by( TCGA_Class ) %>%
    summarize( nums = length(unique(Patient)) )

    dat_plot_tmp <- merge( dat_plot_tmp , sample_num , by = "TCGA_Class" )
    dat_plot_tmp$TCGA_Class_num <- paste0(dat_plot_tmp$TCGA_Class , "\n" , "(" , dat_plot_tmp$nums , ")" )

    my_comparisons_1 <- list( c(1,2) , c(1,3) , c(2,3))


    dat_plot_tmp$TCGA_Class <- factor( dat_plot_tmp$TCGA_Class , levels = c( "GS" , "CIN" , "MSI/POLE" ) , order = T )

    col_use <- c(rgb(red=179,green=34,blue=35,alpha=255,max=255) ,
        rgb(red=247,green=184,blue=71,alpha=255,max=255) ,
        rgb(red=2,green=100,blue=190,alpha=255,max=255) 
    )

    names(col_use) <- c("GS" , "MSI/POLE" , "CIN")

    plot <- ggplot( dat_plot_tmp , aes( x = TCGA_Class , y = BurdenExon , color = TCGA_Class ) ) +
        geom_line( aes( group = Patient ) , size = 0.4 , color = "gray" ) +  
        geom_boxplot(alpha =1 , size = 0.9 , width = 0.6 , outlier.shape = NA) +
        geom_jitter(position = position_jitter(0.2) , size = 1 , alpha = 1) +
        scale_color_manual(values = col_use) +
        facet_grid(.~title) +
        scale_y_continuous(
                trans = sqrt_trans(),
                breaks = c(1:5)
                ) + 
        xlab(NULL) +
        ylab("Mutation rate per MB")+
        theme_bw() +
        stat_compare_means(comparisons = my_comparisons_1,method = "wilcox.test") +
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

type <- "IM"
dat_plot_tmp <- data.frame(dat_plot2)
p1 <- plotBurden(dat_plot_tmp = dat_plot_tmp , type = type )
out_name <- paste0( images_path , "/Fig3C1.pdf" )  
ggsave(file=out_name,plot=p1,width=3.8/1.2,height=4.4/1.2)

##############################################################################

info <- data.frame(fread(input_file))
dat_share <- data.frame(fread(mutshare_file))
dat_share_msi <- data.frame(fread(mutshare_msi_file))
dat_ccf <- data.frame(fread(ccf_file))
dat_ccf_msi <- data.frame(fread(ccf_msi_file))

dat_ccf_use <- rbind( dat_ccf , dat_ccf_msi )
dat_share_use <- rbind( dat_share , dat_share_msi )

dat_share_rate <- c() 
for( sample in unique(dat_share_use$Tumor) ){

    tmp <- subset( dat_share_use , Tumor == sample )
    tmp_ccf <- subset( dat_ccf_use , Sample == sample )

    tmp$vid <- paste( tmp$Chromosome , tmp$Start_Position , tmp$End_Position , sep = ":" )
    tmp_ccf$vid <- paste( tmp_ccf$Chr , tmp_ccf$Start_Position , tmp_ccf$End_Position , sep = ":" )

    tmp_use <- merge( tmp , tmp_ccf[,c("vid" , "VAF" , "CCF_adj" , "CLS")] , by = "vid" )

    share_nums <- median(subset( tmp_use , Share=="TRUE" )$CCF_adj)
    private_nums <- median(subset( tmp_use , Share=="FALSE" )$CCF_adj)
        
    dat_tmp <- data.frame( Tumor = sample , share_ccf = share_nums , private_ccf = private_nums )
    dat_share_rate <- rbind( dat_share_rate , dat_tmp )
}

dat <- merge( info , dat_share_rate , by = "Tumor" )

msi_sample <- unique(subset( dat , TCGA_Class %in% c("POLE" , "MSI") )$Patient)
gs_sample <- unique(subset( dat , TCGA_Class %in% c("GS") )$Patient)
cin_sample <- unique(subset( dat , TCGA_Class %in% c("CIN") )$Patient)

dat_use <- subset( dat , TCGA_Class == "IM" )
dat_use$Molecular <- ifelse( dat_use$Patient %in% msi_sample , "MSI" , "" )
dat_use$Molecular <- ifelse( dat_use$Patient %in% gs_sample , "GS" , dat_use$Molecular )
dat_use$Molecular <- ifelse( dat_use$Patient %in% cin_sample , "CIN" , dat_use$Molecular )

plotBurden <- function(dat_plot_tmp = dat_plot_tmp , type = type , share_class = share_class ){

    dat_plot_tmp$title <- type
    dat_plot_tmp$TCGA_Class <- ifelse( dat_plot_tmp$Patient %in% msi_sample , "MSI\nPOLE" , dat_plot_tmp$TCGA_Class)
    dat_plot_tmp$TCGA_Class <- ifelse( dat_plot_tmp$Patient %in% gs_sample , "GS" , dat_plot_tmp$TCGA_Class)
    dat_plot_tmp$TCGA_Class <- ifelse( dat_plot_tmp$Patient %in% cin_sample , "CIN" , dat_plot_tmp$TCGA_Class)
    dat_plot_tmp$TCGA_Class <- factor( dat_plot_tmp$TCGA_Class , levels = c("GS" , "CIN" , "MSI\nPOLE") , order = T )

    sample_num <- dat_plot_tmp %>%
    group_by( TCGA_Class ) %>%
    summarize( nums = length(unique(Patient)) )

    dat_plot_tmp <- merge( dat_plot_tmp , sample_num , by = "TCGA_Class" )
    dat_plot_tmp$TCGA_Class_num <- paste0(dat_plot_tmp$TCGA_Class , "\n" , "(" , dat_plot_tmp$nums , ")" )

    dat_plot_tmp_use <- dat_plot_tmp %>%
    group_by( TCGA_Class_num , Patient , title , TCGA_Class ) %>%
    summarize( share_ccf = median(share_ccf) , private_ccf = median(private_ccf) )

    dat_plot_tmp_use$useCol <- dat_plot_tmp_use$share_ccf
    ylab <- "Median CCF of share mutations"
    ymax <- 0.11

    my_comparisons_1 <- list( c(1,2) , c(2,3) , c(1,3))
    col_use <- c(rgb(red=179,green=34,blue=35,alpha=255,max=255) ,
        rgb(red=247,green=184,blue=71,alpha=255,max=255) ,
        rgb(red=2,green=100,blue=190,alpha=255,max=255) 
    )

    names(col_use) <- c("GS" , "MSI\nPOLE" , "CIN")

   
    ylab <- ""
    plot <- ggplot( dat_plot_tmp_use , aes( x = TCGA_Class , y = useCol , color = TCGA_Class ) ) +
        geom_line( aes( group = Patient ) , size = 0.4 , color = "gray" ) +  
        geom_boxplot(alpha =1 , size = 0.9 , width = 0.6 , outlier.shape = NA) +
        geom_jitter(position = position_jitter(0.2) , size = 1 , alpha = 1) +
        scale_color_manual(values = col_use) +
        ylim(0,ymax) +
        facet_grid(.~title) +
        xlab(NULL) +
        ylab(ylab)+
        theme_bw() +
        stat_compare_means(comparisons = my_comparisons_1,method = "wilcox.test") +
        theme(
            legend.position = 'none',
            legend.title = element_blank() ,
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            plot.title = element_text(size = 12,color="black",face='bold'),
            legend.text = element_text(size = 12,color="black",face='bold'),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 12,color="black",face='bold'),
            axis.title.y = element_text(size = 12,color="black",face='bold'),
            axis.text.x = element_text(size = 12,color="black",face='bold') ,
            axis.ticks.length = unit(0.2, "cm") ,
            strip.text.x = element_text(size = 15, colour = "black",face='bold') ,
            axis.line = element_line(size = 0.5)
        )
    
    return(plot)
}

share_class <- "Share"
type <- "All"
dat_plot_tmp <- data.frame(dat_use)
p1 <- plotBurden(dat_plot_tmp = dat_plot_tmp , type = type , share_class )
out_name <- paste0( images_path , "/Fig3C2.pdf" )  
ggsave(file=out_name,plot=p1,width=6/1.2,height=4.3/1.2)

type <- "Drinking"
dat_plot_tmp <- data.frame(dat_use)
dat_plot_tmp <- subset( dat_plot_tmp , Class =="IM" & Alcohol == "Drink")
p2 <- plotBurden(dat_plot_tmp = dat_plot_tmp , type = type , share_class )
out_name <- paste0( images_path , "/Fig3D1.pdf" )  
ggsave(file=out_name,plot=p2,width=6/1.2,height=4.3/1.2)


##############################################################################


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

dat_input <- data.frame(fread( input_file ))
dat_plot2 <- dat_input %>%
group_by( Patient , Alcohol , Class , Type , TCGA_Class ) %>%
summarize( BurdenAll = median(BurdenAll) ,BurdenExon = median(BurdenExon) )
dat_plot2$Class <- factor( dat_plot2$Class , levels = class_order , order = T )

gs_sample <- subset(dat_plot2 , TCGA_Class== "GS")$Patient
cin_sample <- subset(dat_plot2 , TCGA_Class== "CIN")$Patient
msi_sample <- subset(dat_plot2 , TCGA_Class== "MSI")$Patient

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

    y_max <- max(dat_tmp$MutBurden_use) + 2

    col_tmp <- c(
        rgb(red=247,green=184,blue=71,alpha=255,max=255) ,
        rgb(red=2,green=100,blue=190,alpha=255,max=255) 
    )

    plot <- ggplot( dat_tmp , aes( x = useCol_nums , y = MutBurden_use , color = useCol ) ) +
        geom_boxplot(alpha =1 , outlier.size=0 , size = 1.5 , width = 0.6) +
        geom_jitter(position = position_jitter(0.17) , size = 1 , alpha = 1) +
        facet_grid(.~Class,space='free_x',scales='free_x') +
        scale_color_manual(values=col_tmp) +
        scale_y_sqrt() + 
        xlab(NULL) +
        ylab("Mutation rate per MB")+
        geom_text(aes(label=p_text , y = y_max , x = 1.5),parse = TRUE,size=4.5 , color = "black") +
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
            axis.title.y = element_text(size = 12,color="black",face='bold'),
            axis.text.x = element_text(size = 12,color="black",face='bold') ,
            axis.ticks.length = unit(0.2, "cm") ,
            strip.text.x = element_text(size = 15, colour = "black",face='bold') ,
            axis.line = element_line(size = 0.5)) 
    ggsave(file=out_name,plot=plot,width=width,height=height)

}

baseuse <- "Alcohol"
dat_plot_tmp <- data.frame(dat_plot2)
dat_plot_tmp$useCol <- dat_plot2[[baseuse]]
dat_plot_tmp <- subset( dat_plot_tmp , !is.na(useCol) )
dat_plot_tmp <- subset( dat_plot_tmp , Class == "IM" )
dat_plot_tmp$Class <- ifelse( dat_plot_tmp$Patient %in% gs_sample , "GS" , dat_plot_tmp$Class )
dat_plot_tmp$Class <- ifelse( dat_plot_tmp$Patient %in% cin_sample , "CIN" , dat_plot_tmp$Class )
dat_plot_tmp$Class <- ifelse( dat_plot_tmp$Patient %in% msi_sample , "MSI/POLE" , dat_plot_tmp$Class )

dat_plot_tmp[dat_plot_tmp$useCol=="Drink","useCol"] <- "Ever"
dat_plot_tmp[dat_plot_tmp$useCol=="No","useCol"] <- "Never"
dat_plot_tmp$useCol <- factor( dat_plot_tmp$useCol , levels = c("Ever" , "Never") , order = T )

dat_plot_tmp$id <- paste0( dat_plot_tmp$Class , "_" , dat_plot_tmp$useCol )
sample_num <- dat_plot_tmp %>%
group_by( id ) %>%
summarize( nums = length(unique( Patient )) )

dat_plot_tmp <- merge( dat_plot_tmp , sample_num , by = "id" )
dat_plot_tmp$useCol_nums <- paste0( dat_plot_tmp$useCol , "\n(" , dat_plot_tmp$nums , ")" )
dat_plot_tmp$MutBurden_use <- dat_plot_tmp$BurdenExon
out_name <- paste0( images_path , "/Fig3D2.pdf" )  

width <- 6.5/1.2
height <- 4.0/1.2
plotBurden(dat_plot_tmp = dat_plot_tmp , out_name = out_name , width = width , height = height)

