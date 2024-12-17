##########################################################################################

library(data.table)
library(optparse)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)
library("scales")

##########################################################################################

info <- data.frame(fread(sample_list_file))
info_public <- data.frame(fread(sample_list_public_file))
dat_expression <- data.frame(fread(rsem_file))
dat_maf_public <- data.frame(fread( maf_public_file ))
use_gene <- gene

##########################################################################################

Variant_Types <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation")
check_gene <- c("GKN1" , "GKN2")

dat_expression <- subset( dat_expression , gene_id %in% check_gene)
dat_maf_public <- subset( dat_maf_public , Hugo_Symbol == gene & Variant_Classification %in% Variant_Types )
colnames(dat_expression) <- gsub( "[.]" , "-" , colnames(dat_expression) )

##########################################################################################

info_public <- subset( info_public , From != "NJMU" )
info_public$ID <- info_public$Tumor
info_public$Class_sub <- info_public$Class
info_use <- rbind( info_public[,c( "ID" , "Tumor" , "Class" , "Class_sub" , "From")] , info[,c( "ID" , "Tumor" , "Class" , "Class_sub" , "From")] )

mutTumor <- unique(dat_maf_public$Tumor_Sample_Barcode)
info_mut <- subset( info_use , Tumor %in% mutTumor )
info_mut <- paste0(info_mut$ID , "_" , info_mut$Class_sub)
info_mut <- info_mut[info_mut %in% colnames(dat_expression)]

info_wild <- subset( info_use , !(Tumor %in% mutTumor) )
info_wild <- paste0(info_wild$ID , "_" , info_wild$Class_sub)
info_wild <- info_wild[info_wild %in% colnames(dat_expression)]

##########################################################################################
tmm_combine_final <- c()

for( gene in unique(dat_expression$gene_id) ) {
    tmp_exp <- subset( dat_expression , gene_id == gene )
    
    type <- paste0(use_gene,"_Mut")
    use_sample <- info_mut
    exp_use <- tmp_exp[use_sample]
    tmm_mut <- getTMM(exp_use = exp_use , type = type)

    type <- paste0(use_gene,"_Wild")
    use_sample <- info_wild
    exp_use <- tmp_exp[use_sample]
    tmm_wild <- getTMM(exp_use = exp_use , type = type)

    tmm_combine <- rbind( tmm_wild , tmm_mut )
    tmm_combine$Class <- factor( tmm_combine$Class , levels = c( "IM" , "IGC" , "DGC" ) , order = T )
    tmm_combine$gene <- gene
    tmm_combine_final <- rbind( tmm_combine_final , tmm_combine )
}

tmm_combine_final <- subset( tmm_combine_final , Class == "IM" )
tmm_combine <- tmm_combine_final

y_tmm <- max(tmm_combine$TMM) * 1.1
tmp_dat_use <- tmm_combine
width <- 5

my_comparisons_1 <- list( c(1, 2) )

y_max <- y_tmm
if(y_max > 2000){

    if(y_max > 10000){
        y_breaks <- 10000
    }else if(y_max > 5000){
        y_breaks <- 5000
    }else{
        y_breaks <- 2000
    }
    y_lab_lim <- scale_y_continuous( limits = c(-0.1 ,y_tmm*1.1) , breaks = c( 0 , 500 , 1000 , 2000 , 5000 , y_breaks , 20000 ) , trans = sqrt_trans() )
}else{
    y_lab_lim <- scale_y_continuous( limits = c(-0.1 ,y_tmm*1.1))
}

col <- c(
        rgb(red=179,green=60,blue=59,alpha=255,max=255) ,
        rgb(red=14,green=90,blue=170,alpha=255,max=255)
    )

tmm_combine$Type <- as.character(tmm_combine$Type)
tmm_combine$Type <- gsub( "[_]" , "\n" ,tmm_combine$Type )

trans <- function(num){
    up <- floor(log10(num))
    down <- round(num*10^(-up),2)
    text <- paste("p == ",down," %*% 10","^",up)
    return(text)
}


dat_tmp <- c()
for(geneN in unique(tmm_combine$gene) ){

    dat_plot_tmp <- subset( tmm_combine , gene == geneN )

    a <- dat_plot_tmp[dat_plot_tmp$Type==unique(dat_plot_tmp$Type)[1],"TMM"]
    b <- dat_plot_tmp[dat_plot_tmp$Type==unique(dat_plot_tmp$Type)[2],"TMM"]
    p <- wilcox.test( a , b )$p.value

    if( p < 0.01 ){
        p_text <- trans(p)
    }else{
        p_text <- paste0( "p == " , round(as.numeric(p) , 3) ) 
    }
    dat_plot_tmp$p_text <- ""
    dat_plot_tmp$p_text[1] <- p_text
    dat_tmp <- rbind( dat_plot_tmp , dat_tmp )
}


plot <- ggplot( dat_tmp , aes( x = Type , y = TMM , color = Type ) ) +
    geom_boxplot(lwd=1.5, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.2)) +
    scale_color_manual(values=col) +
    facet_grid(.~gene)+
    xlab(NULL) +
    ylab("TMM")+
    theme_bw() +
    y_lab_lim +
    geom_text(aes(label=p_text , y = y_tmm , x = 1.5),parse = TRUE,size=6 , color = "black" , face = "bold") +
    theme(
      legend.position = 'none',
      legend.title = element_blank() ,
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 12,color="black",face='bold'),
      legend.text = element_text(size = 12,color="black",face='bold'),
      axis.text.y = element_text(size = 12,color="black",face='bold'),
      axis.title.x = element_text(size = 12,color="black",face='bold'),
      axis.title.y = element_text(size = 12,color="black",face='bold'),
      axis.text.x = element_text(size = 12,color="black",face='bold') ,
      axis.ticks.length = unit(0.2, "cm") ,
      strip.text.x = element_text(size = 15, colour = "black",face='bold') ,
      axis.line = element_line(size = 0.5)
        ) 

out_name <- paste0( image_path , "/Fig7I.pdf" )
ggsave(file=out_name,plot=plot,width=3.5,height=5)
