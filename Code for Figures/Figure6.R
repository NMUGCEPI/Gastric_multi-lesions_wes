##########################################################################################

library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(optparse)
library(ggpubr)
library(ggsci)
library(ggrepel)

###########################################################################################

dat_ccf <- fread( ccf_file )
dat_input <- data.frame(fread(input_file))
dat_info <- data.frame(fread(sample_info))
dat_base_info <- data.frame(fread(base_info_file))
dat_gene <- data.frame(fread(gene_list , header = T))
colnames(dat_gene) <- "Gene_Symbol"

###########################################################################################

Variant_Type <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation")

dat_ccf$Location <- paste( dat_ccf$Chr , dat_ccf$Start_Position , 
    dat_ccf$REF , dat_ccf$ALT , sep=":"  )

###########################################################################################

dat_base_info <- dat_base_info[,c("Patient" , "CNV_Type")]
dat_base_info <- subset(dat_base_info , CNV_Type!="")
dat_input <- unique(merge( dat_input , dat_base_info[,c("Patient" , "CNV_Type")] , by.x = "ID" , by.y = "Patient" ))

###########################################################################################
result_driver <- c()

for( Sample in unique(dat_input$ID) ){

    tumors <- subset( dat_info , ID == Sample )$Tumor

    tmp <- subset( dat_ccf , Sample %in% tumors & Variant_Classification %in% Variant_Type & Hugo_Symbol %in% dat_gene$Gene_Symbol )

    tmp_driver <- tmp %>% 
    group_by(Location) %>%
    summarize( MutTumor = length(Sample) )

    if( nrow(tmp_driver) == 0){
        class <- "NoDriver"
        tmp_res <- data.frame( Normal = Sample , DriverClass = class , 
        Share_gene = "" , Priver_gene = "" , Share_gene_variant = "" , Private_gene_variant = "") 

    }else if( nrow(tmp_driver) > 0 ){
        share_driver <- length(which(tmp_driver$MutTumor == length(tumors)))
        tmp1 <- tmp
        tmp_driver1 <- tmp1 %>% 
        group_by(Location ) %>%
        summarize( MutClass = paste0(Class , collapse = "_") )
        share_driver1 <- length(which(tmp_driver1$MutClass %in% c("IM_IGC" , "IM_DGC" , "IGC_IM" , "DGC_IM")))

        if(share_driver > 0 ){
            class <- "ShareDriver"
        }else if(share_driver1 > 0){
            class <- "PartDriver"            
        }else{
            class <- "PrivateDriver"
        }

        share_gene <- unique(subset( tmp , Location %in% subset( tmp_driver , MutTumor == length(tumors) )$Location )$Hugo_Symbol)
        private_gene <- unique(subset( tmp , Location %in% subset( tmp_driver , MutTumor < length(tumors) )$Location )$Hugo_Symbol)

        share_gene_class <- subset( tmp , Location %in% subset( tmp_driver , MutTumor == length(tumors) )$Location )$Variant_Classification
        private_gene_class <- subset( tmp , Location %in% subset( tmp_driver , MutTumor < length(tumors) )$Location )$Variant_Classification

        tmp_res <- data.frame( Normal = Sample , DriverClass = class , 
        Share_gene = paste0(share_gene , collapse = ",") , Priver_gene = paste0(private_gene , collapse = ",") ,
        Share_gene_variant = paste0(share_gene_class , collapse = ",") , Private_gene_variant = paste0(private_gene_class , collapse = ",") 
        )
    }

   
    result_driver <- rbind( result_driver , tmp_res )
}

###########################################################################################

dat_input2 <- merge( dat_input , result_driver , by.x = "ID" , by.y = "Normal" )
dat_input2$plotType <- ifelse( dat_input2$Share_gene != "" , "TrunkDriver" , "Other" )
dat_input2$plotType <- factor( dat_input2$plotType , levels = c("TrunkDriver" , "Other") , order = T )
dat_input2$Class <- factor( dat_input2$Class , levels = c("IGC" , "DGC") , order = T )

trans <- function(num){
    up <- floor(log10(num))
    down <- round(num*10^(-up),2)
    text <- paste("P == ",down," %*% 10","^",up)
    return(text)
}

col_tmp <- c(
    rgb(red=179,green=60,blue=59,alpha=255,max=255) ,
    rgb(red=14,green=90,blue=170,alpha=255,max=255)
)


###########################################################################################
p <- wilcox.test(subset(dat_input2 , Class=="IGC")$mean , subset(dat_input2 , Class=="DGC")$mean)$p.value

if( p < 0.01 ){
    p_text <- trans(p)
}else{
    p_text <- paste0( "P == " , round(as.numeric(p) , 3) ) 
}

dat_input2$p_text <- ""
dat_input2$p_text[1] <- p_text

y_max <- 20 + 1
y_lab <- "Number of years for \nIM to GC progression"

plot <- ggplot( dat_input2 , aes( x = Class , y = mean , color = Class ) ) +
    geom_boxplot(size = 1.2 , outlier.shape = NA ) + 
    scale_color_manual(values=col_tmp) +
    scale_fill_manual(values =col_tmp) +
    geom_jitter(position=position_jitter(0.2),aes(color=Class)) +
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

width <- 3/1
height <- 4.47/1
out_name <- paste0(out_path , "/Fig6A.pdf")
ggsave(file=out_name,plot=plot,width=width,height=height)


tmp <- dat_input2 %>%
group_by(Class) %>%
summarize( median = median(mean) )
print(tmp)

###########################################################################################

for( class in c("All" , "IGC" , "DGC") ){

    if(class == "All"){
        dat_input3 <- dat_input2
    }else{
        dat_input3 <- subset( dat_input2 , Class == class )
    }

    p <- wilcox.test(subset(dat_input3 , plotType=="TrunkDriver")$mean , subset(dat_input3 , plotType=="Other")$mean)$p.value

    if( p < 0.01 ){
        p_text <- trans(p)
    }else{
        p_text <- paste0( "P == " , round(as.numeric(p) , 3) ) 
    }

    dat_input3$p_text <- ""
    dat_input3$p_text[1] <- p_text
    dat_input3$label <- ifelse( dat_input3$Share_gene != "" , paste0( dat_input3$Share_gene , "," , dat_input3$Priver_gene ) , "" )

    plot <- ggplot( dat_input3 , aes( x = plotType , y = mean , color = plotType ) ) +
        geom_boxplot(size = 1.2 , outlier.shape = NA ) + 
        scale_color_manual(values=col_tmp) +
        scale_fill_manual(values =col_tmp) +
        geom_jitter(position=position_jitter(0.2 , seed = 1),aes(color=plotType)) +
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

    width <- 3/1
    height <- 4.47/1
    out_name <- paste0(out_path , "/Fig6B." , class , ".pdf")
    ggsave(file=out_name,plot=plot,width=width,height=height)
}


###########################################################################################

for( class in c("IGC" , "DGC") ){

    dat_input3 <- subset( dat_input2 , Class == class ) 
    dat_input3 <- subset( dat_input3 , plotType=="TrunkDriver" )
    dat_input3$x_id <- ifelse( 
        dat_input3$Priver_gene != "" , 
        paste0( dat_input3$PlotID_Divide , "(" , dat_input3$Share_gene , "->" , dat_input3$Priver_gene , ")") ,
        paste0( dat_input3$PlotID_Divide , "(" , dat_input3$Share_gene , ")") 
    )

    dat_input3 <- dat_input3[order(dat_input3$mean),]
    dat_input3$x_id <- factor( dat_input3$x_id , levels = dat_input3$x_id )

    median_time <- median(subset( dat_input2 , Class==class )$mean)

    q2 <- ggplot(dat_input3) + 
      theme(
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5),
        text = element_text(family = "mono"), 
        plot.title = element_text(hjust = 0.5,size = 12,color="black",face='bold'),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 10,color="black",face='bold'),
        axis.text.y= element_text(size = 8,color="black",face='bold'),
        axis.text.x = element_text(size = 8,color="black",face='bold',angle = 45,hjust = 1) ,
        legend.position="none"
      )
    q2 <- q2 + aes(x = x_id, y = mean, ymax = lower, ymin = upper , color = x_id )

    q2 <- q2 + 
      geom_linerange(position = position_dodge(width = 0.8), colour = "#A5A4A4") + 
      geom_hline( yintercept = median_time , size = 0.1 , linetype="dashed" , colour="gray56") + 
      geom_point(position = position_dodge(width = 0.5), size = 3) 

    q2 <- q2 + ylim(0 , 20)
    q2 <- q2 + scale_color_manual(values = c(brewer.pal(8,"Dark2")) )
    q2 <- q2 + ylab(y_lab) + ggtitle(class)

    if(class == "IGC"){
        width=5.9/1.5
        height=4.3/1.5
    }else{
        width=4/1.5
        height=4.3/1.5
    }
    
    out_name <- paste0(out_path , "/Fig6C." , class , ".pdf")
    ggsave(q2, file=out_name, width=width, height=height)
}