##########################################################################################

library(data.table)
library(optparse)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(ggsci)

###########################################################################################

dir.create(out_path , recursive = T)
col <- c( "#006699","#DDA520"  )

###########################################################################################

dat_mutiCancer <- fread(muti_cancer  ,sep = "\t",  quote="" ,header = T)
dat_mutiPre <- fread(muti_pre  ,sep = "\t",  quote="" ,header = T)
dat_mutiNodule <- rbind(dat_mutiCancer , dat_mutiPre) 

dat_info <- data.frame(fread(sample_info))
dat_gene <- data.frame(fread(gene_list))
dat_ccf <- data.frame(fread(ccf_file))

###########################################################################################

Variant_Type <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation")

###########################################################################################
dat_mutiNodule <- subset(dat_mutiNodule , t_alt_count > 0)

dat_mutiNodule <- data.frame(Hugo_Symbol = dat_mutiNodule$Hugo_Symbol,
    Chromosome = dat_mutiNodule$Chromosome , Start_Position =  dat_mutiNodule$Start_Position , End_Position = dat_mutiNodule$End_Position ,
    Reference_Allele = dat_mutiNodule$Reference_Allele , Tumor_Seq_Allele2 = dat_mutiNodule$Tumor_Seq_Allele2 , 
    Variant_Classification = dat_mutiNodule$Variant_Classification,
    t_ref_count = dat_mutiNodule$t_ref_count , t_alt_count = dat_mutiNodule$t_alt_count ,
    Tumor_Sample_Barcode = dat_mutiNodule$Tumor_Sample_Barcode)

dat_mutiNodule$Location <- paste( dat_mutiNodule$Hugo_Symbol , dat_mutiNodule$Chromosome , dat_mutiNodule$Start_Position , 
    dat_mutiNodule$Reference_Allele , dat_mutiNodule$Tumor_Seq_Allele2 , sep=":"  )

dat_mutiNodule <- subset( dat_mutiNodule , Variant_Classification %in% Variant_Type)

###########################################################################################

dat_mutiNodule$mergeCCF <- paste( dat_mutiNodule$Tumor_Sample_Barcode , dat_mutiNodule$Location , sep = ":" )

dat_ccf$mergeCCF <- paste( dat_ccf$Sample , dat_ccf$Hugo_Symbol , dat_ccf$Chr , dat_ccf$Start_Position , 
    dat_ccf$REF , dat_ccf$ALT , sep=":"  )

dat_mutiNodule <- merge( dat_mutiNodule , dat_ccf[,c("mergeCCF" , "CCF_adj")] , by = "mergeCCF" )
dat_mutiNodule <- merge( dat_mutiNodule , dat_info[,c("Tumor" , "Class")] , by.x = "Tumor_Sample_Barcode" , by.y = "Tumor" )

###########################################################################################

if(type == "IGC"){
    dat_info <- subset( dat_info , Type == "IM + IGC" )
}else if(type == "DGC"){
    dat_info <- subset( dat_info , Type == "IM + DGC" )
}

###########################################################################################

result <- c()
gene_list <- dat_gene$Gene_Symbol

for( gene in gene_list ){

    trunk_driver_num <- 0
    trunk_driver_clone_num <- 0
    trunk_driver_subclone_num <- 0
    share_driver_num <- 0
    evolutionChoose_driver_num <- 0
    private_driver_num <- 0
    private_driver_pre_num <- 0
    private_driver_inv_num <- 0
    
    trunk_driver_sample <- c()
    trunk_driver_clone_sample <- c()
    trunk_driver_subclone_sample <- c()
    share_driver_sample <- c()
    evolutionChoose_driver_sample <- c()
    private_driver_sample <- c()
    private_driver_pre_sample <- c()
    private_driver_inv_sample <- c()

    for( Sample in unique(dat_info$ID) ){

        tumors <- subset( dat_info , ID == Sample )$Tumor

        tmp <- subset( dat_mutiNodule , Tumor_Sample_Barcode %in% tumors & Variant_Classification %in% Variant_Type & Hugo_Symbol == gene )

        tmp_driver <- tmp %>% 
        group_by(Location) %>%
        summarize( MutTumor = length(Tumor_Sample_Barcode) )
        tmp_driver <- data.frame(tmp_driver)

        combine <- unique(subset( dat_info , Tumor %in% tmp$Tumor_Sample_Barcode)$Class)

        if( nrow(tmp_driver) > 0 ){

            tmp$Share <- FALSE
            for( pos in unique(tmp$Location) ){
                share <- length(which(subset( tmp , Location==pos)$Class %in% c("IM"))) > 0  & 
                length(which(subset( tmp , Location==pos)$Class %in% c("IGC" , "DGC"))) > 0
                tmp[tmp$Location==pos,"Share"] <- share
            }

            for( pos in unique(tmp$Location) ){
                share <- nrow(subset( tmp , Location==pos)) == length(tumors)
                if(share){
                    tmp[tmp$Location==pos,"Share"] <- "Trunk"
                }
            }

            if( length(which(tmp$Share!='FALSE')) > 0 ){
                if( length(which(tmp$Share=='Trunk')) > 0 ){
                    class <- "TrunkDriver"
                    trunk_driver_num <- trunk_driver_num + 1
                    trunk_driver_sample <- c(trunk_driver_sample , Sample )

                    tmp_share <- subset( tmp , Share == 'Trunk' )
                    evolution_choose <- max(tmp_share[tmp_share$Class %in% c("IGC" , "DGC"),"CCF_adj"]) >= clone_t
                    
                    if(evolution_choose){
                        evolutionChoose_driver_num <- evolutionChoose_driver_num + 1
                        evolutionChoose_driver_sample <- c(evolutionChoose_driver_sample , Sample )
                        trunk_driver_clone_num <- trunk_driver_clone_num + 1
                        trunk_driver_clone_sample <- c(trunk_driver_clone_sample , Sample )
                    }else{
                        trunk_driver_subclone_num <- trunk_driver_subclone_num + 1
                        trunk_driver_subclone_sample <- c(trunk_driver_subclone_sample , Sample )
                    }

                }else{
                    class <- "ShareDriver"
                    share_driver_num <- share_driver_num + 1
                    share_driver_sample <- c(share_driver_sample , Sample )
                    tmp_share <- subset( tmp , Share == TRUE )
                    evolution_choose <- max(tmp_share[tmp_share$Class %in% c("IGC" , "DGC"),"CCF_adj"]) >= clone_t
                    if(evolution_choose){
                        evolutionChoose_driver_num <- evolutionChoose_driver_num + 1
                        evolutionChoose_driver_sample <- c(evolutionChoose_driver_sample , Sample )
                    }
                }

            }else{
                class <- "PrivateDriver"
                private_driver_num <- private_driver_num + 1
                private_driver_sample <- c(private_driver_sample , Sample )

                if( length(which(combine %in% c("IM") )) > 0 ){
                    private_driver_pre_num <- private_driver_pre_num + 1
                    private_driver_pre_sample <- c(private_driver_pre_sample , Sample )
                }
                if( length(which(combine %in% c("IGC" , "DGC") )) > 0 ){
                    private_driver_inv_num <- private_driver_inv_num + 1
                    private_driver_inv_sample <- c(private_driver_inv_sample , Sample )
                }
                
            }
        }
    }

    tmp <- data.frame( gene = gene , 
        trunk_driver_clone_num = trunk_driver_clone_num ,
        trunk_driver_subclone_num = trunk_driver_subclone_num ,
        trunk_driver_num = trunk_driver_num ,
        share_driver_num = share_driver_num , 
        evolutionChoose_driver_num = evolutionChoose_driver_num ,
        private_driver_num = private_driver_num ,  
        private_driver_pre_num = private_driver_pre_num , private_driver_inv_num = private_driver_inv_num ,
        share_driver_sample = paste0( share_driver_sample , collapse = "," ) , 
        trunk_driver_sample = paste0( trunk_driver_sample , collapse = "," ) , 
        evolutionChoose_driver_sample = paste0( evolutionChoose_driver_sample , collapse = "," ) , 
        private_driver_sample = paste0( private_driver_sample , collapse = "," ) ,
        private_driver_pre_sample = paste0( private_driver_pre_sample , collapse = "," ) ,
        private_driver_inv_sample = paste0( private_driver_inv_sample , collapse = "," )
    )

    result <- rbind( result , tmp )

}

###########################################################################################
dat <- result

gene_order <- dat[
    order( (dat$trunk_driver_clone_num + dat$trunk_driver_subclone_num + dat$share_driver_num + dat$private_driver_inv_num + dat$private_driver_pre_num) , 
        dat$trunk_driver_clone_num , dat$trunk_driver_num , dat$share_driver_num , dat$private_driver_inv_num  , decreasing=T),"gene"]

dat <- dat[,c("gene" , "trunk_driver_clone_num" , "trunk_driver_subclone_num" , "share_driver_num" , "private_driver_inv_num" , "private_driver_pre_num")]
colnames(dat) <- c("gene" , "Trunk Clone" , "Trunk SubClone" , "Share" , "GC branch" , "IM branch")
dat <- melt(data.table(dat),id=c("gene"))

dat$variable <- as.character(dat$variable)
dat$variable <- ifelse( dat$variable %in% c("Trunk SubClone" , "Trunk Clone") , "Trunk" , dat$variable )
dat$variable <- ifelse( dat$variable %in% c("Trunk") , "Share" , dat$variable )
dat$variable <- factor( dat$variable , levels = c("IM branch" , "GC branch" , "Share") , order = T )
dat$gene <- factor( dat$gene , levels = gene_order , order = T )

col_use <- c(rgb(red=179,green=34,blue=35,alpha=255,max=255) ,
    rgb(red=247,green=184,blue=71,alpha=255,max=255) ,
    rgb(red=2,green=100,blue=190,alpha=255,max=255) ,
    "#4DAF4A"
    )
col_use <- col_use[c(4,3,1)]
names(col_use) <- c("IM branch" , "GC branch" , "Share")

if(type == "IGC"){
    width = 5
    height= 4 
}else if(type == "DGC"){
    width = 5
    height= 4 
}

dat$gene_type <- "GC\nfavored"
dat$gene_type <- ifelse( dat$gene %in% c("TP53" , "APC" , "PIK3CA") , "Maintained" , dat$gene_type )
dat$gene_type <- ifelse( dat$gene %in% c("CDH1") & type == "DGC" , "Maintained" , dat$gene_type )
dat$gene_type <- ifelse( dat$gene %in% c("MUC6" , "BMP6" , "CFTR" , "MTRR") , "IM\nfavored" , dat$gene_type )
dat$gene_type <- factor( dat$gene_type , levels = c("Maintained" , "IM\nfavored" , "GC\nfavored"))

dat <- subset( dat , !(gene %in% c("KRAS" , "CDKN2A")) )

p1 <- ggplot(dat , aes(gene, weight= value , fill=variable))+
    geom_bar(position="stack")+labs(,y="Number of Patients")+
    facet_grid(cols = vars(gene_type) , scales = "free_x" ) +
    theme(panel.background = element_rect(fill = NA, colour = "black", size = 1),
        legend.position ='top',
        panel.grid.major = element_line(colour=NA),
        legend.text = element_text(size = 8,color="black",face='bold'),
        axis.text.y = element_text(size = 8,color="black",face='bold'),
        axis.title.y = element_text(size = 12,color="black",face='bold'),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 12,color="black",face='bold' , angle = 90 , vjust=.5, hjust=1),
        strip.text.x = element_text(size = 8 , face = 'bold'),
        axis.ticks.length = unit(0.2, "cm") ,
        axis.line = element_line(size = 0.5)) +
    guides(fill = guide_legend(reverse=TRUE , title = NULL)) +
    scale_fill_manual(values = col_use)

p <- p1
gp <- ggplotGrob(p)
facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
x.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
                function(l) length(l$range$range))
gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var

image_name <- paste0(out_path , "/GeneTrunk.evolutionChoose.ratio.",type,".pdf")
pdf(image_name , width=width ,height=height)
grid::grid.draw(gp)
dev.off()
