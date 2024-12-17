##########################################################################################

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(optparse)
library(circlize)

###########################################################################################
smg <- data.frame(fread(smg_list , header = T))
dat_ccf <- data.frame(fread(ccf_file))

class_order <- data.frame(fread(class_order_file , header = T))
class_order_sub <- data.frame(fread(class_order_sub_file , header = T))

info <- data.frame(fread(info_file))

dat_tp53_pre <- fread(tp53_pre_file)
dat_tp53_can <- fread(tp53_cancer_file)

###########################################################################################

info$Class <- factor(info$Class,levels= unique(class_order$Class), ordered=TRUE)
info$Class_sub <- factor(info$Class_sub,levels= unique(class_order_sub$Class), ordered=TRUE)

info$ID_order <- paste0(info$ID,"_", as.numeric(info$Class_sub),"_",info$Class)
info <- info[order(info$ID_order),]

###########################################################################################
show_gene <- smg$Gene_Symbol
Variant_Types <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation")

dat_ccf <- subset( dat_ccf , Hugo_Symbol %in% show_gene & Variant_Classification %in% Variant_Types  )
use_sample <- unique(subset( dat_ccf , Share=="TRUE" )$ID)
dat_ccf <- subset( dat_ccf , ID %in% use_sample  )
info <- subset( info , ID %in% use_sample )

###########################################################################################
normal_order1 <- unique(info[info$Normal %in% dat_tp53_pre$Normal,"ID"])
normal_order2 <- unique(info[info$Normal %in% dat_tp53_can$Normal,"ID"])
normal_order2 <-  normal_order2[!(normal_order2 %in% normal_order1)]
normal_order_tp53 <- c( normal_order1 , normal_order2 )
normal_order <- normal_order_tp53

info <- subset( info , ID %in% normal_order )
info$ID <- factor(info$ID,levels= normal_order, ordered=TRUE)
info <- info[order(info$ID , info$ID_order),]

###########################################################################################

col <- c(
  brewer.pal(9,"YlGnBu")[6],
  rgb(234,106,79,alpha=255,maxColorValue=255),
  rgb(203,24,30,alpha=255,maxColorValue=255),
  rgb(255,0,0,alpha=255,maxColorValue=255)
  )

names(col) <- c("IM" , "IGC" , "DGC" , "GC")


col_class <- col[1:3]

CreateMutMatrix <- function(dat = dat , Variant_Type = Variant_Type ){

	pre <- c("IM")
	can <- c("IGC","DGC")

	mut <- dat
	mut$Variant_Classification <- as.character(mut$Variant_Classification)
	mut <- mut[which(mut$Variant_Classification %in% Variant_Type),]
	mut$use_share <- ifelse( mut$Share == "TRUE" , "Share" , "Private" )
	mut$clonal_status <- ifelse( mut$CCF_adj >= 0.6 , "Clonal" , "Subclonal" )
	mut$cnv_status <- ifelse( mut$total_cn==2 & mut$minor_cn==0 , "LOH" , "" )
	mut$cnv_status <- ifelse( mut$total_cn > 2 , "AMP" , mut$cnv_status)
	mut$cnv_status <- ifelse( mut$total_cn == 1 & mut$minor_cn==0 , "LOSS" , mut$cnv_status )

	mut[grep("In_Frame",mut$Variant_Classification),'Variant_Classification'] = "In_Frame"
	mut[grep("Frame_Shift",mut$Variant_Classification),'Variant_Classification'] = "Frame_Shift"

	Gene <- show_gene

	Sample <- info[,"Tumor"]
	maf_matrix <- matrix("" , ncol = length(Sample) , nrow = length(Gene) , dimnames = list(Gene,Sample))

	for(gene in rownames(maf_matrix)){
		print(gene)

		for(tumor in colnames(maf_matrix)){
			index <- which(mut$Hugo_Symbol==gene & mut$Tumor==tumor)
			if(length(index)==0){
				var=""
			}else if(length(index)==1){
				var <- mut[index,'clonal_status']
				var <- paste(mut[index,'clonal_status'] , mut[index,'cnv_status'] , sep = ";")

			}else if(length(index)>1){
				if(length( which(mut[index,'clonal_status'] == "Clonal") != 0 )){
					var <- paste( unique( mut[index,'clonal_status']) , mut[index,'cnv_status'] , sep = ";")
				}else{
					var <- paste( "Subclonal", mut[index,'cnv_status'] , sep = ";")
				}
			}
			print(var)
			maf_matrix[ which(rownames(maf_matrix)==gene) , which(colnames(maf_matrix)==tumor) ] <-  var
		}
	}

	return(maf_matrix)
}

MutMatrixOrder <- function( mut = mut , info = info){

	IGC_sample <- info[info$Class=="IGC","Tumor"]
	DGC_sample <- info[info$Class=="DGC","Tumor"]
	IM_sample <- info[info$Class=="IM","Tumor"]

	mut_tumor <- mut[,colnames(mut) %in% c(IGC_sample,DGC_sample)]
	mut_pre <- mut[,colnames(mut) %in% c(IM_sample)]

	NumMut_Tumor <- apply(mut_tumor , 1 ,function(x){
		length(unique(info[info$Tumor %in% names(which(x!="")) ,"Normal"]))}
	)

	NumMut_Pre <- apply(mut_pre , 1 ,function(x){
		length(unique(info[info$Tumor %in% names(which(x!="")) ,"Normal"]))}
	)

	MutNumMatrix <- data.frame(cbind(NumMut_Pre,NumMut_Tumor))

	mut <- mut[order(MutNumMatrix$NumMut_Tumor , MutNumMatrix$NumMut_Pre , decreasing = T ),]

	return(mut)
}

plotMutWaterFull <- function( MutMatrix = MutMatrix , col_class = col_class , Variant_Type_Combine = Variant_Type_Combine , images_name = images_name ){

	annotation_name <- c("Class")

	mut <- MutMatrix[apply(MutMatrix,1,function(x){length(which(x!=""))>=1}),]

	mut <- MutMatrixOrder( mut = mut , info = info )
	
	class_order <- as.character(info$Class)
	sample_order <- colnames(mut)
	type_order <- as.character(info$Type)
	type_order <- factor(type_order , levels = c("IM + IGC" , "IM + DGC") , order = T)

	################################################################################################
	col = c(
		"orange" , rgb(red=90,green=147,blue=189,alpha=255,max=255) ,
		"white")

	names(col) = c(
	  'Clonal',
	  'Subclonal',
	  'Private'
	)

	col_cnv <-c( "#00A087FF" , "#E64B35FF" , "#3C5488FF" )
	names(col_cnv) <- c("LOH" ,
	  "AMP" ,
	  "LOSS")

	alter_fun = list(
	    background = function(x, y, w, h) {
	        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
	            gp = gpar(fill = "#F2F2F2", col = NA))
	    },
	    Clonal = function(x, y, w, h) {
	        grid.rect(x, y + h*0.45 , w-unit(0.5, "mm"), 0.44* h, gp = gpar(fill = col["Clonal"], col = NA), just = "top")
	    },
	    Subclonal = function(x, y, w, h) {
	        grid.rect(x, y + h*0.45 , w-unit(0.5, "mm"), 0.44* h, gp = gpar(fill = col["Subclonal"], col = NA), just = "top")
	    },
	    LOH = function(x, y, w, h) {
			grid.rect(x, y - h*0.45 , w-unit(0.5, "mm"), 0.44* h, gp = gpar(fill = col_cnv["LOH"], col = NA), just = "bottom")
	    },
	    AMP = function(x, y, w, h) {
			grid.rect(x, y - h*0.45 , w-unit(0.5, "mm"), 0.44* h, gp = gpar(fill = col_cnv["AMP"], col = NA), just = "bottom")
	    },
	    LOSS = function(x, y, w, h) {
			grid.rect(x, y - h*0.45 , w-unit(0.5, "mm"), 0.44* h, gp = gpar(fill = col_cnv["LOSS"], col = NA), just = "bottom")
	    },
	    show_legend = FALSE 
	)

	top_annotation <- HeatmapAnnotation(
		foo = anno_empty(height = unit(2, "cm") , border =F) ,  
		Class = class_order ,
		col = list( 
			Class = col_class 
		) ,
		annotation_name_side = "left" , 
	  	border = T ,
	  	gap = unit(1, "mm") ,
	  	show_annotation_name = c(Mut_Num = FALSE) , 
	  	annotation_name_gp = gpar(fontsize = 12),
	  	show_legend = FALSE  
	)

	################################################################################################

	ldg_Class = Legend(labels = names(col_class) , title = "Pathogenic" ,
	 legend_gp = gpar(fill = col_class , bar_width = 1 , fontsize = 12) , ncol = 3 , 
	 gap = unit(1, "cm") 
	)

	col_Mut <- col[1:2]
	ldg_Variant = Legend(labels = names(col_Mut) , title = "Clonal status" , border = "black" , 
	 	legend_gp = gpar(fill = col_Mut , bar_width = 1 ,fontsize = 12) , ncol = 3 , 
	 	gap = unit(1, "cm") 
	)

	col_Mut <- col 
	ldg_CNV = Legend(labels = names(col_cnv) , title = "CNV status" , border = "black" , 
	 	legend_gp = gpar(fill = col_cnv , bar_width = 1 ,fontsize = 12) , ncol = 3 , 
	 	gap = unit(1, "cm")
	)	

	lgd_all <- packLegend(ldg_Class , ldg_Variant , ldg_CNV , column_gap = unit(1, "cm") , row_gap = unit(10, "mm") , direction = "horizontal"  )

	################################################################################################
	p <- oncoPrint(mut, name = "cases", 
	    alter_fun = alter_fun, col = col, 
	    top_annotation = top_annotation  ,
	    row_names_side = "left", row_names_gp = gpar(fontsize = 16) ,  
	    left_annotation = NULL , right_annotation = NULL ,  
	    show_pct = FALSE , 
	    border = TRUE,
	    row_order = 1:nrow(mut) ,
	    column_order = sample_order,
	    column_split = type_order , column_title = NULL ,
	    show_heatmap_legend = FALSE
	)

	width = 15
	height = 4

	pdf(images_name , width = width , height = height)
	draw(p )
	draw(lgd_all, x = unit(0.3, "npc"), y = unit(0.99, "npc"), just = c("left", "top"))
	dev.off()
}


###########################################################################################
dat <- dat_ccf
dat$Tumor <- dat$Sample

Variant_Type <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation")
Variant_Type_Combine <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift","In_Frame","Splice_Site","Nonstop_Mutation","Multiple_Hits")

MutMatrix <- CreateMutMatrix( dat = dat , Variant_Type = Variant_Type )

images_name <- paste0(images_path , "/Fig5A.pdf")
plotMutWaterFull( MutMatrix = MutMatrix , col_class = col_class , Variant_Type_Combine = Variant_Type_Combine , images_name = images_name )


