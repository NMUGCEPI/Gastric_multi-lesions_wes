##########################################################################################

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(optparse)
library(circlize)

###########################################################################################

dat_im_smg <- data.frame(fread(im_list , header = T))
dat_igc_smg <- data.frame(fread(igc_list , header = T))
dat_dgc_smg <- data.frame(fread(dgc_list , header = T))

class_order <- data.frame(fread(class_order_file , header = T))
info <- data.frame(fread(info_public_file))

dat_im <- data.frame(fread(im_maf_file))
dat_gc <- data.frame(fread(gc_maf_file))

###########################################################################################

smg <- data.frame(Gene_Symbol = unique(c(dat_im_smg$Gene_Symbol , dat_igc_smg$Gene_Symbol , dat_dgc_smg$Gene_Symbol)))

###########################################################################################

info <- subset( info , !(Molecular.subtype %in% c("POLE" , "MSI") ))
info$Molecular.subtype <- ifelse( info$Molecular.subtype %in% c("EBV" , "unknown") , "Other" , info$Molecular.subtype )

###########################################################################################
dat_im$Tumor_Sample_Barcode <- paste0( dat_im$Tumor_Sample_Barcode , "_IM" )

info_im <- subset( info , From == "NJMU"  )
info_im$Tumor <- paste0( info_im$Tumor , "_IM" )
info_im$Class <- "IM"
info_im$Stage <- "unknown"

info <- rbind(info_im , info)

info$Class <- factor(info$Class , levels = c("IM", "IGC" , "DGC" ) , ordered=T)
info$Normal <- gsub( "_IM" , "" , info$Tumor )
rownames(info) <- info$Tumor

info$Gender[info$Gender=="male"] <- "Male"
info$Gender[info$Gender=="female"] <- "Female"

info$HP <- ifelse( info$HP == "unknown" , "Unknown" , info$HP )
info$HP <- factor(info$HP , levels = c("Positive", "Negative" , "Unknown" ) , ordered=T)

info$Alcohol <- ifelse( info$Alcohol == "unknown" , "Unknown" , info$Alcohol )
info$Alcohol <- ifelse( info$Alcohol == "No" , "Non-drinker" , info$Alcohol )
info$Alcohol <- ifelse( info$Alcohol == "Drink" , "Drinker" , info$Alcohol )
info$Alcohol <- factor( info$Alcohol , levels = c("Drinker", "Non-drinker" , "Unknown" ) , ordered=T)

info$Tobacco <- ifelse( info$Tobacco == "unknown" , "Unknown" , info$Tobacco )
info$Tobacco <- ifelse( info$Tobacco == "No" , "Non-smoker" , info$Tobacco )
info$Tobacco <- ifelse( info$Tobacco == "Smoke" , "Smoker" , info$Tobacco )
info$Tobacco <- factor( info$Tobacco , levels = c("Smoker", "Non-smoker" , "Unknown" ) , ordered=T)

###########################################################################################

col <- c(
  brewer.pal(9,"YlGnBu")[6],
  rgb(234,106,79,alpha=255,maxColorValue=255),
  rgb(203,24,30,alpha=255,maxColorValue=255),
  rgb(255,0,0,alpha=255,maxColorValue=255)
  )

names(col) <- c("IM" , "IGC" , "DGC" , "GC")

Variant_Types <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation")

col_class <- col[1:3]

###########################################################################################

CreateMutMatrix <- function(dat = dat , smg = smg , Variant_Type = Variant_Type ){

	mut <- dat
	mut$Variant_Classification <- as.character(mut$Variant_Classification)
	mut <- mut[which(mut$Variant_Classification %in% Variant_Type),]

	mut[grep("In_Frame",mut$Variant_Classification),'Variant_Classification'] = "In_Frame"
	mut[grep("Frame_Shift",mut$Variant_Classification),'Variant_Classification'] = "Frame_Shift"

	Gene <- smg[smg$Gene_Symbol %in% mut$Hugo_Symbol,]

	Sample <- info[,"Tumor"]
	maf_matrix <- matrix("" , ncol = length(Sample) , nrow = length(Gene) , dimnames = list(Gene,Sample))

	for(gene in rownames(maf_matrix)){
		print(gene)

		for(tumor in colnames(maf_matrix)){
			index <- which(mut$Hugo_Symbol==gene & mut$Tumor==tumor)
			if(length(index)==0){
				var=""
			}else if(length(index)==1){
				var <- mut[index,'Variant_Classification']
			}else if(length(index)>1){
				var <- "Multiple_Hits"
			}
			print(var)
			maf_matrix[ which(rownames(maf_matrix)==gene) , which(colnames(maf_matrix)==tumor) ] <-  var
		}
	}

	return(maf_matrix)
}

MutMatrixOrder <- function( mut = mut , info = info){

	mut_per <- paste0( 100 * round(
			apply(mut , 1 , function(x){length(unique(info[info$Tumor %in% names(which(x!="")) ,"Normal"]))} ) / 
			length(unique(info$Normal)) ,
		2)) 
	mut_per <- as.numeric(mut_per)

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

	mut <- rbind( mut , Class = info$Class )
	mut <- mut[ , order(mut["Class",] , 
		mut[1,] , mut[2,] , mut[3,] , mut[4,] , mut[5,] , mut[6,] , 
		mut[7,] , mut[8,] , mut[9,] , mut[10,] , 
		mut[11,] , mut[12,] , mut[13,] , mut[14,] , mut[15,] , mut[16,] , 
		mut[17,] , mut[18,] , mut[19,] , mut[20,] , mut[21,] ,
		colnames(mut) , decreasing = T )]

	mut <- mut[rownames(mut)!=c("Class"),]

	return(mut)
}

plotMutWaterFull <- function( MutMatrix = MutMatrix , col_class = col_class , Variant_Type_Combine = Variant_Type_Combine , images_name = images_name ){

	annotation_name <- c("Class","Mut_Num")

	mut <- MutMatrix[apply(MutMatrix,1,function(x){length(which(x!=""))>1}),]
	mut <- MutMatrixOrder( mut = mut , info = info )
	
	class_order <- info[order(info$Class , decreasing=T),"Class"]
	sample_order <- colnames(mut)

	from_order <- info[sample_order,"From"]
	molecular_order <- info[sample_order,"Molecular.subtype"]
	hp_order <- info[sample_order,"HP"]
	alcohol_order <- info[sample_order,"Alcohol"]
	tobacco_order <- info[sample_order,"Tobacco"]

	################################################################################################
	col = c(rgb(red=48,green=115,blue=186,alpha=255,max=255),
		rgb(red=236,green=27,blue=35,alpha=255,max=255),
		rgb(red=236,green=179,blue=33,alpha=255,max=255),
		rgb(red=235,green=230,blue=26,alpha=255,max=255),
		rgb(red=150,green=131,blue=189,alpha=255,max=255),
		rgb(red=65,green=174,blue=119,alpha=255,max=255))

	names(col) = c(
	  'Missense_Mutation',
	  'Nonsense_Mutation',
	  'Frame_Shift',
	  'In_Frame',
	  'Splice_Site',
	  'Multiple_Hits'
	)

	alter_fun = list(
	    background = function(x, y, w, h) {
	        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
	            gp = gpar(fill = "#F2F2F2", col = NA))
	    },
	    Missense_Mutation = function(x, y, w, h) {
	        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
	            gp = gpar(fill = col["Missense_Mutation"], col = NA))
	    },
	    Nonsense_Mutation = function(x, y, w, h) {
	        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
	            gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
	    },
	  
	    Nonstop_Mutation = function(x, y, w, h) {
	        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
	            gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
	    },
	    Frame_Shift = function(x, y, w, h) {
	        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
	            gp = gpar(fill = col["Frame_Shift"], col = NA))
	    },
	    In_Frame = function(x, y, w, h) {
	        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
	            gp = gpar(fill = col["In_Frame"], col = NA))
	    },
	    Splice_Site = function(x, y, w, h) {
	        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
	            gp = gpar(fill = col["Splice_Site"], col = NA))
	    },
	    Multiple_Hits = function(x, y, w, h) {
	        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
	            gp = gpar(fill = col["Multiple_Hits"], col = NA))
	    },
	    show_legend = FALSE 
	)

	################################################################################################

	col_from <- c(
		rgb(red=244, green=241,blue=222,alpha=255,max=255) ,
		rgb(red=223, green=122,blue=94,alpha=255,max=255) ,
		rgb(red=60, green=64,blue=91,alpha=255,max=255) ,
		rgb(red=130, green=178,blue=154,alpha=255,max=255) ,
		rgb(red=242, green=204,blue=142,alpha=255,max=255) 
		)
	names(col_from) <- c("NJMU" , "OncoSG" ,  "TCGA" , "TMUCIH" , "Utokyo")

	col_molecular <- c(col_from[c(2,3)] , "grey")
	names(col_molecular) <- c("CIN" , "GS" , "Other")

	col_hp <- c(col_from[c(2,3)] , "grey")
	names(col_hp) <- c("Positive" , "Negative" , "Unknown")

	col_drink <- c(col_from[c(2,3)] , "grey")
	names(col_drink) <- c("Drinker" , "Non-drinker" , "Unknown")

	col_smoke <- c(col_from[c(2,3)] , "grey")
	names(col_smoke) <- c("Smoker" , "Non-smoker" , "Unknown")

	top_annotation <- HeatmapAnnotation(
		foo = anno_empty(height = unit(2, "cm") , border =F) ,  
		Molecular_subtype = molecular_order ,
		HP_status = hp_order ,
		Drinking_status = alcohol_order ,
		Smoking_status = tobacco_order ,
		Study = from_order ,
		col = list( 
			Molecular_subtype = col_molecular ,
			HP_status = col_hp , 
			Drinking_status = col_drink ,
			Smoking_status = col_smoke ,
			Study = col_from 
		) ,
		annotation_name_side = "left" , 
	  	border = T ,
	  	gap = unit(1, "mm") ,
	  	annotation_name_gp = gpar(fontsize = 12),
	  	show_legend = FALSE  
	)

	bottom_annotation <- HeatmapAnnotation(
		Class = class_order ,
		col = list( 
			Class = col_class 
		) ,
		annotation_name_side = "left" , 
	  	border = T ,
	  	gap = unit(1, "mm") ,
	  	show_annotation_name = FALSE , 
	  	annotation_name_gp = gpar(fontsize = 12),
	  	show_legend = FALSE  
	)

	################################################################################################
	IGC_sample <- info[info$Class=="IGC","Tumor"]
	DGC_sample <- info[info$Class=="DGC","Tumor"]
	IM_sample <- info[info$Class=="IM","Tumor"]

	mut_igc <- mut[,colnames(mut) %in% c(IGC_sample)]
	mut_dgc <- mut[,colnames(mut) %in% c(DGC_sample)]
	mut_im <- mut[,colnames(mut) %in% c(IM_sample)]

	NumMut_igc <- apply(mut_igc , 1 ,function(x){
			length(unique(info[info$Tumor %in% names(which(x!="")) ,"Normal"]))/length(IGC_sample)
		}
	)

	NumMut_dgc <- apply(mut_dgc , 1 ,function(x){
			length(unique(info[info$Tumor %in% names(which(x!="")) ,"Normal"]))/length(DGC_sample)
		}
	)

	NumMut_im <- apply(mut_im , 1 ,function(x){
			length(unique(info[info$Tumor %in% names(which(x!="")) ,"Normal"]))/length(IM_sample)
		}
	)
	NumMut_im <- round(NumMut_im * 100)
	NumMut_igc <- round(NumMut_igc * 100)
	NumMut_dgc <- round(NumMut_dgc * 100)
	grey_color <- rgb(red=243, green=243,blue=243,alpha=255,max=255)

	num_size <- 14
	col_fun = colorRamp2(c(0 , 20), c("white", rgb(203,24,30,alpha=255,maxColorValue=255)))
	right_annotation <- rowAnnotation(
		IM = anno_text(paste0( " " , NumMut_im , "% "), gp = gpar( fill = col_fun(NumMut_im) , col = "black", border = grey_color , fontsize = num_size) , show_name = TRUE ) ,
		IGC = anno_text(paste0( " " , NumMut_igc , "% "), gp = gpar( fill = col_fun(NumMut_igc) , col = "black", border = grey_color , fontsize = num_size) , show_name = TRUE ) ,
		DGC = anno_text(paste0( " " , NumMut_dgc , "% "), gp = gpar( fill = col_fun(NumMut_dgc) , col = "black", border = grey_color , fontsize = num_size) , show_name = TRUE ) ,
		annotation_name_rot = 90 ,
		annotation_name_side = "top" ,
	    show_legend = FALSE  
	)

	################################################################################################
	im_smg <- ifelse(rownames(mut) %in% dat_im_smg$Gene_Symbol , 1 , 0)
	igc_smg <- ifelse(rownames(mut) %in% dat_igc_smg$Gene_Symbol , 1 , 0)
	dgc_smg <- ifelse(rownames(mut) %in% dat_dgc_smg$Gene_Symbol , 1 , 0)
	

	col_im = colorRamp2(c(0 , 1) , c(grey_color, col_class["IM"]) )
	col_igc = colorRamp2(c(0 , 1) , c(grey_color, col_class["IGC"]) )
	col_dgc = colorRamp2(c(0 , 1) , c(grey_color, col_class["DGC"]) )

	mutclass_df <- data.frame( IM = im_smg , IGC = igc_smg , DGC = dgc_smg )

	left_annotation = rowAnnotation(
		df = mutclass_df,
		col = list(
			IM = col_im , IGC = col_igc , DGC = col_dgc
			) ,
	    annotation_name_rot = 90 ,
	    annotation_name_side = "top" ,
	    show_legend = FALSE  
	)
	
	################################################################################################
	mut_per <- as.numeric(apply(mut , 1 , function(x){length(unique(info[info$Tumor %in% names(which(x!="")) ,"Normal"]))} ))
			
	font_size <- 20

	ldg_mol = Legend(labels = names(col_molecular) , title = "Molecular subtype" ,
	 legend_gp = gpar(fill = col_molecular , bar_width = 1 , fontsize = font_size) , ncol = 1 , 
	 gap = unit(1, "cm") 
	)

	ldg_hp = Legend(labels = names(col_hp) , title = "HP status" ,
	 legend_gp = gpar(fill = col_hp , bar_width = 1 , fontsize = font_size) , ncol = 1 , 
	 gap = unit(1, "cm") 
	)

	ldg_drink = Legend(labels = names(col_drink) , title = "Drinking status" ,
	 legend_gp = gpar(fill = col_drink , bar_width = 1 , fontsize = font_size) , ncol = 1 , 
	 gap = unit(1, "cm") 
	)

	ldg_smoke = Legend(labels = names(col_smoke) , title = "Smoking status" ,
	 legend_gp = gpar(fill = col_smoke , bar_width = 1 , fontsize = font_size) , ncol = 1 , 
	 gap = unit(1, "cm") 
	)

	ldg_from = Legend(labels = names(col_from) , title = "Study" ,
	 legend_gp = gpar(fill = col_from , bar_width = 1 , fontsize = font_size) , ncol = 2 , 
	 gap = unit(1, "cm") 
	)

	col_Mut <- col 

	ldg_Variant = Legend(labels = names(col_Mut) , title = "Mutations" , border = "black" , 
	 	legend_gp = gpar(fill = col_Mut , bar_width = 1 ,fontsize = font_size) , ncol = 2 , 
	 	gap = unit(1, "cm") 
	)	

	lgd_all <- packLegend( ldg_Variant , ldg_mol , ldg_hp , ldg_drink , ldg_smoke , ldg_from , 
		column_gap = unit(1, "cm") , row_gap = unit(15, "mm") , direction = "horizontal"  )

	################################################################################################
	p <- oncoPrint(mut, name = "cases",
	    alter_fun = alter_fun, col = col, 
	    top_annotation = bottom_annotation ,
	    bottom_annotation = NULL ,
	    row_names_side = "left", row_names_gp = gpar(fontsize = 16) ,  
	    right_annotation = right_annotation, 
	    left_annotation = left_annotation, 
	    show_pct = FALSE , 
	    border = TRUE,
	    row_order = 1:nrow(mut) ,
	    column_order = sample_order,
	    show_heatmap_legend = FALSE ,
	    row_title_rot = 0 , row_title_gp = gpar(fontsize = 15) , 
	    column_split = class_order , column_title = NULL ,
	    row_gap = unit(0, "points")
	)

	##
	pdf(images_name , width = 8.95/0.8 , height = 4.59/0.8)
	draw(p)
	dev.off()

}


###########################################################################################
dat <- rbind( dat_gc , dat_im )
dat <- dat[,colnames(dat) != "Variant_Type" ]
dat$Tumor <- dat$Tumor_Sample_Barcode

Variant_Type <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation")
Variant_Type_Combine <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift","In_Frame","Splice_Site","Nonstop_Mutation","Multiple_Hits")

MutMatrix <- CreateMutMatrix( dat = dat , smg = smg , Variant_Type = Variant_Type )

images_name <- paste0(images_path , "/Fig1.waterfall.pdf")
plotMutWaterFull( MutMatrix = MutMatrix , col_class = col_class , Variant_Type_Combine = Variant_Type_Combine , images_name = images_name )


