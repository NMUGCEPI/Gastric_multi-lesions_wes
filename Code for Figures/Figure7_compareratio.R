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

dat_maf_mss <- data.frame(fread( maf_file ))
dat_maf_msi <- data.frame(fread( maf_msi_file ))

sc_dataset <- load(single_cell_file, verbose = F)
sc_dataset_all <- epiall_nor_PCA_50_RE0.5
Idents(sc_dataset_all) <- "sample"   

##########################################################################################

dat_maf_public <- rbind( dat_maf_msi , dat_maf_mss )
dat_maf_public <- merge( dat_maf_public , info_singlecell[,c("Tumor" , "singlecell_ID" , "Class")] , by.x = "Tumor_Sample_Barcode" , by.y = "Tumor" )
dat_maf_public <- subset( dat_maf_public , t_alt_count > 0 )
Variant_Types <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation")
dat_maf_public <- subset( dat_maf_public , Hugo_Symbol == gene & Variant_Classification %in% Variant_Types )
dat_maf_public <- subset( dat_maf_public , Class == "IM" )

##########################################################################################
mutTumor <- unique(dat_maf_public$singlecell_ID)
wildTumor <- subset( info_singlecell , !(singlecell_ID %in% mutTumor) )$singlecell_ID

trans <- function(num){
    up <- floor(log10(num))
    down <- round(num*10^(-up),2)
    text <- paste("P == ",down," %*% 10","^",up)
    return(text)
}

cell_order <- c("Enterocytes" , "Neck" , "Chief" , "Pit" , "Endocrine" , "Goblet" , "Parietal" , "Tumor")

##########################################################################################
result_final <- c()

sc_dataset <- sc_dataset_all
wild_cells <- sc_dataset$celltype[grep( paste0(wildTumor , collapse = "|") , names(sc_dataset$celltype) )]
mut_cells <- sc_dataset$celltype[grep( paste0(mutTumor , collapse = "|") , names(sc_dataset$celltype) )]
wild_cells <- data.frame(table(wild_cells))
mut_cells <- data.frame(table(mut_cells))

wild_cells$Ratio <- wild_cells$Freq/sum(wild_cells$Freq)
mut_cells$Ratio <- mut_cells$Freq/sum(mut_cells$Freq)

wild_cells$Type <- "Wild"
mut_cells$Type <- "Mut"
colnames(wild_cells)[1] <- "Cell_Type"
colnames(mut_cells)[1] <- "Cell_Type"

dat_plot <- rbind(wild_cells , mut_cells )

result_plot <- c()
for( cell_Type in unique(dat_plot$Cell_Type) ){
	dat_tmp <- subset( dat_plot , Cell_Type == cell_Type )

	a <- round(100 * subset( dat_plot , Cell_Type == cell_Type & Type == "Mut"  )$Ratio)
	b <- 100 - a

	c <- round(100 * subset( dat_plot , Cell_Type == cell_Type & Type == "Wild"  )$Ratio)
	d <- 100 - c

	p <- fisher.test( matrix(c(a,b,c,d) , ncol = 2) )$p.value

	dat_tmp <- data.frame( 
		Cell_Type = c("Cells of interest" , "Other Cells" , "Cells of interest" , "Other Cells") ,
		Ratio = c(
			subset( dat_plot , Cell_Type == cell_Type & Type == "Mut"  )$Ratio ,
			1 - subset( dat_plot , Cell_Type == cell_Type & Type == "Mut"  )$Ratio ,
			subset( dat_plot , Cell_Type == cell_Type & Type == "Wild"  )$Ratio ,
			1 - subset( dat_plot , Cell_Type == cell_Type & Type == "Wild"  )$Ratio
			),
		Type = c("Mut" , "Mut" , "Wild" , "Wild") ,
		Cell_Type_use = cell_Type
	 )
	if( p < 0.01 ){
        p_text <- trans(p)
    }else{
        p_text <- paste0( "P == " , round(as.numeric(p) , 3) ) 
    }
	dat_tmp$p <- ""
	dat_tmp$p_text <- ""
	dat_tmp$p[1] <- p
	dat_tmp$p_text[1] <- p_text
	result_plot <- rbind( result_plot , dat_tmp )
}

result_plot$value_text <- paste0( round(result_plot$Ratio , 2) * 100 , "%") 
result_plot$value_text <- ifelse( round(result_plot$Ratio , 2) == 0 , "" , result_plot$value_text  )
result_plot$Class <- class
result_final <- rbind( result_final , result_plot )

col <- c(
    rgb(red=179,green=34,blue=35,alpha=255,max=255), 
    rgb(red=2,green=100,blue=190,alpha=255,max=255) 
)
col <- col[2:1]
names(col) <- c("Other Cells" , "Cells of interest")

result_plot <- data.frame(result_plot)
result_plot <- subset( result_plot , Cell_Type_use != "Tumor" )
result_plot$Cell_Type_use <- factor( result_plot$Cell_Type_use , levels = cell_order , order = T )

p1 <- ggplot(result_plot,aes(x=Type,y=Ratio,fill=factor(Cell_Type))) +
	geom_bar(stat="identity") +
	ylab("Proportion") +
	geom_text(aes(label=p_text , y = 1.05 , x = 1.5),parse = TRUE,size=3.5 , color = "black") +
	geom_text(aes(label=value_text) , position=position_stack(vjust = 0.5) , size=4 , color="white")+
	xlab(NULL) +
	facet_grid( .~Cell_Type_use , scales = "free_x") +
	theme_bw() +
  	theme(
  		panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position ='top',
        legend.title = element_blank() ,
        legend.text = element_text(size = 12,color="black",face='bold'),
        axis.text.x = element_text(size = 12,color="black",face='bold'),
        axis.text.y = element_text(size = 10,color="black",face='bold'),
        axis.title.x = element_text(size = 10,color="black",face='bold'),
        axis.title.y = element_text(size = 14,color="black",face='bold'),
        strip.text.x = element_text(size = 12,color="black",face='bold'),
        axis.ticks.length = unit(0.2, "cm") ,
        axis.line = element_line(size = 0.5))  +
	scale_fill_manual(values=c(col))


out_name <- paste0( out_path , "/Fig7C1.pdf"  )
ggsave(file=out_name,plot=p1,width=10,height=4)



##########################################################################################

result_final <- c()
sc_dataset <- sc_dataset_all

dat_wild_cells <- c()
for( tumor in wildTumor ){
	mut_cells <- sc_dataset$celltype[ grep( tumor , names(sc_dataset$celltype)) ]
	if(length(mut_cells) > 0){
		mut_cells <- data.frame(table(mut_cells))
		mut_cells$Ratio <- mut_cells$Freq/sum(mut_cells$Freq)
		colnames(mut_cells)[1] <- "Cell_Type"
		mut_cells$Sample <- tumor
		dat_wild_cells <- rbind( dat_wild_cells , mut_cells )
	}
}
dat_wild_cells$Type <- "Wild"

dat_mut_cells <- c()
for( tumor in mutTumor ){
	mut_cells <- sc_dataset$celltype[ grep( tumor , names(sc_dataset$celltype)) ]
	if(length(mut_cells) > 0){
		mut_cells <- data.frame(table(mut_cells))
		mut_cells$Ratio <- mut_cells$Freq/sum(mut_cells$Freq)
		colnames(mut_cells)[1] <- "Cell_Type"
		mut_cells$Sample <- tumor
		dat_mut_cells <- rbind( dat_mut_cells , mut_cells )
	}
}

dat_mut_cells$Type <- "Mut"
dat_plot <- rbind(dat_wild_cells , dat_mut_cells )

result_plot <- c()
for( cell_Type in unique(dat_plot$Cell_Type) ){
	dat_tmp <- subset( dat_plot , Cell_Type == cell_Type )

    a <- dat_tmp[dat_tmp$Type=="Mut","Ratio"]
    b <- dat_tmp[dat_tmp$Type=="Wild","Ratio"]
    
    p <- wilcox.test( a , b )$p.value
	if( p < 0.01 ){
        p_text <- trans(p)
    }else{
        p_text <- paste0( "p == " , round(as.numeric(p) , 3) ) 
    }
	dat_tmp$p <- ""
	dat_tmp$p_text <- ""
	dat_tmp$p[1] <- p
	dat_tmp$p_text[1] <- p_text
	result_plot <- rbind( result_plot , dat_tmp )
}

result_plot$Class <- class
result_final <- rbind( result_final , result_plot )

col <- c(
    rgb(red=179,green=34,blue=35,alpha=255,max=255), 
    rgb(red=2,green=100,blue=190,alpha=255,max=255) 
)

names(col) <- c("Mut" , "Wild")
result_plot$Cell_Type <- factor( result_plot$Cell_Type , levels = cell_order , order = T )
result_plot <- subset( result_plot , Cell_Type != "Tumor" )
result_plot$size_dot <- ifelse( result_plot$Type == "Wild" , 1.2 , 1.3 )
result_plot$Type <- factor( result_plot$Type , levels = c("Wild" , "Mut") , order = T )

p1 <- ggplot(result_plot, aes(x = Cell_Type , y = Ratio )) +
    geom_boxplot(size = 1.2 , outlier.alpha=0 , color = col['Wild']) + 
    geom_jitter(position = position_jitter(0.2) , aes(color = Type , shape = Type )  , size = 3) + 
    scale_color_manual(values=c(col)) +
    xlab(NULL) +
    ylab("Proportion")+
    theme_bw() +
    theme(
        legend.position = 'top',
        legend.title = element_blank() ,
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 8,color="black",face='bold'),
        legend.text = element_text(size = 8,color="black",face='bold'),
        axis.text.y = element_text(size = 8,color="black",face='bold'),
        axis.title.x = element_text(size = 8,color="black",face='bold'),
        axis.title.y = element_text(size = 8,color="black",face='bold'),
        axis.text.x = element_text(size = 8,color="black",face='bold') ,
        axis.ticks.length = unit(0.2, "cm") ,
        axis.line = element_line(size = 0.5)) 

out_name <- paste0( out_path , "/Fig7C2.pdf"  )
ggsave(file=out_name,plot=p1,width=4,height=4)


