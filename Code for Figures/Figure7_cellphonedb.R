##########################################################################################

library(ktplots)
library(optparse)

##########################################################################################

pvals <- read.delim(pvals_file, check.names = FALSE)
means <- read.delim(means_file, check.names = FALSE)

##########################################################################################

cell_order <- c("Enterocytes" , "Pit_Mut" , "Pit_Other" , "Neck" , "Goblet" , "Chief" , "Endocrine" , "Parietal" , "Tumor")

pvalue.threshold = 0.05
order.of.celltype = cell_order[length(cell_order):1]
ccc.number.max = 10
size.of.text = c(18,14,14,12)
color.palette = c("#4393C3","#ffdbba","#B2182B")

pvalues=read.table(pvals_file,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
pvalues=pvalues[,12:dim(pvalues)[2]] 
statdf=as.data.frame(colSums(pvalues < pvalue.threshold)) 
colnames(statdf)=c("number")

statdf$indexb=stringr::str_replace(rownames(statdf),"^.*\\|","")
statdf$indexa=stringr::str_replace(rownames(statdf),"\\|.*$","")
statdf$total_number=0

for (i in 1:dim(statdf)[1]) {
  tmp_indexb=statdf[i,"indexb"]
  tmp_indexa=statdf[i,"indexa"]
  if (tmp_indexa == tmp_indexb) {
    statdf[i,"total_number"] = statdf[i,"number"]
  } else {
    statdf[i,"total_number"] = statdf[statdf$indexb==tmp_indexb & statdf$indexa==tmp_indexa,"number"]+
      statdf[statdf$indexa==tmp_indexb & statdf$indexb==tmp_indexa,"number"]
  }
}

if(is.null(order.of.celltype)){
	rankname=sort(unique(statdf$indexa))
} else {
	rankname=order.of.celltype
}

statdf$indexa=factor(statdf$indexa,levels = rankname[length(rankname):1])
statdf$indexb=factor(statdf$indexb,levels = rankname)

if(is.null(ccc.number.max)){
	limits.max = ceiling(max(statdf$total_number) / 10) * 10
}else{
	statdf$total_number[statdf$total_number > ccc.number.max] = ccc.number.max
	limits.max = ccc.number.max
}

p <-ggplot(statdf,aes(x=indexa,y=indexb,fill=total_number))+
	geom_tile(color="white")+
	scale_fill_gradientn(colours = color.palette,limits=c(0,limits.max) , name="number of interactions")+
	scale_x_discrete("")+
	scale_y_discrete("")+
	theme_minimal()+
	theme(
	  axis.title = element_text(size = size.of.text[1]),
	  axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45,size = 14,color = "black",face='bold'),
	  axis.text.y.left = element_text(size = 14,color = "black",face='bold'),
	  legend.title = element_text(size = size.of.text[3]),
	  legend.position = "right" ,
	  panel.grid = element_blank()
	)

out_name <- paste0( out_path , "/Fig7F1.pdf" )
ggsave( out_name , p , width = 8 , height = 6 )

##########################################################################################

a_nocomplex <- grep( "complex" , pvals$partner_a , invert = T)
b_nocomplex <- grep( "complex" , pvals$partner_b , invert = T)

use_index <- a_nocomplex[a_nocomplex %in% b_nocomplex]

pvals_use <- pvals[use_index,]
means_use <- means[use_index,]

neg_log10_th= -log10(0.05) 
means_exp_log2_th=0 
notused.cell=NULL 
used.cell=NULL 
neg_log10_th2=3 
means_exp_log2_th2=c(-2,0.75)
cell.pair=NULL 
cell.pair=c(
  "Enterocytes|Enterocytes","Enterocytes|Pit_Mut","Enterocytes|Pit_Other",
	"Pit_Mut|Enterocytes","Pit_Mut|Pit_Mut" , "Pit_Mut|Pit_Other" ,
	"Pit_Other|Enterocytes" , "Pit_Other|Pit_Mut" , "Pit_Other|Pit_Other"
	)

gene.pair=NULL 
color_palette = c("#313695", "#4575B4", "#ABD9E9", "#FFFFB3", "#FDAE61", "#F46D43", "#D73027", "#A50026")
text_size = 12

pvalues=pvals_use
pvalues=pvalues[,c(2,12:dim(pvalues)[2])]
RMpairs=names(sort(table(pvalues$interacting_pair))[sort(table(pvalues$interacting_pair)) > 1])
pvalues=pvalues[!(pvalues$interacting_pair %in% RMpairs),]
pvalues.df1=reshape2::melt(pvalues,id="interacting_pair")
colnames(pvalues.df1)=c("geneA_geneB","cellA_cellB","pvalue")
pvalues.df1$neg_log10=-log10(pvalues.df1$pvalue)
pvalues.df1$geneA_geneB_cellA_cellB=paste(pvalues.df1$geneA_geneB,pvalues.df1$cellA_cellB,sep = ",")

means=means_use
means=means[,c(2,12:dim(means)[2])]
rmpairs=names(sort(table(means$interacting_pair))[sort(table(means$interacting_pair)) > 1])
means=means[!(means$interacting_pair %in% rmpairs),]
means.df1=reshape2::melt(means,id="interacting_pair")
colnames(means.df1)=c("geneA_geneB","cellA_cellB","means_exp")
means.df1$geneA_geneB_cellA_cellB=paste(means.df1$geneA_geneB,means.df1$cellA_cellB,sep = ",")
means.df1=means.df1[,c("geneA_geneB_cellA_cellB","means_exp")]

raw.df=merge(pvalues.df1,means.df1,by="geneA_geneB_cellA_cellB")
raw.df$means_exp_log2=log2(raw.df$means_exp)
raw.df$means_exp_log2=raw.df$means_exp

raw.df <- subset( raw.df , cellA_cellB %in% cell.pair )
final.df=dplyr::filter(raw.df,neg_log10 > neg_log10_th & means_exp_log2 > means_exp_log2_th)

final.df$geneA=stringr::str_replace(final.df$geneA_geneB,"_.*$","")
final.df$geneB=stringr::str_replace(final.df$geneA_geneB,"^.*_","") 
final.df$cellA=stringr::str_replace(final.df$cellA_cellB,"\\|.*$","") 
final.df$cellB=stringr::str_replace(final.df$cellA_cellB,"^.*\\|","")

 if (!is.null(notused.cell)) {
  final.df=final.df[!(final.df$cellA %in% notused.cell),]
  final.df=final.df[!(final.df$cellB %in% notused.cell),]
}

final.df=final.df[!(final.df$cellA==final.df$cellB),]

final.df.gene=unique(final.df$geneA_geneB)
final.df.cell=unique(final.df$cellA_cellB)
if (!is.null(used.cell)){
  tmp_cell=c()
  for (i in used.cell) {
    tmp_cell=union(tmp_cell,final.df.cell[str_detect(final.df.cell,i)])
  }
  final.df.cell=tmp_cell
}
raw.df=raw.df[raw.df$geneA_geneB %in% final.df.gene, ]
raw.df=raw.df[raw.df$cellA_cellB %in% final.df.cell, ]

raw.df$neg_log10=ifelse(raw.df$neg_log10 > neg_log10_th2,neg_log10_th2,raw.df$neg_log10)
raw.df$means_exp_log2=ifelse(
  raw.df$means_exp_log2 > means_exp_log2_th2[2],
  means_exp_log2_th2[2],
  ifelse(
    raw.df$means_exp_log2 < means_exp_log2_th2[1],
    means_exp_log2_th2[1],
    raw.df$means_exp_log2
  )
)

raw.df$cellA_cellB=as.character(raw.df$cellA_cellB)
if (!is.null(cell.pair)) {
  tmp_pair=intersect(cell.pair,unique(raw.df$cellA_cellB))
  raw.df=raw.df[raw.df$cellA_cellB %in% tmp_pair,]
  raw.df$cellA_cellB=factor(raw.df$cellA_cellB,levels = tmp_pair)
} else {
  tmp_pair=sort(unique(raw.df$cellA_cellB))
  raw.df$cellA_cellB=factor(raw.df$cellA_cellB,levels = tmp_pair)
}
raw.df$geneA_geneB=as.character(raw.df$geneA_geneB)
if (!is.null(gene.pair)) {
  tmp_pair=intersect(gene.pair,unique(raw.df$geneA_geneB))
  raw.df=raw.df[raw.df$geneA_geneB %in% tmp_pair,]
  raw.df$geneA_geneB=factor(raw.df$geneA_geneB,levels = tmp_pair)
} else {
  tmp_pair=sort(unique(raw.df$geneA_geneB))
  raw.df$geneA_geneB=factor(raw.df$geneA_geneB,levels = tmp_pair)
}

p <- ggplot(raw.df,aes(cellA_cellB,geneA_geneB))+
    geom_point(aes(size=neg_log10,color=means_exp_log2))+
    scale_color_gradientn("scaled_means",colors = color_palette)+
    scale_size_continuous("-log10(p value)")+
    theme_bw()+
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text.x.bottom = element_text(hjust = 1, angle = 45, size=14, color = "black" , face = "bold"),
      axis.text.y.left = element_text(size = 14,color = "black" , face = "bold"),
      axis.line = element_line(size = 0.5),
      axis.ticks.length = unit(0.15,"cm") 
    )

out_name <- paste0( out_path , "/Fig7F2.pdf" )
ggsave( out_name , p , width = 8 , height = 6 )

