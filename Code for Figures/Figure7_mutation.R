##########################################################################################

library(data.table)
library(optparse)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(ggsci)

###########################################################################################

dir.create(out_path , recursive = T)
col <- c(
        rgb(red=179,green=34,blue=35,alpha=255,max=255), 
        rgb(red=2,green=100,blue=190,alpha=255,max=255) 
    )

###########################################################################################

dat_info <- data.frame(fread(sample_info))
dat_ccf <- data.frame(fread(ccf_file))
dat_gene <- data.frame(fread(gene_list_file))

###########################################################################################

Variant_Type <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Splice_Site","Nonstop_Mutation")
gene_list <- c("MUC6")

###########################################################################################

dat_ccf$vid <- paste0( dat_ccf$Hugo_Symbol, "\n" , dat_ccf$Chr , ":" , dat_ccf$Start_Position , ":" , dat_ccf$REF , ":" , dat_ccf$ALT  , "\n" , dat_ccf$Variant_Classification )
dat_ccf_all <- subset( dat_ccf , Variant_Classification %in% Variant_Type & Hugo_Symbol %in% gene_list )
dat_ccf_driver <- subset( dat_ccf , Variant_Classification %in% Variant_Type & Hugo_Symbol %in% dat_gene$Gene_Symbol )
recurrent_mut <- names(which(table(dat_ccf_all$vid) > 5))
dat_ccf <- subset( dat_ccf_all , vid %in% recurrent_mut )

###########################################################################################

result <- c()
for( geneN in unique(dat_ccf$Hugo_Symbol) ){
    tmp <- subset( dat_ccf , Hugo_Symbol == geneN )
    geneN <- unique(tmp$Hugo_Symbol)
    tmp$Share <- ifelse( tmp$CLS %in% c("clonal [share]" , "subclonal [share]") , "Share" , "Private"  )

    share_sample <- unique(subset( tmp , Share=="Share")$ID)
    im_privatesample <- unique(subset( tmp , Share!="Share" & Class == "IM" )$ID)
    im_privatesample <- im_privatesample[ !im_privatesample %in% share_sample ]
    gc_privatesample <- unique(subset( tmp , Share!="Share" & Class != "IM" )$ID)
    gc_privatesample <- gc_privatesample[ !gc_privatesample %in% share_sample ]

    shareNum <- length(unique(subset( tmp , Share=="Share" )$ID))
    im_privateNum <- length(im_privatesample)
    gc_privateNum <- length(gc_privatesample)

    share_cloneNum <- length(unique(subset( tmp , CLS=="clonal [share]" & Class != "IM" )$ID))
    share_subcloneNum <- length(unique(subset( tmp , CLS=="subclonal [share]" & Class != "IM" )$ID))

    tmp_res <- data.frame( geneN = geneN , vid = geneN , 
        shareNum = shareNum , im_privateNum = im_privateNum , gc_privateNum = gc_privateNum ,
        share_cloneNum = share_cloneNum , share_subcloneNum = share_subcloneNum
        )
    result <- rbind( result , tmp_res )
}

result$vid <- ifelse( result$vid == "MUC6" , "MUC6\nS2130del\nS2266del" , result$vid )

result1 <- data.frame(
    vid = result$vid , Num = result$shareNum , type = "Share"
    )
result2 <- data.frame(
    vid = result$vid , Num = result$im_privateNum , type = "IM branch"
    )
result3 <- data.frame(
    vid = result$vid , Num = result$gc_privateNum , type = "GC branch"
    )
result_plot <- rbind( result1 , result2 , result3 )
totol_num <- result_plot %>%
group_by(vid) %>%
summarize( totol_num = sum(Num) )

result_plot <- merge( result_plot , totol_num , by="vid" )
result_plot$Ratio <- result_plot$Num/result_plot$totol_num
result_plot$value_text <- paste0( round(result_plot$Ratio , 2) * 100 , "%")
result_plot$vid <- factor( result_plot$vid , levels = c("MUC6\nS2130del\nS2266del") )
result_plot <- subset( result_plot , Ratio != 0 )

cols <- c(
    rgb(red=179,green=34,blue=35,alpha=255,max=255), 
    rgb(red=2,green=100,blue=190,alpha=255,max=255) 
)

names(cols) <- c("Share" , "IM branch" )


print("Share IM branch")
print( result_plot )

plot <- ggplot( data = result_plot , aes( x = vid , y = Ratio , fill = type ))+
    geom_bar(position = "stack", stat = "identity") + 
    theme_bw()+
    labs(x="",y="Proportion")+
    theme(panel.grid = element_blank())+
    scale_fill_manual(values = cols[c(1,2)]) +
    ylim(0,1.05)+
    geom_text(aes(label=paste0(type, "\n" , value_text)) , position=position_stack(vjust = 0.5) , size=3 , color="white")+
    theme(panel.background = element_blank(),
        legend.position ='none',
        legend.title = element_blank() ,
        panel.grid.major=element_line(colour=NA),
        plot.title = element_text(size = 12,color="black",face='bold'),
        legend.text = element_text(size = 12,color="black",face='bold'),
        axis.text.y = element_text(size = 12,color="black",face='bold'),
        axis.title.x = element_text(size = 12,color="black",face='bold'),
        axis.title.y = element_text(size = 12,color="black",face='bold'),
        axis.text.x = element_text(size = 12,color="black",face='bold') ,
        axis.ticks.length = unit(0.2, "cm") ,
        axis.line = element_line(size = 0.5)
        )

out_name <- paste0(out_path , "/Fig7B.Share_Private.pdf")
ggsave( out_name , plot , width = 1.5 , height = 4 )

result1 <- data.frame(
    vid = result$vid , Num = result$share_cloneNum , type = "Clone"
    )
result2 <- data.frame(
    vid = result$vid , Num = result$share_subcloneNum , type = "Subclone"
    )
result_plot <- rbind( result1 , result2)
totol_num <- result_plot %>%
group_by(vid) %>%
summarize( totol_num = sum(Num) )

result_plot <- merge( result_plot , totol_num , by="vid" )
result_plot$Ratio <- result_plot$Num/result_plot$totol_num
result_plot$value_text <- paste0( round(result_plot$Ratio , 2) * 100 , "%")

names(cols) <- c("Clone" , "Subclone")
result_plot$vid <- factor( result_plot$vid , levels = c("MUC6\nS2130del\nS2266del") )
result_plot$type <- factor( result_plot$type , levels = c("Subclone" , "Clone") )

print("Subclone  Clone")
print( result_plot )

plot <- ggplot( data = result_plot , aes( x = vid , y = Ratio , fill = type ))+
    geom_bar(position = "stack", stat = "identity") + 
    theme_bw()+
    labs(x="",y="Proportion")+
    theme(panel.grid = element_blank())+
    scale_fill_manual(values = cols[1:2]) +
    ylim(0,1.05)+
    geom_text(aes(label=paste0(type, "\n" , value_text)) , position=position_stack(vjust = 0.5) , size=3 , color="white")+
    theme(panel.background = element_blank(),
        legend.position ='none',
        legend.title = element_blank() ,
        panel.grid.major=element_line(colour=NA),
        plot.title = element_text(size = 12,color="black",face='bold'),
        legend.text = element_text(size = 12,color="black",face='bold'),
        axis.text.y = element_text(size = 12,color="black",face='bold'),
        axis.title.x = element_text(size = 12,color="black",face='bold'),
        axis.title.y = element_text(size = 12,color="black",face='bold'),
        axis.text.x = element_text(size = 12,color="black",face='bold') ,
        axis.ticks.length = unit(0.2, "cm") ,
        axis.line = element_line(size = 0.5)
        )

out_name <- paste0(out_path , "/Fig7B.Share_Clone.pdf")
ggsave( out_name , plot , width = 1.5 , height = 4 )

all_mut_all <- dat_ccf_all %>%
group_by( Hugo_Symbol , Class ) %>%
summarize( sample_num_all = length(unique(ID)) )
all_mut_all <- subset( all_mut_all , Class == "IM" )
all_mut_all <- all_mut_all[,-2]

all_mut <- dat_ccf %>%
group_by( Hugo_Symbol , Class ) %>%
summarize( sample_num_rec = length(unique(ID)) )
all_mut <- subset( all_mut , Class == "IM" )
all_mut <- all_mut[,-2]

dat_plot <- merge( all_mut , all_mut_all )
dat_plot$sample_num_nonrec <-  dat_plot$sample_num_all - dat_plot$sample_num_rec

dat_plot1 <- dat_plot[,c(1:3)]
dat_plot2 <- dat_plot[,c(1,4,3)]
colnames(dat_plot1) <- c("Hugo_Symbol" , "sample_num" , "all_num")
colnames(dat_plot2) <- c("Hugo_Symbol" , "sample_num" , "all_num")
dat_plot1$type <- "Recurrent"
dat_plot2$type <- "Other"

dat_plot <- rbind( dat_plot1 , dat_plot2 )
dat_plot$Ratio <- dat_plot$sample_num/dat_plot$all_num
dat_plot$value_text <- paste0( round(dat_plot$Ratio , 2) * 100 , "%")

names(cols) <- c("Recurrent" , "Other")
dat_plot$Hugo_Symbol <- factor( dat_plot$Hugo_Symbol , levels = c("MUC6") )

print("Recurrent  Other")
print( dat_plot )

plot <- ggplot( data = dat_plot , aes( x = Hugo_Symbol , y = Ratio , fill = type ))+
    geom_bar(position = "stack", stat = "identity") + 
    theme_bw()+
    labs(x="",y="Proportion")+
    theme(panel.grid = element_blank())+
    scale_fill_manual(values = cols[1:2]) +
    ylim(0,1.05)+
    geom_text(aes(label=paste0(type, "\n" , value_text)) , position=position_stack(vjust = 0.5) , size=3 , color="white")+
    theme(panel.background = element_blank(),
        legend.position ='none',
        legend.title = element_blank() ,
        panel.grid.major=element_line(colour=NA),
        plot.title = element_text(size = 12,color="black",face='bold'),
        legend.text = element_text(size = 12,color="black",face='bold'),
        axis.text.y = element_text(size = 12,color="black",face='bold'),
        axis.title.x = element_text(size = 12,color="black",face='bold'),
        axis.title.y = element_text(size = 12,color="black",face='bold'),
        axis.text.x = element_text(size = 12,color="black",face='bold') ,
        axis.ticks.length = unit(0.2, "cm") ,
        axis.line = element_line(size = 0.5)
        )


out_name <- paste0(out_path , "/Fig7B.Recurrent_Ratio.pdf")
ggsave( out_name , plot , width = 1.5 , height = 4 )

###########################################################################################

result <- c()

for( geneN in unique(dat_ccf$Hugo_Symbol) ){
    tmp <- subset( dat_ccf , Hugo_Symbol == geneN )
    geneN <- unique(tmp$Hugo_Symbol)
    tmp$Share <- ifelse( tmp$CLS %in% c("clonal [share]" , "subclonal [share]") , "Share" , "Private"  )

    tmp <- subset(tmp , Share == "Share" & Class != "IM")

    tmp_driver <- subset( dat_ccf_driver , Sample %in% tmp$Sample )

    shareNum <- length(unique(subset( tmp , Share=="Share" )$ID))
    gcDriverNum <- length(unique(tmp_driver$ID))

    tmp_res <- data.frame( geneN = geneN , vid = geneN , 
        shareNum = shareNum , gcDriverNum = gcDriverNum
        )
    result <- rbind( result , tmp_res )
}

result$vid <- ifelse( result$vid == "MUC6" , "MUC6\nS2130del\nS2266del" , result$vid )

result1 <- data.frame(
    vid = result$vid , Num = result$gcDriverNum , type = "Driver\nCon-occurence" , totol_num = result$shareNum
    )
result2 <- data.frame(
    vid = result$vid , Num = result$shareNum - result$gcDriverNum , type = "Other" , totol_num = result$shareNum
    )
result_plot <- rbind( result1 , result2)

result_plot$Ratio <- result_plot$Num/result_plot$totol_num
result_plot$value_text <- paste0( round(result_plot$Ratio , 2) * 100 , "%")
result_plot$vid <- factor( result_plot$vid , levels = c("MUC6\nS2130del\nS2266del") )
result_plot$type <- factor( result_plot$type , levels = c("Other" , "Driver\nCon-occurence") )

col <- c(
    rgb(red=179,green=34,blue=35,alpha=255,max=255), 
    rgb(red=2,green=100,blue=190,alpha=255,max=255) 
)
names(cols) <- c("Driver\nCon-occurence" , "Other")


print("Con-occurence")
print( result_plot )

plot <- ggplot( data = result_plot , aes( x = vid , y = Ratio , fill = type ))+
    geom_bar(position = "stack", stat = "identity") + 
    theme_bw()+
    labs(x="",y="")+
    theme(panel.grid = element_blank())+
    scale_fill_manual(values = cols[1:2]) +
    ylim(0,1.05)+
    geom_text(aes(label=paste0(type, "\n" , value_text)) , position=position_stack(vjust = 0.5) , size=3 , color="white")+
    theme(panel.background = element_blank(),#设置背影为白色#清除网格线
        legend.position ='none',
        legend.title = element_blank() ,
        panel.grid.major=element_line(colour=NA),
        plot.title = element_text(size = 12,color="black",face='bold'),
        legend.text = element_text(size = 12,color="black",face='bold'),
        axis.text.y = element_text(size = 12,color="black",face='bold'),
        axis.title.x = element_text(size = 12,color="black",face='bold'),
        axis.title.y = element_text(size = 12,color="black",face='bold'),
        axis.text.x = element_text(size = 12,color="black",face='bold') ,
        axis.ticks.length = unit(0.2, "cm") ,
        axis.line = element_line(size = 0.5)
        )

out_name <- paste0(out_path , "/Fig7B.Driver_Con-occurence.pdf")
ggsave( out_name , plot , width = 1.5 , height = 4 )

