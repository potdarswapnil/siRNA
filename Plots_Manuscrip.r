options(java.parameters = "-Xmx6g")
library(gdata)
library(MASS)
library(raster)
library(stats)
library(gplots)
library(graphics)
require(gridExtra)
library(reshape2)
library(xlsx)
library(XLConnect)
library(reshape)
library(plyr)

source("C:\\Swapnil\\R\\qc2pdf.r")
setwd("C:\\Swapnil\\siRNA\\Results_Secondary")


hmcols<-colorRampPalette(c("blue","white","red"))(256)
mycol=c("sample"="gray50","OM"="green1","cells only"="darkgreen","cells"="darkgreen","cellsTR"="green1","antibody block"="darkorchid1","internal1"="darkorchid1","internal2"="darkorchid2","internal3"="darkorchid3","internal4"="darkorchid4","internal5"="magenta1","internal6"="magenta2","internal7"="magenta3","internal8"="magenta4","internal9"="mediumorchid","A_NEG_2 (TR+OM)"="red1","neg1"="red1","neg2"="red4","neg3"="firebrick1","neg4"="firebrick2","mneg1"="tomato1","mneg2"="tomato3","Q_DEATH_1 (TR+OM)"="blue1","pos1"="blue1","pos2"="slateblue1","mpos1"="steelblue1","mpos2"="steelblue2","blanks"="cyan","Low_Bac_Nuc_Ratio"="mediumvioletred","High_Bac_Nuc_Ratio"="slateblue1")
cellcol=c("darkgreen","brown","magenta1","red1","blue1","cyan","lightblue","orange","lightgreen","tomato","dimgray","mediumorchid4","darkolivegreen4","palevioletred")
theme_set(theme_bw())
old.par <- par(no.readonly = TRUE)

#######################################
#######################################
load(file="..\\Annotations\\library.rda")
load(file="data_table.rda")
load(file="total_set.rda")

data_table$content=data_table$Gene_Symbol
data_table$content<-ifelse(data_table$content %in% c("A_NEG_2 (TR+OM)","OM","Q_DEATH_1 (TR+OM)","antibody block","cells only"),as.character(data_table$content),"sample")

jpeg("bac_nuc_ratio_ab_ratio.jpeg", width = 24, height =12, units = 'in', res = 300)

dt=data_table
dt=dt[order(dt$bac_nuc_ratio, decreasing=F),][1:400,]
dt_ab=data_table
dt_ab=dt_ab[order(dt_ab$bac_nuc_ratio, decreasing=F),][1:400,]
dt_ab=dt_ab[,c("Gene_Symbol","kgXref_gene_symbol","content","bac_nuc_ratio","bac_nuc_ratio_neg_ratio","bac_nuc_ratio_ab_ratio")]

mean_dt=subset(dt_ab, FALSE)
mean_dt$sirnas=numeric()
for(gs in unique(dt_ab$Gene_Symbol))
{
  gs_table <- subset(dt_ab,Gene_Symbol==gs)
  if(nrow(gs_table) > 1)
  {
    sirna_numbers=ifelse(gs_table$content %in% c("A_NEG_2 (TR+OM)","OM","Q_DEATH_1 (TR+OM)","antibody block","cells only"),5,nrow(gs_table))
    mean_gs=data.frame("Gene_Symbol"=gs,"kgXref_gene_symbol"=head(gs_table$kgXref_gene_symbol,1),"content"=unique(gs_table$content),"bac_nuc_ratio"=mean(gs_table$bac_nuc_ratio),"bac_nuc_ratio_neg_ratio"=mean(gs_table$bac_nuc_ratio_neg_ratio),"bac_nuc_ratio_ab_ratio"=mean(gs_table$bac_nuc_ratio_ab_ratio),"sirnas"=max(sirna_numbers))
    mean_dt=rbind(mean_dt,mean_gs)
  }
}
p <- ggplot(data=mean_dt, aes(Gene_Symbol,bac_nuc_ratio_ab_ratio)) + geom_point(aes(color = factor(content),size=sirnas)) + scale_size(range = c(2,6)) + geom_text(aes(label=kgXref_gene_symbol),hjust=0, vjust=0,size=4) + scale_colour_manual(values = mycol) + ggtitle(paste("bac_nuc_ratio_ab_ratio With Latest Annotations")) + theme(axis.text.x = element_text(angle = 90, hjust=0,vjust=0)) + labs(x="Gene Symbol",y="bac_nuc_ratio_ab_ratio")
print(p)
dev.off()



