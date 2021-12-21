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
#library(Biobase)
library(plyr)
source("C:\\Swapnil\\R\\qc2pdf.r")
setwd("C:\\Swapnil\\siRNA\\Results")

load(file="..\\Annotations\\library.rda")
load(file="data_table.rda")
load(file="total_set.rda")

jpeg("Dapi_Objects_vs_Baku_green_objects.jpeg", width = 24, height =12, units = 'in', res = 300)

hmcols<-colorRampPalette(c("blue","white","red"))(256)
mycol=c("sample"="gray50","cells"="darkgreen","cellsTR"="green1","internal1"="darkorchid1","internal2"="darkorchid2","internal3"="darkorchid3","internal4"="darkorchid4","internal5"="magenta1","internal6"="magenta2","internal7"="magenta3","internal8"="magenta4","internal9"="mediumorchid","neg1"="red1","neg2"="red4","neg3"="firebrick1","neg4"="firebrick2","mneg1"="tomato1","mneg2"="tomato3","pos1"="blue1","pos2"="slateblue1","mpos1"="steelblue1","mpos2"="steelblue2","blanks"="cyan","Low_Bac_Nuc_Ratio"="mediumvioletred","High_Bac_Nuc_Ratio"="slateblue1")
mycol2=c("sample"="gray50","cells"="darkgreen","cellsTR"="green1","internal1"="grey50","internal2"="turquoise1","internal3"="khaki1","internal4"="hotpink","internal5"="coral4","internal6"="navy","internal7"="MediumSpringGreen","internal8"="magenta4", "internal9"="darkviolet","neg1"="red1","neg2"="red4","neg3"="firebrick1","neg4"="firebrick2","mneg1"="tomato1","mneg2"="tomato3","pos1"="blue1","pos2"="slateblue1","mpos1"="steelblue1","mpos2"="steelblue2","blanks"="cyan")
cellcol=c("darkgreen","brown","magenta1","red1","blue1","cyan","lightblue","orange","lightgreen","tomato","dimgray","mediumorchid4","darkolivegreen4","palevioletred")
theme_set(theme_bw())
old.par <- par(no.readonly = TRUE)

data_table$object_cutoff <- NA
data_table$ratio_cutoff <- ""
data_table$Low_Bac_Nuc_Ratio <- NA
data_table$High_Bac_Nuc_Ratio <- NA
data_table$label_list <- NA

###Top 500 hits cutoff for bac_nuc_ratio
dt=data_table[data_table$content=="sample",]
dt=dt[order(dt$bac_nuc_ratio, decreasing=F),][1:500,]
o_lower_cutoff=max(dt$bac_nuc_ratio)

dt=data_table[data_table$content=="sample",]
dt=dt[order(dt$bac_nuc_ratio, decreasing=T),][1:300,]
o_upper_cutoff=min(dt$bac_nuc_ratio)

###Top 500 hits cutoff for bac_nuc_ratio_neg_ratio
dt=data_table[data_table$content=="sample",]
dt=dt[order(dt$bac_nuc_ratio_neg_ratio, decreasing=F),][1:500,]
r_lower_cutoff=max(dt$bac_nuc_ratio_neg_ratio)
dt_Low_Bac_Nuc_Ratio=dt

dt=data_table[data_table$content=="sample",]
dt=dt[order(dt$bac_nuc_ratio_neg_ratio, decreasing=T),][1:300,]
r_upper_cutoff=min(dt$bac_nuc_ratio_neg_ratio)
dt_bottom300=dt

data_table$ratio_cutoff <- ifelse(data_table$bac_nuc_ratio_neg_ratio < r_lower_cutoff | data_table$bac_nuc_ratio_neg_ratio > r_upper_cutoff,as.character(data_table$Gene.Symbol_1),"")

data_table$object_cutoff <- ifelse(data_table$bac_nuc_ratio < o_lower_cutoff | data_table$bac_nuc_ratio > o_upper_cutoff,as.character(data_table$Gene.Symbol_1),"")

data_table$Low_Bac_Nuc_Ratio <- ifelse(data_table$bac_nuc_ratio_neg_ratio < r_lower_cutoff,"Low_Bac_Nuc_Ratio",data_table$Low_Bac_Nuc_Ratio)
data_table$High_Bac_Nuc_Ratio <- ifelse(data_table$bac_nuc_ratio_neg_ratio > r_upper_cutoff,"High_Bac_Nuc_Ratio",data_table$High_Bac_Nuc_Ratio)

data_table$label_list <- ifelse(data_table$bac_nuc_ratio_neg_ratio < r_lower_cutoff,"Low_Bac_Nuc_Ratio",data_table$label_list)
data_table$label_list <- ifelse(data_table$bac_nuc_ratio_neg_ratio > r_upper_cutoff,"High_Bac_Nuc_Ratio",data_table$label_list)

p <- ggplot(data=data_table, aes(x=Dapi_Objects, y=baku_green_objects)) + geom_point(shape=19,alpha=0.75,size=1,aes(colour=factor(content)),show_guide = TRUE) + scale_colour_manual(values = mycol) + geom_text(size=2,aes(x=Dapi_Objects, y=baku_green_objects,label=ratio_cutoff,colour=factor(label_list)),show_guide = F)+ ggtitle("Dapi_Objects vs baku_green_objects")
print(p)

dev.off()



