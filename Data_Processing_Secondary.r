options(java.parameters = "-Xmx8g")
options(java.parameters ="-XX:-UseConcMarkSweepGC")
#options(error = browser)
library(gdata)
library(MASS)
library(raster)
library(stats)
library(gplots)
library(graphics)
require(graphics)
require(biwt)
require(gridExtra)
library(reshape2)
library(xlsx)
library(XLConnect)
library(reshape)
library(Biobase)
library(plyr)
source("C:\\Swapnil\\R\\qc2pdf.r")
setwd("..\\siRNA\\Results_Secondary")


load(file="C:\\Swapnil\\Annotations\\siRNA\\library.rda")
load(file="C:\\Swapnil\\Annotations\\siRNA\\new_ambion_annoframe.rda")
newannoframe=new_ambion_annoframe

plateframe= read.xlsx("..\\Annotations\\Controls_siRNA_Secondary.xlsx", 1, rowIndex=2:17,colIndex=2:25, as.data.frame=TRUE, header=FALSE, colClasses="character")
platemat=as.matrix(plateframe)
file_names <- list.files("..\\Data_Secondary",pattern = ".xlsx", all.files = TRUE,full.names = TRUE, recursive = TRUE,ignore.case = TRUE, include.dirs = TRUE)

ecpframe= read.xlsx2("..\\Annotations\\ECP_Secondary.xlsx", 1, rowIndex=1:1128,colIndex=1:8, as.data.frame=TRUE, header=TRUE, colClasses="character")
#load(file="ecpframe.rda")

#Main/All/Mean Intensity DAPI/Objects     Main/All/Mean Intensity DAPI/Mean	    baku_green/All/Mean Intensity Alexa488-QUAD/Objects	baku_green/All/Mean Intensity Alexa488-QUAD/Mean	bac/nuc
#Main/All/Mean Intensity Hoechst/Objects  Main/All/Mean Intensity Hoechst/Mean	baku_green/All/Mean Intensity Alexa488/Objects	    baku_green/All/Mean Intensity Alexa488/Mean	bac/nyc

data_table <- do.call("rbind",lapply(file_names,function(file_name,skip=1){
  dataframe= read.xlsx(file_name, 1, rowIndex=1:385, colIndex= 1:15,as.data.frame=TRUE, header=TRUE, colClasses=NA,stringsAsFactors=FALSE)
  
  #Index  Row	Column	Main/All/Mean Intensity DAPI/Objects	Main/All/Mean Intensity DAPI/Mean	baku_green/All/Mean Intensity Alexa488-QUAD/Objects	baku_green/All/Mean Intensity Alexa488-QUAD/Mean	bac/nuc	Swell	Plate.ID_1	Plate.ID_2	Plate.ID_3	RefSeq.Accession.Number_1	Gene.Symbol_1	Full.Gene.Name_1
  #Index  Row	Column	Main/All/Mean Intensity DAPI/Objects	Main/All/Mean Intensity DAPI/Mean	baku_green/All/Mean Intensity Alexa488-QUAD/Objects	baku_green/All/Mean Intensity Alexa488-QUAD/Mean	bac/nuc	Swell	Plate.ID_1	Plate.ID_2	Plate.ID_3	RefSeq.Accession.Number_1	Gene.Symbol_1	Full.Gene.Name_1
  #Controls###A_NEG_2 (TR+OM) #OM #Q_DEATH_1 (TR+OM)#antibody block
  ##Index  Description	Main/All/Mean Intensity Hoechst/Objects	Main/All/Mean Intensity Hoechst/Mean	baku_green/All/Mean Intensity Alexa488/Objects	      baku_green/All/Mean Intensity Alexa488/Mean	      bac/nyc
  ##Index  Row	Column	Main/All/Mean Intensity DAPI/Objects	  Main/All/Mean Intensity DAPI/Mean	    baku_green/All/Mean Intensity Alexa488-QUAD/Objects	  baku_green/All/Mean Intensity Alexa488-QUAD/Mean	bac/nuc	  Swell	Plate.ID_1	Plate.ID_2	Plate.ID_3	RefSeq.Accession.Number_1	Gene.Symbol_1	Full.Gene.Name_1
  #bac/nuc    Swell	Plate.ID_1	Plate.ID_2	Plate.ID_3	RefSeq.Accession.Number_1	Gene.Symbol_1	Full.Gene.Name_1
  colnames(dataframe)=c("Index","Row","Column","Dapi_Objects","Dapi_objects_intensity_mean","baku_green_objects","baku_green_intensity_mean","bac_nuc","SWell","Plate.ID_1","Plate.ID_2","Plate.ID_3","RefSeq","Gene_Symbol","Gene_Name")
    
  dataframe$bac_nuc_ratio=round(dataframe$baku_green_objects/dataframe$Dapi_Objects,1)
  dataframe$DWell=paste(dataframe$Row,dataframe$Column,sep="")
  fhead=unlist(strsplit(file_name,"/",fixed=TRUE))
  plate_names <- fhead[length(fhead)]
  plate_names <- gsub(".xlsx","",plate_names)
  dataframe$Plate_Id <- plate_names
  rownames(dataframe) <- paste(dataframe$Plate_Id,dataframe$DWell,sep="")
  median_neg=median(dataframe$bac_nuc_ratio[which(dataframe$Gene_Symbol=="A_NEG_2 (TR+OM)")])
  dataframe$bac_nuc_ratio_neg_ratio=round(dataframe$bac_nuc_ratio/median_neg,3)
 
  median_ab=median(dataframe$bac_nuc_ratio[which(dataframe$Gene_Symbol=="antibody block")])
  dataframe$bac_nuc_ratio_ab_ratio=round(dataframe$bac_nuc_ratio/median_ab,3)
  
  
  ecp=ecpframe[ecpframe$Destination.Plate.Barcode==plate_names,]
  dataframe=merge(dataframe,ecp,by.x="DWell",by.y="Destination.Well",all.x=TRUE)
  
  dataframe$Source.Plate.Barcode=gsub("Ambion-A","X",dataframe$Source.Plate.Barcode)
  dataframe$Source.Plate.Barcode=gsub("Ambion-B","X",dataframe$Source.Plate.Barcode)
  dataframe$Source.Plate.Barcode=gsub("Ambion-C","X",dataframe$Source.Plate.Barcode)
  dataframe=merge(dataframe,newannoframe,by.x=c("Source.Plate.Barcode","Source.Well"),by.y=c("PP.plate","SWell"),all.x=TRUE)
  dataframe <- with(dataframe, dataframe[!(DWell == "" | is.na(DWell)), ])
  dataframe 
}))

fname <- paste("Combined_Data",".txt",sep="")
write.table(data_table,file=fname,append=FALSE,quote=FALSE,sep="\t",eol="\n",na="NA",dec=".",row.names=TRUE,col.names=NA,fileEncoding="")

fname <- paste("Combined_Data",".xlsx",sep="")
write.xlsx2(data_table,fname,sheetName="Data",row.names=TRUE, append=FALSE)

total_set=data_table
total_set=data_table[,c("SWell","DWell","Plate_Id","Dapi_Objects","Dapi_objects_intensity_mean","baku_green_objects","baku_green_intensity_mean","bac_nuc_ratio","bac_nuc_ratio_neg_ratio","bac_nuc_ratio_ab_ratio","Gene_Symbol","Gene_Name","RefSeq","RefSeq_67","refGene_gene_symbol","kgXref_description")]

fname <- paste("total_data",".xlsx",sep="")
write.xlsx2(total_set,fname,sheetName="Data",row.names=TRUE, append=FALSE)

# fname <- paste("total_data",".txt",sep="")
# write.table(total_set,file=fname,append=FALSE,quote=FALSE,sep="\t",eol="\n",na="NA",dec=".",row.names=TRUE,col.names=NA,fileEncoding="")

save(total_set,file ="total_set.rda")
save(data_table,file ="data_table.rda")