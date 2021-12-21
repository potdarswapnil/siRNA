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
setwd("..\\siRNA\\Results")

plateframe= read.xlsx("..\\Annotations\\Controls_siRNA.xlsx", 1, rowIndex=2:17,colIndex=2:25, as.data.frame=TRUE, header=FALSE, colClasses="character")
platemat=as.matrix(plateframe)
file_names <- list.files("..\\Data",pattern = ".xlsx", all.files = TRUE,full.names = TRUE, recursive = TRUE,ignore.case = TRUE, include.dirs = TRUE)

#annoframe= read.xlsx2("..\\Annotations\\111031 AMBION-Silencer Select Human-384 well-ALL.xlsx", 1, rowIndex=1:71425,colIndex=1:26, as.data.frame=TRUE, header=TRUE, colClasses="character")
#annoframe=read.csv("..\\Annotations\\111031 AMBION-Silencer Select Human-384 well-ALL.txt", sep="\t",header = TRUE, skip=0,nrows = 71425)
#save(annoframe,file="library.rda")

load(file="..\\Annotations\\library.rda")

barcodeframe= read.xlsx("..\\Annotations\\Barcode_131204_34.xlsx", 1, rowIndex=1:70,colIndex=1:16, as.data.frame=TRUE, header=TRUE, colClasses="character")
#Plate Name  PP plate	LDV	Lot Number	Plate ID	Sample ID	Location (Row-Col)	Row	Col	RefSeq Accession Number	Gene Symbol	Full Gene Name	Gene ID	siRNA ID	Amount	Exon(s) Targeted	Exon Image	Sense siRNA Sequence	Antisense siRNA Sequence	Validated	Mean RNA Levels Remaining	Plus Error Bar	Minus Error Bar	Validation Cell Line	Matrix 2D Code	Customer Barcode
# colnames(annoframe)
# [1] "Plate.Name"                "PP.plate"                  "LDV"                       "Lot.Number"               
# [5] "Plate.ID"                  "Sample.ID"                 "Location..Row.Col."        "Row"                      
# [9] "Col"                       "RefSeq.Accession.Number"   "Gene.Symbol"               "Full.Gene.Name"           
# [13] "Gene.ID"                   "siRNA.ID"                  "Amount"                    "Exon.s..Targeted"         
# [17] "Exon.Image"                "Sense.siRNA.Sequence"      "Antisense.siRNA.Sequence"  "Validated"                
# [21] "Mean.RNA.Levels.Remaining" "Plus.Error.Bar"            "Minus.Error.Bar"           "Validation.Cell.Line"     
# [25] "Matrix.2D.Code"            "Customer.Barcode"   
  data_table <- do.call("rbind",lapply(file_names,function(file_name,skip=1){
  dataframe= read.xlsx(file_name, 1, rowIndex=1:386, colIndex= 1:6,as.data.frame=TRUE, header=TRUE, colClasses=NA,stringsAsFactors=FALSE)
  
  #Index  Description	Main/All/Mean Intensity Hoechst/Objects	Main/All/Mean Intensity Hoechst/Mean	baku_green/All/Mean Intensity Alexa488/Objects	baku_green/All/Mean Intensity Alexa488/Mean	bac/nuc
  colnames(dataframe)=c("Index","Description","Dapi_Objects","Dapi_objects_intensity_mean","baku_green_objects","baku_green_intensity_mean")
  
  dataframe$bac_nuc_ratio=round(dataframe$baku_green_objects/dataframe$Dapi_Objects,1)
  

  rawdatamat = dataframe[1:384,1:7]
  #fhead=unlist(strsplit(file_name,"\\",fixed=TRUE))
  fhead=unlist(strsplit(file_name,"/",fixed=TRUE))
  plate_names <- fhead[length(fhead)]
  plate_names <- gsub(".xlsx","",plate_names)
  rawdatamat$Plate_Id <- plate_names
  rownames(rawdatamat) <- paste(rawdatamat$Plate_Id,rawdatamat$Description,sep="_")
  barcodeframe$Dest.BC <- gsub("-","_",barcodeframe$Dest.BC)
  
   #barcode=subset(barcodeframe, Dest.BC==plate_names)
  barcode <- as.character(as.matrix(barcodeframe[which(barcodeframe$Dest.BC==plate_names),]))
  srcbc1=barcode[2]
  srcbc2=barcode[5]
  srcbc3=barcode[8]

  annoplate1=subset(annoframe,annoframe$Plate.ID==srcbc1)
  annoplate2=subset(annoframe,annoframe$Plate.ID==srcbc2)
  annoplate3=subset(annoframe,annoframe$Plate.ID==srcbc3)
  colnames(annoplate1) <- paste(colnames(annoplate1),"1",sep="_")
  colnames(annoplate2) <- paste(colnames(annoplate2),"2",sep="_")
  colnames(annoplate3) <- paste(colnames(annoplate3),"3",sep="_")
  annoplate=merge(annoplate1,annoplate2,by.x="Location..Row.Col._1", by.y = "Location..Row.Col._2",all=TRUE)
  annoplate=merge(annoplate,annoplate3,by.x="Location..Row.Col._1", by.y = "Location..Row.Col._3",all=TRUE)
 
  rownames(platemat) <- LETTERS[1:16]
  colnames(platemat) <- 1:24
  plate_anot <- melt(platemat)
  well_names <- paste(plate_anot[,1],plate_anot[,2],sep="")
  
  srcplatemat <- matrix(nrow =16, ncol = 24)
  srcplatemat[,1:23]=platemat[,2:24]
  srcplatemat[,24]=platemat[,1]
  rownames(srcplatemat)=rownames(platemat)
  colnames(srcplatemat)=c(24,1:23)
  srcplate_anot <-melt(srcplatemat)
  src_well_names <- paste(srcplate_anot[,1],srcplate_anot[,2],sep="")
      
  plate_anot_data <-  plate_anot[,3]
  names(plate_anot_data) <- well_names
  blanks <- ifelse(plate_anot_data == "NA","yes","no")
  cells <- ifelse(plate_anot_data == "cells","yes","no")
  cellsoptimem <- ifelse(plate_anot_data == "cells+Optimem","yes","no")
  cellsTR <- ifelse(plate_anot_data == "cells+TR","yes","no")
  internal1 <- ifelse(plate_anot_data == "internal1","yes","no")
  internal2 <- ifelse(plate_anot_data == "internal2","yes","no")
  internal3 <- ifelse(plate_anot_data == "internal3","yes","no")
  internal4 <- ifelse(plate_anot_data == "internal4","yes","no")
  internal5 <- ifelse(plate_anot_data == "internal5","yes","no")
  internal6 <- ifelse(plate_anot_data == "internal6","yes","no")
  internal7 <- ifelse(plate_anot_data == "internal7","yes","no")
  internal8 <- ifelse(plate_anot_data == "internal8","yes","no")
  internal9 <- ifelse(plate_anot_data == "internal9","yes","no")
  neg1 <- ifelse(plate_anot_data == "neg1","yes","no")
  neg2 <- ifelse(plate_anot_data == "neg2","yes","no")
  neg3 <- ifelse(plate_anot_data == "neg3","yes","no")
  neg4 <- ifelse(plate_anot_data == "neg4","yes","no")
  mneg1 <- ifelse(plate_anot_data == "mneg1","yes","no")
  mneg2 <- ifelse(plate_anot_data == "mneg2","yes","no")
  pos1 <- ifelse(plate_anot_data == "pos","yes","no")
  pos2 <- ifelse(plate_anot_data == "pos2","yes","no")
  mpos1 <- ifelse(plate_anot_data == "mpos","yes","no")
  mpos2 <- ifelse(plate_anot_data == "mpos2","yes","no")
  sample <- ifelse(plate_anot_data == "sample","yes","no")
  
  well_annot_tbl <- data.frame(well_names,src_well_names,blanks,cells,cellsoptimem,cellsTR,internal1,internal2,internal3,internal4,internal5,internal6,internal7,internal8,internal9,neg1,neg2,neg3,neg4,mneg1,mneg2,pos1,pos2,mpos1,mpos2,sample,stringsAsFactors=FALSE)
  
  well_annot_tbl$content <- apply(well_annot_tbl,1,function(row){
    i=grep("yes",row)
    content=colnames(well_annot_tbl)[i]
    content
  })
  well_annot_tbl=well_annot_tbl[,c("well_names","src_well_names","content")]
  #by= "Plate ID",all = TRUE
  colnames(annoplate)[1] <- "Well"

  well_annot_tbl_combined <- merge(rawdatamat,well_annot_tbl,by.x="Description",by.y="well_names",all=TRUE)
  
  median_neg= median(well_annot_tbl_combined[as.character(well_annot_tbl_combined$content)=="neg1","bac_nuc_ratio"])
  well_annot_tbl_combined$bac_nuc_ratio_neg_ratio=round(well_annot_tbl_combined$bac_nuc_ratio/median_neg,3)
  
  well_annot_tbl_combined <- merge(well_annot_tbl_combined,annoplate,by.x="src_well_names",by.y="Well",all.y=TRUE)
  well_annot_tbl_combined <- with(well_annot_tbl_combined, well_annot_tbl_combined[!(well_names == "" | is.na(well_names)), ])
  names(well_annot_tbl_combined)[names(well_annot_tbl_combined) == "Description"] <- "DWell"
  well_annot_tbl_combined 
}))

# ###########Robust Zscore calculation for Objects
# median_wells=median(data_table$MT.R01.Elongation.Factor.Objects,na.rm=TRUE)
# mad_wells=mad(data_table$MT.R01.Elongation.Factor.Objects,na.rm=TRUE)
# data_table$MT.R01.Elongation.Factor.Objects_zscore= (data_table$MT.R01.Elongation.Factor.Objects - median_wells)/mad_wells
###########Robust Zscore calculation for Ratio
#median_wells=median(data_table$bac_nuc_ratio,na.rm=TRUE)
#mad_wells=mad(data_table$bac_nuc_ratio,na.rm=TRUE)
#data_table$bac_nuc_ratio_zscore= round((data_table$bac_nuc_ratio - median_wells)/mad_wells,1)
###########Robust Zscore calculation for DAPI
# median_wells=median(data_table$Main.All.Mean.Intensity.Hoechst.Objects,na.rm=TRUE)
# mad_wells=mad(data_table$Main.All.Mean.Intensity.Hoechst.Objects,na.rm=TRUE)
# data_table$Main.All.Mean.Intensity.Hoechst.Objects_zscore= (data_table$Main.All.Mean.Intensity.Hoechst.Objects - median_wells)/mad_wells

#fname <- paste("Combined_Data",".xlsx",sep="")
#write.xlsx2(data_table,fname,sheetName="Data",row.names=TRUE, append=FALSE)

fname <- paste("Combined_Data",".txt",sep="")
write.table(data_table,file=fname,append=FALSE,quote=FALSE,sep="\t",eol="\n",na="NA",dec=".",row.names=TRUE,col.names=NA,fileEncoding="")

total_set=data_table[,c("src_well_names","DWell","Dapi_Objects","Dapi_objects_intensity_mean","baku_green_objects","baku_green_intensity_mean","bac_nuc_ratio","bac_nuc_ratio_neg_ratio","Plate_Id","content",grep("Plate.Name_",colnames(data_table),value=TRUE),grep("Plate.ID_",colnames(data_table),value=TRUE),grep("RefSeq.Accession.Number_",colnames(data_table),value=TRUE),grep("Gene.Symbol_",colnames(data_table),value=TRUE),grep("Full.Gene.Name_",colnames(data_table),value=TRUE),grep("Gene.ID_",colnames(data_table),value=TRUE),grep("siRNA.ID_",colnames(data_table),value=TRUE))]

fname <- paste("total_data",".xlsx",sep="")
write.xlsx2(total_set,fname,sheetName="Data",row.names=TRUE, append=FALSE)

fname <- paste("total_data",".txt",sep="")
write.table(total_set,file=fname,append=FALSE,quote=FALSE,sep="\t",eol="\n",na="NA",dec=".",row.names=TRUE,col.names=NA,fileEncoding="")

save(total_set,file ="total_set.rda")
save(data_table,file ="data_table.rda")