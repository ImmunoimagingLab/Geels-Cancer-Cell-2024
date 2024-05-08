###############
## Libraries ##
###############
library(dplyr)
library(patchwork)

###########################
## Read in scRNAseq data ##
###########################
df <- read.table("C:/Users/FM Lab/Desktop/PD1BioInformaticsData/SadeFeldman2018/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz", sep="\t", header=T,fill=TRUE)
#df <- read.table("C:/Users/rache/Downloads/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz", sep="\t", header=T,fill=TRUE)
df[1:5,1:5]

#find unique patient IDs
apply(df,1,function(x) unique(x))[1]

############################################################
## Adjusting the column names to fix Sade-Feldman dataset ##
############################################################
#get the column names of the original dataset and then remove the first name (i.e. "x")
actualcolnames <- colnames(df)[-1]
#remove the last column from the dataset
newdf <- df[-length(df[1,])]
#make the col names of the dataset the new names without the "x" as a col name
colnames(newdf) <- actualcolnames

#remove T_enriched and myeloid_enriched that is at end of some pre/post patient ID
newdf[1,] <- sub('_T_enriched', '', newdf[1,])
newdf[1,] <- sub('_myeloid_enriched', '', newdf[1,])

#add row one pre/post patient ID to column name
colnames(newdf) <- paste(newdf[1,], colnames(newdf), sep="_")
#remove first row with pre/post patient ID from dataset
finaldf <- newdf[-1,]


##################################################
## Removing patients not treated only with aPD1 ##
##################################################
#Read in clinical data
dfclinical <- read.csv("C:/Users/FM Lab/Desktop/PD1BioInformaticsData/SadeFeldman2018/GSE120575_patient_ID_single_cells.txt", sep="\t", header=T, skip=19)
#dfclinical <- read.table("C:/Users/rache/Downloads/GSE120575_patient_ID_single_cells.txt.gz", sep="\t", header=T, skip=19)
#Remove Protocols at bottom of document
dfclinical <- head(dfclinical, -38)
#Subset clinical data to have the therapy be anti-PD1
dfclinicalPD1 <- dfclinical[dfclinical$"characteristics..therapy" == "anti-PD1",]

#These patients were previously treated with aCTLA4 so we will remove them
dfclinicalPD1 <- dfclinicalPD1[- grep("Post_P1$", dfclinicalPD1$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.),]
dfclinicalPD1 <- dfclinicalPD1[- grep("Post_P1_2$", dfclinicalPD1$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.),]
dfclinicalPD1 <- dfclinicalPD1[- grep("Post_P6$", dfclinicalPD1$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.),]


###################################################################################
## Removing patients/cells from RNAseq data that were not treated only with aPD1 ##
###################################################################################
#add column of pre/post patient ID to title of cell in patient dataset
patientinfo <- dfclinicalPD1$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.
dfclinicalPD1$title <- paste(patientinfo, dfclinicalPD1$title, sep="_")

#Replace "-" with "." so the title nomenclature matches the finaldf
dfclinicalPD1$title <- gsub("-", ".", dfclinicalPD1$title)

#Filter cells of patients who were treated with only aPD1
filtered.df <- select(finaldf, dfclinicalPD1$title)

#Remove large datasets that we don't need
rm(actualcolnames)
rm(patientinfo)
rm(df)
rm(dfclinical)
rm(finaldf)
rm(newdf)

#Save filtered.df and then start a new file for Seurat
write.csv(filtered.df,"C:\\Users\\FM Lab\\Desktop\\PD1BioInformaticsData\\SadeFeldman2018\\filtered.SF.df.csv", row.names = TRUE)

write.csv(dfclinicalPD1,"C:\\Users\\FM Lab\\Desktop\\PD1BioInformaticsData\\SadeFeldman2018\\filtered.dfclinicalPD1.csv", row.names = TRUE)

