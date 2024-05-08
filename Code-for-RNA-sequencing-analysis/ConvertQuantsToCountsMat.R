#####################################################
## Code to convert quant.genes.sf to RData matrix ##
## Code to get gene symbols in tpm matrix ##
#####################################################

library(stringr)
library(dplyr)
library(data.table)

options(stringsAsFactors=FALSE)

samples=dir()
samples ## make sure to take all the folders and not log files etc

thisMat=read.table(paste0(samples[1],"/quant.genes.sf"),sep='\t',head=TRUE)
tpm=thisMat[,c("Name","TPM")]

for(i in 2:length(samples)){
	thisMat=read.table(paste0(samples[i],"/quant.genes.sf"),sep='\t',head=TRUE)
	thisMat=thisMat[na.omit(match(tpm$Name,thisMat$Name)),]
	tpm=as.data.frame(cbind(tpm,thisMat[,"TPM"]))
}

#make column names p#
colnames(tpm)[2:ncol(tpm)]=str_extract(samples, "[^_]+")

# load gene symbol data
gencode <- read.delim('~/MarangoniLab/PD1RNAseq/EGAData/gencode.v32_gene_annotation_table.txt', header=TRUE, sep='\t', dec='.')
gencode <- gencode[,1:2]

#match column 1 of tpm with column 1 of gencode
ind <- match(tpm$Name, gencode$Geneid)
table(is.na(ind)) #how many matches did we get?

#remove indices that did not match with tpm
gencode.subset <- gencode[na.omit(ind),]

#merge the data sets together so we have the gene symbols
data.Merged <- merge(tpm, gencode.subset, by.x='Name', by.y='Geneid')

#Move the gene symbol column to second column
data.Merged.rearrange <- data.Merged %>% select(Name, GeneSymbol, everything())

tpm <- data.Merged.rearrange

save(tpm,file="~/nb227/practice/output/GeneExpressionAnalysis/TPM_Matrix_Raw.RData")

#Get row that has FOXP3 as gene symbol
#tpm[tpm$GeneSymbol %like% "FOXP3", ]
