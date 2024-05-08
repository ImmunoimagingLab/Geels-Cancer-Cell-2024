# Libraries ---------------------------------------------------------------

knitr::opts_chunk$set(echo = TRUE)
library("Seurat")
library("dplyr")

#this is useful function to replicate the default color scheme for clustering in Seruat for other plots
#say you have 7 clusters - to get the color of cluster 3, you would put gg_color_hue(7)[3]
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Original Seurat Pipeline ------------------------------------------------

#load Sade-Feldman filtered dataset
df_aPD1scRNAseq <- read.table("C:\\Users\\FM Lab\\Desktop\\PD1BioInformaticsData\\SadeFeldman2018\\filtered.SF.df.csv", sep=",", header=T,fill=TRUE, row.names = 1)
#df_aPD1scRNAseq <- read.table("C:\\Users\\rache\\Documents\\filtered.SF.df.csv", sep=",", header=T,fill=TRUE, row.names = 1)

#setup Seurat object, filtering any cells with fewer than 200 unique genes and genes expressed in fewer than 3 cells
#currently this will be delim by "pre" and "post" -- not by which patient it is
data_seurat <- CreateSeuratObject(counts = df_aPD1scRNAseq, min.cells = 3, min.features = 200, names.field = 1, names.delim = "_", row.names = rownames(df_aPD1scRNAseq)) 

#Add patient ID to meta data
patientIDs <- sapply(strsplit(colnames(df_aPD1scRNAseq), "_"), function(x) paste(x[1:2], collapse = "_"))
data_seurat$patientID <- patientIDs

#add percent of genes coming from mitochondrial chromosome as an additional metadata column 
data_seurat[["percent.mito"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")

#look at the metadata to see what we've added
#you can use this at any point to remember the name of a column
head(data_seurat@meta.data)
tail(data_seurat@meta.data)

#what is the maxium percent mito -- less than 4 so Sade-Feldman already filtered this
max(data_seurat@meta.data$percent.mito)

#investigate data distributions 
VlnPlot(data_seurat, c("nFeature_RNA","nCount_RNA","percent.mito"), pt.size = 0.01)

#Don't need to do this filtering since it was already filtered by Sade-Feldman
#filter visual outlier cells - most will be doublets, empty, or dying
# data_seurat <- subset(data_seurat, subset = nFeature_RNA > 300 & nFeature_RNA < 4000)
# #investigate data distributions again to ensure more filtering isn't necessary
# VlnPlot(data_seurat, c("nFeature_RNA","nCount_RNA"), pt.size = 0.01)

#run general Seurat pipeline - default parameters work very well for most cases
data_seurat <- NormalizeData(data_seurat)
data_seurat <- FindVariableFeatures(data_seurat)
data_seurat <- ScaleData(object = data_seurat, verbose = TRUE)
data_seurat <- RunPCA(object = data_seurat, npcs = 30, verbose = TRUE)

#look for "elbow" to decide how many PCs to include
ElbowPlot(object = data_seurat)

#perform clustering and dimensionality reduction based on PCs chosen above 
#clustering resolution will make more clusters if the number is larger, and less clusters if it is smaller. 
#make sure to use the same dims argument for FindNeighbors and RunUMAP/RunTSNE or clusters will look "odd"
data_seurat <- FindNeighbors(object = data_seurat, dims = 1:13) # dims= PCs
data_seurat <- FindClusters(object = data_seurat, resolution = 0.7) #initially 0.1 then one may tweak to not higher than 1
data_seurat <- RunUMAP(object = data_seurat, dims = 1:13) #keep same dimensions (=PCs) as for FindNeighbors

#look at clustering - does the coloring match what you would expect visually?
DimPlot(data_seurat, pt.size=1)

#do your conditions seperate in space? 
DimPlot(data_seurat, group.by = "orig.ident", pt.size=1)
DimPlot(data_seurat, group.by = "orig.ident", pt.size=1)
DimPlot(data_seurat, group.by = "patientID", pt.size=1)

FeaturePlot(data_seurat,c("CD8B","CD4","FOXP3","ICOS","CD3E","CD14"))
FeaturePlot(data_seurat, c("FOXP3","ICOS", "CCR8", "CXCR6",  "TK1", "IL1R2", "NRBF2", "DCUN1D5", "ZBTB32", "C15orf57", "EBI3"))
FeaturePlot(data_seurat, c("FOXP3","ICOS", "CCR8", "CXCR6",  "TK1", "IL1R2"))
FeaturePlot(data_seurat, c("FOXP3","ICOS", "NRBF2", "DCUN1D5", "ZBTB32",  "EBI3"))
FeaturePlot(data_seurat, c("FOXP3","ICOS", "MTMR6", "IL21R", "CTSC",  "NPTN"))
FeaturePlot(data_seurat, c("FOXP3","ICOS", "EZH2", "SDC4","NT5C3A","NFE2L2","ERH","TNIP1"))
FeaturePlot(data_seurat, c("FOXP3","ICOS", "STMN1", "FOSL2", "KSR1", "CLTC","HIVEP3","NF2"))
FeaturePlot(data_seurat, c("FOXP3","ICOS", "IDH3B", "PELI1", "C14orf2", "ATP1B1", "ZC3H15" ))

FeaturePlot(data_seurat, c("FOXP3","ICOS", "SYNE2", "DGKA", "ALOX5AP", "CCNL2", "FBXW4"))
FeaturePlot(data_seurat, c("FOXP3","ICOS", "HACL1", "SELM", "GCC2", "FKBP5","PLA2G4B"))
FeaturePlot(data_seurat, c("FOXP3","ICOS", "FCER1A", "CST3" , "MS4A1"))
FeaturePlot(data_seurat, c("MS4A1", 	"GNLY", "NKG7", "PPBP", "FCGR3A", "MS4A7"))


#Rename clusters 
new.cluster.ids <- c("CD8 T", "CD4 T", "B", "NK", "Monocytes", "5")
names(new.cluster.ids) <- levels(data_seurat)
data_seurat <- RenameIdents(data_seurat, new.cluster.ids)
DimPlot(data_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


#what are your clusters? find the genes that define them.
cluster.markers.res.0.1 <- FindAllMarkers(data_seurat, only.pos = T)  # if one wants also downregulated genes, eliminate only.pos=T

#let's look at the top 10 genes (by log fold change) for each cluster
View(cluster.markers.res.0.1 %>% group_by(cluster) %>% top_n(10, avg_log2FC))


# Subset to only FOXP3 cells ----------------------------------------------

#now, let's subset to only Tregs (based on expression level) and reanalyze the data
# data_CD4T_cells_only <- subset(data_seurat, ident = c("1","4", "7"))
# data_CD4T_cells_only <- subset(data_seurat, ident = c("1"))
data_Tregs_cells_only <- subset(data_seurat, FOXP3 > 0)

#we can now rerun the seurat pipeline on only these cells
data_Tregs_cells_only <- NormalizeData(data_Tregs_cells_only)
data_Tregs_cells_only <- FindVariableFeatures(data_Tregs_cells_only)
data_Tregs_cells_only <- ScaleData(object = data_Tregs_cells_only, verbose = TRUE)
data_Tregs_cells_only <- RunPCA(object = data_Tregs_cells_only, npcs = 30, verbose = TRUE)

#look for "elbow" to decide how many PCs to include
ElbowPlot(object = data_Tregs_cells_only)

#perform clustering and dimensionality reduction based on PCs chosen above 
#clustering resolution will make more clusters if the number is larger, and less clusters if it is smaller. 
#make sure to use the same dims argument for FindNeighbors and RunUMAP/RunTSNE or clusters will look "odd"
data_Tregs_cells_only <- FindNeighbors(object = data_Tregs_cells_only, dims = 1:12)
data_Tregs_cells_only <- FindClusters(object = data_Tregs_cells_only, resolution = 0.3)
data_Tregs_cells_only <- RunUMAP(object = data_Tregs_cells_only, dims = 1:12)

#again, let's look at our clusters and condition effects
DimPlot(data_Tregs_cells_only, pt.size=1)
DimPlot(data_Tregs_cells_only, group.by = "orig.ident", pt.size=1)
DimPlot(data_Tregs_cells_only, split.by = "orig.ident", pt.size=1)
#FeaturePlot(data_Tregs_cells_only, c("ICOS","CD4","CD8A","FOXP3", "CD3E", "PTPRC"))
FeaturePlot(data_Tregs_cells_only, c("ICOS","CD4","CD8A","CD8B","FOXP3", "CD3E", "PTPRC"))

# Subsetting again to get FOXP3 cells only NO CD8+ dudes ------------------
data_Tregs_cells_only <- subset(data_Tregs_cells_only, ident = c("0","1","2"))

#we can now rerun the seurat pipeline on only these cells
data_Tregs_cells_only <- NormalizeData(data_Tregs_cells_only)
data_Tregs_cells_only <- FindVariableFeatures(data_Tregs_cells_only)
data_Tregs_cells_only <- ScaleData(object = data_Tregs_cells_only, verbose = TRUE)
data_Tregs_cells_only <- RunPCA(object = data_Tregs_cells_only, npcs = 30, verbose = TRUE)

#look for "elbow" to decide how many PCs to include
ElbowPlot(object = data_Tregs_cells_only)

#perform clustering and dimensionality reduction based on PCs chosen above 
#clustering resolution will make more clusters if the number is larger, and less clusters if it is smaller. 
#make sure to use the same dims argument for FindNeighbors and RunUMAP/RunTSNE or clusters will look "odd"
data_Tregs_cells_only <- FindNeighbors(object = data_Tregs_cells_only, dims = 1:13)
data_Tregs_cells_only <- FindClusters(object = data_Tregs_cells_only, resolution = 0.3)
data_Tregs_cells_only <- RunUMAP(object = data_Tregs_cells_only, dims = 1:13)

#again, let's look at our clusters and condition effects
DimPlot(data_Tregs_cells_only, pt.size=1)
DimPlot(data_Tregs_cells_only, group.by = "orig.ident", pt.size=1)
DimPlot(data_Tregs_cells_only, group.by = "patientID", pt.size=1)
DimPlot(data_Tregs_cells_only, split.by = "orig.ident", pt.size=1)
FeaturePlot(data_Tregs_cells_only, c("ICOS","CD8A","FOXP3"))
FeaturePlot(data_Tregs_cells_only, c("ICOS","CD4","CD8A","FOXP3", "CD3E", "PTPRC","PER1","IFI44L"))
FeaturePlot(data_Tregs_cells_only, c("PER1","IFI44L"))

FeaturePlot(data_Tregs_cells_only, c("ICOS","CCR8","FOXP3","CXCR6"))


# Differential gene expression for FOXP3+ ICOS+ and FOXP3+ ICOS-  --------

#Get ICOS positive cells from the FOXP3 positive cells
data_Tregs_cells_only_ICOSpos <- subset(data_Tregs_cells_only, ICOS > 0)
#Add meta data column for ICOS positive or negative 
data_Tregs_cells_only@meta.data$ICOScol <-
  ifelse(
    rownames(data_Tregs_cells_only@meta.data) %in% colnames(data_Tregs_cells_only_ICOSpos),
    "ICOSpos",
    "ICOSneg"
  )

#we can also do a supervised marker analysis between conditions - a simple option is ICOS positive vs ICOS negative across all of the Tregs
Idents(data_Tregs_cells_only) <- "ICOScol"  # gives ICOS positive vs ICOS negative identities back on the Treg filtered matrix
condition.markers.res.0.3 <- FindAllMarkers(data_Tregs_cells_only, only.pos = T)  # finds positively regulated genes across all Tregs cells, ICOS negative or positive.  One may want to do this clusted by cluster to find differntially expressed genes.


#let's look at the top 10 genes (by log fold change) for each cluster
# In the output, pct.1 is the percentage of cells in the cluster where the gene is detected, while pct.2 is the percentage of cells on average in all the other clusters where the gene is detected. 
View(condition.markers.res.0.3 %>% group_by(cluster) %>% top_n(100, avg_log2FC))

# writes down the table on Desktop (delimiters are spacers)
write.table(condition.markers.res.0.3 %>% group_by(cluster) %>% top_n(10, avg_log2FC), "C:/Users/FM Lab/Rachel/Lowengrub-Marangoni-Research/PD1_Bioinformatics/FOXP3ICOSposnegDEgenesUPDATED.csv", sep = ",")

# writes down the table on Desktop (delimiters are spacers)
write.table(condition.markers.res.0.3[condition.markers.res.0.3$cluster=="ICOSpos",], "C:/Users/FM Lab/Rachel/Lowengrub-Marangoni-Research/PD1_Bioinformatics/FOXP3ICOSposDEG.csv", sep = ",")

#Feature plot to look at ICOS positive and negative differentially expressed genes
FeaturePlot(data_Tregs_cells_only, c("ICOS","CCR8","CXCR6","SYNE2","DGKA"))

# Differential gene expression for FOXP3+ CCR8+ and FOXP3+ CCR8-  --------

#Get CCR8 positive cells from the FOXP3 positive cells
data_Tregs_cells_only_CCR8pos <- subset(data_Tregs_cells_only, CCR8 > 0)
#Add meta data column for CCR8 positive or negative 
data_Tregs_cells_only@meta.data$CCR8col <-
  ifelse(
    rownames(data_Tregs_cells_only@meta.data) %in% colnames(data_Tregs_cells_only_CCR8pos),
    "CCR8pos",
    "CCR8neg"
  )

#we can also do a supervised marker analysis between conditions - a simple option is CCR8 positive vs CCR8 negative across all of the Tregs
Idents(data_Tregs_cells_only) <- "CCR8col"  # gives CCR8 positive vs CCR8 negative identities back on the Treg filtered matrix
condition.markers.res.0.3 <- FindAllMarkers(data_Tregs_cells_only, only.pos = T)  # finds positively regulated genes across all Tregs cells, CCR8 negative or positive.  One may want to do this clusted by cluster to find differntially expressed genes.


#let's look at the top 10 genes (by log fold change) for each cluster
# In the output, pct.1 is the percentage of cells in the cluster where the gene is detected, while pct.2 is the percentage of cells on average in all the other clusters where the gene is detected. 
View(condition.markers.res.0.3 %>% group_by(cluster) %>% top_n(20, avg_log2FC))

# writes down the table on Desktop (delimiters are spacers)
write.table(condition.markers.res.0.3 %>% group_by(cluster) %>% top_n(10, avg_log2FC), "C:/Users/FM Lab/Rachel/Lowengrub-Marangoni-Research/PD1_Bioinformatics/FOXP3CCR8posnegDEgenesUPDATED.csv", sep = ",")




# Differential gene expression for PRE vs POST ----------------------------

data_Tregs_cells_only_PrevPost <- subset(data_seurat, FOXP3 > 0)

#we can now rerun the seurat pipeline on only these cells
data_Tregs_cells_only_PrevPost <- NormalizeData(data_Tregs_cells_only_PrevPost)
data_Tregs_cells_only_PrevPost <- FindVariableFeatures(data_Tregs_cells_only_PrevPost)
data_Tregs_cells_only_PrevPost <- ScaleData(object = data_Tregs_cells_only_PrevPost, verbose = TRUE)
data_Tregs_cells_only_PrevPost <- RunPCA(object = data_Tregs_cells_only_PrevPost, npcs = 30, verbose = TRUE)

#look for "elbow" to decide how many PCs to include
ElbowPlot(object = data_Tregs_cells_only_PrevPost)

#perform clustering and dimensionality reduction based on PCs chosen above 
#clustering resolution will make more clusters if the number is larger, and less clusters if it is smaller. 
#make sure to use the same dims argument for FindNeighbors and RunUMAP/RunTSNE or clusters will look "odd"
data_Tregs_cells_only_PrevPost <- FindNeighbors(object = data_Tregs_cells_only_PrevPost, dims = 1:13)
data_Tregs_cells_only_PrevPost <- FindClusters(object = data_Tregs_cells_only_PrevPost, resolution = 0.3)
data_Tregs_cells_only_PrevPost <- RunUMAP(object = data_Tregs_cells_only_PrevPost, dims = 1:13)

#Subsetting again to get FOXP3 cells only NO CD8+ dudes 
data_Tregs_cells_only_PrevPost <- subset(data_Tregs_cells_only_PrevPost, ident = c("0","2"))

#we can now rerun the seurat pipeline on only these cells
data_Tregs_cells_only_PrevPost <- NormalizeData(data_Tregs_cells_only_PrevPost)
data_Tregs_cells_only_PrevPost <- FindVariableFeatures(data_Tregs_cells_only_PrevPost)
data_Tregs_cells_only_PrevPost <- ScaleData(object = data_Tregs_cells_only_PrevPost, verbose = TRUE)
data_Tregs_cells_only_PrevPost <- RunPCA(object = data_Tregs_cells_only_PrevPost, npcs = 30, verbose = TRUE)

#look for "elbow" to decide how many PCs to include
ElbowPlot(object = data_Tregs_cells_only_PrevPost)

#perform clustering and dimensionality reduction based on PCs chosen above 
#clustering resolution will make more clusters if the number is larger, and less clusters if it is smaller. 
#make sure to use the same dims argument for FindNeighbors and RunUMAP/RunTSNE or clusters will look "odd"
data_Tregs_cells_only_PrevPost <- FindNeighbors(object = data_Tregs_cells_only_PrevPost, dims = 1:13)
data_Tregs_cells_only_PrevPost <- FindClusters(object = data_Tregs_cells_only_PrevPost, resolution = 0.3)
data_Tregs_cells_only_PrevPost <- RunUMAP(object = data_Tregs_cells_only_PrevPost, dims = 1:13)

#we can also do a supervised marker analysis between conditions - a simple option is pre vs post across all of the Tregs
Idents(data_Tregs_cells_only_PrevPost) <- "orig.ident"  # gives pre vs post identities back on the Treg filtered matrix
condition.markers.res.0.3.PREPOST <- FindAllMarkers(data_Tregs_cells_only_PrevPost, only.pos = T)  # finds positively regulated genes across all Tregs cells, ICOS negative or positive.  One may want to do this clusted by cluster to find differntially expressed genes.

View(condition.markers.res.0.3.PREPOST %>% group_by(cluster) %>% top_n(10, avg_log2FC))

# writes down the table on Desktop (delimiters are spacers)
write.table(condition.markers.res.0.3.PREPOST %>% group_by(cluster) %>% top_n(10, avg_log2FC), "C:/Users/FM Lab/Rachel/Lowengrub-Marangoni-Research/PD1_Bioinformatics/PREvPOST_TREG_DEgenes.csv", sep = ",")

#Feature plot to look at ICOS positive and negative differentially expressed genes
FeaturePlot(data_Tregs_cells_only_PrevPost, c("ICOS", "IFI44L", "MIR29B1"))
DimPlot(data_Tregs_cells_only_PrevPost, pt.size=1)
DimPlot(data_Tregs_cells_only_PrevPost, group.by = "orig.ident", pt.size=1)
DimPlot(data_Tregs_cells_only_PrevPost, split.by = "orig.ident", pt.size=1)


# PRE Treg Subset ---------------------------------------------------------

PREdata_Tregs_cells_only <- subset(data_Tregs_cells_only, orig.ident=="Pre")

#we can now rerun the seurat pipeline on only these cells
PREdata_Tregs_cells_only <- NormalizeData(PREdata_Tregs_cells_only)
PREdata_Tregs_cells_only <- FindVariableFeatures(PREdata_Tregs_cells_only)
PREdata_Tregs_cells_only <- ScaleData(object = PREdata_Tregs_cells_only, verbose = TRUE)
PREdata_Tregs_cells_only <- RunPCA(object = PREdata_Tregs_cells_only, npcs = 30, verbose = TRUE)

#look for "elbow" to decide how many PCs to include
ElbowPlot(object = PREdata_Tregs_cells_only)

#perform clustering and dimensionality reduction based on PCs chosen above 
#clustering resolution will make more clusters if the number is larger, and less clusters if it is smaller. 
#make sure to use the same dims argument for FindNeighbors and RunUMAP/RunTSNE or clusters will look "odd"
PREdata_Tregs_cells_only <- FindNeighbors(object = PREdata_Tregs_cells_only, dims = 1:13)
PREdata_Tregs_cells_only <- FindClusters(object = PREdata_Tregs_cells_only, resolution = 0.5)
PREdata_Tregs_cells_only <- RunUMAP(object = PREdata_Tregs_cells_only, dims = 1:13)

#again, let's look at our clusters and condition effects
DimPlot(PREdata_Tregs_cells_only,pt.size=1)
DimPlot(PREdata_Tregs_cells_only, group.by = "orig.ident",pt.size=1)
DimPlot(PREdata_Tregs_cells_only, split.by = "orig.ident",pt.size=1)
FeaturePlot(PREdata_Tregs_cells_only, c("ICOS","CD4","CD8A","FOXP3", "CTLA4","PDCD1"))
FeaturePlot(PREdata_Tregs_cells_only, c("PER1","AC004453.8"))

#what are your clusters? find the genes that define them.
PREcluster.markers.res.0.1 <- FindAllMarkers(PREdata_Tregs_cells_only, only.pos = T)  # if one wants also downregulated genes, eliminate only.pos=T

#let's look at the top 10 genes (by log fold change) for each cluster
View(PREcluster.markers.res.0.1 %>% group_by(cluster) %>% top_n(10, avg_log2FC))


# POST Treg Subset --------------------------------------------------------

POSTdata_Tregs_cells_only <- subset(data_Tregs_cells_only, orig.ident=="Post")

#we can now rerun the seurat pipeline on only these cells
POSTdata_Tregs_cells_only <- NormalizeData(POSTdata_Tregs_cells_only)
POSTdata_Tregs_cells_only <- FindVariableFeatures(POSTdata_Tregs_cells_only)
POSTdata_Tregs_cells_only <- ScaleData(object = POSTdata_Tregs_cells_only, verbose = TRUE)
POSTdata_Tregs_cells_only <- RunPCA(object = POSTdata_Tregs_cells_only, npcs = 30, verbose = TRUE)

#look for "elbow" to decide how many PCs to include
ElbowPlot(object = POSTdata_Tregs_cells_only)

#perform clustering and dimensionality reduction based on PCs chosen above 
#clustering resolution will make more clusters if the number is larger, and less clusters if it is smaller. 
#make sure to use the same dims argument for FindNeighbors and RunUMAP/RunTSNE or clusters will look "odd"
POSTdata_Tregs_cells_only <- FindNeighbors(object = POSTdata_Tregs_cells_only, dims = 1:14)
POSTdata_Tregs_cells_only <- FindClusters(object = POSTdata_Tregs_cells_only, resolution = 0.3)
POSTdata_Tregs_cells_only <- RunUMAP(object = POSTdata_Tregs_cells_only, dims = 1:14)

#again, let's look at our clusters and condition effects
DimPlot(POSTdata_Tregs_cells_only,pt.size=1)
DimPlot(POSTdata_Tregs_cells_only, group.by = "orig.ident",pt.size=1)
DimPlot(POSTdata_Tregs_cells_only, split.by = "orig.ident",pt.size=1)
FeaturePlot(POSTdata_Tregs_cells_only, c("ICOS","CD4","CD8A","FOXP3", "CD3E", "PTPRC"))
FeaturePlot(POSTdata_Tregs_cells_only, c("PER1","AC004453.8"))
