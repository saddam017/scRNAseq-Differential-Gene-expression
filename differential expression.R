#seting working directory
setwd("C:/Users/sadda/Downloads/gene expression paper/RNAseq single cell wrok reproduce/")
getwd()
#Lad Library Seurat packages
library(Seurat)

#install SeuratData packages
#install.packages("SeuratData")
#library(SeuratData)

# Load the PBMC dataset
projDir <- "C:/Users/sadda/Downloads/gene expression paper/RNAseq single cell wrok reproduce"
srt.data <- Read10X(data.dir = "./single cell/")
srt <- CreateSeuratObject(counts = srt.data, project = "NSCLC")
srt


colnames(srt)
head(colnames(srt))
tail(colnames(srt))
table(srt$orig.ident)

srt@meta.data



##########A. Data Processing: Standard pre-processing workflow###
##The step below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. These represent the selection and filtration of cells based on QC metrics, data normalization and scaling, and the detection of highly variable features.

# 1.QC and selecting cells for further analysis

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^MT-") #calculate parcentage of MT gene and store the resuult into new column (percent.mt)
# Visualize QC metrics as a violin plot
VlnPlot(srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

############The number of unique genes and total molecules are automatically calculated during CreateSeuratObject###
# You can find them stored in the object meta data, Show QC metrics for the first 5 cells
head(srt@meta.data, 5)
tail(srt@meta.data, 5)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(srt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(srt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#here we subset the bad cell by using subset function and value 
srt <- subset(srt, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
srt

plot3 <- FeatureScatter(srt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(srt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 + plot3 + plot4


#2. Normalizing the data:  After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method "LogNormalize" that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in pbmc[["RNA"]]@data.

srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
srt <- NormalizeData(srt)


###let us check the normalized log expression values###
head(srt[["RNA"]]@counts)
head(srt[["RNA"]]@data)

#3. Identification of highly variable features (feature selection)
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(srt), 15)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(srt)
plot2 <- LabelPoints(plot = plot1, points = top10, )
plot1 + plot2

############DATA-SCALING#####################
#Why Data-Scaling ?? Data scaling  a linear transformation ('scaling') that is a standard pre-processing step prior to dimensional reduction techniques like PCA. 
#The ScaleData function: Shifts the expression of each gene, so that the mean expression across cells is 0, Scales the expression of each gene, so that the variance across cells is 1. 
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate. The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(srt)
head(all.genes)
srt <- ScaleData(srt, features = all.genes)

#4. perform linear dimensional reduction

srt <- RunPCA(srt, features = VariableFeatures(object = srt))

#### Examine and visualize PCA results a few diifferent ways
#PCA is to represent a multivariate data table as smaller set of variables (summary indices) in order to observe trends, jumps, clusters and outliers.

print(srt[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(srt,dims = 1:2, reduction = "pca")
DimPlot(srt, reduction = "pca")
DimHeatmap(srt, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(srt, dims = 1:15, cells = 500, balanced = TRUE)

#5. Cell clustering

####5a- Manual way: one resoluation at a time
srt <- FindNeighbors(srt, dims = 1:10)
srt <- FindClusters(srt, resolution = 0.5)

#Look at cluster ID's of the first 5 cells
head(Idents(srt),5)

#5b in loop way:

res = c(0.2,0.8,1,2,5) #Deovo cluster resulation

for (i in seq_along(res))
{
  print(res[i])
  pbmc.combined = FindClusters(srt, resolution = res[i], verbase = T, n.start = 100)
}


#6. Run non-linear dimensional reduction (UMAP/tSEN)

# if you haven't installed UMAP, you can do so reticulate::py_install(packages = 'umap-learn')

srt <- RunUMAP(srt, dims = 1:10)

#individual cluster
DimPlot(srt, reduction = "umap")
srt@reductions$pca
srt@reductions$umap
saveRDS(srt, file = "./rds_file.rds")



##### Find Cell types Differentially Expressed markers ###

# step 1: find cell-type markers

#seuurat uses Limma package i Bioconductor to use Wilcoxon Rank Sum Test, for efficient identification of cell types marker
#install.packages('BiocManager')
#BiocManager::install('limma')

# find all markers of cluster 1
cluster1.markers <- FindMarkers(srt, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

##find markers for every cluster compared to all remaining cells, report only the positive ones
srt.markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(srt.markers,paste0(projDir,'NSCLC_deg_markers.csv'))
write.csv(srt.markers, "NSCLC_deg_markers.csv")


## Top 10 genes#

library(dplyr)
srt.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10


DoHeatmap(srt, features = top10$gene)


## Top 5 genes 
srt.markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt = avg_log2FC) -> top5

# write.csv (marker. files)
Idents(srt) ###ATTTT_1_3K: 1K : T
"pbmc4k"

srt$seurat_clusters
Idents(srt) <- "seurat_clusters"







library(DESeq2)

##find markers for every cluster compared to all remaining cells, report only the positive ones
head(srt.markers, n=5)

resSig <- subset(srt.markers, avg_log2FC >2 & p_val_adj <0.01 | avg_log2FC < -2 & p_val_adj < 0.01)
write.csv(resSig , "differential_gene.csv")





upregulated_gene <- subset(srt.markers,avg_log2FC >2 & p_val_adj <0.01)
write.csv(upregulated_gene, "upregulated_gene.csv")

down_regulated_gene <- subset(srt.markers,avg_log2FC < -2 & p_val_adj <0.01)
write.csv(down_regulated_gene, "down_regulated_gene.csv")



library(EnhancedVolcano)
EnhancedVolcano(resSig,
                lab = "",
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = NULL,
                title = NULL,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                xlim = c(-8,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylim = c(0,12),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #transcriptPointSize = 0.5,
                #transcriptLabSize = 4.0,
                colAlpha = 1,
                shape = 19,
                subtitle = NULL,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                colConnectors = 'grey50',
                border = 'full' )






