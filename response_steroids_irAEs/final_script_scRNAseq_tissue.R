library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
reticulate::py_install(packages = 'umap-learn')
library(ggplot2)
library(readxl)
library(scales)
library(sctransform)
BiocManager::install("glmGamPoi")
library(RColorBrewer)
library(clustree)
library(future)
library(Matrix)
library(tidyr)
library(reshape2)
library(rstatix)
library(ggsignif)
library(GSVA)
library(rstatix)
library(DESeq2)
library(ggrepel)

set.seed(1234)
options(Seurat.object.assay.version = "v3")

########################################
####### LOADING Thomas 2024 data #######
########################################

setwd("~")
setwd(".../scRNAseq_reanalysis/Thomas2024/")

#load data (from .h5ad file)
sce <- zellkonverter::readH5AD("GSE206299_ircolitis-tissue-cd4.h5ad")

#Create count matrix and create Seurat object
matrix <- as.data.frame(sce@assays@data$X)
new_rownames <- sapply(strsplit(rownames(matrix), "\\|"), function(x) x[2])
rownames(matrix) <- make.unique(new_rownames)
colnames(matrix) <- gsub("\\|", ".", colnames(matrix))
data <- CreateSeuratObject(counts=as.matrix(matrix))

metadata <- data.frame(Patient = sapply(strsplit(colnames(data), "\\."), function(x) x[1]))
rownames(metadata) <- colnames(data)
data <- AddMetaData(data,metadata)

#Only select patient samples of Colitis treated with steroids for which Response was reported.
patients <- c("SIC_100","SIC_121","SIC_126","SIC_134","SIC_140","SIC_141","SIC_32","SIC_36",
              "SIC_40","SIC_43","SIC_71","SIC_89","SIC_97")

matching_cells <- sapply(data@meta.data$Patient, function(patient) any(sapply(patients, function(p) grepl(p, patient))))
data <- subset(data, cells = rownames(metadata)[matching_cells])
unique(data@meta.data$Patient)

#Remove cells from uninflamed loci, or from non-CD45pos samples
remove <- c("SIC_36_Colon_140_3p_GEX","SIC_43_LeftColon-uninflamed_157_CD45pos_3p_GEX",
            "SIC_32_Colon_128_CAP_unselect_5p_GEX","SIC_40_Colon_155_3p_GEX",
            "SIC_43_LeftColon-uninflamed_157_3p_GEX")
`%notin%` <- Negate(`%in%`)
data <- subset(data,subset=Patient %notin% remove)
unique(data@meta.data$Patient)

#Clean up data
data <- PercentageFeatureSet(data, pattern = "^MT-",col.name="percent.mt")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

#Remove TCR genes
data <- data[!grepl("^TR[ABDG][VJC]",rownames(data))]

#Pull apart different samples for SCTransform, then merge again to correct for patient batch effect
seurats <- list()
for(i in unique(data@meta.data$Patient)){
  seurats[[i]] <- subset(data,subset=Patient==i)
}

for(i in seq_along(seurats)){
  seurats[[i]] <- SCTransform(seurats[[i]], method = "glmGamPoi", vars.to.regress = c("percent.mt"), verbose = TRUE)
}

seurat.merged <- merge(seurats[[1]], y = seurats[2:length(seurats)], 
                       add.cell.ids = c("SIC_100","SIC_121","SIC_126","SIC_134","SIC_140","SIC_141_A",
                                        "SIC_141_B","SIC_36_A","SIC_36_B","SIC_40","SIC_43","SIC_71",
                                        "SIC_89","SIC_97"), 
                       project = "response_pred")

Thomas_data <- seurat.merged

#Clear up memory
rm(data,matrix,metadata,sce,seurat.merged,seurats)
gc()

########################################
####### LOADING Gupta 2024 data #######
########################################
setwd("~")
setwd(".../scRNAseq_reanalysis/Gupta2024/")
data <- readRDS("CD3_scRNAseq_data.RDS")

#Select only the 11 donors treated with steroids and for whom steroid response was reported
donors <- c("SA","SC","GI3692","GI4685","GI4740","GI6274","GI6275","GI6852","PRISE14","PRISE50",
            "GI2935","GI3152","GI4059","GI4069","GI4158","GI4738","PRISE5")
data <- subset(data,subset=Donor %in% donors)

seurats <- list()
for(i in unique(data@meta.data$Donor)){
  seurats[[i]] <- subset(data,subset=Donor==i)
}

for(i in seq_along(seurats)){
  DefaultAssay(seurats[[i]]) <- "RNA"
  seurats[[i]] <- SCTransform(seurats[[i]], method = "glmGamPoi", verbose = TRUE)
}

for(i in seq_along(seurats)){
  print(seurats[[i]]$Donor[1])
}

seurat.merged <- merge(seurats[[1]], y = seurats[2:length(seurats)], 
                       add.cell.ids = c("PRISE50","GI6852","GI3692","GI4740","GI4685","SC","SA","GI4059","GI3152",
                                        "GI4158","PRISE5","GI4069","GI2935","GI4738","GI6275","GI6274","PRISE14"), 
                       project = "response_pred")

avg_exp <- AggregateExpression(seurat.merged)
avg_exp <- as.data.frame(avg_exp$RNA)
CD8A_expr <- avg_exp[rownames(avg_exp)=="CD8A",]/(rowSums(avg_exp[rownames(avg_exp)=="CD8A",])/ncol(avg_exp)) # Divide by mean as normalization measure
CD4_expr <- avg_exp[rownames(avg_exp)=="CD4",]/(rowSums(avg_exp[rownames(avg_exp)=="CD4",])/ncol(avg_exp)) # Divide by mean as normalization measure
CD4_expr-CD8A_expr

CD4_cells <- c("Th17 PD1+","Naive CD4+","Tph","Th1","Tfh","Tregs","Th17 PD1-","IFN Response")
subset_CD4 <- subset(seurat.merged,idents=CD4_cells)

Gupta_data <- subset_CD4

#Clear up memory
rm(avg_exp,CD4_expr,CD8A_expr,data,subset_CD4,seurat.merged,seurats)
gc()

########################################
####### LOADING Luoma 2020 data #######
########################################

options(Seurat.object.assay.version = "v3")
setwd("~")
setwd(".../scRNAseq_reanalysis/Luoma2020/")
path1 <- "C1-CD3-IFX_no/"
path2 <- "C2-CD3-IFX_no/" 
path3 <- "C3-CD3-IFX_yes/"
path4 <- "C4-CD3-IFX_yes/"
path5 <- "C5-CD3-IFX_yes/"
path6 <- "C6-CD3-IFX_no/"
path7 <- "C7-CD3-IFX_yes/"
path8 <- "C8-CD3-IFX_no/" 

#Load scRNAseq CD3+ T cells dataset
for(i in 1:8){
  path <- get(paste0("path",i))
  sc.data <- ReadMtx(mtx=paste0(path,"matrix.mtx"),
                     features=paste0(path,"features.tsv"),
                     cells=paste0(path,"barcodes.tsv"))
  seurat <- paste0("seurat",i)
  assign(seurat,CreateSeuratObject(counts=sc.data))
}

# store mitochondrial percentage in object meta data
#ELISE: seurat1 <- SetIdent(seurat1, value="orig.ident") --> hier logischere var namen voor maken die patient beschrijven.
seurat1 <- PercentageFeatureSet(seurat1, pattern = "^MT-",col.name="percent.mt")
seurat2 <- PercentageFeatureSet(seurat2, pattern = "^MT-",col.name="percent.mt")
seurat3 <- PercentageFeatureSet(seurat3, pattern = "^MT-",col.name="percent.mt")
seurat4 <- PercentageFeatureSet(seurat4, pattern = "^MT-",col.name="percent.mt")
seurat5 <- PercentageFeatureSet(seurat5, pattern = "^MT-",col.name="percent.mt")
seurat6 <- PercentageFeatureSet(seurat6, pattern = "^MT-",col.name="percent.mt")
seurat7 <- PercentageFeatureSet(seurat7, pattern = "^MT-",col.name="percent.mt")
seurat8 <- PercentageFeatureSet(seurat8, pattern = "^MT-",col.name="percent.mt")

# Visualize QC metrics as a violin plot
for(i in 1:8){
  plot(VlnPlot(get(paste0("seurat",i)), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),main=paste0("Plots for seurat ",i))
}

#QC pass (based on Violin plots; values set slightly less conservative than standard nFeature_RNA<2500 and percent.mt<5
seurat1 <- subset(seurat1, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
seurat2 <- subset(seurat2, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
seurat3 <- subset(seurat3, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
seurat4 <- subset(seurat4, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
seurat5 <- subset(seurat5, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
seurat6 <- subset(seurat6, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
seurat7 <- subset(seurat7, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
seurat8 <- subset(seurat8, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

#Remove TCR genes (script/method based on Sundell et al. Briefings in Functional Genomics 22(3):2023;263-73.)
seurat1 <- seurat1[!grepl("^TR[ABDG][VJC]",rownames(seurat1))]
seurat2 <- seurat2[!grepl("^TR[ABDG][VJC]",rownames(seurat2))]
seurat3 <- seurat3[!grepl("^TR[ABDG][VJC]",rownames(seurat3))]
seurat4 <- seurat4[!grepl("^TR[ABDG][VJC]",rownames(seurat4))]
seurat5 <- seurat5[!grepl("^TR[ABDG][VJC]",rownames(seurat5))]
seurat6 <- seurat6[!grepl("^TR[ABDG][VJC]",rownames(seurat6))]
seurat7 <- seurat7[!grepl("^TR[ABDG][VJC]",rownames(seurat7))]
seurat8 <- seurat8[!grepl("^TR[ABDG][VJC]",rownames(seurat8))]

#Perform SCTransform on each sample individually (this will correct for batch effects; only possible if cells are approximately the same)
# The latest version of sctransform also supports using glmGamPoi package which substantially improves the speed of the learning procedure.
# Single command SCTransform replaces commands NormalizeData(), ScaleData() and FindVariableFeatures()

#First increase globals.maxSize to ~2.7 GiB
options(future.globals.maxSize = 3000 * 1024^2)

#Then perform SCTransform()
seurat1 <- SCTransform(seurat1, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
seurat2 <- SCTransform(seurat2, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
seurat3 <- SCTransform(seurat3, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
seurat4 <- SCTransform(seurat4, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
seurat5 <- SCTransform(seurat5, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
seurat6 <- SCTransform(seurat6, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
seurat7 <- SCTransform(seurat7, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
seurat8 <- SCTransform(seurat8, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

# Merge seurat objects
seurat.list <- c(seurat1,seurat2,seurat3,seurat4,seurat5,seurat6,seurat7,seurat8)

rm(seurat1,seurat2,seurat3,seurat4,seurat5,seurat6,seurat7,seurat8)
gc()

seurat.merged <- merge(seurat.list[[1]], y = seurat.list[2:length(seurat.list)], 
                       add.cell.ids = c("NR1","NR2","R1","R2","R3","NR3","R4","NR4"), 
                       project = "response_pred")

# Intersecting variable features must be selected for downstream PCA and UMAP.
seurat.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
VariableFeatures(seurat.merged[["SCT"]]) <- seurat.features

# Remove separate seurat objects that will not be used any more (to free up memory)
rm(seurat.list)
gc() # run garbage collect to free up memory

#Dimensionality reduction and clustering
seurat.merged <- RunPCA(seurat.merged)
seurat.merged <- RunUMAP(seurat.merged,dims=1:30)
seurat.merged <- FindNeighbors(seurat.merged,dims=1:30)
seurat.merged <- FindClusters(seurat.merged,resolution=1.0)
DimPlot(seurat.merged, reduction = "umap",label=T,label.size=5,repel=T,label.box=T)
FeaturePlot(seurat.merged,features=c("CD4","CD8A"))

setwd("~")
setwd(".../scRNAseq_reanalysis/Luoma2020/data/CD3_FACS/")
saveRDS(seurat.merged, file = "seurat_Luoma.rds")

setwd("~")
setwd(".../scRNAseq_reanalysis/Luoma2020/data/CD3_FACS/")
seurat.merged <- readRDS("seurat_Luoma.rds")
setwd("~")
setwd(".../scRNAseq_reanalysis/")

# Select CD4+ T cells only
avg_exp <- AggregateExpression(seurat.merged)
avg_exp <- as.data.frame(avg_exp$SCT)
CD8A_expr <- avg_exp[rownames(avg_exp)=="CD8A",]/(rowSums(avg_exp[rownames(avg_exp)=="CD8A",])/ncol(avg_exp)) # Divide by mean as normalization measure
CD4_expr <- avg_exp[rownames(avg_exp)=="CD4",]/(rowSums(avg_exp[rownames(avg_exp)=="CD4",])/ncol(avg_exp)) # Divide by mean as normalization measure
CD4_expr-CD8A_expr

#Select all clusters with CD4 (mean-normalized) > CD8A (mean-normalized)
subset_CD4 <- subset(seurat.merged,idents=c(2,3,4,5,7,8,10,12,13,16))
subset_CD4@meta.data$group <- as.factor(ifelse(grepl("^NR",substr(Cells(subset_CD4),1,2)),"NR","R"))
subset_CD4@meta.data$orig.cluster <- Idents(subset_CD4) #Store cluster for each cell
subset_CD4@meta.data$group_cluster <- paste(Idents(subset_CD4),subset_CD4@meta.data$group,sep="_")
cell_ids <- colnames(subset_CD4)
subset_CD4$Patient <- sapply(strsplit(cell_ids, "_"), `[`, 1)

Luoma_data <- subset_CD4

#Clear up memory
rm(avg_exp,CD4_expr,CD8A_expr,seurat.merged,subset_CD4)
gc()

#Combine all datasets
Gupta_data$Source <- "Gupta_2024"
Gupta_data$Patient <- Gupta_data$Donor
Assays(Gupta_data)
Luoma_data$Source <- "Luoma_2020"
# Luoma_data$Patient already exists
Assays(Luoma_data)
Thomas_data$Source <- "Thomas_2024"
Thomas_data$Patient <- sub("^(([^_]*_){2}[^_]*).*", "\\1", Thomas_data$Patient)
Assays(Thomas_data)

Gupta_data@meta.data <- Gupta_data@meta.data %>% dplyr::select(Patient,Source)
Luoma_data@meta.data <- Luoma_data@meta.data %>% dplyr::select(Patient,Source)
Thomas_data@meta.data <- Thomas_data@meta.data %>% dplyr::select(Patient,Source)

#Only keep features that are detected in all datasets
common_features_RNA <- Reduce(intersect, list(rownames(Gupta_data[["RNA"]]), rownames(Luoma_data[["RNA"]]), rownames(Thomas_data[["RNA"]])))
common_features_SCT <- Reduce(intersect, list(rownames(Gupta_data[["SCT"]]), rownames(Luoma_data[["SCT"]]), rownames(Thomas_data[["SCT"]])))

Gupta_data[["RNA"]] <- subset(Gupta_data[["RNA"]], features = common_features_RNA)
Luoma_data[["RNA"]] <- subset(Luoma_data[["RNA"]], features = common_features_RNA)
Thomas_data[["RNA"]] <- subset(Thomas_data[["RNA"]], features = common_features_RNA)

Gupta_data[["SCT"]] <- subset(Gupta_data[["SCT"]], features = common_features_SCT)
Luoma_data[["SCT"]] <- subset(Luoma_data[["SCT"]], features = common_features_SCT)
Thomas_data[["SCT"]] <- subset(Thomas_data[["SCT"]], features = common_features_SCT)

# Ensure unique cell names for each object
Gupta_data <- RenameCells(Gupta_data, add.cell.id = "Gupta")
Luoma_data <- RenameCells(Luoma_data, add.cell.id = "Luoma")
Thomas_data <- RenameCells(Thomas_data, add.cell.id = "Thomas")

#Now w e want to drop all Seurat elements (like reductions) that are not present for each object
Gupta_data@graphs <- list()
Gupta_data@reductions <- list()
Gupta_data@commands <- list()
Gupta_data@active.ident <- factor()
Gupta_data@assays$HTO <- NULL
Gupta_data@assays$ADT <- NULL
Gupta_data@active.ident <- factor(rep("Gupta",ncol(Gupta_data)))
Gupta_data <- SetIdent(object = Gupta_data, value = "Gupta")
Gupta_data@assays$SCT@key <- "sct_"

Luoma_data@graphs <- list()
Luoma_data@reductions <- list()
Luoma_data@commands <- list()
Luoma_data@active.ident <- factor()
Luoma_data@active.ident <- factor(rep("Luoma",ncol(Luoma_data)))
Luoma_data <- SetIdent(object = Luoma_data, value = "Luoma")
Luoma_data@assays$SCT@key <- "sct_"

Thomas_data@graphs <- list()
Thomas_data@reductions <- list()
Thomas_data@commands <- list()
Thomas_data@active.ident <- factor()
Thomas_data@active.ident <- factor(rep("Thomas",ncol(Thomas_data)))
Thomas_data <- SetIdent(object = Thomas_data, value = "Thomas")
Thomas_data@assays$SCT@key <- "sct_"

merged <- merge(x=Gupta_data,y=list(Luoma_data,Thomas_data))
merged[["RNA"]] <- split(merged[["RNA"]],f=merged$Source)

#Create UMAP of unintegrated datasets
options(future.globals.maxSize = 4000 * 1024^2)
merged2 <- merged
merged2 <- SCTransform(merged2)
merged2 <- RunPCA(merged2, npcs = 30, verbose = T)
merged2 <- FindNeighbors(merged2, dims = 1:30, reduction = "pca")
merged2 <- FindClusters(merged2, resolution = 2)
avg_exp <- AggregateExpression(merged2)
avg_exp <- as.data.frame(avg_exp$SCT)
CD8A_expr <- avg_exp[rownames(avg_exp)=="CD8A",]/(rowSums(avg_exp[rownames(avg_exp)=="CD8A",])/ncol(avg_exp)) # Divide by mean as normalization measure
CD4_expr <- avg_exp[rownames(avg_exp)=="CD4",]/(rowSums(avg_exp[rownames(avg_exp)=="CD4",])/ncol(avg_exp)) # Divide by mean as normalization measure
CD4_expr-CD8A_expr
clusters_to_drop <- c("10","25","29","33") #all clusters with CD4-CD8A < -1
CD4_subset2 <- subset(merged2, idents = setdiff(Idents(merged2), clusters_to_drop))
CD4_subset2 <- RunPCA(CD4_subset2, npcs = 30, verbose = T)
CD4_subset2 <- FindNeighbors(CD4_subset2, dims = 1:30, reduction = "pca")
CD4_subset2 <- FindClusters(CD4_subset2, resolution = 1)
CD4_subset2 <- RunUMAP(CD4_subset2, reduction = "pca", dims = 1:30, reduction.name = "umap")
colors <- c("#F98E09","#BC3754","#440154")

setwd("~")
setwd(".../")

pdf("dimplot_not_integrated.pdf",width=6,height=5)
DimPlot(CD4_subset2,reduction = "umap",group.by = c("Source"),cols=colors, label.size = 2)
dev.off()

rm(CD4_subset2,merged2)
gc()

setwd("~")
setwd(".../")
saveRDS(merged, file = "merged.rds")
gc()

merged <- readRDS("merged.rds")
#Integrate data and create new UMAP
options(future.globals.maxSize = 6000 * 1024^2)
merged <- SCTransform(merged)
merged <- RunPCA(merged, npcs = 30, verbose = T)

merged <- IntegrateLayers(
  object = merged,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = T
)
merged <- FindNeighbors(merged, dims = 1:30, reduction = "integrated.dr")
merged <- FindClusters(merged, resolution = 2)
avg_exp <- AggregateExpression(merged)
avg_exp <- as.data.frame(avg_exp$SCT)
CD8A_expr <- avg_exp[rownames(avg_exp)=="CD8A",]/(rowSums(avg_exp[rownames(avg_exp)=="CD8A",])/ncol(avg_exp)) # Divide by mean as normalization measure
CD4_expr <- avg_exp[rownames(avg_exp)=="CD4",]/(rowSums(avg_exp[rownames(avg_exp)=="CD4",])/ncol(avg_exp)) # Divide by mean as normalization measure
CD4_expr-CD8A_expr

#Drop CD8A-expressing clusters
clusters_to_drop <- c("13","16","19","23","25","33") #drop all with CD4 - CD8A < -0.1
CD4_subset <- subset(merged, idents = setdiff(Idents(merged), clusters_to_drop))

#Run PCA and clustering again
CD4_subset <- RunPCA(CD4_subset, npcs = 30, verbose = T)
CD4_subset <- FindNeighbors(CD4_subset, dims = 1:30, reduction = "integrated.dr")
CD4_subset <- FindClusters(CD4_subset, resolution = 1)
CD4_subset <- RunUMAP(CD4_subset, reduction = "integrated.dr", dims = 1:30, reduction.name = "umap.integrated")

#saveRDS(CD4_subset,"CD4_subset.RDS")
setwd("~")
setwd(".../")
CD4_subset <- readRDS("CD4_subset.RDS")

mycol1 <- colorRampPalette(brewer.pal(8,"Set1"))(25)

pdf("dimplot_integrated.pdf",width=6,height=5)
DimPlot(CD4_subset,reduction = "umap.integrated",group.by = c("Source"), cols=colors,label.size = 2)
dev.off()

pdf("dimplot_25_clusters.pdf",width=6,height=5)
DimPlot(CD4_subset,reduction = "umap.integrated", cols=mycol1, label.size = 2)
dev.off()

colors2 <- c("grey","#F98E09","#440154")
pdf("featurePlot_IL17_IL22.pdf",width=15,height=5)
FeaturePlot(CD4_subset,features = c("IL17A","IL22"),blend=TRUE,cols = colors2)
dev.off()

FeaturePlot(CD4_subset,features = c("CD3D","CD4"),blend=TRUE)

#Find markers per cluster
interest.markers.all <- FindAllMarkers(CD4_subset,only.pos=T,min.pct=0.1,logfc.threshold = 0.1,recorrect_umi=FALSE)
check <- interest.markers.all %>% group_by(cluster) %>% slice_max(n=40,order_by=avg_log2FC)

#Specifically check Th17-associated genes
interest.markers.Th17 <- FindAllMarkers(CD4_subset,only.pos=T,features=c("IL17A"),recorrect_umi=FALSE)

#Cell counts per patient
cell_counts <- as.data.frame(table(CD4_subset@meta.data$Patient,Idents(CD4_subset)))
colnames(cell_counts) <- c("Patient","Cluster","Value")
cell_counts <- pivot_wider(cell_counts,names_from = "Cluster",values_from = "Value")
cell_counts$Patient <- sub("([^_]*_[^_]*)_.*", "\\1", cell_counts$Patient)
cell_counts <- cell_counts %>%
  group_by(Patient) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

cell_counts$Total <- rowSums(cell_counts[,c(2:ncol(cell_counts))])

for(i in 2:(ncol(cell_counts)-1)){
  cell_counts[,i] <- cell_counts[,i]/cell_counts[,ncol(cell_counts)]
}

cell_counts$Total <- cell_counts$Total/cell_counts$Total
cell_counts <- melt(cell_counts,id.vars = "Patient")
metadata <- read_excel("metadata_complete.xlsx")
cell_counts <- merge(cell_counts,metadata,by="Patient",all=T)

cell_counts <- cell_counts %>% drop_na(value)

#Th17 cluster = 18
#Th22 cluster = 19

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

data <- cell_counts
data <- data %>% dplyr::filter(variable!="Total")
data$Response <- factor(data$Response,levels=c("R","NR"),labels=c("Responder","Non-responder"))

plot <- ggplot(data,aes(variable,100*value,fill=Response))+
  geom_split_violin(trim=T,position=position_dodge(width=0.5),scale="width")+
  scale_fill_manual(values=c("#5ec962","#440154"))+
  ylab("% cells of all CD4+ T cells")+
  xlab("Cluster")+
  theme_classic()+
  theme(legend.title=element_blank())
pdf("scRNAseq_all_subsets.pdf",width=9,height=6)
plot(plot)
dev.off()

test_table <- data.frame(cluster=NA,pval=NA)
r <- 1
for(i in c(0:24)){
  data_sub <- data %>% dplyr::filter(variable==i)
  test <- wilcox_test(data_sub,value~Response)
  test_table[r,1] <- i
  test_table[r,2] <- test$p
  r <- r+1
}

data <- cell_counts %>% dplyr::filter(variable==18 | variable==19)
data$Response <- factor(data$Response,levels=c("R","NR"),labels=c("Responder","Non-responder"))
summed_Th17_Th22 <- data %>%
  group_by(Patient, Response) %>%
  summarise(total_value = sum(value))
summed_Th17_Th22$variable <- "1"
plot <- ggplot(summed_Th17_Th22,aes(variable,100*total_value,fill=Response))+
  geom_boxplot(outlier.shape=NA,alpha=0.7)+
  scale_fill_manual(values=c("#5ec962","#440154"))+
  geom_point(shape=21,position=position_jitterdodge())+
  theme_classic()+
  geom_signif(annotation=paste0(expression(italic(P) == 0.95)),parse=T,y_position=10,xmin=0.8,xmax=1.2,tip_length=c(0,0))+
  ylab("% IL17A and/or IL22 expressing cells of all CD4+ T cells")+
  xlab("")+
  theme(axis.text.x = element_blank(),
        legend.title=element_blank())
pdf("scRNAseq_plot_Th17_Th22.pdf",width=4,height=5)
plot(plot)
dev.off()
wilcox.test(summed_Th17_Th22$total_value~summed_Th17_Th22$Response)

#Pseudo-bulk analysis gene set
CD4_prep <- CD4_subset
Idents(CD4_prep) <- CD4_subset@meta.data$Patient
CD4_prep_exp <- AggregateExpression(CD4_prep)
CD4_prep_exp_SCT <- as.data.frame(CD4_prep_exp$SCT)
CD4_prep_exp_RNA <- as.data.frame(CD4_prep_exp$RNA)

Th17.1_selective <- c("ME1","ENPP1","ADAM23","COLQ","SLC4A10","CA2","SCRN1",
                      "ABCB1","CLEC2B","CCR2","DPP4","KLRB1")
Th17_core <- c("IL17RE","PIK3R6","HLF","CCR6","RORC","CTSH","LPR12","CHN1",
               "HPGD","ADAM12","PTPN13")
Th17_selective <- c("SHF","IL17A","IL17F","ATP1B1","S100A4","NTRK2","PLD1",
                    "PPARG","ITGAE","PAK3","TMOD1","IL26","TIMP4","HRH4","TIMP1")

gene.sets <- list(Th17.1_selective=Th17.1_selective,
                  Th17_selective=Th17_selective,
                  Th17_core=Th17_core)

pseudobulk_GSVA_param <- gsvaParam(as.matrix(CD4_prep_exp_SCT),
                               geneSets = gene.sets,
                               kcdf="Gaussian")
pseudobulk_GSVA <- gsva(pseudobulk_GSVA_param)

results <- as.data.frame(t(pseudobulk_GSVA))
results$Patient <- rownames(results)
results$Patient <- sub("-Colon", "", results$Patient)
results$Patient <- sub("-RectosigmoidColon", "", results$Patient)
results$Patient <- sub("-", "_", results$Patient)
results <- merge(results,metadata,by="Patient")

plot <- ggplot(results,aes(Response,Th17.1_selective,fill=Response))+
  geom_boxplot(outlier.shape=NA,alpha=0.7)+
  scale_fill_manual(values=c("#5ec962","#440154"))+
  geom_point(shape=21,position=position_jitterdodge())+
  theme_classic()+
  geom_signif(annotation=paste0(expression(italic(P) == 0.99)),parse=T,y_position=0.8,xmin=1,xmax=2,tip_length=c(0,0))+
  ylab("GSVA score Th17.1 signature")+
  xlab("")+
  theme(axis.text.x = element_blank(),
        legend.title=element_blank())
pdf("scRNAseq_plot_Th17_1_signature.pdf",width=4,height=5)
plot(plot)
dev.off()
wilcox.test(results$Th17.1_selective~results$Response)

#Differential gene expression non-responder versus responder
CD4_prep <- CD4_subset
Idents(CD4_prep) <- CD4_subset@meta.data$Patient
CD4_prep_exp <- AggregateExpression(CD4_prep)
CD4_prep_exp_RNA <- as.data.frame(CD4_prep_exp$RNA)
bulk_data <- CD4_prep_exp_RNA
colnames(bulk_data) <- sub("-Colon", "", colnames(bulk_data))
colnames(bulk_data) <- sub("-RectosigmoidColon", "", colnames(bulk_data))
colnames(bulk_data) <- sub("-", "_", colnames(bulk_data))
metadata <- as.data.frame(read_excel("metadata_complete.xlsx"))
rownames(metadata) <- metadata$Patient
metadata <- metadata[order(match(rownames(metadata),colnames(bulk_data))),]
all(rownames(metadata) == colnames(bulk_data))

dds <- DESeqDataSetFromMatrix(countData = bulk_data,
                       colData = metadata,
                       design = ~ Dataset + Response)
rownames(bulk_data)
dds <- DESeq(dds)
resultsNames(dds)
res1 <- results(dds,name="Response_R_vs_NR")
volcano_results <- as.data.frame(res1@listData)
rownames(volcano_results) <- res1@rownames
volcano_results$log2FoldChange <- -volcano_results$log2FoldChange

volcano_results$label <- ifelse(volcano_results$pvalue<0.05 & abs(volcano_results$log2FoldChange)>1,rownames(volcano_results),"")
volcano_results$color <- ifelse(volcano_results$pvalue<0.05 & volcano_results$log2FoldChange>1,"NR_up_sign",
                                ifelse(volcano_results$pvalue<0.05 & volcano_results$log2FoldChange<(-1),"NR_dn_sign",
                                       ifelse(volcano_results$pvalue>=0.05 & volcano_results$log2FoldChange<(-1),"NR_dn",
                                              ifelse(volcano_results$pvalue>=0.05 & volcano_results$log2FoldChange>1,"NR_up",""))))
volcano_results$color <- factor(volcano_results$color,levels=c("","NR_up","NR_up_sign","NR_dn","NR_dn_sign"))
labels_to_show <- rownames(volcano_results)#c("SELL","MT1F","MT1M","PDPN","SIGLEC9","WARS1","ALPL","CXCL8",
                    #"LILRB3","CXCL11","CXCL9")
volcano_NR_v_R <- ggplot(volcano_results,aes(log2FoldChange,-log10(pvalue),color=color))+
  geom_point()+
  geom_hline(yintercept = 1.3,linetype="dashed")+
  geom_vline(xintercept = 1,linetype="dashed")+
  geom_vline(xintercept = -1,linetype="dashed")+
  scale_color_manual(values=c("grey","#8E6698","#440154","#9edea0","#5ec962"))+
  geom_text_repel(aes(label=ifelse(label %in% labels_to_show,label,"")),max.overlaps=Inf)+
  theme_minimal()+
  theme(legend.position="none")+
  #xlim(-3.5,3)+
  #ylim(0,1.8)+
  ylab("-log10(raw P value)")+
  xlab("log2FC Non-responder versus Responder")
pdf("volcano_NR_v_R_RAW_pval.pdf",width=6,height=5)
plot(volcano_NR_v_R)
dev.off()

#Adjusted Pvals
res1 <- results(dds,name="Response_R_vs_NR")
volcano_results <- as.data.frame(res1@listData)
rownames(volcano_results) <- res1@rownames
volcano_results$log2FoldChange <- -volcano_results$log2FoldChange
volcano_results$label <- ifelse(volcano_results$padj<0.05 & abs(volcano_results$log2FoldChange)>2,rownames(volcano_results),"")
volcano_results$color <- ifelse(volcano_results$padj<0.05 & volcano_results$log2FoldChange>1,"NR_up_sign",
                                ifelse(volcano_results$padj<0.05 & volcano_results$log2FoldChange<(-1),"NR_dn_sign",
                                       ifelse(volcano_results$padj>=0.05 & volcano_results$log2FoldChange<(-1),"NR_dn",
                                              ifelse(volcano_results$padj>=0.05 & volcano_results$log2FoldChange>1,"NR_up",""))))
volcano_results$color <- factor(volcano_results$color,levels=c("","NR_up","NR_up_sign","NR_dn","NR_dn_sign"))
labels_to_show <- rownames(volcano_results)#c("SELL","MT1F","MT1M","PDPN","SIGLEC9","WARS1","ALPL","CXCL8",
#"LILRB3","CXCL11","CXCL9")
volcano_NR_v_R <- ggplot(volcano_results,aes(log2FoldChange,-log10(padj),color=color))+
  geom_point()+
  geom_hline(yintercept = 1.3,linetype="dashed")+
  geom_vline(xintercept = 1,linetype="dashed")+
  geom_vline(xintercept = -1,linetype="dashed")+
  scale_color_manual(values=c("grey","#8E6698","#9edea0"))+
  geom_text_repel(aes(label=ifelse(label %in% labels_to_show,label,"")),max.overlaps=Inf)+
  theme_minimal()+
  theme(legend.position="none")+
  #xlim(-3.5,3)+
  #ylim(0,1.8)+
  ylab("-log10(Padj)")+
  xlab("log2FC Non-responder versus Responder")
pdf("volcano_NR_v_R_Padj.pdf",width=4,height=5)
plot(volcano_NR_v_R)
dev.off()


volcano_results$Gene <- rownames(volcano_results)
select <- volcano_results[volcano_results$Gene=="MS4A1",]

