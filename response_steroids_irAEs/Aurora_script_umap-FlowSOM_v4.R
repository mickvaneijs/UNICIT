library(flowCore)
library(flowViz)
library(flowVS)
library(flowAI)
library(PeacoQC)
library(devtools)
#install_github("ssmpsn2/flowAssist",force=TRUE)
#install_github("sararselitsky/FastPG",force=TRUE)
library(FastPG)
library(CytoNorm)
library(SingleCellExperiment)
library(writexl)
library(readxl)
library(ggplot2)
library(dplyr)
library(uwot)
library(Seurat)
library(stringr)
library(tidyr)
library(writexl)
library(ggridges)
library(viridis)
library(RColorBrewer)
library(clustree)
library(FlowSOM)
library(flowAssist)
library(Polychrome)
library(pheatmap)
#library(scCustomize)

set.seed(1234)
setwd("~")
setwd(".../Analysis/")

dataset <- readRDS(file="full_dataset.rds")
gc()

dataset$timepoint <- ifelse(grepl("T1",dataset$sample),1,
                            ifelse(grepl("T2",dataset$sample),2,
                                   ifelse(grepl("T3",dataset$sample),3,
                                          ifelse(grepl("T4",dataset$sample),4,
                                                 ifelse(grepl("HD",dataset$sample),1,"Check")))))
dataset$subject <- str_extract(dataset$sample,"[^_]+")

#Defines parameters to be used for dimensionality reduction and clustering
fsom_markers <- c("CD57 APC","CD56 APC-Fire750","FoxP3 AF488","CD161 AF700",
                  "CD39 BB700","CD3 BUV496","CD16 BV570","CCR7 BV605",
                  "HLA.DR PerCP","CD4 PerCP-Fire806",
                  "CD8 SparkBlue550","CD14 SparkBlue574","CD45RA SparkUV387")

umap_markers <- c("FSC-A","SSC-A","CD57 APC","CD56 APC-Fire750","FoxP3 AF488","CD161 AF700",
                  "CD39 BB700","CD3 BUV496","CD16 BV570","CCR7 BV605",
                  "HLA.DR PerCP","CD4 PerCP-Fire806",
                  "CD8 SparkBlue550","CD14 SparkBlue574","CD45RA SparkUV387")

meta_columns <- c("FSC-A","FSC-H","SSC-A","SSC-B-A","SSC-B-H",              
                  "SSC-H","Autofluorescence","CD57 APC","CD56 APC-Fire750","FoxP3 AF488",          
                  "CD161 AF700","CD39 BB700","CD3 BUV496","CD38 BUV563","CD134 BUV661",         
                  "TNFa BUV737","TCF1 BV421","IntB7 BV480","CD16 BV570","CCR7 BV605",           
                  "CD25 BV650","IL17A BV711","PD1 BV750","CXCR3 BV785","IL10 PE",              
                  "TIM3 PE-Cy5","Ki67 PE-Cy7","LAG3 PE-Dazzle594","CD69 PE-Fire640","CCR4 PE-Fire810",      
                  "GzmB PB","HLA.DR PerCP","IL13 PerCP-Cy5.5","CD4 PerCP-Fire806","CTLA4 PerCP-eFluor710",
                  "IL5 RY586","CD8 SparkBlue550","CD14 SparkBlue574","IFNy SparkNIR695","CD45RA SparkUV387",    
                  "IntB1 SuperBright436","dye","Time","Original_ID","sample","timepoint","subject",
                  "clr_clusters","clr_metaclusters")    

#Create FlowSOM for entire dataset
numeric_cols <- colnames(dataset[sapply(dataset,is.numeric)])
dataset_ff <- DFtoFF(dataset[,numeric_cols])
names(dataset_ff@exprs[,fsom_markers]) <- gsub(x=names(dataset_ff@exprs[,fsom_markers]),pattern=" .*",replacement = "")
for(i in colnames(dataset_ff@exprs[,fsom_markers])){
  print(paste0("Mean is ",mean(dataset_ff@exprs[,i])," and minimum is ",min(dataset_ff@exprs[,i])))
  minimum <- min(dataset_ff@exprs[,i])
  dataset_ff@exprs[,i] <- dataset_ff@exprs[,i]-minimum
  print(paste0("Mean is ",mean(dataset_ff@exprs[,i])," and minimum is ",min(dataset_ff@exprs[,i])))
  print(paste0("Done with column ",i))
  print("")
}

#Perform CLR normalization from Seurat (https://divingintogeneticsandgenomics.com/post/details-in-centered-log-ratio-clr-normalization-for-cite-seq-protein-count-data/#:~:text=CLR%20normalization%20in%20detail&text=Seurat%20CLR%20removes%200%20counts,Seurat%20CLR%20normalization.&text=We%20do%20not%20see%20negative%20values.)
clr_function <- function(x) {
  return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
}
dataset_ff@exprs[,fsom_markers] <- apply(dataset_ff@exprs[,fsom_markers],MARGIN=2,FUN=clr_function)

#Run flowSOM
fSOM_CLR <- FlowSOM(dataset_ff,
                       compensate=F,
                       transform=F,
                       scale=F,
                       colsToUse = fsom_markers,
                       nClus=50)

plot <- PlotStars(fSOM_CLR,backgroundValues=fSOM_CLR$metaclustering)
pdf("PlotStars_CLR.pdf")
plot(plot)
dev.off()

fsom_data_clr_clusters <- as.data.frame(GetClusters(fSOM_CLR))
fsom_data_clr_metaclusters <- as.data.frame(GetMetaclusters(fSOM_CLR))
colnames(fsom_data_clr_clusters) <- "clr_clusters"
colnames(fsom_data_clr_metaclusters) <- "clr_metaclusters"
rownames(fsom_data_clr_clusters) <- rownames(dataset)
rownames(fsom_data_clr_metaclusters) <- rownames(dataset)
dataset$clr_clusters <- fsom_data_clr_clusters$clr_clusters
dataset$clr_metaclusters <- fsom_data_clr_metaclusters$clr_metaclusters

#Subsample dataset for UMAP and Seurat clustering
sub_data <- dataset[sample(nrow(dataset), 0.01*nrow(dataset)), ]
cluster_data <- sub_data[,umap_markers]
names(cluster_data) <- gsub(x=names(cluster_data),pattern=" .*",replacement = "")
for(i in 1:ncol(cluster_data)){
 print(paste0("Mean is ",mean(cluster_data[,i])," and minimum is ",min(cluster_data[,i])))
 minimum <- min(cluster_data[,i])
 cluster_data[,i] <- cluster_data[,i]-minimum
 print(paste0("Mean is ",mean(cluster_data[,i])," and minimum is ",min(cluster_data[,i])))
 print(paste0("Done with column ",i))
 print("")
}

cluster_data <- as.data.frame(t(cluster_data))
cluster_data <- CreateSeuratObject(cluster_data)
meta_data <- sub_data[,meta_columns]

meta_add <- meta_data[rownames(cluster_data@meta.data),meta_columns]
cluster_data@meta.data[,meta_columns] <- meta_add

cluster_data <- NormalizeData(cluster_data,normalization.method = "CLR")
cluster_data <- ScaleData(cluster_data)
VariableFeatures(cluster_data)<-rownames(cluster_data)
cluster_data <- RunPCA(cluster_data)
cluster_data <- FindNeighbors(cluster_data,dims=1:14)

#Find optimal number of clusters
# cluster_data <- FindClusters(cluster_data,resolution=0.0)
# cluster_data <- FindClusters(cluster_data,resolution=0.2)
# cluster_data <- FindClusters(cluster_data,resolution=0.4)
# cluster_data <- FindClusters(cluster_data,resolution=0.6)
# cluster_data <- FindClusters(cluster_data,resolution=0.8)
# cluster_data <- FindClusters(cluster_data,resolution=1.0)
# cluster_data <- FindClusters(cluster_data,resolution=1.2)
# cluster_data <- FindClusters(cluster_data,resolution=1.4)
# cluster_data <- FindClusters(cluster_data,resolution=1.6)
# cluster_data <- FindClusters(cluster_data,resolution=1.8)
# cluster_data <- FindClusters(cluster_data,resolution=2.0)
# cluster_data <- FindClusters(cluster_data,resolution=2.2)
# cluster_data <- FindClusters(cluster_data,resolution=2.4)
# cluster_data <- FindClusters(cluster_data,resolution=2.6)
# cluster_data <- FindClusters(cluster_data,resolution=2.8)
# cluster_data <- FindClusters(cluster_data,resolution=3.0)

# pdf("Clustree_1_percent.pdf",width=10,height=10)
# plot(clustree(cluster_data))
# dev.off()

cluster_data <- FindClusters(cluster_data,resolution=0.8)
cluster_data <- RunUMAP(cluster_data,dims=1:14,return.model = TRUE)

DimPlot(cluster_data,reduction = "umap",raster=FALSE,label = TRUE,group.by="clr_clusters")
DimPlot(cluster_data,reduction = "umap",raster=FALSE,label = TRUE,group.by="clr_metaclusters")
DimPlot(cluster_data,reduction = "umap",raster=FALSE,label = TRUE)

cluster_relationships <- as.data.frame.matrix(table(cluster_data$clr_metaclusters,cluster_data$seurat_clusters))
cluster_relationships$row_sum <- rowSums(cluster_relationships)
for(i in 1:(ncol(cluster_relationships)-1)){
  cluster_relationships[,i] <- (cluster_relationships[,i]/cluster_relationships[,ncol(cluster_relationships)])*100
}
cluster_relationships$row_max <- apply(cluster_relationships[,1:(ncol(cluster_relationships)-1)],1,max,na.rm=T)

second_max <- cluster_relationships
for(i in 1:(ncol(second_max)-1)){
  second_max[,i] <- (second_max[,i]-second_max[,ncol(second_max)])
  second_max[,i] <- ifelse(second_max[,i]==0,NA,second_max[,i])
  second_max[,i] <- second_max[,i]+second_max[,ncol(second_max)]
}
second_max$second_max <- apply(second_max[,1:(ncol(second_max)-2)],1,max,na.rm=T)
second_max$first_and_second <- second_max$row_max+second_max$second_max

message(paste0("Percentage of cells belonging with >95% certainty to one Seurat cluster: ",
               round((sum(cluster_relationships$row_sum[cluster_relationships$row_max>95])/ncol(cluster_data))*100)),"%")
message(paste0("Percentage of cells belonging with >90% certainty to one Seurat cluster: ",
               round((sum(cluster_relationships$row_sum[cluster_relationships$row_max>90])/ncol(cluster_data))*100)),"%")
message(paste0("Percentage of cells belonging with >85% certainty to one Seurat cluster: ",
               round((sum(cluster_relationships$row_sum[cluster_relationships$row_max>85])/ncol(cluster_data))*100)),"%")
message(paste0("Percentage of cells belonging with >80% certainty to one Seurat cluster: ",
               round((sum(cluster_relationships$row_sum[cluster_relationships$row_max>80])/ncol(cluster_data))*100)),"%")
message(paste0("Percentage of cells belonging with >75% certainty to one Seurat cluster: ",
               round((sum(cluster_relationships$row_sum[cluster_relationships$row_max>75])/ncol(cluster_data))*100)),"%")
message(paste0("Percentage of cells belonging with >50% certainty to one Seurat cluster: ",
               round((sum(cluster_relationships$row_sum[cluster_relationships$row_max>50])/ncol(cluster_data))*100)),"%")

cluster_relationships_simple <- cluster_relationships
for(i in 1:ncol(cluster_relationships_simple)){
  for(j in 1:nrow(cluster_relationships_simple)){
    cluster_relationships_simple[j,i] <- ifelse(cluster_relationships_simple[j,i]<5,".",cluster_relationships_simple[j,i])
  }
}

cluster_data_exprs <- as.data.frame(t(cluster_data@assays$RNA@scale.data))
cluster_data_meta <- cluster_data@meta.data
cluster_data_meta <- as.data.frame(cluster_data_meta[,c("clr_metaclusters")])
colnames(cluster_data_meta) <- "clr_metaclusters"
cluster_data_merged <- cbind(cluster_data_exprs,cluster_data_meta)
cluster_data_heatmap <- cluster_data_merged %>% group_by(clr_metaclusters) %>% summarise_all(mean)
cluster_data_heatmap <- as.data.frame(cluster_data_heatmap)
rownames(cluster_data_heatmap) <- cluster_data_heatmap$clr_metaclusters
cluster_data_heatmap <- cluster_data_heatmap[,-1]
  
pheatmap(cluster_data_heatmap,
         scale="row",
         color=colorRampPalette(c("navy","white","darkred"))(25))

for(i in rownames(cluster_data@assays$RNA)){
  plot_feature <- FeaturePlot(object=cluster_data,features=paste0(i),min.cutoff='q10',max.cutoff='q99',order=T,raster=F) & scale_color_gradientn(colors=plasma(n=10))
  pdf(paste0("Featureplot_",i,".pdf"))
  plot(plot_feature)
  dev.off()
}

meta_data_features <- c("FSC-H","SSC-B-A","SSC-B-H","SSC-H","Autofluorescence","CD38 BUV563",
                        "CD134 BUV661","TNFa BUV737","TCF1 BV421","IntB7 BV480","CD25 BV650",
                        "IL17A BV711","PD1 BV750","CXCR3 BV785","IL10 PE","TIM3 PE-Cy5","Ki67 PE-Cy7",          
                        "LAG3 PE-Dazzle594","CD69 PE-Fire640","CCR4 PE-Fire810","GzmB PB","HLA.DR PerCP",
                        "IL13 PerCP-Cy5.5","CTLA4 PerCP-eFluor710","IL5 RY586","IFNy SparkNIR695",
                        "CD45RA SparkUV387","IntB1 SuperBright436","dye")      

for(i in meta_data_features){
  plot_feature <- FeaturePlot(object=cluster_data,features=paste0(i),min.cutoff='q10',max.cutoff='q99',order=T,raster=F) & scale_color_gradientn(colors=plasma(n=10))
  pdf(paste0("Featureplot_",i,".pdf"))
  plot(plot_feature)
  dev.off()
}

#Save dataset:
write_xlsx(cluster_relationships,"cluster_relationships.xlsx")
write_xlsx(cluster_relationships_simple,"cluster_relationships_simple.xlsx")
saveRDS(dataset,file="full_mapped_dataset.rds")
saveRDS(sub_data,file="mapped_subset.rds")
saveRDS(cluster_data,file="mapped_seurat_object.RObj")
saveRDS(fSOM,file="FlowSOM.RObj")

rm(fsom_data)
gc()
