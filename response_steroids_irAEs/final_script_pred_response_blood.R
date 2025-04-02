library(writexl)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(ggridges)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(patchwork)
library(ggbiplot)
library(ggpubr)
library(ggrepel)
library(Seurat)
library(flowCore)
library(clustree)
library(FlowSOM)
library(flowAssist)
library(clustree)
library(stats)
library(coin) #for exact wilcox test with ties
library(car) #for vif() function
library(pROC)
library(caret)
library(MASS)
library(randomForest)
library(glmnet)
library(table1)
library(presto)
library(alphahull)

#Custom functions
#Get density from data points
get_density<-function(x,y,...){
  dens<-MASS::kde2d(x,y,...)
  ix<-findInterval(x,dens$x)
  iy<-findInterval(y,dens$y)
  ii<-cbind(ix,iy)
  return(dens$z[ii])
}

#Save pheatmap
save_pheatmap_pdf<-function(x,filename,width=5,height=9){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename,width=width,height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

setwd("~")
setwd(".../")
dataset<-readRDS(file="full_mapped_dataset.rds")

#Create heatmap of marker expression
cluster_relationships<-read_xlsx("cluster_relationships.xlsx")
cluster_relationships<-as.data.frame(cluster_relationships)
colnames(cluster_relationships)<-paste0("Seurat_",gsub("X",",",colnames(cluster_relationships)))
cluster_relationships<-as.data.frame(lapply(cluster_relationships,as.numeric))
heatmap_data<-cluster_relationships[,c(1:(ncol(cluster_relationships)-2))]
heatmap_data<-heatmap_data[c(43,22,46,31,2,5,26,41,45,39,8,9,11,19,25,29,
                             18,36,44,47,3,20,49,35,13,16,17,14,15,23,27,
                             28,32,34,48,24,33,30,10,4,42,38,40,50,6,12,
                             7,37,1),
                           c(1:26)]
rownames(heatmap_data)<-paste0("FlowSOM_",rownames(heatmap_data))

plot1<-pheatmap(heatmap_data,
                cluster_cols=F,
                cluster_rows=F,
                breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                color=colorRampPalette(c("navy","white","darkred"))(10))

save_pheatmap_pdf<-function(x,filename,width=5,height=9){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename,width=width,height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(plot1,"Heatmap_seurat_fSOM_clusters.pdf")

#FlowSOM cluster 1 (= smallest cluster, Seurat cluster 25)is dropped because
#biologically unknown which cell type this cluster represents and features
#extreme expression of several unrelated markers (artefacts?)

#Percentage of cells in metacluster 1:
(nrow(dataset[dataset$clr_metaclusters==1,])/nrow(dataset))*100

#Discard metacluster 1:
dataset<-dataset%>%dplyr::filter(clr_metaclusters!=1)

cluster_data_exprs<-as.data.frame(dataset[,c(8:41)])

clr_function<-function(x){
  return(log1p(x=x/(exp(x=sum(log1p(x=x[x>0]),na.rm=TRUE)/length(x=x)))))
}

for(i in 1:ncol(cluster_data_exprs)){
  print(paste0("Mean is ",mean(cluster_data_exprs[,i])," and minimum is ",min(cluster_data_exprs[,i])))
  minimum<-min(cluster_data_exprs[,i])
  cluster_data_exprs[,i]<-cluster_data_exprs[,i]-minimum
  print(paste0("Mean is ",mean(cluster_data_exprs[,i])," and minimum is ",min(cluster_data_exprs[,i])))
  print(paste0("Done with column ",i))
  print(",")
}

cluster_data_exprs<-apply(cluster_data_exprs,MARGIN=2,FUN=clr_function)

cluster_data_meta<-as.data.frame(dataset$clr_metaclusters)
colnames(cluster_data_meta)<-"clr_metaclusters"
cluster_data_merged<-cbind(cluster_data_exprs,cluster_data_meta)
cluster_data_heatmap<-cluster_data_merged%>%dplyr::group_by(clr_metaclusters)%>%dplyr::summarise_all(median)
cluster_data_heatmap<-as.data.frame(cluster_data_heatmap)
rownames(cluster_data_heatmap)<-cluster_data_heatmap$clr_metaclusters
cluster_data_heatmap<-cluster_data_heatmap[,-1]
cluster_data_heatmap_core<-cluster_data_heatmap

#And create expression heatmap
cluster_data_heatmap_core<-cluster_data_heatmap_core[,c("CD3 BUV496","CD4 PerCP-Fire806","CD8 SparkBlue550","CD45RA SparkUV387",
                                                        "CCR7 BV605","FoxP3 AF488","HLA.DR PerCP","CD39 BB700","CD14 SparkBlue574",
                                                        "CD56 APC-Fire750","CD16 BV570")]

cluster_data_heatmap_core<-cluster_data_heatmap[c(40,44,#CD4 Tn
                                                  39,49,#CD4 Tcm
                                                  21,25,29,30,32,41,42,45,#CD4 Tem
                                                  22,26,27,31,#CD4 Tregs
                                                  37,#CD4CD8 DP
                                                  35,#CD8 Tn
                                                  36,#CD8 Tcm
                                                  23,24,28,33,#CD8 Tem
                                                  17,43,46,47,#CD8 Temra
                                                  34,38,48,#NKT-like cells
                                                  2,3,5,9,10,13,18,19,#NK cells
                                                  1,4,#B cells
                                                  7,8,12,14,#Monocytes
                                                  15,16,#Dendritic-like cells
                                                  6,11,20#Undefined
),
c("CD3 BUV496","CD4 PerCP-Fire806","CD8 SparkBlue550","CD45RA SparkUV387",
  "CCR7 BV605","FoxP3 AF488","HLA.DR PerCP","CD39 BB700","CD14 SparkBlue574",
  "CD56 APC-Fire750","CD16 BV570")]

rownames(cluster_data_heatmap_core)<-paste0("FlowSOM_",rownames(cluster_data_heatmap_core))

save_pheatmap_pdf<-function(x,filename,width=4,height=11){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename,width=width,height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#Figure S1A
plot2<-pheatmap(cluster_data_heatmap_core,
                scale="column",
                cluster_cols=F,
                cluster_rows=F,
                gaps_row=c(2,4,12,16,17,18,19,23,27,30,38,40,44,46),
                color=colorRampPalette(c("navy","white","darkred"))(25))

save_pheatmap_pdf(plot2,"Heatmap_key_markers_clr_metaclusters.pdf")

#Combine highly similar metaclusters as main_cluster

#backup <- dataset
dataset <- backup

main_clusters<-read_xlsx("main_clusters.xlsx")
dataset$clr_metaclusters<-as.numeric(dataset$clr_metaclusters)
dataset <- left_join(dataset,main_clusters,by=join_by("clr_metaclusters"=="clr_metaclusters"))

#colnames(dataset)<-sub(".*",",",colnames(dataset))
cluster_data_exprs<-as.data.frame(dataset[,c(8:41)])
cluster_data_meta<-as.data.frame(dataset$main_cluster)
colnames(cluster_data_meta)<-"main_cluster"
cluster_data_merged<-cbind(cluster_data_exprs,cluster_data_meta)
cluster_data_heatmap<-cluster_data_merged%>%group_by(main_cluster)%>%summarise_all(median)
cluster_data_heatmap<-as.data.frame(cluster_data_heatmap)
rownames(cluster_data_heatmap)<-cluster_data_heatmap$main_cluster
cluster_data_heatmap<-cluster_data_heatmap[,-1]
colnames(cluster_data_heatmap) <- sub(" .*","",colnames(cluster_data_heatmap))

save_pheatmap_pdf<-function(x,filename,width=6,height=4){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename,width=width,height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

annotation_row <- data.frame(cell_type = rownames(cluster_data_heatmap))
rownames(annotation_row) <- rownames(cluster_data_heatmap)

# Custom colors for annotation
annotation_colors <- list(Condition = c(A = "#1b9e77", B = "#d95f02"))

annot_col <- list(cell_type=c("CD4CD8_DP_T"="#05605E",
                              "CD4_Tn"="#072E67",
                              "CD4_Tcm"="#B5D9F1",
                              "CD4_Tem"="#1F6EB0",
                              "CD4_Treg"="#DB3395",
                              "CD8_Tn"="#FFDB7C",
                              "CD8_Tcm"="#FB4A22",
                              "CD8_Tem"="#FFC229",
                              "CD8_TEMRA"="#7D0024",
                              "NKT_like"="#4ADC61",
                              "NK_cells"="#156D24",
                              "B_cells"="#F5669E",
                              "Monocytes"="#770074",
                              "DC_like"="#34958C",
                              "Undefined"="grey"))

#Figure 1B
plot3<-pheatmap(cluster_data_heatmap[c("CD4CD8_DP_T","CD4_Tn","CD4_Tcm","CD4_Tem",
                                       "CD4_Treg","CD8_Tn","CD8_Tcm","CD8_Tem","CD8_TEMRA",
                                       "NKT_like","NK_cells","B_cells","Monocytes","DC_like",
                                       "Undefined"),
                                     c("CD3","CD4","CD8","CD45RA",
                                       "CCR7","FoxP3","HLA.DR","CD39","CD14",
                                       "CD56","CD16")],
                cluster_cols=F,
                cluster_rows=F,
                scale="column",
                annotation_row = annotation_row,
                annotation_colors = annot_col,
                gaps_row=c(1,5,9,11,12,14),
                color=colorRampPalette(c("navy","white","darkred"))(25))

save_pheatmap_pdf(plot3,"Heatmap_key_markers_main_clusters.pdf")

#Now recreate umap plot, with metacluster 1 discarded and relabeled with biological labels
cluster_data<-readRDS("mapped_seurat_object.RObj")
cluster_data_excl1<-subset(cluster_data,clr_metaclusters==1,invert=TRUE)

cluster_data_excl1@meta.data$main_cluster<-main_clusters$main_cluster[match(cluster_data_excl1@meta.data$clr_metaclusters,main_clusters$clr_metaclusters)]
cluster_data_excl1@meta.data$detail_cluster<-main_clusters$detail_cluster[match(cluster_data_excl1@meta.data$clr_metaclusters,main_clusters$clr_metaclusters)]

#Figure 1A
pdf("umap_main_clusters.pdf",width=8,height=5)
DimPlot(cluster_data_excl1,reduction="umap",raster=FALSE,label=FALSE,group.by="main_cluster",cols=c(
  "CD4CD8_DP_T"="#05605E",
  "CD4_Tn"="#072E67",
  "CD4_Tcm"="#B5D9F1",
  "CD4_Tem"="#1F6EB0",
  "CD4_Treg"="#DB3395",
  "CD8_Tn"="#FFDB7C",
  "CD8_Tcm"="#FB4A22",
  "CD8_Tem"="#FFC229",
  "CD8_TEMRA"="#7D0024",
  "NKT_like"="#4ADC61",
  "NK_cells"="#156D24",
  "B_cells"="#F5669E",
  "Monocytes"="#770074",
  "DC_like"="#34958C",
  "Undefined"="grey"))
dev.off()

#Bin cells in large cell types (CD4, CD8 and other)
CD4_cells<-c("CD4_Tn","CD4_Tcm","CD4_Tem","CD4_Treg","CD4CD8_DP_T")
CD8_cells<-c("CD8_Tn","CD8_Tcm","CD8_Tem","CD8_TEMRA","CD4CD8_DP_T")
dataset$large_cluster<-ifelse(dataset$main_cluster%in%CD4_cells,"CD4_cell",
                              ifelse(dataset$main_cluster%in%CD8_cells,"CD8_cell","Other"))

#Savedataset
#saveRDS(dataset,file="full_annotated_dataset.rds")

####RESTARTPOINT#####
setwd("~")
setwd(".../")
dataset<-readRDS("full_annotated_dataset.rds")

#Create count tables of main clusters
cell_counts<-dataset%>%group_by(subject,timepoint,main_cluster)%>%tally()
cell_counts<-pivot_wider(cell_counts,names_from="main_cluster",values_from="n")
cell_counts<-as.data.frame(cell_counts)
cell_counts$total<-rowSums(cell_counts[3:ncol(cell_counts)],na.rm=T)
cell_counts[is.na(cell_counts)]<-0

cell_counts_relative<-cell_counts
for(r in 1:nrow(cell_counts_relative)){
  for(c in 3:ncol(cell_counts_relative)){
    cell_counts_relative[r,c]<-cell_counts_relative[r,c]/cell_counts_relative[r,ncol(cell_counts_relative)]
  }
}

cell_counts_longitudinal<-cell_counts_relative
attach(cell_counts_longitudinal)
for(i in unique(subject)){
  if(length(timepoint[subject==i])>1){
    for(j in 2:length(timepoint[subject==i])){
      for(k in 3:ncol(cell_counts_longitudinal)){
        cell_counts_longitudinal[cell_counts_longitudinal$timepoint==j&cell_counts_longitudinal$subject==i,k]<-cell_counts_longitudinal[cell_counts_longitudinal$timepoint==j&cell_counts_longitudinal$subject==i,k]/cell_counts_longitudinal[cell_counts_longitudinal$timepoint==1&cell_counts_longitudinal$subject==i,k]
      }
    }
  }
}

for(i in unique(subject)){
  for(k in 3:ncol(cell_counts_longitudinal)){
    cell_counts_longitudinal[cell_counts_longitudinal$timepoint==1&cell_counts_longitudinal$subject==i,k]<-1
  }
}

#Load meta data and Fortessa data
fortessa<-read_xlsx("Data_incl_metadata_all_batches_relative.xlsx")
fortessa<-fortessa%>%dplyr::rename(
  subject=Sample,
  timepoint=Timepoint
)
newdata<-merge(cell_counts_relative,fortessa,by=c("subject","timepoint"),all=T)
newdata$main_group<-ifelse(newdata$Treatment=="anti-PD-1","anti-PD-(L)1",
                           ifelse(newdata$Treatment=="anti-PD-L1","anti-PD-(L)1",
                                  ifelse(newdata$Treatment=="ipi-nivo","ipi-nivo",
                                         ifelse(newdata$Treatment=="ipi_mono","ipi_mono",
                                                "healthy"))))

#Check baseline abundance of main cell populations by response status(main outcome)
baseline_data<-newdata %>% dplyr::filter(timepoint==1 & Response != "HD")
baseline_long<-baseline_data %>% dplyr::select(c("subject","B_cells","CD4_Tcm","CD4_Tem","CD4_Tn",
                                                 "CD4_Treg","CD8_Tcm","CD8_Tem","CD8_TEMRA","CD8_Tn",
                                                 "CD4CD8_DP_T","Monocytes","DC_like","NK_cells","NKT_like","Undefined",
                                                 "Response","main_group"))%>%
  pivot_longer(!c(subject,Response,main_group),names_to="main_cluster",values_to="fraction")

baseline_long<-baseline_long%>%arrange(main_cluster,fraction)

baseline_long$Response <- factor(baseline_long$Response,levels=c("R","NR"))

plot2<-ggplot(baseline_long,aes(factor(main_cluster),fraction*100,fill=factor(Response)))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(pch=21,position=position_jitterdodge())+
  scale_fill_manual(labels=c("Responder","Non-responder"),
                    values=c("#5ec962","#440154"))+
  labs(title="Baseline cell abundance",fill=",")+
  ylab("% of all PBMCs")+
  xlab(",")+
  theme_classic()+
  #geom_signif(annotations=rep("ns",15),
  #            y_position=rep(50,15),
  #            xmin=c(0.8,1.8,2.8,3.8,4.8,5.8,6.8,7.8,8.8,9.8,10.8,11.8,12.8,13.8,14.8),
  #            xmax=c(1.2,2.2,3.2,4.2,5.2,6.2,7.2,8.2,9.2,10.2,11.2,12.2,13.2,14.2,15.2),
  #            tip_length=c(0,0))+
  ylim(0,45)+
  theme(axis.text.x=element_text(angle=45,vjust=1,
                                 hjust=1))
#Figure 1C
pdf("Baseline_abundance_by_response.pdf",width=9,height=4)
plot(plot2)
dev.off()

test_table<-data.frame(Group=NA,Pvalue=NA)
j<-1
for(i in unique(baseline_long$main_cluster)){
  test_set<-baseline_long%>%dplyr::filter(main_cluster==i)
  test_table[j,1]<-paste0(i)
  result<-kruskal.test(test_set$fraction,test_set$Response)
  test_table[j,2]<-paste0(as.numeric(result$p.value))
  j<-j+1
}

#Check baseline abundance of main cell populations by treatment group
baseline_data<-newdata%>%dplyr::filter(timepoint==1)
baseline_data$main_group<-ifelse(baseline_data$main_group=="ipi_mono"|baseline_data$main_group=="ipi-nivo","ipi+/-nivo",baseline_data$main_group)
baseline_data$main_group<-factor(baseline_data$main_group,levels=c('healthy','anti-PD-(L)1','ipi+/-nivo'))

baseline_long<-baseline_data %>% dplyr::select(c("subject","B_cells","CD4_Tcm","CD4_Tem","CD4_Tn",
                                                 "CD4_Treg","CD8_Tcm","CD8_Tem","CD8_TEMRA","CD8_Tn",
                                                 "CD4CD8_DP_T","Monocytes","DC_like","NK_cells","NKT_like","Undefined",
                                                 "Response","main_group","Colitis"))%>%
  pivot_longer(!c(subject,Response,main_group,Colitis),names_to="main_cluster",values_to="fraction")

baseline_long<-baseline_long%>%arrange(main_cluster,fraction)

plot3<-ggplot(baseline_long,aes(factor(main_cluster),fraction*100,fill=factor(main_group)))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(pch=21,position=position_jitterdodge())+
  scale_fill_manual(labels=c("Healthy donor","Anti-PD-(L)1 mono","Ipilimumab +/- nivolumab"),
                    values=c("#FCFFA4","#F98E09","#BC3754"))+
  labs(title="Baseline cell abundance",fill="")+
  ylab("% of all PBMCs")+
  xlab("")+
  theme_classic()+
  geom_signif(annotations=paste0(expression(italic(P) == 0.067)),parse=T,
              y_position=45,
              xmin=2.8,
              xmax=3.2,
              tip_length=c(0,0))+
  theme(axis.text.x=element_text(angle=45,vjust=1,
                                 hjust=1))+
  ylim(0,50)
#Figure S1B
pdf("Baseline_abundance_by_treatment.pdf",width=10,height=3.5)
plot(plot3)
dev.off()

test_table<-data.frame(Group=NA,Pvalue=NA)
j<-1
for(i in unique(baseline_long$main_cluster)){
  test_set<-baseline_long%>%dplyr::filter(main_cluster==i)
  test_table[j,1]<-paste0(i)
  result<-kruskal.test(test_set$fraction,test_set$main_group)
  test_table[j,2]<-paste0(as.numeric(result$p.value))
  j<-j+1
}

#By response adjusted for ICI treatment
baseline_long_adj <- baseline_long %>% dplyr::filter(Response!="HD")
df_adjusted <- data.frame(matrix(ncol=4))
colnames(df_adjusted) <- c("cell_type","coef_response","ci_lower","ci_upper")
r <- 1
for(c in unique(baseline_long_adj$main_cluster)){
  data <- baseline_long_adj %>% dplyr::filter(main_cluster==c)
  data$Response <- factor(data$Response,levels=c("R","NR"))
  model <- lm(fraction~Response+main_group,data=data)
  df_adjusted[r,1] <- c
  df_adjusted[r,2] <- as.numeric(coefficients(model)[2])
  df_adjusted[r,3] <- as.data.frame(confint(model))[2,1]
  df_adjusted[r,4] <- as.data.frame(confint(model))[2,2]
  r <- r +1
}

summary(model)

df_adjusted <- df_adjusted %>% 
  arrange(coef_response) %>% 
  mutate(cell_type = factor(cell_type, levels = cell_type))
df_adjusted$direction <- ifelse(df_adjusted$coef_response>0,"NR","R")

adjusted_plot <- ggplot(df_adjusted, aes(x = cell_type, y = 100*coef_response, fill = direction)) +
  geom_hline(yintercept=0,color="grey",linetype="dashed")+
  geom_bar(stat = "identity", position = "dodge", width = 0.7,alpha=0.6) +
  geom_errorbar(aes(ymin = 100*ci_lower, ymax = 100*ci_upper), width = 0.2) +
  scale_fill_manual(labels=c("Non-responder","Responder"),values = c("NR" = "#440154", "R" = "#5ec962")) +
  coord_flip() +
  labs(x = "", y = "percentage point diff. in baseline abundance", fill = "Adjusted for ICI-\nregimen, up in:") +
  theme_classic()
#Figure 1D
pdf("Adjusted_plot_main_group_response.pdf",width=5,height=6)
plot(adjusted_plot)
dev.off()

# Functional analyses / T cell subsets
functional_data<-dataset%>%dplyr::filter(large_cluster!="Other")
functional_markers<-c("TNFa","GzmB","IFNy","IL5","IL13","IL17A","IL10","Ki67","CD134","HLA.DR","CD38","CD25","CD69","PD1","TIM3","LAG3","TCF1","CD57","IntB1","IntB7","CXCR3","CCR4")

#Calculate which functional markers are statistically up in subclusters
calc_data <- functional_data[,c(functional_markers,"detail_cluster")]

#Subsample to 10% to increase speed of comparisons and scale data.
sample_size <- round(nrow(calc_data) * 0.01)
calc_data <- calc_data[sample(nrow(calc_data), sample_size), ]
colnames(calc_data)[ncol(calc_data)] <- "detail_cluster"
calc_data[, c(1:(ncol(calc_data)-1))] <- lapply(calc_data[, c(1:(ncol(calc_data)-1))], as.numeric)

subsets <- unique(calc_data[, ncol(calc_data)])
results <- list()
FC_list <- list()
for (subset in subsets) {
  subset_data <- calc_data[calc_data[, ncol(calc_data)] == subset, ]
  other_data <- calc_data[calc_data[, ncol(calc_data)] != subset, ]
  
  p_values <- numeric(ncol(calc_data)-1)
  fold_change <- numeric(ncol(calc_data)-1)
  
  for (i in 1:(ncol(calc_data)-1)) {
    test <- wilcox.test(subset_data[, i], other_data[, i])
    p_values[i] <- test$p.value
    minimum <- abs(min(min(subset_data[, i]),min(other_data[, i])))
    FC <- median(subset_data[, i]+minimum)/median(other_data[, i]+minimum)
    fold_change[i] <- FC
  }
  
  results[[subset]] <- p_values
  FC_list[[subset]] <- fold_change
  message(paste0("Comparisons for ",subset," are ready."))
}
results_df_heatmap <- do.call(rbind, results)
results_df <- results_df_heatmap
FC_df_heatmap <- do.call(rbind, FC_list)
FC_df <- FC_df_heatmap
colnames(results_df) <- colnames(calc_data)[1:(ncol(calc_data)-1)]
colnames(FC_df) <- colnames(calc_data)[1:(ncol(calc_data)-1)]
results_df <- cbind(Subset = rownames(results_df), results_df)
results_df <- as.data.frame(results_df)
results_df <- melt(results_df,id.vars="Subset")
results_df$p.adj <- p.adjust(results_df$value,method="BH")
colnames(results_df) <- c("subset","protein","p.val","p.adj")
FC_df <- cbind(Subset = rownames(FC_df), FC_df)
FC_df <- as.data.frame(FC_df)
FC_df <- melt(FC_df,id.vars="Subset")
colnames(FC_df) <- c("subset","protein","FoldChange")
results_df <- merge(results_df,FC_df,by=c("subset","protein"))
results_df <- results_df[order(results_df$subset, results_df$p.adj), ]
results_df_sign <- results_df %>% dplyr::filter(p.adj<0.05)
results_df_up <- results_df_sign %>% dplyr::filter(FoldChange>1)
results_df_dn <- results_df_sign %>% dplyr::filter(FoldChange<1)

#Create heatmap for expression
heatmap_data<-functional_data[,c(functional_markers,"detail_cluster")]
heatmap_data<-heatmap_data%>%group_by(detail_cluster)%>%dplyr::summarise_all(median)
heatmap_data<-as.data.frame(heatmap_data)
rownames(heatmap_data)<-heatmap_data$detail_cluster
heatmap_data<-heatmap_data[,-1]
heatmap_data<-as.data.frame(t(heatmap_data))

heatmap_data <- heatmap_data[,c("CD4_Tn_1","CD4_Tn_2","CD4_Tcm_1","CD4_Tcm_2","CD4_Tem_1",
                                    "CD4_Tem_2","CD4_Tem_3","CD4_Tem_4","CD4_Tem_5","CD4_Tem_6",
                                    "CD4_Tem_7","CD4_Tem_8","CD4_Treg_1","CD4_Treg_2","CD4_Treg_3",
                                    "CD4_Treg_4","CD4CD8_DP_T","CD8_Tn_1","CD8_Tcm","CD8_Tem_1",
                                    "CD8_Tem_2","CD8_Tem_3","CD8_Tem_4","CD8_TEMRA_1","CD8_TEMRA_2",
                                    "CD8_TEMRA_3","CD8_TEMRA_4")]

heatmap_data_CD4 <- heatmap_data[,c("CD4_Tn_1","CD4_Tn_2","CD4_Tcm_1","CD4_Tcm_2","CD4_Tem_1",
                              "CD4_Tem_2","CD4_Tem_3","CD4_Tem_4","CD4_Tem_5","CD4_Tem_6",
                              "CD4_Tem_7","CD4_Tem_8","CD4_Treg_1","CD4_Treg_2","CD4_Treg_3",
                              "CD4_Treg_4","CD4CD8_DP_T")]
heatmap_data_CD8 <- heatmap_data[,c("CD4CD8_DP_T","CD8_Tn_1","CD8_Tcm","CD8_Tem_1",
                              "CD8_Tem_2","CD8_Tem_3","CD8_Tem_4","CD8_TEMRA_1","CD8_TEMRA_2",
                              "CD8_TEMRA_3","CD8_TEMRA_4")]

results_df_heatmap <- as.data.frame(results_df_heatmap)
colnames(results_df_heatmap) <- colnames(calc_data)[1:(ncol(calc_data)-1)]
results_df_heatmap <- as.data.frame(t(results_df_heatmap))
results_df_heatmap <- results_df_heatmap[,colnames(heatmap_data)]
results_df_heatmap <- results_df_heatmap[rownames(heatmap_data),]
significance_symbols <- apply(results_df_heatmap, 2, function(x) {
  ifelse(x < 0.05, "*", "")
})
FC_df_heatmap <- as.data.frame(FC_df_heatmap)
colnames(FC_df_heatmap) <- colnames(calc_data)[1:(ncol(calc_data)-1)]
FC_df_heatmap <- as.data.frame(t(FC_df_heatmap))
significance_symbols <- as.data.frame(significance_symbols)
for(c in colnames(significance_symbols)){
  for(r in rownames(significance_symbols)){
    significance_symbols[r,c] <- ifelse(FC_df_heatmap[r,c]>1,significance_symbols[r,c],"")
  }
}

#Figure S1C
plot_CD4<-pheatmap(heatmap_data[,c(1:17)],
                display_numbers = significance_symbols[,c(1:17)],
                scale="row",
                cluster_cols=F,
                cluster_rows=F,
                gaps_col=c(2,4,12,16),
                gaps_row=c(3,5,6,7,13,18),
                color=colorRampPalette(c("navy","white","darkred"))(35),
                angle_col = 45,
                cellwidth=12,
                cellheight=11)

plot_CD8<-pheatmap(heatmap_data[,c(17:27)],
                   display_numbers = significance_symbols[,c(17:27)],
                   scale="row",
                   cluster_cols=F,
                   cluster_rows=F,
                   gaps_col=c(1,2,3,7),
                   gaps_row=c(3,5,6,7,13,18),
                   color=colorRampPalette(c("navy","white","darkred"))(35),
                   angle_col = 45,
                   cellwidth=12,
                   cellheight=11)

save_pheatmap_pdf<-function(x,filename,width=5,height=4.8){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename,width=width,height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(plot_CD4,"Heatmap_functional_CD4.pdf")

save_pheatmap_pdf<-function(x,filename,width=4.3,height=4.8){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename,width=width,height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(plot_CD8,"Heatmap_functional_CD8.pdf")

#Check baseline abundance of detailed CD4 T cell populations by response status (main outcome)
cell_counts<-dataset%>%dplyr::filter(large_cluster=="CD4_cell")%>%group_by(subject,timepoint,detail_cluster)%>%tally()
cell_counts<-pivot_wider(cell_counts,names_from="detail_cluster",values_from="n")
cell_counts<-as.data.frame(cell_counts)
cell_counts$total<-rowSums(cell_counts[3:ncol(cell_counts)],na.rm=T)
cell_counts[is.na(cell_counts)]<-0

cell_counts_relative<-cell_counts
for(r in 1:nrow(cell_counts_relative)){
  for(c in 3:ncol(cell_counts_relative)){
    cell_counts_relative[r,c]<-cell_counts_relative[r,c]/cell_counts_relative[r,ncol(cell_counts_relative)]
  }
}

#Load metadata and Fortessa data
fortessa<-read_xlsx("Data_incl_metadata_all_batches_relative.xlsx")
fortessa<-fortessa%>%dplyr::rename(
  subject=Sample,
  timepoint=Timepoint
)
newdata<-merge(cell_counts_relative,fortessa,by=c("subject","timepoint"),all=T)
newdata$main_group<-ifelse(newdata$Treatment=="anti-PD-1","anti-PD-(L)1",
                           ifelse(newdata$Treatment=="anti-PD-L1","anti-PD-(L)1",
                                  ifelse(newdata$Treatment=="ipi-nivo","ipi-nivo",
                                         ifelse(newdata$Treatment=="ipi_mono","ipi_mono",
                                                "healthy"))))

baseline_data<-newdata%>%dplyr::filter(timepoint==1 & Response != "HD")
baseline_long<-baseline_data%>%dplyr::select(c("subject","CD4_Tn_1","CD4_Tn_2","CD4_Tcm_1","CD4_Tcm_2","CD4_Tem_1",
                                               "CD4_Tem_2","CD4_Tem_3","CD4_Tem_4","CD4_Tem_5","CD4_Tem_6",
                                               "CD4_Tem_7","CD4_Tem_8","CD4_Treg_1","CD4_Treg_2","CD4_Treg_3",
                                               "CD4_Treg_4","CD4CD8_DP_T","Response","main_group"))%>%
  pivot_longer(!c(subject,Response,main_group),names_to="main_cluster",values_to="fraction")

baseline_long<-baseline_long%>%arrange(main_cluster,fraction)

baseline_long$Response <- factor(baseline_long$Response,levels=c("R","NR"))
plot2<-ggplot(baseline_long,aes(factor(main_cluster),log10(fraction*100),fill=factor(Response)))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(pch=21,position=position_jitterdodge())+
  scale_fill_manual(labels=c("Responder","Non-responder"),
                    values=c("#5ec962","#440154"))+
  labs(title="Baseline cell abundance",fill=",")+
  ylab("% of all CD4+ T cells")+
  annotation_logticks(sides="l")+
  theme(axis.ticks.x=element_blank())+
  scale_y_continuous(breaks=c(-2,-1,0,1,2),labels=c(0.010,0.10,1.0,10,100))+
  xlab("")+
  theme_classic()+
  geom_signif(annotations=paste0(expression(italic(P) == 0.087)),parse=T,
              y_position=1.6,
              xmin=9.8,
              xmax=10.2,
              tip_length=c(0,0))+
  theme(axis.text.x=element_text(angle=45,vjust=1,
                                 hjust=1))
#Figure S1D
pdf("Baseline_abundance_CD4_by_response.pdf",width=12,height=4)
plot(plot2)
dev.off()

test_table<-data.frame(Group=NA,Pvalue=NA)
j<-1
for(i in unique(baseline_long$main_cluster)){
  test_set<-baseline_long%>%dplyr::filter(main_cluster==i)
  test_table[j,1]<-paste0(i)
  result<-kruskal.test(test_set$fraction,test_set$Response)
  test_table[j,2]<-paste0(as.numeric(result$p.value))
  j<-j+1
}

#Check baseline abundance of detailed CD4 T cell populations by treatment
baseline_data<-newdata%>%dplyr::filter(timepoint==1)
baseline_long<-baseline_data%>%dplyr::select(c("subject","CD4_Tn_1","CD4_Tn_2","CD4_Tcm_1","CD4_Tcm_2","CD4_Tem_1",
                                               "CD4_Tem_2","CD4_Tem_3","CD4_Tem_4","CD4_Tem_5","CD4_Tem_6",
                                               "CD4_Tem_7","CD4_Tem_8","CD4_Treg_1","CD4_Treg_2","CD4_Treg_3",
                                               "CD4_Treg_4","CD4CD8_DP_T","Response","main_group"))%>%
  pivot_longer(!c(subject,Response,main_group),names_to="main_cluster",values_to="fraction")

baseline_long<-baseline_long%>%arrange(main_cluster,fraction)
baseline_long$main_group<-ifelse(baseline_long$main_group=="ipi_mono"|baseline_long$main_group=="ipi-nivo","ipi+/-nivo",baseline_long$main_group)
baseline_long$main_group<-factor(baseline_long$main_group,levels=c("healthy","anti-PD-(L)1","ipi+/-nivo"))
plot3<-ggplot(baseline_long,aes(factor(main_cluster),log10(fraction*100),fill=factor(main_group)))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(pch=21,position=position_jitterdodge())+
  scale_fill_manual(labels=c("Healthy donor","Anti-PD-(L)1 mono","Ipilimumab +/- nivolumab"),
                    values=c("#FCFFA4","#F98E09","#BC3754"))+
  labs(title="Baseline cell abundance",fill="")+
  ylab("% of all CD4+ T cells")+
  annotation_logticks(sides="l")+
  theme(axis.ticks.x=element_blank())+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),labels=c(0.010,0.10,1.0,10,100,1000))+
  xlab("")+
  theme_classic()+
  geom_signif(annotations=paste0(expression(italic(P) == "6.1e-4")),parse=T,
              y_position=1.6,
              xmin=4.8,
              xmax=5.2,
              tip_length=c(0,0))+
  geom_signif(annotations=paste0(expression(italic(P) == "1.4e-3")),parse=T,
              y_position=1.6,
              xmin=6.8,
              xmax=7.2,
              tip_length=c(0,0))+
  geom_signif(annotations=paste0(expression(italic(P) == 0.063)),parse=T,
              y_position=1.6,
              xmin=14.8,
              xmax=15.2,
              tip_length=c(0,0))+
  theme(axis.text.x=element_text(angle=45,vjust=1,
                                 hjust=1))
#Figure S1E
pdf("Baseline_abundance_CD4_by_treatment.pdf",width=13,height=4)
plot(plot3)
dev.off()

test_table<-data.frame(Group=NA,Pvalue=NA,Padj=NA)
j<-1
for(i in unique(baseline_long$main_cluster)){
  test_set<-baseline_long%>%dplyr::filter(main_cluster==i)
  test_table[j,1]<-paste0(i)
  result<-kruskal.test(test_set$fraction,test_set$main_group)
  test_table[j,2]<-paste0(as.numeric(result$p.value))
  j<-j+1
}
test_table$Padj <- p.adjust(test_table$Pvalue,method = "BH",n=nrow(test_table))

#By response adjusted for ICI treatment
baseline_long_adj <- baseline_long %>% dplyr::filter(Response!="HD")
df_adjusted <- data.frame(matrix(ncol=4))
colnames(df_adjusted) <- c("cell_type","coef_response","ci_lower","ci_upper")
r <- 1
for(c in unique(baseline_long_adj$main_cluster)){
  data <- baseline_long_adj %>% dplyr::filter(main_cluster==c)
  data$Response <- factor(data$Response,levels=c("R","NR"))
  model <- lm(fraction~Response+main_group,data=data)
  df_adjusted[r,1] <- c
  df_adjusted[r,2] <- as.numeric(coefficients(model)[2])
  df_adjusted[r,3] <- as.data.frame(confint(model))[2,1]
  df_adjusted[r,4] <- as.data.frame(confint(model))[2,2]
  r <- r +1
}

df_adjusted <- df_adjusted %>% 
  arrange(coef_response) %>% 
  mutate(cell_type = factor(cell_type, levels = cell_type))
df_adjusted$direction <- ifelse(df_adjusted$coef_response>0,"NR","R")

adjusted_plot <- ggplot(df_adjusted, aes(x = cell_type, y = 100*coef_response, fill = direction)) +
  geom_hline(yintercept=0,color="grey",linetype="dashed")+
  geom_bar(stat = "identity", position = "dodge", width = 0.7,alpha=0.6) +
  geom_errorbar(aes(ymin = 100*ci_lower, ymax = 100*ci_upper), width = 0.2) +
  scale_fill_manual(labels=c("Non-responder","Responder"),values = c("NR" = "#440154", "R" = "#5ec962")) +
  coord_flip() +
  labs(x = "", y = "percentage point diff. in baseline abundance", fill = "Adjusted for ICI-\nregimen, up in:") +
  theme_classic()
#Figure S1F
pdf("Adjusted_plot_CD4_subsets_response.pdf",width=5,height=6)
plot(adjusted_plot)
dev.off()

#By treatment adjusted for response
baseline_long_adj <- baseline_long %>% dplyr::filter(Response!="HD")
df_adjusted <- data.frame(matrix(ncol=4))
colnames(df_adjusted) <- c("cell_type","coef_treatment","ci_lower","ci_upper")
r <- 1
for(c in unique(baseline_long_adj$main_cluster)){
  data <- baseline_long_adj %>% dplyr::filter(main_cluster==c)
  model <- lm(fraction~main_group+Response,data=data)
  df_adjusted[r,1] <- c
  df_adjusted[r,2] <- as.numeric(coefficients(model)[2])
  df_adjusted[r,3] <- as.data.frame(confint(model))[2,1]
  df_adjusted[r,4] <- as.data.frame(confint(model))[2,2]
  r <- r +1
}

df_adjusted <- df_adjusted %>% 
  arrange(coef_treatment) %>% 
  mutate(cell_type = factor(cell_type, levels = cell_type))
df_adjusted$direction <- ifelse(df_adjusted$coef_treatment>0,"ipi+/-nivo","anti-PD-(L)1")

adjusted_plot <- ggplot(df_adjusted, aes(x = cell_type, y = 100*coef_treatment, fill = direction)) +
  geom_hline(yintercept=0,color="grey",linetype="dashed")+
  geom_bar(stat = "identity", position = "dodge", width = 0.7,alpha=0.6) +
  geom_errorbar(aes(ymin = 100*ci_lower, ymax = 100*ci_upper), width = 0.2) +
  scale_fill_manual(labels=c("Anti-PD-(L)1 mono","Ipilimumab+/-nivolumab"),values = c("ipi+/-nivo" = "#BC3754", "anti-PD-(L)1" = "#F98E09")) +
  coord_flip() +
  labs(x = "", y = "percentage point diff. in baseline abundance", fill = "Adjusted for\nresponse, up in:") +
  theme_classic()
#Figure S1G
pdf("Adjusted_plot_CD4_subsets_treatment.pdf",width=5,height=6)
plot(adjusted_plot)
dev.off()

#Same, now check baseline abundance of detailed CD8 T cell populations by response status (main outcome)
cell_counts<-dataset%>%dplyr::filter(large_cluster=="CD8_cell")%>%group_by(subject,timepoint,detail_cluster)%>%tally()
cell_counts<-pivot_wider(cell_counts,names_from="detail_cluster",values_from="n")
cell_counts<-as.data.frame(cell_counts)
cell_counts$total<-rowSums(cell_counts[3:ncol(cell_counts)],na.rm=T)
cell_counts[is.na(cell_counts)]<-0

cell_counts_relative<-cell_counts
for(r in 1:nrow(cell_counts_relative)){
  for(c in 3:ncol(cell_counts_relative)){
    cell_counts_relative[r,c]<-cell_counts_relative[r,c]/cell_counts_relative[r,ncol(cell_counts_relative)]
  }
}

#Load metadata and Fortessa data
fortessa<-read_xlsx("Data_incl_metadata_all_batches_relative.xlsx")
fortessa<-fortessa%>%dplyr::rename(
  subject=Sample,
  timepoint=Timepoint
)
newdata<-merge(cell_counts_relative,fortessa,by=c("subject","timepoint"),all=T)
newdata$main_group<-ifelse(newdata$Treatment=="anti-PD-1","anti-PD-(L)1",
                           ifelse(newdata$Treatment=="anti-PD-L1","anti-PD-(L)1",
                                  ifelse(newdata$Treatment=="ipi-nivo","ipi-nivo",
                                         ifelse(newdata$Treatment=="ipi_mono","ipi_mono",
                                                "healthy"))))

baseline_data<-newdata%>%dplyr::filter(timepoint==1 & Response != "HD")
baseline_long<-baseline_data%>%dplyr::select(c("subject","CD8_TEMRA_1",                      
                                               "CD8_TEMRA_2","CD8_TEMRA_3","CD8_TEMRA_4",                      
                                               "CD8_Tcm","CD8_Tn_1","CD8_Tem_1","CD8_Tem_2",                        
                                               "CD8_Tem_3","CD8_Tem_4",
                                               "Response","main_group"))%>%
  pivot_longer(!c(subject,Response,main_group),names_to="main_cluster",values_to="fraction")

baseline_long<-baseline_long%>%arrange(main_cluster,fraction)

baseline_long$Response <- factor(baseline_long$Response,levels=c("R","NR"))
plot2<-ggplot(baseline_long,aes(factor(main_cluster),log10(fraction*100),fill=factor(Response)))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(pch=21,position=position_jitterdodge())+
  scale_fill_manual(labels=c("Responder","Non-responder"),
                    values=c("#5ec962","#440154"))+
  labs(title="Baseline cell abundance",fill=",")+
  ylab("% of all CD8+ T cells")+
  annotation_logticks(sides="l")+
  theme(axis.ticks.x=element_blank())+
  scale_y_continuous(breaks=c(-2,-1,0,1,2),labels=c(0.010,0.10,1.0,10,100))+
  xlab("")+
  theme_classic()+
  geom_signif(annotations=paste0(expression(italic(P) == "0.060")),parse=T,
              y_position=1.3,
              xmin=1.8,
              xmax=2.2,
              tip_length=c(0,0))+
  geom_signif(annotations=paste0(expression(italic(P) == 0.068)),parse=T,
              y_position=1.3,
              xmin=7.8,
              xmax=8.2,
              tip_length=c(0,0))+
  geom_signif(annotations=paste0(expression(italic(P) == 0.068)),parse=T,
              y_position=1.3,
              xmin=3.8,
              xmax=4.2,
              tip_length=c(0,0))+
  theme(axis.text.x=element_text(angle=45,vjust=1,
                                 hjust=1))
#Figure S1D
pdf("Baseline_abundance_CD8_by_response.pdf",width=10,height=4)
plot(plot2)
dev.off()

test_table<-data.frame(Group=NA,Pvalue=NA)
j<-1
for(i in unique(baseline_long$main_cluster)){
  test_set<-baseline_long%>%dplyr::filter(main_cluster==i)
  test_table[j,1]<-paste0(i)
  result<-kruskal.test(test_set$fraction,test_set$Response)
  test_table[j,2]<-paste0(as.numeric(result$p.value))
  j<-j+1
}

#Check baseline abundance of detailed CD8 T cell populations by treatment
baseline_data<-newdata%>%dplyr::filter(timepoint==1)
baseline_long<-baseline_data%>%dplyr::select(c("subject","CD8_TEMRA_1",                      
                                               "CD8_TEMRA_2","CD8_TEMRA_3","CD8_TEMRA_4",                      
                                               "CD8_Tcm","CD8_Tn_1","CD8_Tem_1","CD8_Tem_2",                        
                                               "CD8_Tem_3","CD8_Tem_4",
                                               "Response","main_group"))%>%
  pivot_longer(!c(subject,Response,main_group),names_to="main_cluster",values_to="fraction")

baseline_long<-baseline_long%>%arrange(main_cluster,fraction)
baseline_long$main_group<-ifelse(baseline_long$main_group=="ipi_mono"|baseline_long$main_group=="ipi-nivo","ipi+/-nivo",baseline_long$main_group)
baseline_long$main_group<-factor(baseline_long$main_group,levels=c("healthy","anti-PD-(L)1","ipi+/-nivo"))
plot3<-ggplot(baseline_long,aes(factor(main_cluster),log10(fraction*100),fill=factor(main_group)))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(pch=21,position=position_jitterdodge())+
  scale_fill_manual(labels=c("Healthy donor","Anti-PD-(L)1 mono","Ipilimumab +/- nivolumab"),
                    values=c("#FCFFA4","#F98E09","#BC3754"))+
  labs(title="Baseline cell abundance",fill="")+
  ylab("% of all CD8+ T cells")+
  annotation_logticks(sides="l")+
  theme(axis.ticks.x=element_blank())+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),labels=c(0.010,0.10,1.0,10,100,1000))+
  xlab("")+
  theme_classic()+
  geom_signif(annotations=paste0(expression(italic(P) == "9.4e-3")),parse=T,
              y_position=1.6,
              xmin=3.8,
              xmax=4.2,
              tip_length=c(0,0))+
  geom_signif(annotations=paste0(expression(italic(P) == "0.093")),parse=T,
              y_position=1.6,
              xmin=6.8,
              xmax=7.2,
              tip_length=c(0,0))+
  theme(axis.text.x=element_text(angle=45,vjust=1,
                                 hjust=1))
#Figure S1E
pdf("Baseline_abundance_CD8_by_treatment.pdf",width=10,height=4)
plot(plot3)
dev.off()

test_table<-data.frame(Group=NA,Pvalue=NA,Padj=NA)
j<-1
for(i in unique(baseline_long$main_cluster)){
  test_set<-baseline_long%>%dplyr::filter(main_cluster==i)
  test_table[j,1]<-paste0(i)
  result<-kruskal.test(test_set$fraction,test_set$main_group)
  test_table[j,2]<-paste0(as.numeric(result$p.value))
  j<-j+1
}
test_table$Padj <- p.adjust(test_table$Pvalue,method = "BH",n=nrow(test_table))

#By response adjusted for ICI treatment
baseline_long_adj <- baseline_long %>% dplyr::filter(Response!="HD")
df_adjusted <- data.frame(matrix(ncol=4))
colnames(df_adjusted) <- c("cell_type","coef_response","ci_lower","ci_upper")
r <- 1
for(c in unique(baseline_long_adj$main_cluster)){
  data <- baseline_long_adj %>% dplyr::filter(main_cluster==c)
  data$Response <- factor(data$Response,levels=c("R","NR"))
  model <- lm(fraction~Response+main_group,data=data)
  df_adjusted[r,1] <- c
  df_adjusted[r,2] <- as.numeric(coefficients(model)[2])
  df_adjusted[r,3] <- as.data.frame(confint(model))[2,1]
  df_adjusted[r,4] <- as.data.frame(confint(model))[2,2]
  r <- r +1
}

df_adjusted <- df_adjusted %>% 
  arrange(coef_response) %>% 
  mutate(cell_type = factor(cell_type, levels = cell_type))
df_adjusted$direction <- ifelse(df_adjusted$coef_response>0,"NR","R")

adjusted_plot <- ggplot(df_adjusted, aes(x = cell_type, y = 100*coef_response, fill = direction)) +
  geom_hline(yintercept=0,color="grey",linetype="dashed")+
  geom_bar(stat = "identity", position = "dodge", width = 0.7,alpha=0.6) +
  geom_errorbar(aes(ymin = 100*ci_lower, ymax = 100*ci_upper), width = 0.2) +
  scale_fill_manual(labels=c("Non-responder","Responder"),values = c("NR" = "#440154", "R" = "#5ec962")) +
  coord_flip() +
  labs(x = "", y = "percentage point diff. in baseline abundance", fill = "Adjusted for ICI-\nregimen, up in:") +
  theme_classic()
#Figure
pdf("Adjusted_plot_CD8_subsets_response.pdf",width=5,height=6)
plot(adjusted_plot)
dev.off()

#By treatment adjusted for response
baseline_long_adj <- baseline_long %>% dplyr::filter(Response!="HD")
df_adjusted <- data.frame(matrix(ncol=4))
colnames(df_adjusted) <- c("cell_type","coef_treatment","ci_lower","ci_upper")
r <- 1
for(c in unique(baseline_long_adj$main_cluster)){
  data <- baseline_long_adj %>% dplyr::filter(main_cluster==c)
  model <- lm(fraction~main_group+Response,data=data)
  df_adjusted[r,1] <- c
  df_adjusted[r,2] <- as.numeric(coefficients(model)[2])
  df_adjusted[r,3] <- as.data.frame(confint(model))[2,1]
  df_adjusted[r,4] <- as.data.frame(confint(model))[2,2]
  r <- r +1
}

df_adjusted <- df_adjusted %>% 
  arrange(coef_treatment) %>% 
  mutate(cell_type = factor(cell_type, levels = cell_type))
df_adjusted$direction <- ifelse(df_adjusted$coef_treatment>0,"ipi+/-nivo","anti-PD-(L)1")

adjusted_plot <- ggplot(df_adjusted, aes(x = cell_type, y = 100*coef_treatment, fill = direction)) +
  geom_hline(yintercept=0,color="grey",linetype="dashed")+
  geom_bar(stat = "identity", position = "dodge", width = 0.7,alpha=0.6) +
  geom_errorbar(aes(ymin = 100*ci_lower, ymax = 100*ci_upper), width = 0.2) +
  scale_fill_manual(labels=c("Anti-PD-(L)1 mono","Ipilimumab+/-nivolumab"),values = c("ipi+/-nivo" = "#BC3754", "anti-PD-(L)1" = "#F98E09")) +
  coord_flip() +
  labs(x = "", y = "percentage point diff. in baseline abundance", fill = "Adjusted for\nresponse, up in:") +
  theme_classic()
#Figure
pdf("Adjusted_plot_CD8_subsets_treatment.pdf",width=5,height=6)
plot(adjusted_plot)
dev.off()

#Luminex analyses (in extended cohort)
#Load metadata
luminex_spectral<-as.data.frame(read_xlsx("luminex_spectral.xlsx"))
luminex_validation<-as.data.frame(read_xlsx("luminex_validation.xlsx"))

luminex_merged<-merge(luminex_spectral,
                      luminex_validation,
                      by=intersect(colnames(luminex_spectral),colnames(luminex_validation)),
                      all=T)

luminex_merged <- luminex_merged %>% dplyr::filter(Project=="Validation" | Project =="Spectral")

newdata<-luminex_merged
newdata$main_group<-ifelse(newdata$ICI_type=="anti-PD-1","anti-PD-(L)1",
                           ifelse(newdata$ICI_type=="anti-PD-L1","anti-PD-(L)1",
                                  ifelse(newdata$ICI_type=="ipinivo","ipi-nivo",
                                         ifelse(newdata$ICI_type=="ipi_mono","ipi-nivo",
                                                "anti-PD-(L)1"))))


baseline_data <- newdata %>% dplyr::filter(Timepoint==1 & (Response_steroids=="R" | Response_steroids=="NR"))

#Filter out patients with too many days between blood collection and start of steroids, or blood collection after starts of steroids
baseline_data$Time_blood_start_pred <- as.numeric(baseline_data$Time_blood_start_pred)
baseline_data$Time_blood_start_pred[is.na(baseline_data$Time_blood_start_pred)] <- 999
baseline_data <- baseline_data %>% dplyr::filter((Time_blood_start_pred>(-1) & Time_blood_start_pred<5) | Time_blood_start_pred==999)

baseline_long <- baseline_data %>% dplyr::select(c("UNI","IL5","IL6",
                                                   "IL10","IL12","IL13",
                                                   "IL17","IL21","IL23",
                                                   "TNFa","IFNg","APRIL",
                                                   "CCL2","CCL4","CCL17",
                                                   "CXCL9","CXCL10","CXCL13",
                                                   "CD40L","MMP10","sIL2Ra",
                                                   "GranB","TGFb1","TACI",
                                                   "Gal9","MMP3","MMP9",
                                                   "Response_steroids","main_group"))%>%
  pivot_longer(!c(UNI,Response_steroids,main_group),names_to="Analyte",values_to="concentration")

baseline_long$Analyte <- factor(baseline_long$Analyte,
                                levels=c("TGFb1","IL10",
                                         "APRIL","CXCL13","IL21","TACI",
                                         "Gal9","CXCL9","CXCL10","GranB","IL12","IFNg","TNFa",
                                         "CCL17","IL5","IL13",
                                         "IL6","IL23","IL17",
                                         "CCL2","CCL4","sIL2Ra","CD40L","MMP10","MMP3","MMP9"))
baseline_long$main_group <- factor(baseline_long$main_group,
                                   levels=c("anti-PD-(L)1","ipi-nivo"))
baseline_long<-baseline_long%>%arrange(Analyte,concentration)

#Volcano plots of Luminex analytes at baseline 
#No MA plots on purpose, because expression of different analytes has wide range without biological meaning
#Create dataframe with average expression, FC R vs NR, and wilcox_test unadj Pval for every analyte
volcano_frame <- data.frame()
j <- 1
for(i in unique(baseline_long$Analyte)){
  data <- baseline_long %>% dplyr::filter(Analyte==i)
  log2FC <- log2(mean(data$concentration[data$Response_steroids=="NR"])/mean(data$concentration[data$Response_steroids=="R"]))
  pval <- coin::wilcox_test(data$concentration~factor(data$Response_steroids))
  volcano_frame[j,1] <- i
  volcano_frame[j,2] <- log2FC
  volcano_frame[j,3] <- paste0(as.numeric(coin::pvalue(pval)))
  j <- j+1
}
colnames(volcano_frame)<-c("Analyte","log2FC","pval_raw")
volcano_frame$pval_raw <- as.numeric(volcano_frame$pval_raw)
q <- p.adjust(volcano_frame$pval_raw,method="BH",n=length(volcano_frame$pval_raw))
volcano_frame$padj <- q
volcano_frame$P.value.adj <- ifelse(volcano_frame$padj<0.05,"*","NS")

volcano_frame$color <- ifelse(volcano_frame$log2FC>0.58,"NR",
                              ifelse(volcano_frame$log2FC<(-0.58),"R","Neutral"))
volcano_frame$color <- factor(volcano_frame$color,levels=c("NR","R","Neutral"))
volcano_plot <- ggplot(volcano_frame,aes(log2FC,-log10(padj),color=color))+
  geom_point()+
  geom_text_repel(aes(label=Analyte),
                  max.overlaps = 20)+
  scale_color_manual(values=c(NR="#440154",R="#5ec962",Neutral="grey"))+
  labs(title="Baseline Luminex")+
  ylab("-log10(Padj)")+
  xlab("log2FC Non-responder versus Responder")+
  xlim(-2,2.5)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  geom_vline(xintercept = -log2(1.5),linetype="dashed")+
  geom_vline(xintercept = log2(1.5),linetype="dashed")+
  theme_classic()+
  theme(legend.position = "none")
pdf("Volcano_Luminex_baseline_by_response.pdf",height=3.5,width=4)
plot(volcano_plot)
dev.off()

volcano_frame <- data.frame()
j <- 1
for(i in unique(baseline_long$Analyte)){
  data <- baseline_long %>% dplyr::filter(Analyte==i)
  log2FC <- log2(mean(data$concentration[data$main_group=="ipi-nivo"])/mean(data$concentration[data$main_group=="anti-PD-(L)1"]))
  pval <- coin::wilcox_test(data$concentration~factor(data$main_group))
  volcano_frame[j,1] <- i
  volcano_frame[j,2] <- log2FC
  volcano_frame[j,3] <- paste0(as.numeric(coin::pvalue(pval)))
  j <- j+1
}
colnames(volcano_frame)<-c("Analyte","log2FC","pval_raw")
volcano_frame$pval_raw <- as.numeric(volcano_frame$pval_raw)
q <- p.adjust(volcano_frame$pval_raw,method="BH",n=length(volcano_frame$pval_raw))
volcano_frame$padj <- q
volcano_frame$P.value.adj <- ifelse(volcano_frame$padj<0.05,"*","NS")

volcano_frame$color <- ifelse(volcano_frame$log2FC>0.58,"ipi-nivo",
                              ifelse(volcano_frame$log2FC<(-0.58),"anti-PD-(L)1","Neutral"))
volcano_frame$color <- factor(volcano_frame$color,levels=c("ipi-nivo","anti-PD-(L)1","Neutral"))
volcano_plot <- ggplot(volcano_frame,aes(log2FC,-log10(padj),color=color))+
  geom_point()+
  geom_text_repel(aes(label=Analyte),
                  max.overlaps = 20)+
  scale_color_manual(values=c(`ipi-nivo`="#BC3754",`anti-PD-(L)`="#F98E09",Neutral="grey"))+
  labs(title="")+
  ylab("-log10(Padj)")+
  xlab("log2FC [ipi+/-nivo] versus [anti-PD-(L)1+/-targeted/chemo]")+
  xlim(-2,2.5)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  geom_vline(xintercept = -log2(1.5),linetype="dashed")+
  geom_vline(xintercept = log2(1.5),linetype="dashed")+
  theme_classic()+
  theme(legend.position = "none")
pdf("Volcano_Luminex_baseline_by_treatment.pdf",height=4,width=5)
plot(volcano_plot)
dev.off()

newdata$combined <- paste0(newdata$main_group,"_",newdata$Response_steroids)
subdata <- newdata %>% dplyr::filter(Response_steroids=="NR" | Response_steroids=="R")
subdata$combined <- factor(subdata$combined,levels=c("anti-PD-(L)1_R","anti-PD-(L)1_NR",
                                                     "ipi-nivo_R","ipi-nivo_NR"))
plot <- ggplot(subdata,aes(combined,log10(IL17)))+
  geom_boxplot(outlier.shape=NA,aes(fill=main_group),alpha=0.5)+
  geom_point(shape=21,position=position_jitter(height=0),aes(color=Response_steroids))+
  theme_classic()+
  scale_color_manual(labels=c("Steroid non-responder","Steroid responder"),values=c("#440154","#5ec962"))+
  scale_fill_manual(labels=c("anti-PD-(L)1 +/- chemo/targeted","ipilimumab +/- nivolumab"),values=c("#F98E09","#BC3754"))+
  xlab("ICI group x Steroid response")+
  ylab("Log10 (serum IL-17 [pg/mL])")+
  theme(axis.text.x = element_blank())+
  labs(color=NULL,fill=NULL)+
  geom_signif(annotations="",
              y_position=3,
              xmin=1,
              xmax=2,
              tip_length=c(0,0))+
  geom_signif(annotations="",
              y_position=3,
              xmin=3,
              xmax=4,
              tip_length=c(0,0))+
  geom_signif(annotations=paste0(expression(italic(P) == "0.047")),parse=T,
              y_position=3.3,
              xmin=1.5,
              xmax=3.5,
              tip_length=c(0.14,0.14))
pdf("IL17_by_response_treatment.pdf",height=3.5,width=5.5)
plot(plot)
dev.off()
wilcox_test(subdata$IL17~factor(subdata$main_group))

#Baseline blood cell counts / CRP
luminex_spectral<-as.data.frame(read_xlsx("luminex_spectral.xlsx"))
luminex_validation<-as.data.frame(read_xlsx("luminex_validation.xlsx"))

luminex_merged<-merge(luminex_spectral,
                      luminex_validation,
                      by=intersect(colnames(luminex_spectral),colnames(luminex_validation)),
                      all=T)

luminex_merged$Response_steroids <- ifelse(luminex_merged$Project=="Healthy","HD",
                                           ifelse(luminex_merged$Project=="ICI_controle_melanoom","ICI_no_irAE",
                                                  luminex_merged$Response_steroids))

newdata <- luminex_merged
newdata$main_group<-ifelse(newdata$ICI_type=="ipinivo","ipi-nivo",
                           ifelse(newdata$ICI_type=="ipi_mono","ipi-nivo",
                                  ifelse(is.na(newdata$ICI_type),"HD","anti-PD-(L)1_based")))
newdata$main_group[is.na(newdata$ICI_type)] <- "HD" 

baseline_data <- newdata %>% dplyr::filter(Timepoint==1 & Response_steroids!="noIS" & Response_steroids!="R_infection" & Response_steroids!="ICI_no_irAE" & Response_steroids!="HD")

#Filter out patients with too many days between blood collection and start of steroids, or blood collection after starts of steroids
baseline_data$Time_blood_start_pred <- as.numeric(baseline_data$Time_blood_start_pred)
baseline_data$Time_blood_start_pred[is.na(baseline_data$Time_blood_start_pred)] <- 999
baseline_data <- baseline_data %>% dplyr::filter((Time_blood_start_pred>(-1) & Time_blood_start_pred<5) | Time_blood_start_pred==999)

colors <- c("R"="#5ec962","NR"="#440154")
baseline_data$Response_steroids <- factor(baseline_data$Response_steroids,levels=c("R","NR"))
CRP_response <- ggplot(baseline_data,aes(Response_steroids,CRP,fill=Response_steroids))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(shape=21,position=position_jitter(height=0))+
  geom_signif(annotation=paste0(expression(italic(P) == "0.55")),parse=T,y_position=350,xmin=1,xmax=2,tip_length=c(0,0))+
  scale_fill_manual(labels=c("Responder","Non-responder"),values=colors)+
  theme_classic()+
  theme(legend.title=element_blank(),
        axis.text.x=element_blank())+
  xlab("")+
  ylab("Blood CRP (mg/L)")
pdf("CRP_response_boxplot.pdf",width=4,height=4)
plot(CRP_response)
dev.off()
wilcox_test(baseline_data$CRP~baseline_data$Response_steroids)

baseline_data$leukos <- baseline_data$neutros+baseline_data$lymphos+baseline_data$monos+baseline_data$eos
baseline_data$NLR <- baseline_data$neutros/baseline_data$lymphos
Leuko_response <- ggplot(baseline_data,aes(Response_steroids,leukos,fill=Response_steroids))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(shape=21,position=position_jitter(height=0))+
  geom_signif(annotation=paste0(expression(italic(P) == "6.6e-3")),parse=T,y_position=20,xmin=1,xmax=2,tip_length=c(0,0))+
  scale_fill_manual(labels=c("Responder","Non-responder"),values=colors)+
  theme_classic()+
  theme(legend.title=element_blank(),
        axis.text.x=element_blank())+
  xlab("")+
  ylab("Blood WBC (x10^9/L)")
pdf("Leuko_response_boxplot.pdf",width=4,height=4)
plot(Leuko_response)
dev.off()
wilcox_test(baseline_data$leukos~baseline_data$Response_steroids)

Eos_response <- ggplot(baseline_data,aes(Response_steroids,eos,fill=Response_steroids))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(shape=21,position=position_jitter(height=0))+
  geom_signif(annotation=paste0(expression(italic(P) == 0.047)),parse=T,y_position=3,xmin=1,xmax=2,tip_length=c(0,0))+
  scale_fill_manual(labels=c("Responder","Non-responder"),values=colors)+
  theme_classic()+
  theme(legend.title=element_blank(),
        axis.text.x=element_blank())+
  xlab("")+
  ylab("Blood eosinophil count (x10^9/L)")+
  ylim(0,3.2)
pdf("Eos_response_boxplot.pdf",width=4,height=3)
plot(Eos_response)
dev.off()
wilcox_test(baseline_data$eos~baseline_data$Response_steroids)

Neutros_response <- ggplot(baseline_data,aes(Response_steroids,neutros,fill=Response_steroids))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(shape=21,position=position_jitter(height=0))+
  geom_signif(annotation=paste0(expression(italic(P) == "0.030")),parse=T,y_position=15,xmin=1,xmax=2,tip_length=c(0,0))+
  scale_fill_manual(labels=c("Responder","Non-responder"),values=colors)+
  theme_classic()+
  theme(legend.title=element_blank(),
        axis.text.x=element_blank())+
  xlab("")+
  ylab("Blood neutrophil count (x10^9/L)")+
  ylim(0,16)
pdf("Neutros_response_boxplot.pdf",width=4,height=3)
plot(Neutros_response)
dev.off()
wilcox_test(baseline_data$neutros~baseline_data$Response_steroids)

Lymphos_response <- ggplot(baseline_data,aes(Response_steroids,lymphos,fill=Response_steroids))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(shape=21,position=position_jitter(height=0))+
  geom_signif(annotation=paste0(expression(italic(P) == "0.020")),parse=T,y_position=5,xmin=1,xmax=2,tip_length=c(0,0))+
  scale_fill_manual(labels=c("Responder","Non-responder"),values=colors)+
  theme_classic()+
  theme(legend.title=element_blank(),
        axis.text.x=element_blank())+
  xlab("")+
  ylab("Blood lymphocyte count (x10^9/L)")+
  ylim(0,6)
pdf("Lymphos_response_boxplot.pdf",width=4,height=3)
plot(Lymphos_response)
dev.off()
wilcox_test(baseline_data$lymphos~baseline_data$Response_steroids)

baseline_data$NLR <- baseline_data$neutros/baseline_data$lymphos
NLR_response <- ggplot(baseline_data,aes(Response_steroids,NLR,fill=Response_steroids))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(shape=21,position=position_jitter(height=0))+
  geom_signif(annotation=paste0(expression(italic(P) == "0.91")),parse=T,y_position=15,xmin=1,xmax=2,tip_length=c(0,0))+
  scale_fill_manual(labels=c("Responder","Non-responder"),values=colors)+
  theme_classic()+
  theme(legend.title=element_blank(),
        axis.text.x=element_blank())+
  xlab("")+
  ylab("Neutrophil-to-lymphocyte ratio")+
  ylim(0,17)
pdf("NLR_response_boxplot.pdf",width=4,height=3)
plot(NLR_response)
dev.off()
wilcox_test(baseline_data$NLR~baseline_data$Response_steroids)

#Figure 2 change over time
#Change of cytokines after relative to before steroids
change <- newdata %>% dplyr::filter(Timepoint==1 | Timepoint==2)
change_long <- change %>% dplyr::select(c("UNI","IL5","IL6",
                                   "IL10","IL12","IL13",
                                   "IL17","IL21","IL23",
                                   "TNFa","IFNg","APRIL",
                                   "CCL2","CCL4","CCL17",
                                   "CXCL9","CXCL10","CXCL13",
                                   "CD40L","MMP10","sIL2Ra",
                                   "GranB","TGFb1","TACI",
                                   "Gal9","MMP3","MMP9",
                                   "Response_steroids","main_group","Timepoint"))%>%
  pivot_longer(!c(UNI,Response_steroids,main_group,Timepoint),names_to="Analyte",values_to="concentration")

change_long$Response_steroids <- factor(change_long$Response_steroids,
                                        levels=c("HD","NR","R"))

change_long$Analyte <- factor(change_long$Analyte,
                              levels=c("TGFb1","IL10",
                                       "APRIL","CXCL13","IL21","TACI",
                                       "Gal9","CXCL9","CXCL10","GranB","IL12","IFNg","TNFa",
                                       "CCL17","IL5","IL13",
                                       "IL6","IL23","IL17",
                                       "CCL2","CCL4","sIL2Ra","CD40L","MMP10","MMP3","MMP9"))

change_long$FC <- 1
for(i in unique(change_long$Analyte)){
  for(j in unique(change_long$UNI)){
    change_long$FC[change_long$Analyte==i&change_long$UNI==j&change_long$Timepoint==2] <- change_long$concentration[change_long$Analyte==i&change_long$UNI==j&change_long$Timepoint==2]/change_long$concentration[change_long$Analyte==i&change_long$UNI==j&change_long$Timepoint==1] 
  }
}

change_long <- change_long %>% dplyr::filter(Timepoint==2)
change_long<-change_long%>%arrange(Analyte,concentration)

test_table<-data.frame(Group=NA,Pvalue_kw=NA,Pvalue_t=NA)
j<-1
for(i in unique(change_long$Analyte)){
  test_set<-change_long%>%dplyr::filter(Analyte==i)
  test_table[j,1]<-paste0(i)
  result<-wilcox_test(FC~Response_steroids,data=test_set)
  test_table[j,2]<-paste0(as.numeric(coin::pvalue(result)))
  result<-t.test(test_set$FC~test_set$Response_steroids)
  test_table[j,3]<-paste0(as.numeric(result$p.value))
  j<-j+1
}

test_table$Pvalue_kw <- as.numeric(test_table$Pvalue_kw)
q <- p.adjust(test_table$Pvalue_kw,method="hochberg",n=length(test_table$Pvalue_kw))
test_table$padj_kw <- q

plot_table <- data.frame()
j<-1
for(i in unique(change_long$Analyte)){
  data <- change_long %>% dplyr::filter(Analyte==i)
  plot_table[j,1]<-paste0(i)
  plot_table[j,2]<-log2(median(data$FC[data$Response_steroids=="R"]))
  plot_table[j,3]<-log2(median(data$FC[data$Response_steroids=="NR"]))
  plot_table[j,4]<-test_table$Pvalue_kw[test_table$Group==i]
  plot_table[j,5]<-test_table$padj_kw[test_table$Group==i]
  j<-j+1
}
colnames(plot_table)<-c("Analyte","log2FC_in_R","log2FC_in_NR","Pval_kw","Padj_kw")
plot_table$P.value.raw <- ifelse(plot_table$Pval_kw<0.05,"*","NS")
plot_table$direction <- plot_table$log2FC_in_NR-plot_table$log2FC_in_R
plot_table$color <- ifelse(plot_table$P.value.raw<0.05&plot_table$direction>0,"NR",
                           ifelse(plot_table$P.value.raw<0.05&plot_table$direction<0,"R","neutral"))
  
plotje<-ggplot(plot_table,aes(log2FC_in_R,log2FC_in_NR))+
  geom_point(aes(color=color))+
  scale_color_manual(labels=c("NR"="Skewed towards NR at P<0.05",
                              "neutral"="Equal change in R and NR",
                              "R"="Skewed towards R at P<0.05"),
                     values=c("NR"="#440154","neutral"="grey","R"="#5ec962"))+
  geom_text_repel(
    data = subset(plot_table, color %in% c("NR", "R")),
    aes(label = Analyte),
    max.overlaps = Inf)+
  labs(title="",fill="",color="")+
  geom_hline(yintercept=0,linetype="dashed",color="grey")+
  geom_vline(xintercept=0,linetype="dashed",color="grey")+
  ylab("Median log2FC in NR after vs. before steroids")+
  xlab("Median log2FC in R after vs. before steroids")+
  xlim(-2,2)+
  ylim(-2,2)+
  theme_classic()
pdf("Change_T2_Luminex_by_response_FCx2.pdf",width=7,height=5)
plot(plotje)
dev.off()

boxplot_data <- newdata
boxplot_data <- boxplot_data %>% dplyr::filter(Project=="Spectral" & (Timepoint==1 | Timepoint==2))

boxplot_data$Response_steroids <- factor(boxplot_data$Response_steroids,levels=c("NR","R"))
median_df <- boxplot_data %>% 
  group_by(Timepoint,Response_steroids) %>% 
  summarize(median_value = median(log10(CXCL10)),.groups="drop")
boxplot1 <- ggplot(boxplot_data,aes(Timepoint,log10(CXCL10),group=UNI,color=Response_steroids,fill=Response_steroids))+
  geom_point(shape=21)+
  geom_line(alpha=0.3)+
  geom_line(data=median_df,aes(y=median_value,group=Response_steroids),linewidth=1,linetype="dashed")+
  scale_color_manual(labels=c("Non-responder","Responder"),
                     values=c("#440154","#5ec962"))+
  scale_fill_manual(labels=c("Non-responder","Responder"),
                    values=c("#440154","#5ec962"))+
  labs(fill=NULL,group=NULL,color=NULL)+
  ylab("Serum CXCL10 (pg/mL)")+
  xlab("")+
  annotation_logticks(sides="l")+
  scale_y_continuous(breaks=c(0,1,2,3,4),labels=c(1.0,10,100,1000,10000))+
  scale_x_discrete(labels=c("1" = "Pre steroids", "2" = "On steroids"))+
  theme_classic()
pdf("Boxplot_change_CXCL10.pdf",width=3.7,height=2.1)
plot(boxplot1)
dev.off()

median_df <- boxplot_data %>% 
  group_by(Timepoint,Response_steroids) %>% 
  summarize(median_value = median(log10(CCL4)),.groups="drop")
boxplot2 <- ggplot(boxplot_data,aes(Timepoint,log10(CCL4),group=UNI,color=Response_steroids,fill=Response_steroids))+
  geom_point(shape=21)+
  geom_line(alpha=0.3)+
  geom_line(data=median_df,aes(y=median_value,group=Response_steroids),linewidth=1,linetype="dashed")+
  scale_color_manual(labels=c("Non-responder","Responder"),
                     values=c("#440154","#5ec962"))+
  scale_fill_manual(labels=c("Non-responder","Responder"),
                    values=c("#440154","#5ec962"))+
  labs(fill=NULL,group=NULL,color=NULL)+
  ylab("Serum CCL4 (pg/mL)")+
  xlab("")+
  annotation_logticks(sides="l")+
  scale_y_continuous(breaks=c(0,1,2,3,4),labels=c(1.0,10,100,1000,10000))+
  scale_x_discrete(labels=c("1" = "Pre steroids", "2" = "On steroids"))+
  theme_classic()
pdf("Boxplot_change_CCL4.pdf",width=3.7,height=2.1)
plot(boxplot2)
dev.off()

median_df <- boxplot_data %>% 
  group_by(Timepoint,Response_steroids) %>% 
  summarize(median_value = median(log10(MMP9)),.groups="drop")
boxplot2 <- ggplot(boxplot_data,aes(Timepoint,log10(MMP9),group=UNI,color=Response_steroids,fill=Response_steroids))+
  geom_point(shape=21)+
  geom_line(alpha=0.3)+
  geom_line(data=median_df,aes(y=median_value,group=Response_steroids),linewidth=1,linetype="dashed")+
  scale_color_manual(labels=c("Non-responder","Responder"),
                     values=c("#440154","#5ec962"))+
  scale_fill_manual(labels=c("Non-responder","Responder"),
                    values=c("#440154","#5ec962"))+
  labs(fill=NULL,group=NULL,color=NULL)+
  ylab("Serum MMP9 (pg/mL)")+
  xlab("")+
  annotation_logticks(sides="l")+
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9),labels=c(1.0,10,100,1000,10000,1e5,1e6,1e7,1e8,1e9))+
  scale_x_discrete(labels=c("1" = "Pre steroids", "2" = "On steroids"))+
  theme_classic()
pdf("Boxplot_change_MMP9.pdf",width=3.7,height=2.1)
plot(boxplot2)
dev.off()

median_df <- boxplot_data %>% 
  group_by(Timepoint,Response_steroids) %>% 
  summarize(median_value = median(log10(APRIL)),.groups="drop")
boxplot2 <- ggplot(boxplot_data,aes(Timepoint,log10(APRIL),group=UNI,color=Response_steroids,fill=Response_steroids))+
  geom_point(shape=21)+
  geom_line(alpha=0.3)+
  geom_line(data=median_df,aes(y=median_value,group=Response_steroids),linewidth=1,linetype="dashed")+
  scale_color_manual(labels=c("Non-responder","Responder"),
                     values=c("#440154","#5ec962"))+
  scale_fill_manual(labels=c("Non-responder","Responder"),
                    values=c("#440154","#5ec962"))+
  labs(fill=NULL,group=NULL,color=NULL)+
  ylab("Serum APRIL (pg/mL)")+
  xlab("")+
  annotation_logticks(sides="l")+
  scale_y_continuous(breaks=c(0,1,2,3,4),labels=c(1.0,10,100,1000,10000))+
  scale_x_discrete(labels=c("1" = "Pre steroids", "2" = "On steroids"))+
  theme_classic()
pdf("Boxplot_change_APRIL.pdf",width=3.7,height=2.1)
plot(boxplot2)
dev.off()

median_df <- boxplot_data %>% 
  group_by(Timepoint,Response_steroids) %>% 
  summarize(median_value = median(log10(CXCL13)),.groups="drop")
boxplot2 <- ggplot(boxplot_data,aes(Timepoint,log10(CXCL13),group=UNI,color=Response_steroids,fill=Response_steroids))+
  geom_point(shape=21)+
  geom_line(alpha=0.3)+
  geom_line(data=median_df,aes(y=median_value,group=Response_steroids),linewidth=1,linetype="dashed")+
  scale_color_manual(labels=c("Non-responder","Responder"),
                     values=c("#440154","#5ec962"))+
  scale_fill_manual(labels=c("Non-responder","Responder"),
                    values=c("#440154","#5ec962"))+
  labs(fill=NULL,group=NULL,color=NULL)+
  ylab("Serum CXCL13 (pg/mL)")+
  xlab("")+
  annotation_logticks(sides="l")+
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9),labels=c(1.0,10,100,1000,10000,1e5,1e6,1e7,1e8,1e9))+
  scale_x_discrete(labels=c("1" = "Pre steroids", "2" = "On steroids"))+
  theme_classic()
pdf("Boxplot_change_CXCL13.pdf",width=3.7,height=2.1)
plot(boxplot2)
dev.off()

#Now for main_cluster cell types
cell_counts<-dataset%>%group_by(subject,timepoint,main_cluster)%>%tally()
cell_counts<-pivot_wider(cell_counts,names_from="main_cluster",values_from="n")
cell_counts<-as.data.frame(cell_counts)
cell_counts$total<-rowSums(cell_counts[3:ncol(cell_counts)],na.rm=T)
cell_counts[is.na(cell_counts)]<-0

cell_counts_relative<-cell_counts
for(r in 1:nrow(cell_counts_relative)){
  for(c in 3:ncol(cell_counts_relative)){
    cell_counts_relative[r,c]<-cell_counts_relative[r,c]/cell_counts_relative[r,ncol(cell_counts_relative)]
  }
}

cell_counts_longitudinal<-cell_counts_relative
attach(cell_counts_longitudinal)
for(i in unique(subject)){
  if(length(timepoint[subject==i])>1){
    for(j in 2:length(timepoint[subject==i])){
      for(k in 3:ncol(cell_counts_longitudinal)){
        cell_counts_longitudinal[cell_counts_longitudinal$timepoint==j&cell_counts_longitudinal$subject==i,k]<-cell_counts_longitudinal[cell_counts_longitudinal$timepoint==j&cell_counts_longitudinal$subject==i,k]/cell_counts_longitudinal[cell_counts_longitudinal$timepoint==1&cell_counts_longitudinal$subject==i,k]
      }
    }
  }
}

for(i in unique(subject)){
  for(k in 3:ncol(cell_counts_longitudinal)){
    cell_counts_longitudinal[cell_counts_longitudinal$timepoint==1&cell_counts_longitudinal$subject==i,k]<-1
  }
}

#Load metadata and Fortessa data
fortessa<-read_xlsx("Data_incl_metadata_all_batches_relative.xlsx")
fortessa<-fortessa%>%dplyr::rename(
  subject=Sample,
  timepoint=Timepoint
)
newdata<-merge(cell_counts_longitudinal,fortessa,by=c("subject","timepoint"),all=T)
newdata$main_group<-ifelse(newdata$Treatment=="anti-PD-1","anti-PD-(L)1",
                           ifelse(newdata$Treatment=="anti-PD-L1","anti-PD-(L)1",
                                  ifelse(newdata$Treatment=="ipi-nivo","ipi-nivo",
                                         ifelse(newdata$Treatment=="ipi_mono","ipi_mono",
                                                "healthy"))))
newdata<-newdata%>%dplyr::filter(timepoint==2)
newdata_long<-newdata%>%dplyr::select(c("subject","B_cells","CD4CD8_DP_T",                      
                                        "CD4_Tcm","CD4_Tem","CD4_Tn","CD4_Treg",                         
                                        "CD8_TEMRA","CD8_Tcm","CD8_Tem","CD8_Tn",                            
                                        "DC_like","Monocytes","NKT_like","NK_cells",         
                                        "NKT_like","NK_cells","Undefined","Response",
                                        "main_group","Time_start_steroids_rel",
                                        "Time_blood_rel"))%>%
  pivot_longer(!c(subject,Response,main_group,Time_start_steroids_rel,Time_blood_rel),names_to="main_cluster",values_to="fold_change")

test_table<-data.frame(Group=NA,Pvalue_kw=NA,Pvalue_t=NA)
j<-1
for(i in unique(newdata_long$main_cluster)){
  test_set<-newdata_long%>%dplyr::filter(main_cluster==i)
  test_table[j,1]<-paste0(i)
  result<-wilcox_test(fold_change~as.factor(Response),data=test_set)
  test_table[j,2]<-paste0(as.numeric(coin::pvalue(result)))
  result<-t.test(test_set$fold_change~test_set$Response)
  test_table[j,3]<-paste0(as.numeric(result$p.value))
  j<-j+1
}
test_table$padj_kw <- p.adjust(test_table$Pvalue_kw,method="BH")

plot_table <- data.frame()
j<-1
for(i in unique(newdata_long$main_cluster)){
  data <- newdata_long %>% dplyr::filter(main_cluster==i)
  plot_table[j,1]<-paste0(i)
  plot_table[j,2]<-log2(median(data$fold_change[data$Response=="R"]))
  plot_table[j,3]<-log2(median(data$fold_change[data$Response=="NR"]))
  plot_table[j,4]<-test_table$Pvalue_kw[test_table$Group==i]
  plot_table[j,5]<-test_table$padj_kw[test_table$Group==i]
  j<-j+1
}
colnames(plot_table)<-c("main_cluster","log2FC_in_R","log2FC_in_NR","Pval_kw","Padj_kw")
plot_table$P.value.raw <- ifelse(plot_table$Pval_kw<0.05,"*","NS")
plot_table$direction <- plot_table$log2FC_in_NR-plot_table$log2FC_in_R
plot_table$color <- ifelse(plot_table$P.value.raw<0.05&plot_table$direction>0,"NR",
                           ifelse(plot_table$P.value.raw<0.05&plot_table$direction<0,"R","neutral"))

plotje<-ggplot(plot_table,aes(log2FC_in_R,log2FC_in_NR))+
  geom_point(aes(color=color))+
  scale_color_manual(labels=c("NR"="Skewed towards NR at P<0.05",
                              "neutral"="Equal change in R and NR",
                              "R"="Skewed towards R at P<0.05"),
                     values=c("NR"="#440154","neutral"="grey","R"="#5ec962"))+
  geom_text_repel(
    data = subset(plot_table, color %in% c("NR", "R")),
    aes(label = main_cluster),
    max.overlaps = Inf)+
  labs(title="",fill="",color="")+
  geom_hline(yintercept=0,linetype="dashed",color="grey")+
  geom_vline(xintercept=0,linetype="dashed",color="grey")+
  ylab("Median log2FC in NR after vs. before steroids")+
  xlab("Median log2FC in R after vs. before steroids")+
  xlim(-1.5,1.5)+
  ylim(-1.5,1.5)+
  theme_classic()
pdf("Change_T2_main_clusters_by_response_FCx2.pdf",width=7,height=5)
plot(plotje)
dev.off()

#Same analysis, but for CD4 T cells:
cell_counts<-dataset%>%dplyr::filter(large_cluster=="CD4_cell")%>%group_by(subject,timepoint,detail_cluster)%>%tally()
cell_counts<-pivot_wider(cell_counts,names_from="detail_cluster",values_from="n")
cell_counts<-as.data.frame(cell_counts)
cell_counts$total<-rowSums(cell_counts[3:ncol(cell_counts)],na.rm=T)
cell_counts[is.na(cell_counts)]<-0

cell_counts_relative<-cell_counts
for(r in 1:nrow(cell_counts_relative)){
  for(c in 3:ncol(cell_counts_relative)){
    cell_counts_relative[r,c]<-cell_counts_relative[r,c]/cell_counts_relative[r,ncol(cell_counts_relative)]
  }
}

#Because CD4_Tem_3 and CD4_Tem_5 are highly similar, proliferative clusters, we add one variable in which both subsets are combined
#(for more power)
cell_counts_relative$CD4_Tem_3_and_5 <- cell_counts_relative$CD4_Tem_3 + cell_counts_relative$CD4_Tem_5

cell_counts_longitudinal<-cell_counts_relative
attach(cell_counts_longitudinal)
for(i in unique(subject)){
  if(length(timepoint[subject==i])>1){
    for(j in 2:length(timepoint[subject==i])){
      for(k in 3:ncol(cell_counts_longitudinal)){
        cell_counts_longitudinal[cell_counts_longitudinal$timepoint==j&cell_counts_longitudinal$subject==i,k]<-cell_counts_longitudinal[cell_counts_longitudinal$timepoint==j&cell_counts_longitudinal$subject==i,k]/cell_counts_longitudinal[cell_counts_longitudinal$timepoint==1&cell_counts_longitudinal$subject==i,k]
      }
    }
  }
}

for(i in unique(subject)){
  for(k in 3:ncol(cell_counts_longitudinal)){
    cell_counts_longitudinal[cell_counts_longitudinal$timepoint==1&cell_counts_longitudinal$subject==i,k]<-1
  }
}

#Load metadata and Fortessa data
fortessa<-read_xlsx("Data_incl_metadata_all_batches_relative.xlsx")
fortessa<-fortessa%>%dplyr::rename(
  subject=Sample,
  timepoint=Timepoint
)
newdata<-merge(cell_counts_longitudinal,fortessa,by=c("subject","timepoint"),all=T)
newdata$main_group<-ifelse(newdata$Treatment=="anti-PD-1","anti-PD-(L)1",
                           ifelse(newdata$Treatment=="anti-PD-L1","anti-PD-(L)1",
                                  ifelse(newdata$Treatment=="ipi-nivo","ipi-nivo",
                                         ifelse(newdata$Treatment=="ipi_mono","ipi_mono",
                                                "healthy"))))
newdata<-newdata%>%dplyr::filter(timepoint==2)

newdata_long<-newdata%>%dplyr::select(c("subject","CD4_Tn_1","CD4_Tn_2","CD4_Tcm_1","CD4_Tcm_2","CD4_Tem_1",
                                        "CD4_Tem_2","CD4_Tem_3","CD4_Tem_4","CD4_Tem_5","CD4_Tem_6",
                                        "CD4_Tem_7","CD4_Tem_8","CD4_Treg_1","CD4_Treg_2","CD4_Treg_3",
                                        "CD4_Treg_4","CD4CD8_DP_T","CD4_Tem_3_and_5","Response","main_group","Time_start_steroids_rel",
                                        "Time_blood_rel"))%>%
  pivot_longer(!c(subject,Response,main_group,Time_start_steroids_rel,Time_blood_rel),names_to="main_cluster",values_to="fold_change")

test_table<-data.frame(Group=NA,Pvalue_kw=NA,Pvalue_t=NA)
j<-1
for(i in unique(newdata_long$main_cluster)){
  test_set<-newdata_long%>%dplyr::filter(main_cluster==i)
  test_table[j,1]<-paste0(i)
  result<-wilcox_test(fold_change~as.factor(Response),data=test_set)
  test_table[j,2]<-paste0(as.numeric(coin::pvalue(result)))
  result<-t.test(test_set$fold_change~test_set$Response)
  test_table[j,3]<-paste0(as.numeric(result$p.value))
  j<-j+1
}
test_table$padj_kw <- p.adjust(test_table$Pvalue_kw,method="BH")

plot_table <- data.frame()
j<-1
for(i in unique(newdata_long$main_cluster)){
  data <- newdata_long %>% dplyr::filter(main_cluster==i)
  plot_table[j,1]<-paste0(i)
  plot_table[j,2]<-log2(median(data$fold_change[data$Response=="R"]))
  plot_table[j,3]<-log2(median(data$fold_change[data$Response=="NR"]))
  plot_table[j,4]<-test_table$Pvalue_kw[test_table$Group==i]
  plot_table[j,5]<-test_table$padj_kw[test_table$Group==i]
  j<-j+1
}
colnames(plot_table)<-c("main_cluster","log2FC_in_R","log2FC_in_NR","Pval_kw","Padj_kw")
plot_table$P.value.raw <- ifelse(plot_table$Pval_kw<0.05,"*","NS")
plot_table$direction <- plot_table$log2FC_in_NR-plot_table$log2FC_in_R
plot_table$color <- ifelse(plot_table$P.value.raw<0.05&plot_table$direction>0,"NR",
                           ifelse(plot_table$P.value.raw<0.05&plot_table$direction<0,"R","neutral"))

plotje<-ggplot(plot_table,aes(log2FC_in_R,log2FC_in_NR))+
  geom_point(aes(color=color))+
  scale_color_manual(labels=c("NR"="Skewed towards NR at P<0.05",
                              "neutral"="Equal change in R and NR",
                              "R"="Skewed towards R at P<0.05"),
                     values=c("NR"="#440154","neutral"="grey","R"="#5ec962"))+
  geom_text_repel(
    data = subset(plot_table, color %in% c("NR", "R")),
    aes(label = main_cluster),
    max.overlaps = Inf)+
  labs(title="",fill="",color="")+
  geom_hline(yintercept=0,linetype="dashed",color="grey")+
  geom_vline(xintercept=0,linetype="dashed",color="grey")+
  ylab("Median log2FC in NR after vs. before steroids")+
  xlab("Median log2FC in R after vs. before steroids")+
  xlim(-1.5,1.5)+
  ylim(-1.5,1.5)+
  theme_classic()
pdf("Change_T2_CD4_subsets_by_response_FCx2.pdf",width=7,height=5)
plot(plotje)
dev.off()

boxplot_data <- cell_counts_relative
fortessa<-read_xlsx("Data_incl_metadata_all_batches_relative.xlsx")
fortessa<-fortessa%>%dplyr::rename(
  subject=Sample,
  timepoint=Timepoint
)
boxplot_data<-merge(boxplot_data,fortessa,by=c("subject","timepoint"),all=T)
boxplot_data <- boxplot_data %>% dplyr::filter(Response!="HD" & (timepoint==1 | timepoint==2))

boxplot_data$Response <- factor(boxplot_data$Response,levels=c("NR","R"))
boxplot1 <- ggplot(boxplot_data,aes(timepoint,100*(CD4_Tem_5+CD4_Tem_3),group=subject,color=Response,fill=Response))+
  geom_point(shape=21)+
  geom_line()+
  scale_color_manual(labels=c("Non-responder","Responder"),
                     values=c("#440154","#5ec962"))+
  scale_fill_manual(labels=c("Non-responder","Responder"),
                    values=c("#440154","#5ec962"))+
  labs(fill=NULL,group=NULL,color=NULL)+
  ylab("% CD8 Tem 3 of all CD8 T cells")+
  xlab("")+
  scale_x_discrete(labels=c("1" = "Pre steroids", "2" = "After steroids"))+
  theme_classic()
pdf("Boxplot_change_CD8_tem3.pdf",width=4,height=2.8)
plot(boxplot1)
dev.off()
boxplot_data$CD4_prolif <- boxplot_data$CD4_Tem_5 + boxplot_data$CD4_Tem_3
wilcox_test(boxplot_data$CD4_prolif~boxplot_data$Response)

#Specifically for CD8 T cells:
cell_counts<-dataset%>%dplyr::filter(large_cluster=="CD8_cell")%>%group_by(subject,timepoint,detail_cluster)%>%tally()
cell_counts<-pivot_wider(cell_counts,names_from="detail_cluster",values_from="n")
cell_counts<-as.data.frame(cell_counts)
cell_counts$total<-rowSums(cell_counts[3:ncol(cell_counts)],na.rm=T)
cell_counts[is.na(cell_counts)]<-0

cell_counts_relative<-cell_counts
for(r in 1:nrow(cell_counts_relative)){
  for(c in 3:ncol(cell_counts_relative)){
    cell_counts_relative[r,c]<-cell_counts_relative[r,c]/cell_counts_relative[r,ncol(cell_counts_relative)]
  }
}

cell_counts_longitudinal<-cell_counts_relative
attach(cell_counts_longitudinal)
for(i in unique(subject)){
  if(length(timepoint[subject==i])>1){
    for(j in 2:length(timepoint[subject==i])){
      for(k in 3:ncol(cell_counts_longitudinal)){
        cell_counts_longitudinal[cell_counts_longitudinal$timepoint==j&cell_counts_longitudinal$subject==i,k]<-cell_counts_longitudinal[cell_counts_longitudinal$timepoint==j&cell_counts_longitudinal$subject==i,k]/cell_counts_longitudinal[cell_counts_longitudinal$timepoint==1&cell_counts_longitudinal$subject==i,k]
      }
    }
  }
}

for(i in unique(subject)){
  for(k in 3:ncol(cell_counts_longitudinal)){
    cell_counts_longitudinal[cell_counts_longitudinal$timepoint==1&cell_counts_longitudinal$subject==i,k]<-1
  }
}

#Load metadata and Fortessa data
fortessa<-read_xlsx("Data_incl_metadata_all_batches_relative.xlsx")
fortessa<-fortessa%>%dplyr::rename(
  subject=Sample,
  timepoint=Timepoint
)
newdata<-merge(cell_counts_longitudinal,fortessa,by=c("subject","timepoint"),all=T)
newdata$main_group<-ifelse(newdata$Treatment=="anti-PD-1","anti-PD-(L)1",
                           ifelse(newdata$Treatment=="anti-PD-L1","anti-PD-(L)1",
                                  ifelse(newdata$Treatment=="ipi-nivo","ipi-nivo",
                                         ifelse(newdata$Treatment=="ipi_mono","ipi_mono",
                                                "healthy"))))
newdata<-newdata%>%dplyr::filter(timepoint==2)

newdata_long<-newdata%>%dplyr::select(c("subject","CD8_Tn_1","CD8_Tcm","CD8_Tem_1",
                                        "CD8_Tem_2","CD8_Tem_3","CD8_Tem_4","CD8_TEMRA_1","CD8_TEMRA_2",
                                        "CD8_TEMRA_3","CD8_TEMRA_4","Response","main_group"))%>%
  pivot_longer(!c(subject,Response,main_group),names_to="main_cluster",values_to="fold_change")

test_table<-data.frame(Group=NA,Pvalue_kw=NA,Pvalue_t=NA)
j<-1
for(i in unique(newdata_long$main_cluster)){
  test_set<-newdata_long%>%dplyr::filter(main_cluster==i)
  test_table[j,1]<-paste0(i)
  result<-wilcox_test(fold_change~as.factor(Response),data=test_set)
  test_table[j,2]<-paste0(as.numeric(coin::pvalue(result)))
  result<-t.test(test_set$fold_change~test_set$Response)
  test_table[j,3]<-paste0(as.numeric(result$p.value))
  j<-j+1
}
test_table$padj_kw <- p.adjust(test_table$Pvalue_kw,method="BH")

plot_table <- data.frame()
j<-1
for(i in unique(newdata_long$main_cluster)){
  data <- newdata_long %>% dplyr::filter(main_cluster==i)
  plot_table[j,1]<-paste0(i)
  plot_table[j,2]<-log2(median(data$fold_change[data$Response=="R"]))
  plot_table[j,3]<-log2(median(data$fold_change[data$Response=="NR"]))
  plot_table[j,4]<-test_table$Pvalue_kw[test_table$Group==i]
  plot_table[j,5]<-test_table$padj_kw[test_table$Group==i]
  j<-j+1
}
colnames(plot_table)<-c("main_cluster","log2FC_in_R","log2FC_in_NR","Pval_kw","Padj_kw")
plot_table$P.value.raw <- ifelse(plot_table$Pval_kw<0.05,"*","NS")
plot_table$direction <- plot_table$log2FC_in_NR-plot_table$log2FC_in_R
plot_table$color <- ifelse(plot_table$P.value.raw<0.05&plot_table$direction>0,"NR",
                           ifelse(plot_table$P.value.raw<0.05&plot_table$direction<0,"R","neutral"))

plotje<-ggplot(plot_table,aes(log2FC_in_R,log2FC_in_NR))+
  geom_point(aes(color=color))+
  scale_color_manual(labels=c("NR"="Skewed towards NR at P<0.05",
                              "neutral"="Equal change in R and NR",
                              "R"="Skewed towards R at P<0.05"),
                     values=c("NR"="#440154","neutral"="grey","R"="#5ec962"))+
  geom_text_repel(
    data = subset(plot_table, color %in% c("NR", "R")),
    aes(label = main_cluster),
    max.overlaps = Inf)+
  labs(title="",fill="",color="")+
  geom_hline(yintercept=0,linetype="dashed",color="grey")+
  geom_vline(xintercept=0,linetype="dashed",color="grey")+
  ylab("Median log2FC in NR after vs. before steroids")+
  xlab("Median log2FC in R after vs. before steroids")+
  xlim(-1.5,1.5)+
  ylim(-1.5,1.5)+
  theme_classic()
pdf("Change_T2_CD8_subsets_by_response_FCx2.pdf",width=7,height=5)
plot(plotje)
dev.off()

boxplot_data <- cell_counts_relative
fortessa<-read_xlsx("Data_incl_metadata_all_batches_relative.xlsx")
fortessa<-fortessa%>%dplyr::rename(
  subject=Sample,
  timepoint=Timepoint
)
boxplot_data<-merge(boxplot_data,fortessa,by=c("subject","timepoint"),all=T)
boxplot_data <- boxplot_data %>% dplyr::filter(Response!="HD" & (timepoint==1 | timepoint==2))

boxplot_data$Response <- factor(boxplot_data$Response,levels=c("NR","R"))
median_df <- boxplot_data %>% 
  group_by(timepoint,Response) %>% 
  summarize(median_value = median(CD8_Tem_3),.groups="drop")
boxplot1 <- ggplot(boxplot_data,aes(timepoint,100*CD8_Tem_3,group=subject,color=Response,fill=Response))+
  geom_point(shape=21)+
  geom_line(alpha=0.3)+
  geom_line(data=median_df,aes(y=100*median_value,group=Response),linewidth=1,linetype="dashed")+
  scale_color_manual(labels=c("Non-responder","Responder"),
                     values=c("#440154","#5ec962"))+
  scale_fill_manual(labels=c("Non-responder","Responder"),
                     values=c("#440154","#5ec962"))+
  labs(fill=NULL,group=NULL,color=NULL)+
  ylab("% CD8 Tem 3 of all CD8 T cells")+
  xlab("")+
  scale_x_discrete(labels=c("1" = "Pre steroids", "2" = "On steroids"))+
  theme_classic()
pdf("Boxplot_change_CD8_tem3.pdf",width=3.5,height=2.8)
plot(boxplot1)
dev.off()

#Figure 3 and associated supplementary Figure 
setwd("~")
setwd(".../")

#Create gating strategy plots
get_density<-function(x,y,...){
  dens<-MASS::kde2d(x,y,...)
  ix<-findInterval(x,dens$x)
  iy<-findInterval(y,dens$y)
  ii<-cbind(ix,iy)
  return(dens$z[ii])
}

healthy_CD4 <- read.FCS("MDDXXXXXXX_HD_CD4.fcs",truncate_max_range = FALSE)
healthy_CD4_unbound <- read.FCS("MDDXXXXXXX_HD_CD4_unbound.fcs",truncate_max_range = FALSE)
CD4_unbound <- read.FCS("UNI-XXX_T1_CD4_unbound_Specimen_003_C7_C07_018_CD4 T cells_PD-1+, unbound by pembro_nivo.fcs",truncate_max_range = FALSE)
CD4_T <- read.FCS("UNI-XXX_T1_CD4Tcells_Specimen_003_C7_C07_018_CD4 T cells.fcs",truncate_max_range = FALSE)
live_T_2 <- read.FCS("UNI-XXX_T1_live_Tcells_Specimen_003_C7_C07_018_Alive CD3.fcs",truncate_max_range = FALSE)

full <- read.FCS("UNI-YYY_T1.fcs",truncate_max_range = FALSE)
lymphos <- read.FCS("UNI-YYY_T1_lymphocytes_Specimen_001_B10_B10_009_Lymphocytes.fcs",truncate_max_range = FALSE)
singlets <- read.FCS("UNI-YYY_T1_singlets_Specimen_001_B10_B10_009_Single Cells.fcs",truncate_max_range = FALSE)
live_T_1 <- read.FCS("UNI-YYY_T1_live_Tcells_Specimen_001_B10_B10_009_Alive CD3.fcs",truncate_max_range = FALSE)

CD4_unbound <- as.data.frame(CD4_unbound@exprs)
CD4_T <- as.data.frame(CD4_T@exprs)
live_T_2 <- as.data.frame(live_T_2@exprs)

full <- as.data.frame(full@exprs)
lymphos <- as.data.frame(lymphos@exprs)
singlets <- as.data.frame(singlets@exprs)
live_T_1 <- as.data.frame(live_T_1@exprs)
healthy <- as.data.frame(healthy@exprs)
healthy_CD4 <- as.data.frame(healthy_CD4@exprs)
healthy_CD4_unbound <- as.data.frame(healthy_CD4_unbound@exprs)

#GATE LYMPHOS
full$density<-get_density(full$`FSC-A`,full$`SSC-A`,n=200)
lymphos_jittered <- lymphos
lymphos_jittered$`FSC-A` <- jitter(lymphos_jittered$`FSC-A`, amount = 1)
lymphos_jittered$`SSC-A` <- jitter(lymphos_jittered$`SSC-A`, amount = 1)
lymphos_jittered <- lymphos_jittered[sample(nrow(lymphos_jittered), 0.01*nrow(lymphos_jittered)), ]
gate <- ashape(lymphos_jittered$`FSC-A`,lymphos_jittered$`SSC-A`, alpha = 100000)
gate_edges <- as.data.frame(gate$edges)

nrow(lymphos)/nrow(full)
plot <- ggplot(full,aes(`FSC-A`,`SSC-A`,color=density))+
  geom_point(shape=46,size=0.1)+
  geom_segment(data = gate_edges, aes(x = x1, y = y1, xend = x2, yend = y2), color = "#bc3754",size=2)+
  scale_color_viridis(option="inferno")+
  ggtitle("Lymphocytes (61.7%)")+
  theme_classic()+
  theme(legend.position = "none")+
  xlim(0,240000)+
  ylim(0,240000)
pdf("lympho_gate.pdf",width=5,height=5)
plot(plot)
dev.off()

#GATE SINGLETS
lymphos$density<-get_density(lymphos$`FSC-A`,lymphos$`FSC-H`,n=200)
singlets_jittered <- singlets
singlets_jittered$`FSC-A` <- jitter(singlets_jittered$`FSC-A`, amount = 1)
singlets_jittered$`SSC-A` <- jitter(singlets_jittered$`FSC-H`, amount = 1)
singlets_jittered <- singlets_jittered[sample(nrow(singlets_jittered), 0.1*nrow(singlets_jittered)), ]
gate <- ashape(singlets_jittered$`FSC-A`,singlets_jittered$`SSC-A`, alpha = 100000)
gate_edges <- as.data.frame(gate$edges)

nrow(singlets)/nrow(lymphos)
plot <- ggplot(lymphos,aes(`FSC-A`,`FSC-H`,color=density))+
  geom_point(shape=46,size=0.1)+
  geom_segment(data = gate_edges, aes(x = x1, y = y1, xend = x2, yend = y2), color = "#bc3754",size=2)+
  scale_color_viridis(option="inferno")+
  ggtitle("Singlets (98.7%)")+
  theme_classic()+
  theme(legend.position = "none")+
  xlim(0,150000)+
  ylim(0,150000)
pdf("singlet_gate.pdf",width=5,height=5)
plot(plot)
dev.off()

#GATE LIVE CELLS (spectral only)
singlets <- read.FCS("UNI-YYY_T1_singlets_Specimen_001_B10_B10_009_Single Cells.fcs",truncate_max_range = FALSE)
singlets <- as.data.frame(singlets@exprs)
hist(singlets$`FJComp-AmCyan-A`,breaks=100)
min_y <- min(singlets$`FJComp-AmCyan-A`)
singlets$`FJComp-AmCyan-A` <- log10((singlets$`FJComp-AmCyan-A`-min_y)+1)
singlets$`FJComp-AmCyan-A` <- ifelse(singlets$`FJComp-AmCyan-A`<4.2,4.2,singlets$`FJComp-AmCyan-A`)
singlets$density <- get_density(singlets$`FSC-A`,singlets$`FJComp-AmCyan-A`,n=300)
nrow(singlets[singlets$`FJComp-AmCyan-A`<4.305,])/nrow(singlets)
plot <- ggplot(singlets,aes(`FSC-A`,`FJComp-AmCyan-A`,color=density))+
  geom_point(shape=46,size=0.1)+
  geom_rect(aes(xmin = 25100, xmax = 124900, ymin = 4.25, ymax = 4.305),fill=NA, color = "#bc3754",size=2)+
  scale_color_viridis(option="inferno")+
  ggtitle("Live cells (88.6%)")+
  theme_classic()+
  theme(legend.position = "none")+
  xlim(25000,125000)+
  ylim(4.2,4.45)+
  xlab("FSC-A")+
  ylab("eFluor506 live/dead")
pdf("live_all_spectral_gate.pdf",width=5,height=5)
plot(plot)
dev.off()

#GATE LIVE T CELLS
singlets <- read.FCS("UNI-YYY_T1_singlets_Specimen_001_B10_B10_009_Single Cells.fcs",truncate_max_range = FALSE)
singlets <- as.data.frame(singlets@exprs)
hist(singlets$`FJComp-Alexa Fluor 700-A`,breaks=100)
hist(singlets$`FJComp-AmCyan-A`,breaks=100)

min_x <- min(singlets$`FJComp-Alexa Fluor 700-A`)
min_y <- min(singlets$`FJComp-AmCyan-A`)

singlets$`FJComp-Alexa Fluor 700-A` <- log10((singlets$`FJComp-Alexa Fluor 700-A`-min_x)+1)
singlets$`FJComp-AmCyan-A` <- log10((singlets$`FJComp-AmCyan-A`-min_y)+1)

singlets$`FJComp-Alexa Fluor 700-A` <- ifelse(singlets$`FJComp-Alexa Fluor 700-A`<2.9,2.9,singlets$`FJComp-Alexa Fluor 700-A`)
singlets$`FJComp-AmCyan-A` <- ifelse(singlets$`FJComp-AmCyan-A`<4.2,4.2,singlets$`FJComp-AmCyan-A`)
singlets$density <- get_density(singlets$`FJComp-Alexa Fluor 700-A`,singlets$`FJComp-AmCyan-A`,n=300)
nrow(live_T_1)/nrow(singlets)
plot <- ggplot(singlets,aes(`FJComp-Alexa Fluor 700-A`,`FJComp-AmCyan-A`,color=density))+
  geom_point(shape=46,size=0.1)+
  geom_rect(aes(xmin = 3.2, xmax = 4, ymin = 4.25, ymax = 4.31),fill=NA, color = "#bc3754",size=2)+
  scale_color_viridis(option="inferno")+
  ggtitle("Live T cells (46.3%)")+
  theme_classic()+
  theme(legend.position = "none")+
  xlim(2.9,4.1)+
  ylim(4.2,4.45)+
  xlab("AF700-CD3")+
  ylab("eFluor506 live/dead")
pdf("live_T_1_gate.pdf",width=5,height=5)
plot(plot)
dev.off()

#GATE CD4 T CELLS
min_x <- min(live_T_2$`FJComp-PE-Cy7-A`)
min_y <- min(live_T_2$`FJComp-BV-785-A`)

live_T_2$`FJComp-PE-Cy7-A` <- log10((live_T_2$`FJComp-PE-Cy7-A`-min_x)+1)
live_T_2$`FJComp-BV-785-A` <- log10((live_T_2$`FJComp-BV-785-A`-min_y)+1)

live_T_2$density <- get_density(live_T_2$`FJComp-PE-Cy7-A`,live_T_2$`FJComp-BV-785-A`,n=300)
nrow(CD4_T)/nrow(live_T_2)
plot <- ggplot(live_T_2,aes(`FJComp-PE-Cy7-A`,`FJComp-BV-785-A`,color=density))+
  geom_point(shape=46,size=0.1)+
  geom_rect(aes(xmin = 2.3, xmax = 3.3, ymin = 3.9, ymax = 5),fill=NA, color = "#bc3754",size=2)+
  scale_color_viridis(option="inferno")+
  ggtitle("CD4+ T cells (71.0%)")+
  theme_classic()+
  theme(legend.position = "none")+
  xlim(2.2,4.5)+
  ylim(3,5)+
  xlab("PE-Cy7-CD8")+
  ylab("BV785-CD4")
pdf("CD4_gate.pdf",width=5,height=5)
plot(plot)
dev.off()

#GATE PD-1 in patient
min_x <- min(min(CD4_T$`FJComp-BV-711-A`),min(healthy_CD4$`FJComp-BV-711-A`))
min_y <- min(min(CD4_T$`FJComp-PE-A`),min(healthy_CD4$`FJComp-PE-A`))

CD4_T$`FJComp-BV-711-A` <- log10((CD4_T$`FJComp-BV-711-A`-min_x)+1)
CD4_T$`FJComp-PE-A` <- log10((CD4_T$`FJComp-PE-A`-min_y)+1)

CD4_T$density <- get_density(CD4_T$`FJComp-BV-711-A`,CD4_T$`FJComp-PE-A`,n=600)

unbound_pt <- CD4_T %>% dplyr::filter(`FJComp-BV-711-A`>=3.25 & `FJComp-PE-A`<3.13)
part_pt <- CD4_T %>% dplyr::filter(`FJComp-BV-711-A`>=3.25 & `FJComp-PE-A`>=3.13)
comp_pt <- CD4_T %>% dplyr::filter(`FJComp-BV-711-A`<3.25 & `FJComp-PE-A`>=3.13)
negative_pt <- CD4_T %>% dplyr::filter(`FJComp-BV-711-A`<3.25 & `FJComp-PE-A`<3.13)

nrow(unbound_pt)/nrow(CD4_T)*100
nrow(part_pt)/nrow(CD4_T)*100
nrow(comp_pt)/nrow(CD4_T)*100
nrow(negative_pt)/nrow(CD4_T)*100

plot <- ggplot(CD4_T,aes(`FJComp-BV-711-A`,`FJComp-PE-A`,color=density))+
  geom_point(shape=46,size=0.1)+
  geom_hline(yintercept=3.13,color = "#bc3754",size=2)+
  geom_vline(xintercept=3.25,color = "#bc3754",size=2)+
  scale_color_viridis(option="inferno")+
  ggtitle("Unbound PD-1+ CD4+ T cells (4.4%)")+
  theme_classic()+
  theme(legend.position = "none")+
  xlim(2.7,3.8)+
  ylim(2.75,3.5)+
  xlab("BV711-PD-1")+
  ylab("PE-IgG4")
pdf("PD_1_gate_pt.pdf",width=5,height=5)
plot(plot)
dev.off()

#GATE PD-1 in healthy
healthy_CD4$`FJComp-BV-711-A` <- log10((healthy_CD4$`FJComp-BV-711-A`-min_x)+1)
healthy_CD4$`FJComp-PE-A` <- log10((healthy_CD4$`FJComp-PE-A`-min_y)+1)

healthy_CD4$density <- get_density(healthy_CD4$`FJComp-BV-711-A`,healthy_CD4$`FJComp-PE-A`,n=600)

unbound_hd <- healthy_CD4 %>% dplyr::filter(`FJComp-BV-711-A`>=3.25 & `FJComp-PE-A`<3.13)
part_hd <- healthy_CD4 %>% dplyr::filter(`FJComp-BV-711-A`>=3.25 & `FJComp-PE-A`>=3.13)
comp_hd <- healthy_CD4 %>% dplyr::filter(`FJComp-BV-711-A`<3.25 & `FJComp-PE-A`>=3.13)
negative_hd <- healthy_CD4 %>% dplyr::filter(`FJComp-BV-711-A`<3.25 & `FJComp-PE-A`<3.13)

nrow(unbound_hd)/nrow(healthy_CD4)*100
nrow(part_hd)/nrow(healthy_CD4)*100
nrow(comp_hd)/nrow(healthy_CD4)*100
nrow(negative_hd)/nrow(healthy_CD4)*100


plot <- ggplot(healthy_CD4,aes(`FJComp-BV-711-A`,`FJComp-PE-A`,color=density))+
  geom_point(shape=46,size=0.1)+
  geom_hline(yintercept=3.13,color = "#bc3754",size=2)+
  geom_vline(xintercept=3.25,color = "#bc3754",size=2)+
  scale_color_viridis(option="inferno")+
  ggtitle("Unbound PD-1+ CD4+ T cells (16.1%)")+
  theme_classic()+
  theme(legend.position = "none")+
  xlim(2.7,3.8)+
  ylim(2.75,3.5)+
  xlab("BV711-PD-1")+
  ylab("PE-IgG4")
pdf("PD_1_gate_healthy.pdf",width=5,height=5)
plot(plot)
dev.off()

#GATE Ki67 in patient
min_x <- min(CD4_unbound$`FJComp-FITC-A`)
CD4_unbound$`FJComp-FITC-A` <- log10((CD4_unbound$`FJComp-FITC-A`-min_x)+1)
CD4_unbound$density <- get_density(CD4_unbound$`FJComp-FITC-A`,CD4_unbound$`FSC-A`,n=600)
CD4_prolif <- CD4_unbound %>% dplyr::filter(`FJComp-FITC-A`>2.95)
nrow(CD4_prolif)/nrow(CD4_unbound)
plot <- ggplot(CD4_unbound,aes(`FJComp-FITC-A`,`FSC-A`,color=density))+
  geom_point(shape=46,size=0.1)+
  geom_rect(aes(xmin = 2.95, xmax = 4.05, ymin = 25000, ymax = 125000),fill=NA, color = "#bc3754",size=2)+
  scale_color_viridis(option="inferno")+
  ggtitle("Ki67+ unbound PD-1+ CD4+ T cells (46.3%)")+
  theme_classic()+
  theme(legend.position = "none")+
  xlim(2.7,4.1)+
  xlab("FITC-Ki67")+
  ylab("FSC-A")
pdf("Ki67_gate_pt.pdf",width=5,height=5)
plot(plot)
dev.off()

#General spectral panel pre-gating
all_events <- read.FCS("spectral_1_all_events_A6 Well_006.fcs",truncate_max_range = FALSE)
noDebris <- read.FCS("spectral_2_noDebris_A6 Well_006_All cells.fcs",truncate_max_range = FALSE)
alive <- read.FCS("spectral_3_alive_A6 Well_006_Alive cells.fcs",truncate_max_range = FALSE)
singlets <- read.FCS("spectral_4_singlets_A6 Well_006_Single Cells.fcs",truncate_max_range = FALSE)

all_events <- as.data.frame(all_events@exprs)
noDebris <- as.data.frame(noDebris@exprs)
alive <- as.data.frame(alive@exprs)
singlets <- as.data.frame(singlets@exprs)

#GATE CELLS
all_events$density<-get_density(all_events$`FSC-A`,all_events$`SSC-A`,n=200)
noDebris_jittered <- noDebris
noDebris_jittered$`FSC-A` <- jitter(noDebris_jittered$`FSC-A`, amount = 1)
noDebris_jittered$`SSC-A` <- jitter(noDebris_jittered$`SSC-A`, amount = 1)
noDebris_jittered <- noDebris_jittered[sample(nrow(noDebris_jittered), 0.01*nrow(noDebris_jittered)), ]
gate <- ashape(noDebris_jittered$`FSC-A`,noDebris_jittered$`SSC-A`, alpha = 10^15)
gate_edges <- as.data.frame(gate$edges)

nrow(noDebris)/nrow(all_events)
plot <- ggplot(all_events,aes(`FSC-A`,`SSC-A`,color=density))+
  geom_point(shape=46,size=0.1)+
  geom_segment(data = gate_edges, aes(x = x1, y = y1, xend = x2, yend = y2), color = "#bc3754",size=2)+
  scale_color_viridis(option="inferno")+
  ggtitle("No debris (71.6%)")+
  theme_classic()+
  theme(legend.position = "none")
#xlim(0,240000)+
#ylim(0,240000)
pdf("spectral_debris_gate.pdf",width=5,height=5)
plot(plot)
dev.off()

#GATE ALIVE cells
nrow(alive)/nrow(noDebris)
noDebris$`eFluor 506-A` <- -(noDebris$`eFluor 506-A`)
noDebris$`eFluor 506-A` <- noDebris$`eFluor 506-A` - min(noDebris$`eFluor 506-A`)
noDebris$`eFluor 506-A` <- ifelse(noDebris$`eFluor 506-A`<quantile(noDebris$`eFluor 506-A`,probs=0.01),quantile(noDebris$`eFluor 506-A`,probs=0.01),
                              ifelse(noDebris$`eFluor 506-A`>quantile(noDebris$`eFluor 506-A`,probs=0.99),quantile(noDebris$`eFluor 506-A`,probs=0.99),noDebris$`eFluor 506-A`))
noDebris$`eFluor 506-A` <- -(noDebris$`eFluor 506-A`)
noDebris$`eFluor 506-A` <- noDebris$`eFluor 506-A` - min(noDebris$`eFluor 506-A`)
noDebris$`eFluor 506-A` <- log10(noDebris$`eFluor 506-A`+0.01)
noDebris <- noDebris %>% dplyr::filter(`eFluor 506-A`>2 & `eFluor 506-A`<4.7) 
noDebris$density<-get_density(noDebris$`FSC-A`,noDebris$`eFluor 506-A`,n=200)

plot <- ggplot(noDebris,aes(`FSC-A`,`eFluor 506-A`,color=density))+
  geom_point(shape=46,size=0.1)+
  geom_rect(aes(xmin = 1, xmax = 4000000, ymin = 2.5, ymax = 4.15),fill=NA, color = "#bc3754",size=2)+
  scale_color_viridis(option="inferno")+
  ggtitle("Alive cells (74.8%)")+
  ylab("eFluor506 Live/dead")+
  theme_classic()+
  theme(legend.position = "none")
pdf("spectral_alive_gate.pdf",width=5,height=5)
plot(plot)
dev.off()

#GATE SINGLETS
alive$density<-get_density(alive$`FSC-A`,alive$`FSC-H`,n=200)
singlets_jittered <- singlets
singlets_jittered$`FSC-A` <- jitter(singlets_jittered$`FSC-A`, amount = 1)
singlets_jittered$`FSC-H` <- jitter(singlets_jittered$`FSC-H`, amount = 1)
singlets_jittered <- singlets_jittered[sample(nrow(singlets_jittered), 0.1*nrow(singlets_jittered)), ]
gate <- ashape(singlets_jittered$`FSC-A`,singlets_jittered$`FSC-H`, alpha = 10^12)
gate_edges <- as.data.frame(gate$edges)

nrow(singlets)/nrow(alive)
plot <- ggplot(alive,aes(`FSC-A`,`FSC-H`,color=density))+
  geom_point(shape=46,size=0.1)+
  geom_segment(data = gate_edges, aes(x = x1, y = y1, xend = x2, yend = y2), color = "#bc3754",size=2)+
  scale_color_viridis(option="inferno")+
  ggtitle("Singlets (97.2%)")+
  theme_classic()+
  theme(legend.position = "none")
pdf("spectral_singlet_gate.pdf",width=5,height=5)
plot(plot)
dev.off()

#Plots granzyme B CD4 / CD8 T cells
#General spectral panel pre-gating
setwd("~")
setwd(".../")
CD4_T2 <- read.FCS("UNI-ZZZ_T2_CD4_F2 Well_013_CD4 T cells.fcs",truncate_max_range = FALSE)
CD8_T2 <- read.FCS("UNI-ZZZ_T2_CD8_F2 Well_013_CD8 T cells.fcs",truncate_max_range = FALSE)
CD4_T3 <- read.FCS("UNI-ZZZ_T3_CD4_F3 Well_014_CD4 T cells.fcs",truncate_max_range = FALSE)
CD8_T3 <- read.FCS("UNI-ZZZ_T3_CD8_F3 Well_014_CD8 T cells.fcs",truncate_max_range = FALSE)

CD4_T2 <- as.data.frame(CD4_T2@exprs)
CD8_T2 <- as.data.frame(CD8_T2@exprs)
CD4_T3 <- as.data.frame(CD4_T3@exprs)
CD8_T3 <- as.data.frame(CD8_T3@exprs)

CD4_T2 <- CD4_T2 %>% dplyr::select(`FSC-A`,`Pacific Blue-A`)
CD8_T2 <- CD8_T2 %>% dplyr::select(`FSC-A`,`Pacific Blue-A`)
CD4_T3 <- CD4_T3 %>% dplyr::select(`FSC-A`,`Pacific Blue-A`)
CD8_T3 <- CD8_T3 %>% dplyr::select(`FSC-A`,`Pacific Blue-A`)

CD4_T2$file <- "CD4_T2" 
CD8_T2$file <- "CD8_T2" 
CD4_T3$file <- "CD4_T3" 
CD8_T3$file <- "CD8_T3"

merged_T2_T3 <- rbind(CD4_T2,CD8_T2)
merged_T2_T3 <- rbind(merged_T2_T3,CD4_T3)
merged_T2_T3 <- rbind(merged_T2_T3,CD8_T3)

#GzmB in CD4 T cells at T2 and T3
hist(merged_T2_T3$`Pacific Blue-A`,breaks=100)
merged_T2_T3$`Pacific Blue-A` <- ifelse(merged_T2_T3$`Pacific Blue-A`<quantile(merged_T2_T3$`Pacific Blue-A`,probs=0.01),quantile(merged_T2_T3$`Pacific Blue-A`,probs=0.01),merged_T2_T3$`Pacific Blue-A`)

merged_T2_T3$`Pacific Blue-A` <- merged_T2_T3$`Pacific Blue-A`-min(merged_T2_T3$`Pacific Blue-A`)
merged_T2_T3$`Pacific Blue-A` <- log10(merged_T2_T3$`Pacific Blue-A`+0.01)
hist(merged_T2_T3$`Pacific Blue-A`,breaks=100)

merged_T2_T3$density<-get_density(merged_T2_T3$`FSC-A`,merged_T2_T3$`Pacific Blue-A`,n=200)
merged_T2_T3 <- merged_T2_T3 %>% dplyr::filter(`Pacific Blue-A`>1.5)
merged_T2_T3$Tcell <- ifelse(grepl("CD4",merged_T2_T3$file),"CD4","CD8")
merged_T2_T3$time <- ifelse(grepl("T2",merged_T2_T3$file),"T1","T2")

plot <- ggplot(merged_T2_T3,aes(`Pacific Blue-A`,`FSC-A`,color=density))+
  geom_point(shape=46,size=0.1)+
  scale_color_viridis(option="inferno")+
  theme_classic()+
  xlab("PB granzyme B")+
  theme(legend.position = "none")+
  facet_wrap(~file,ncol=4,labeller = labeller(file = 
                                         c("CD4_T2" = "Timepoint 1, CD4+ T cells",
                                           "CD4_T3" = "Timepoint 2, CD4+ T cells",
                                           "CD8_T2" = "Timepoint 1, CD8+ T cells",
                                           "CD8_T3" = "Timepoint 2, CD8+ T cells")))
pdf("granzyme_B_example.pdf",width=11,height=3)
plot(plot)
dev.off()

#Analyses 
data <- read_excel("Data_incl_metadata_all_batches_relative.xlsx")

#only select patients treated with anti-PD-1 +/- anti-CTLA-4 and only T1 and T2
subset <- data %>% dplyr::filter((Treatment=="ipi-nivo" | Treatment=="anti-PD-1") & (Timepoint==1 | Timepoint==2))
subset2 <- subset %>% dplyr::filter(Timepoint == 1)
subset3 <- subset2 %>% dplyr::select(c(Sample,Response,Ki67pos_of_PD1neg_CD4,Ki67pos_of_PD1pos_unbound_CD4,Ki67pos_of_PD1pos_bound_CD4))
subset3 <- melt(subset3,id.vars=c("Sample","Response"))
subset4 <- subset2 %>% dplyr::select(c(Sample,Response,Ki67pos_of_PD1neg_CD8,Ki67pos_of_PD1pos_unbound_CD8,Ki67pos_of_PD1pos_bound_CD8))
subset4 <- melt(subset4,id.vars=c("Sample","Response"))

CD4_bound <- ggplot(subset,aes(Time_blood_rel,
                               100*(as.numeric(`CD4_sat_PD-1`)+as.numeric(`CD4_unsat_PD-1`))/(as.numeric(`CD4_sat_PD-1`)+as.numeric(`CD4_unsat_PD-1`)+as.numeric(`CD4_unbound_PD-1`)),
                               colour=factor(Response)))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey")+
  geom_point()+
  geom_line(aes(group=Sample))+
  theme_bw()+
  scale_colour_manual(values=c("#440154","#5ec962"),labels=c("Non-responder","Responder"))+
  labs(title="Relative to last dose anti-PD-1 treatment",
       color=NULL)+
  ylim(0,100)+
  xlim((-2),250)+
  ylab("% PD-1+ CD4+ T cells completely/partially bound")+
  xlab("Time since last anti-PD-1 dose (days)")+
  theme_classic()
pdf("CD4_bound_rel_last_ICI.pdf",width=5,height=4)
plot(CD4_bound)
dev.off()

CD4_bound <- ggplot(subset,aes(Time_blood_rel-Time_start_steroids_rel,
                               100*(as.numeric(`CD4_sat_PD-1`)+as.numeric(`CD4_unsat_PD-1`))/(as.numeric(`CD4_sat_PD-1`)+as.numeric(`CD4_unsat_PD-1`)+as.numeric(`CD4_unbound_PD-1`)),
                               colour=factor(Response)))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey")+
  geom_point()+
  geom_point()+
  geom_line(aes(group=Sample))+
  theme_bw()+
  scale_colour_manual(values=c("#440154","#5ec962"),labels=c("Non-responder","Responder"))+
  labs(title="Fraction of PD-1-expressing CD4 T cells completely\nor partially bound by anti-PD-1",
       color=NULL)+
  ylim(0,100)+
  ylab("% PD-1+ CD4+ T cells completely/partially bound")+
  xlab("Time since prednisone initiation (days)")+
  theme_classic()
pdf("CD4_bound_rel_steroids.pdf",width=5,height=5)
plot(CD4_bound)
dev.off()

CD4_bound_new <- ggplot(subset,aes(Time_blood_rel,
                                   100*(as.numeric(`CD4_sat_PD-1`)+as.numeric(`CD4_unsat_PD-1`))/(as.numeric(`CD4_sat_PD-1`)+as.numeric(`CD4_unsat_PD-1`)+as.numeric(`CD4_unbound_PD-1`)),
                                   colour=factor(Treatment)))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey")+
  geom_point()+
  geom_point()+
  geom_line(aes(group=Sample))+
  theme_bw()+
  scale_colour_manual(values=c("#F98E09","#BC3754"),labels=c("anti-PD-1-based","ipilimumab + nivolumab"))+
  labs(title="Fraction of PD-1-expressing CD4 T cells completely\nor partially bound by anti-PD-1",
       color=NULL)+
  ylim(0,100)+
  xlim((-2),250)+
  ylab("% PD-1+ CD4+ T cells completely/partially bound")+
  xlab("Time since last anti-PD-1 dose (days)")+
  theme_classic()
pdf("CD4_bound_rel_last_ICI_treatment.pdf",width=5,height=5)
plot(CD4_bound_new)
dev.off()

subset3 <- subset2
subset3 <- subset3 %>% dplyr::filter(`CD4_unsat_PD-1`>0.01) 

Ki67_corr_CD4 <- ggplot(subset3,aes(100*Ki67pos_of_CD4,
                                    100*(as.numeric(`CD4_sat_PD-1`)+as.numeric(`CD4_unsat_PD-1`))/(as.numeric(`CD4_sat_PD-1`)+as.numeric(`CD4_unsat_PD-1`)+as.numeric(`CD4_unbound_PD-1`))))+
  geom_point(aes(color=factor(Response)))+
  stat_cor(method="spearman",label.x=25,label.y=90)+
  theme_bw()+
  scale_colour_manual(values=c("#440154","#5ec962"),labels=c("Non-responder","Responder"))+
  labs(title="Fraction of PD-1-expressing CD4 T cells completely\nor partially bound by anti-PD-1",
       color=NULL)+
  ylim(0,100)+
  ylab("% PD-1+ CD4+ T cells completely/partially bound")+
  xlab("% Ki67+ CD4+ T cells")+
  theme_classic()
pdf("Ki67_corr_CD4.pdf",width=5,height=5)
plot(Ki67_corr_CD4)
dev.off()

CD8_bound <- ggplot(subset,aes(Time_blood_rel,
                               100*(as.numeric(`CD8_sat_PD-1`)+as.numeric(`CD8_unsat_PD-1`))/(as.numeric(`CD8_sat_PD-1`)+as.numeric(`CD8_unsat_PD-1`)+as.numeric(`CD8_unbound_PD-1`)),
                               colour=factor(Response)))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey")+
  geom_point()+
  geom_point()+
  geom_line(aes(group=Sample))+
  theme_bw()+
  scale_colour_manual(values=c("#440154","#5ec962"),labels=c("Non-responder","Responder"))+
  labs(title="Fraction of PD-1-expressing CD8 T cells completely\nor partially bound by anti-PD-1",
       color=NULL)+
  ylim(0,100)+
  xlim((-2),250)+
  ylab("% PD-1+ CD8+ T cells completely/partially bound")+
  xlab("Time since last anti-PD-1 dose (days)")+
  theme_classic()
pdf("CD8_bound_rel_last_ICI.pdf",width=5,height=5)
plot(CD8_bound)
dev.off()

CD8_bound <- ggplot(subset,aes(Time_blood_rel,
                               100*(as.numeric(`CD8_sat_PD-1`)+as.numeric(`CD8_unsat_PD-1`))/(as.numeric(`CD8_sat_PD-1`)+as.numeric(`CD8_unsat_PD-1`)+as.numeric(`CD8_unbound_PD-1`)),
                               colour=factor(Treatment)))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey")+
  geom_point()+
  geom_point()+
  geom_line(aes(group=Sample))+
  theme_bw()+
  scale_colour_manual(values=c("#F98E09","#BC3754"),labels=c("anti-PD-1","ipilimumab + nivolumab"))+
  labs(title="Fraction of PD-1-expressing CD8 T cells completely\nor partially bound by anti-PD-1",
       color=NULL)+
  ylim(0,100)+
  xlim((-2),250)+
  ylab("% PD-1+ CD8+ T cells completely/partially bound")+
  xlab("Time since last anti-PD-1 dose (days)")+
  theme_classic()
pdf("CD8_bound_rel_last_ICI_treatment.pdf",width=5,height=5)
plot(CD8_bound)
dev.off()

CD8_bound <- ggplot(subset,aes(Time_blood_rel-Time_start_steroids_rel,
                               100*(as.numeric(`CD8_sat_PD-1`)+as.numeric(`CD8_unsat_PD-1`))/(as.numeric(`CD8_sat_PD-1`)+as.numeric(`CD8_unsat_PD-1`)+as.numeric(`CD8_unbound_PD-1`)),
                               colour=factor(Response)))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey")+
  geom_point()+
  geom_point()+
  geom_line(aes(group=Sample))+
  theme_bw()+
  scale_colour_manual(values=c("#440154","#5ec962"),labels=c("Non-responder","Responder"))+
  labs(title="Fraction of PD-1-expressing CD8 T cells completely\nor partially bound by anti-PD-1",
       color=NULL)+
  ylim(0,100)+
  ylab("% PD-1+ CD8+ T cells completely/partially bound")+
  xlab("Time since prednisone initiation (days)")+
  theme_classic()
pdf("CD8_bound_rel_steroids.pdf",width=5,height=5)
plot(CD8_bound)
dev.off()

Ki67_corr_CD8 <- ggplot(subset3,aes(100*Ki67pos_of_CD8,
                                    100*(as.numeric(`CD8_sat_PD-1`)+as.numeric(`CD8_unsat_PD-1`))/(as.numeric(`CD8_sat_PD-1`)+as.numeric(`CD8_unsat_PD-1`)+as.numeric(`CD8_unbound_PD-1`))))+
  geom_point(aes(color=factor(Response)))+
  stat_cor(method="spearman",label.x=25,label.y=90)+
  scale_colour_manual(values=c("#440154","#5ec962"),labels=c("Non-responder","Responder"))+
  labs(title="Fraction of PD-1-expressing CD8 T cells completely\nor partially bound by anti-PD-1")+
  ylim(0,100)+
  ylab("% PD-1+ CD8+ T cells completely/partially bound")+
  xlab("Fraction Ki67+ CD8 T cells")+
  theme_classic()
pdf("Ki67_corr_CD8.pdf",width=5,height=5)
plot(Ki67_corr_CD8)
dev.off()

#ImageStream data
setwd("~")
setwd(".../")
data <- read_xlsx("ImageStream_data.xlsx")
median_data <- data[,c(1,5,10,13,16,19,22,25)]
sd_data <- data[,c(1,5,12,15,18,21,24,27)]
median_long <- median_data %>%
  pivot_longer(
    cols = -c(Sample, Fraction_bound),  # Columns to keep as keys
    names_to = c("Marker", "Time"),      # Name for new columns
    names_pattern = "(.*)_Median_Internalization_(.*)",  # Regular expression to split names
    values_to = "Median_Internalization"  # Name for values column
  )

median_long <- median_long %>% dplyr::filter(Marker=="IgG4")
median_long$Time <- sub("min","",median_long$Time)

sd_long <- sd_data %>%
  pivot_longer(
    cols = -c(Sample, Fraction_bound),  # Columns to keep as keys
    names_to = c("Marker", "Time"),      # Name for new columns
    names_pattern = "(.*)_SD_Internalization_(.*)",  # Regular expression to split names
    values_to = "SD_Internalization"  # Name for values column
  )

sd_long <- sd_long %>% dplyr::filter(Marker=="IgG4")
sd_long$Time <- sub("min","",median_long$Time)

merged <- merge(median_long,sd_long)

Figure <- ggplot(data = merged, aes(x = Time, y = Median_Internalization, group = Sample)) +
  geom_point(shape = 21,position=position_dodge(width=0.2)) +
  geom_line(position=position_dodge(width=0.2)) +
  geom_ribbon(aes(ymin = Median_Internalization - SD_Internalization, ymax = Median_Internalization + SD_Internalization, 
                  color = "grey"), alpha = 0.1,color=NA,
                position=position_dodge(width=0.2)) +
  facet_wrap(~ as.character(paste0(100*round(Fraction_bound,2),"% of PD-1 bound")),ncol=4) +
  labs(color="Fraction PD-1\n(in)completely bound\nby anti-PD-1 drug")+
  xlab("Time (min) longitudinal ImageStream measurements")+
  ylab("Internalization IgG4 (ICI-bound PD-1)")+
  theme_classic()
pdf("ImageStream.pdf",height=3,width=8)
plot(Figure)
dev.off()