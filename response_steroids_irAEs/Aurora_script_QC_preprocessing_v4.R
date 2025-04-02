#Script based on: https://github.com/HdBraanker/Spectral_Flow_Workflow/blob/main/script.Rmd

library(flowCore)
library(flowViz)
library(flowVS)
library(flowAI)
library(PeacoQC)
library(CATALYST)
library(devtools)
install_github("saeyslab/CytoNorm",force=TRUE)
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
#library(scCustomize)

#Start script
set.seed(123)

## If already done with QC, lookup "START here" to continue script ##
setwd("~")
setwd(".../Analysis/")
fcs.dir <- file.path(getwd())

#Load FCS files batch 1 and add batch number as metadata
fcs.data_batch1 <- read.flowSet(path=paste0(fcs.dir,"/viable_singlets/viable_singlets_batch1/"),pattern="*.fcs",transformation=FALSE,truncate_max_range=FALSE)
for(i in 1:length(fcs.data_batch1@frames)){
  fcs.data_batch1[[i]]@description$batch <- 1
}

#Load FCS files batch 2 and add batch number as metadata
fcs.data_batch2 <- read.flowSet(path=paste0(fcs.dir,"/viable_singlets/viable_singlets_batch2/"),pattern="*.fcs",transformation=FALSE,truncate_max_range=FALSE)
for(i in 1:length(fcs.data_batch2@frames)){
  fcs.data_batch2[[i]]@description$batch <- 2
}

#Load FCS files batch 3 and add batch number as metadata
fcs.data_batch3 <- read.flowSet(path=paste0(fcs.dir,"/viable_singlets/viable_singlets_batch3/"),pattern="*.fcs",transformation=FALSE,truncate_max_range=FALSE)
for(i in 1:length(fcs.data_batch3@frames)){
  fcs.data_batch3[[i]]@description$batch <- 3
}

#Combine all three batches into one FlowSet
fcs.data <- rbind2(fcs.data_batch1,fcs.data_batch2)
fcs.data <- rbind2(fcs.data,fcs.data_batch3)
rm(fcs.data_batch1,fcs.data_batch2,fcs.data_batch3)
gc()

#Cofactors have manually been determined for all channels These are now loaded and channel names are changed accordingly.
cofactors <- read_excel("cofactordata_manual.xlsx")
cofactors$cofactors <- as.numeric(cofactors$cofactors)
colnames(fcs.data) <- cofactors$markerstotransform[match(colnames(fcs.data),cofactors$Colnames)]
markerstotransform <- c("FoxP3 AF488","CD161 AF700","CD57 APC","CD56 APC-Fire750",     
                        "CD39 BB700","CD3 BUV496","CD38 BUV563","CD134 BUV661",         
                        "TNFa BUV737","TCF1 BV421","IntB7 BV480","CD16 BV570",           
                        "CCR7 BV605","CD25 BV650","IL17A BV711","PD1 BV750",            
                        "CXCR3 BV785","GzmB PB","IL10 PE","TIM3 PE-Cy5",          
                        "Ki67 PE-Cy7","LAG3 PE-Dazzle594","CD69 PE-Fire640","CCR4 PE-Fire810",      
                        "HLA.DR PerCP","IL13 PerCP-Cy5.5","CTLA4 PerCP-eFluor710","CD4 PerCP-Fire806",    
                        "IL5 RY586","CD8 SparkBlue550","CD14 SparkBlue574","IFNy SparkNIR695",     
                        "CD45RA SparkUV387","IntB1 SuperBright436")

cofactors <- cofactors[order(cofactors$Colnames),]
cofactors <- cofactors %>% drop_na(cofactors)
cofactor_list <- unlist(cofactors$cofactors)

#Transform data with arcsinh transformation
fcs_transform <- transFlowVS(fcs.data, channels = markerstotransform, cofactor_list)
filenames <- sampleNames(fcs.data)
sampleNames(fcs_transform) <- filenames

#Save transformed dataset
write.flowSet(fcs_transform,outdir=paste0(getwd(),"/Transformed"))

transformed_files <- list.files(paste0(fcs.dir,"/Transformed"),pattern="fcs$")
transformed_files <- as.data.frame(transformed_files)
transformed_files

#Check if transformation went OK for samples from different batch
sampling.ceiling <- 5000
fcs_transform.ds <- fsApply(fcs_transform,function(ff){
  idx <- sample.int(nrow(ff),min(sampling.ceiling,nrow(ff)))
  ff[idx,] # alt. ff[order(idx),]
})
pdf("Manually_adjusted_density_plot_1.pdf")
plot(densityplot(~.,fcs_transform.ds[[1]]))
dev.off()
pdf("Manually_adjusted_density_plot_6.pdf")
plot(densityplot(~.,fcs_transform.ds[[6]]))
dev.off()

#Automatic quality control of flow data with peacoQC package (flowAI didn't work: not enough memory)

#Use first flowFrame to optimize parameters for QC
#ff <- fcs_transform[[1]]
#peacoqc_res <- PeacoQC(ff=ff,channels=markerstotransform,determine_good_cells="all",save_fcs=FALSE,plot=TRUE,output_directory = "PeacoQCresults",IT_limit=0.55,MAD=6)

#Then perform QC using these parameters
for(i in 1:length(fcs_transform@frames)){
  ff <- fcs_transform[[i]]
  channels=markerstotransform
  peacoqc_res <- PeacoQC(ff,channels,determine_good_cells="all",IT_limit=0.55,MAD=6,save_fcs=TRUE,suffix_fcs = paste0(sub(".fcs","",sampleNames(fcs_transform[i])),"_QC"), plot=TRUE,output_directory = "PeacoQCresults")
}

#Create table with sample names and file names so that number of cells removed can be traced back
removed <- data.frame()
j <- 1
for(i in unique(names(fcs_transform@frames))){
  removed[j,1]<-i
  removed[j,2]<-fcs_transform@frames[[i]]@description$FILENAME
  j <- j+1
}
colnames(removed) <- c("Sample_name","File_name")
removed

#### START here ####
#Construct new dataset from cleaned files
getwd()
setwd(".../Analysis/")
fcs.dir <- file.path(getwd(),"PeacoQCresults/PeacoQC_results/fcs_files")
fcs_QC_transform <- read.flowSet(path=fcs.dir,pattern="*.fcs",transformation = FALSE,truncate_max_range = FALSE)

fcs.dir <- file.path(getwd(),"tech_replicates_transf_QC")
tech_replicates <- read.flowSet(path=fcs.dir,pattern="*.fcs",transformation = FALSE,truncate_max_range = FALSE)
#Make sure object with which analysis is continued comprises 1 FlowSet with all QC-adjusted transformed FCS frames!

#Check channel deviation for all channels in technical replicates
#First for patient sample
batch_data <- data.frame()
for(i in c(2,4,6)){
  data <- prepData(tech_replicates[[i]],FACS=TRUE,transform = FALSE)
  assayNames(data)<- "exprs"
  exp_data <- assay(data,"exprs")
  exp_data <- t(exp_data)
  exp_data <- as.data.frame(exp_data)
  exp_data <- exp_data[sample(nrow(exp_data), 0.05*nrow(exp_data)), ]
  exp_data$sample <- ifelse(tech_replicates[[i]]@description$FILENAME==".../Analysis/tech_replicates_transf_QC/Bridge_2_UNI-XXX_T1_viable_QC.fcs",2,
                            ifelse(tech_replicates[[i]]@description$FILENAME==".../Analysis/tech_replicates_transf_QC/Bridge_3_UNI-XXX_T1_viable_QC.fcs",3,1))
  batch_data <- rbind(batch_data,exp_data)
}

for(i in 1:ncol(batch_data)){
histogram <- ggplot(batch_data, aes(x = get(colnames(batch_data)[i]), y = factor(batch_data$sample), fill = stat(x))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE,quantiles = c(0.025,0.25,0.75,0.975)) +
  scale_fill_viridis_c(name = paste0(print(colnames(batch_data)[i])), option = "C")
pdf(paste0("Histogram_patient_",colnames(batch_data)[i],".pdf"))
plot(histogram)
dev.off()
}

#Second for healthy control sample
batch_data <- data.frame()
for(i in c(1,3,5)){
  data <- prepData(tech_replicates[[i]],FACS=TRUE,transform = FALSE)
  assayNames(data)<- "exprs"
  exp_data <- assay(data,"exprs")
  exp_data <- t(exp_data)
  exp_data <- as.data.frame(exp_data)
  exp_data <- exp_data[sample(nrow(exp_data), 0.05*nrow(exp_data)), ]
  exp_data$sample <- ifelse(tech_replicates[[i]]@description$FILENAME==".../Analysis/tech_replicates_transf_QC/Bridge_2_MDDXXXXXXX_HD_viable_QC.fcs",2,
                            ifelse(tech_replicates[[i]]@description$FILENAME==".../Analysis/tech_replicates_transf_QC/Bridge_3_MDDXXXXXXX_HD_viable_QC.fcs",3,1))
  batch_data <- rbind(batch_data,exp_data)
}

for(i in 1:ncol(batch_data)){
  histogram <- ggplot(batch_data, aes(x = get(colnames(batch_data)[i]), y = factor(batch_data$sample), fill = stat(x))) + 
    geom_density_ridges_gradient(quantile_lines = TRUE,quantiles = c(0.025,0.25,0.75,0.975)) +
    scale_fill_viridis_c(name = paste0(print(colnames(batch_data)[i])), option = "C")
  pdf(paste0("Histogram_healthy_",colnames(batch_data)[i],".pdf"))
  plot(histogram)
  dev.off()
}

#Now check granzyme B distribution for all samples
gzmb_data <- data.frame()

for(i in 1:length(fcs_QC_transform)){
  data <- prepData(fcs_QC_transform[[i]],FACS=TRUE,transform = FALSE)
  assayNames(data)<- "exprs"
  exp_data <- assay(data,"exprs")
  exp_data <- t(exp_data)
  exp_data <- as.data.frame(exp_data)
  exp_data <- exp_data[colnames(exp_data) %in% c("GzmB PB")]
  exp_data <- as.data.frame(exp_data)
  exp_data <- exp_data[sample(nrow(exp_data), 0.05*nrow(exp_data)), ]
  exp_data <- as.data.frame(exp_data)
  exp_data$sample <- paste0(basename(fcs_QC_transform[[i]]@description$FILENAME))
  gzmb_data <- rbind(gzmb_data,exp_data)
  print(paste0("ready with sample ",i))
}
colnames(gzmb_data) <- c("GzmB","sample")

attach(gzmb_data)
histogram_gzmb <- ggplot(gzmb_data, aes(x=GzmB,y=factor(sample), fill = stat(x))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE,quantiles = c(0.025,0.25,0.75,0.975)) +
  scale_fill_viridis_c(name = "GzmB", option = "C")+
  xlim(-1,5)  
pdf("Histogram_gzmb.pdf",width=6,height=15)
plot(histogram_gzmb)
dev.off()

#Now derive the downslope inflection point of the density functions: all gzmb data will be aligned with this inflection point as landmark
infliction_points <- data.frame(sample=1,infliction_point=1)
k <- 1
attach(gzmb_data)
for(i in unique(sample)){
  data <- gzmb_data %>% filter(sample==i)
  density_gzmb <- ggplot(data, aes(x=GzmB))
  density_gzmb <- density_gzmb + geom_density()
  p <- ggplot_build(density_gzmb)
  dens <- as.data.frame(p$data)
  row_max_dens <- as.numeric(rownames(dens)[which.max(dens$density)])
  dens <- dens[-c(1:row_max_dens),]
  first <- data.frame()
  first <- diff(dens$y)/diff(dens$x)
  first <- as.data.frame(first)
  name <- paste0("dens_",i)
  name <- gsub("_viable_QC.fcs","",name)
  name1 <- paste0(name,"_first")
  name2 <- paste0(name,"_dens")
  name3 <- paste0(name,"_infliction")
  assign(name1,first)
  assign(name2,dens)
  print(paste0("Done with sample ",i))
  desc <- dens[as.numeric(rownames(first)[which.min(first$first)]),2]
  assign(name3,desc)
  infliction_points[k,1] <- paste0(i)
  infliction_points[k,2] <- desc
  k <- k+1
}

infliction_points$corr_factor <- mean(infliction_points$infliction_point) - infliction_points$infliction_point
gzmb_data_corr <- gzmb_data
gzmb_data_corr$GzmB <- gzmb_data_corr$GzmB + infliction_points$corr_factor[match(gzmb_data_corr$sample,infliction_points$sample)]

attach(gzmb_data_corr)
histogram_gzmb_corr <- ggplot(gzmb_data_corr, aes(x=GzmB,y=factor(sample), fill = stat(x))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE,quantiles = c(0.025,0.25,0.75,0.975)) +
  scale_fill_viridis_c(name = "GzmB", option = "C")+
  xlim(-1,5)  
pdf("Histogram_gzmb_corr.pdf",width=6,height=15)
plot(histogram_gzmb_corr)
dev.off()

#If correct adjustment of GzmB data is confirmed, apply correction to GzmB data in full dataset
#First prepare full dataset
dataset <- data.frame()
for(i in 1:length(fcs_QC_transform)){
  data <- prepData(fcs_QC_transform[[i]],FACS=TRUE,transform = FALSE)
  assayNames(data)<- "exprs"
  exp_data <- assay(data,"exprs")
  exp_data <- t(exp_data)
  exp_data <- as.data.frame(exp_data)
  exp_data$sample <- paste0(basename(fcs_QC_transform[[i]]@description$FILENAME))
  dataset <- rbind(dataset,exp_data)
  print(paste0("ready with sample ",i))
}
head(dataset)

#Then, adjust GzmB values
dataset$`GzmB PB` <- dataset$`GzmB PB` + infliction_points$corr_factor[match(dataset$sample,infliction_points$sample)]

#Check that adjusted GzmB distributions are correct
histogram_gzmb_full_corr <- ggplot(dataset, aes(x=`GzmB PB`,y=factor(sample), fill = stat(x))) + 
  geom_density_ridges_gradient(quantile_lines = FALSE) +
  scale_fill_viridis_c(name = "GzmB", option = "C")+
  xlim(-1,5)  
pdf("histogram_gzmb_full_corr.pdf",width=6,height=15)
plot(histogram_gzmb_full_corr)
dev.off()

#Now check FoxP3 distribution for all samples
foxp3_data <- data.frame()

for(i in 1:length(fcs_QC_transform)){
  data <- prepData(fcs_QC_transform[[i]],FACS=TRUE,transform = FALSE)
  assayNames(data)<- "exprs"
  exp_data <- assay(data,"exprs")
  exp_data <- t(exp_data)
  exp_data <- as.data.frame(exp_data)
  exp_data <- exp_data[colnames(exp_data) %in% c("FoxP3 AF488")]
  exp_data <- as.data.frame(exp_data)
  exp_data <- exp_data[sample(nrow(exp_data), 0.05*nrow(exp_data)), ]
  exp_data <- as.data.frame(exp_data)
  exp_data$sample <- paste0(basename(fcs_QC_transform[[i]]@description$FILENAME))
  foxp3_data <- rbind(foxp3_data,exp_data)
  print(paste0("ready with sample ",i))
}
colnames(foxp3_data) <- c("FoxP3","sample")

attach(foxp3_data)
histogram_foxp3 <- ggplot(foxp3_data, aes(x=FoxP3,y=factor(sample), fill = stat(x))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE,quantiles = c(0.025,0.25,0.75,0.975)) +
  scale_fill_viridis_c(name = "FoxP3", option = "C")+
  xlim(-1,3)  
pdf("Histogram_FoxP3.pdf",width=6,height=15)
plot(histogram_foxp3)
dev.off()

#Now derive the downslope inflection point of the density functions: all foxp3 data will be aligned with this inflection point as landmark
infliction_points <- data.frame(sample=1,infliction_point=1)
k <- 1
attach(foxp3_data)
for(i in unique(sample)){
  data <- foxp3_data %>% filter(sample==i)
  density_foxp3 <- ggplot(data, aes(x=FoxP3))
  density_foxp3 <- density_foxp3 + geom_density()
  p <- ggplot_build(density_foxp3)
  dens <- as.data.frame(p$data)
  row_max_dens <- as.numeric(rownames(dens)[which.max(dens$density)])
  dens <- dens[-c(1:row_max_dens),]
  first <- data.frame()
  first <- diff(dens$y)/diff(dens$x)
  first <- as.data.frame(first)
  name <- paste0("dens_",i)
  name <- gsub("_viable_QC.fcs","",name)
  name1 <- paste0(name,"_first")
  name2 <- paste0(name,"_dens")
  name3 <- paste0(name,"_infliction")
  assign(name1,first)
  assign(name2,dens)
  print(paste0("Done with sample ",i))
  desc <- dens[as.numeric(rownames(first)[which.min(first$first)]),2]
  assign(name3,desc)
  infliction_points[k,1] <- paste0(i)
  infliction_points[k,2] <- desc
  k <- k+1
}

infliction_points$corr_factor <- mean(infliction_points$infliction_point) - infliction_points$infliction_point
foxp3_data_corr <- foxp3_data
foxp3_data_corr$FoxP3 <- foxp3_data_corr$FoxP3 + infliction_points$corr_factor[match(foxp3_data_corr$sample,infliction_points$sample)]

attach(foxp3_data_corr)
histogram_foxp3_corr <- ggplot(foxp3_data_corr, aes(x=FoxP3,y=factor(sample), fill = stat(x))) + 
  geom_density_ridges_gradient(quantile_lines = TRUE,quantiles = c(0.025,0.25,0.75,0.975)) +
  scale_fill_viridis_c(name = "FoxP3", option = "C")+
  xlim(-1,5)  
pdf("Histogram_FoxP3_corr.pdf",width=6,height=15)
plot(histogram_foxp3_corr)
dev.off()

#If correct adjustment of FoxP3 data is confirmed, apply correction to FoxP3 data in full dataset
#Then, adjust FOXP3 values
dataset$`FoxP3 AF488` <- dataset$`FoxP3 AF488` + infliction_points$corr_factor[match(dataset$sample,infliction_points$sample)]

#Check that adjusted FOXP3 distributions are correct
histogram_foxp3_full_corr <- ggplot(dataset, aes(x=`FoxP3 AF488`,y=factor(sample), fill = stat(x))) + 
  geom_density_ridges_gradient(quantile_lines = FALSE) +
  scale_fill_viridis_c(name = "FoxP3", option = "C")+
  xlim(-1,5)  
pdf("histogram_FoxP3_full_corr.pdf",width=6,height=15)
plot(histogram_foxp3_full_corr)
dev.off()

#Save dataset
saveRDS(dataset,file="full_dataset.rds")