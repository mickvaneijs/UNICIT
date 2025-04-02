library(devtools)
library(readr)
library(readxl)
library(DESeq2)
library(apeglm)
library(GSVA)
library(dplyr)
library(tidyr)
library(limma)
library(ggplot2)
library(MASS)
library(caret)
library(pROC)
library(purrr)
library(pheatmap)
library(ggsignif)
library(ggrepel)
library(reshape2)
library(car)
library(CIBERSORT)
library(GSEABase)
library(mastR)
library(immunedeconv)
library(remotes)
library(preprocessCore)
library(ggpubr)
library(corrplot)
library(Hmisc)
library(coin) #for exact wilcox_test() with ties
library(openxlsx)
library(patchwork)

setwd("~")
setwd(".../Counts/")

#########################################
########## BASELINE ANALYSES ############
#########################################

#Create raw count expression matrix
count_files <- list.files(pattern="count$")
read_count_file <- function(file){
  df <- read_delim(file,delim="\t",col_names=FALSE)
  colnames(df) <- c("Ensembl","Gene",sub("\\.count$","",basename(file)))
  return(df)
}
count_list <- lapply(count_files,read_count_file)
count_matrix <- Reduce(function(x,y) merge(x,y,by=c("Ensembl","Gene"),all=TRUE),count_list)

#Filter low-abundant genes
keep <- rowSums(count_matrix[,c(3:33)]) > 5
count_matrix <- count_matrix[keep,]

setwd("~")
setwd(".../")
metadata <- read_excel("metadata.xlsx")
colnames_mapping <- setNames(metadata$pseudo_ID, metadata$GS_ID)
colnames(count_matrix) <- colnames_mapping[colnames(count_matrix)]
colnames(count_matrix)[1:2] <- c("Ensembl","Gene")
count_matrix_input <- count_matrix[,-1]
count_matrix_input <- count_matrix_input %>% dplyr::filter(Gene!="-")
count_matrix_input_new <- count_matrix_input %>% group_by(Gene) %>% summarise(across(everything(), sum))
count_matrix_input_new <- as.data.frame(count_matrix_input_new)
rownames(count_matrix_input_new) <- count_matrix_input_new$Gene
count_matrix_input_new <- count_matrix_input_new[,-1]
rownames(metadata) <- metadata$pseudo_ID
exclude_IDs <- metadata$pseudo_ID[which(metadata$response_pred=="NT")]
count_matrix_R_NR <- count_matrix_input_new[!colnames(count_matrix_input_new) %in% exclude_IDs]
exclude_IDs <- metadata$pseudo_ID[which(metadata$pre_post_pred=="post")]
count_matrix_R_NR <- count_matrix_R_NR[!colnames(count_matrix_R_NR) %in% exclude_IDs]
metadata_R_NR <- metadata %>% dplyr::filter(rownames(metadata) %in% colnames(count_matrix_R_NR))
metadata_R_NR$ICI_main <- ifelse(metadata_R_NR$ICI_type=="combination ICI","ipi_nivo",
                           ifelse(metadata_R_NR$ICI_type=="ipi_mono","ipi_nivo","anti_PD_1"))
metadata_R_NR$response_pred <- factor(metadata_R_NR$response_pred,levels=c("R","NR"))

setwd("~")
setwd(".../")


#Perform DESeq2 on steroid-naive samples of R and NR only (NT are excluded)
dds <- DESeqDataSetFromMatrix(countData=count_matrix_R_NR,
                              colData=metadata_R_NR,
                              design=~response_pred)
dds <- DESeq(dds)
resultsNames(dds)
res1 <- results(dds,name="response_pred_NR_vs_R")

volcano_results <- as.data.frame(res1@listData)
rownames(volcano_results) <- res1@rownames
volcano_results$label <- ifelse(volcano_results$padj<0.05 & abs(volcano_results$log2FoldChange)>2,rownames(volcano_results),"")
volcano_results$color <- ifelse(volcano_results$padj<0.05 & volcano_results$log2FoldChange>1,"NR_up_sign",
                                ifelse(volcano_results$padj<0.05 & volcano_results$log2FoldChange<(-1),"NR_dn_sign",
                                       ifelse(volcano_results$padj>=0.05 & volcano_results$log2FoldChange<(-1),"NR_dn",
                                              ifelse(volcano_results$padj>=0.05 & volcano_results$log2FoldChange>1,"NR_up",""))))
volcano_results$color <- factor(volcano_results$color,levels=c("","NR_up","NR_up_sign","NR_dn","NR_dn_sign"))
labels_to_show <- c("SELL","MT1F","MT1M","PDPN","SIGLEC9","WARS1","ALPL","CXCL8",
                    "LILRB3","CXCL11","CXCL9")
volcano_NR_v_R <- ggplot(volcano_results,aes(log2FoldChange,-log10(padj),color=color))+
  geom_point()+
  geom_hline(yintercept = 1.3,linetype="dashed")+
  geom_vline(xintercept = 1,linetype="dashed")+
  geom_vline(xintercept = -1,linetype="dashed")+
  scale_color_manual(values=c("grey","#8E6698","#9edea0"))+
  geom_text_repel(aes(label=ifelse(label %in% labels_to_show,label,"")),max.overlaps=Inf)+
  theme_minimal()+
  theme(legend.position="none")+
  xlim(-3.5,3)+
  ylim(0,1.8)+
  ylab("-log10(Padj)")+
  xlab("log2FC Non-responder versus Responder")
pdf("volcano_NR_v_R_baseline.pdf",width=7,height=5)
plot(volcano_NR_v_R)
dev.off()

#Next: create logTPM matrix for GSVA / TPM matrix for CIBERSORT
setwd("~")
setwd(".../TPM/")

#Create TPM expression matrix
count_files <- list.files(pattern="*.tpm.tsv$")

read_count_file <- function(file){
  df <- read_delim(file,delim="\t",col_names=TRUE)
  df <- df[,c("gene_id","gene_short_name","FPKM")]
  colnames(df) <- c("Ensembl","Gene",sub("\\.tpm.tsv$","",basename(file)))
  return(df)
}
count_list <- lapply(count_files,read_count_file)

merged <- count_list[[1]]
for(i in 2:length(count_list)){
  count_list[[i]] <- count_list[[i]][count_list[[i]][,3] > 0,]
  Column=colnames(count_list[[i]])[3]
  print(paste0("Duplicate Ensembl IDs for ",i,": ",(nrow(count_list[[i]])-nrow(unique(count_list[[i]][,1])))/nrow(count_list[[i]])," %."))
  print(paste0("Duplicate Gene names for ",i," excluding Ensembl IDs not mapping to a gene name: ",round(100*(nrow(count_list[[i]][count_list[[i]]$Gene!="-",])-nrow(unique(count_list[[i]][,2])))/nrow(count_list[[i]][count_list[[i]]$Gene!="-",]),2)," %."))
  count_list[[i]] = count_list[[i]] %>% group_by(Ensembl,Gene) %>% summarise(sum_value = sum(get(Column)))
  colnames(count_list[[i]]) <- c("Ensembl",
                                 "Gene",
                                 ifelse(i<10,
                                        paste0("106131-001-00",i),
                                        paste0("106131-001-0",i)))
  merged <- merge(merged,count_list[[i]],by=c("Ensembl","Gene"),all=TRUE)
  print(paste0("Merged file ",i," of ",length(count_list)))
}
merged[is.na(merged)] <- 0

#Drop low-abundant genes across patients
keep <- rowSums(merged[,c(3:ncol(merged))]) > 0.01
merged <- merged[keep,]

#First, NON-log normalized data: see this website where use of NON-log-normalized TPM 
#data is advised for CIBERSORT: https://omnideconv.org/immunedeconv/articles/immunedeconv.html
TPM_matrix <- merged

#Merge with metadata
setwd("~")
setwd(".../")
metadata <- read_excel("metadata.xlsx")
colnames_mapping <- setNames(metadata$pseudo_ID, metadata$GS_ID)
colnames(TPM_matrix) <- colnames_mapping[colnames(TPM_matrix)]
colnames(TPM_matrix)[1:2] <- c("Ensembl","Gene")
TPM_matrix_input <- TPM_matrix[,-1]
TPM_matrix_input <- TPM_matrix_input %>% dplyr::filter(Gene!="-")
TPM_matrix_input_new <- TPM_matrix_input %>% group_by(Gene) %>% summarise(across(everything(), sum))
TPM_matrix_input_new <- as.data.frame(TPM_matrix_input_new)
rownames(TPM_matrix_input_new) <- TPM_matrix_input_new$Gene
TPM_matrix_input_new <- TPM_matrix_input_new[,-1]

##############################################################################
################################ CIBERSORT ###################################
##############################################################################
lm22 <- CIBERSORT::LM22
lm22 <- as.data.frame(lm22)
lm22 <- as.matrix(lm22)
CIBERSORT_matrix <- TPM_matrix_input_new
deconv_res_22 <- cibersort(lm22, 
                          as.matrix(CIBERSORT_matrix),
                          perm = 1000)
setwd("~")
setwd(".../")
write.xlsx(as.data.frame(deconv_res_22),"deconv_res_22_1000perm.xlsx",rowNames=TRUE)

#reload
deconv_res_22 <- read_excel("deconv_res_22_1000perm.xlsx")
deconv_res_22 <- as.data.frame(deconv_res_22)
rownames(deconv_res_22) <- deconv_res_22$...1
deconv_res_22 <- deconv_res_22[,-1]

abs_score <- list()
Ymedian <- max(median(unlist(CIBERSORT_matrix)),1)
for(i in 1:ncol(CIBERSORT_matrix)){
  abs_score[[i]] <- median(CIBERSORT_matrix[,i])/Ymedian
}

deconv_res_22_abs <- deconv_res_22
for(i in 1:nrow(deconv_res_22_abs)){
  deconv_res_22_abs[i,c(1:22)] <- deconv_res_22_abs[i,c(1:22)]*abs_score[[i]]
}
write.xlsx(as.data.frame(deconv_res_22_abs),"deconv_res_22_abs_1000perm.xlsx",rowNames=TRUE)

panel_22colors <- c("B cells naive"="#E6194B","B cells memory"="#3CB44B", "Plasma cells"="#FFE119",
                    "T cells CD8"="#4363D8","T cells CD4 naive"="#F58231","T cells CD4 memory resting"="#911EB4",
                    "T cells CD4 memory activated"="#46F0F0","T cells follicular helper"="#F032E6",
                    "T cells regulatory (Tregs)"="#BCF60C","T cells gamma delta"="#FABEBE",
                    "NK cells resting"="#008080","NK cells activated"="#E6BEFF","Monocytes"="#9A6324",
                    "Macrophages M0"="#FFFAC8","Macrophages M1"="#800000","Macrophages M2"="#AAFFC3",
                    "Dendritic cells resting"="#808000","Dendritic cells activated"="#FFD8B1",
                    "Mast cells resting"="#000075","Mast cells activated"="#808080","Eosinophils"="#FFFFFF",
                    "Neutrophils"="#000000")

response_colors <- c("NT"="#FDE725","R"="#5ec962","NR"="#440154")

response_cat <- metadata[,c(2,10,13,16,17)]

#Filter samples with P<0.05
deconv_res_22_pass <- deconv_res_22 %>% dplyr::filter(`P-value`<0.05)
deconv_res_22_abs_pass <- deconv_res_22_abs %>% dplyr::filter(`P-value`<0.05)
deconv_res_22_pass$pseudo_ID <- rownames(deconv_res_22_pass)
deconv_res_22_abs_pass$pseudo_ID <- rownames(deconv_res_22_abs_pass)

data_abs_lm22 <- melt(deconv_res_22_abs_pass[,c(1:22,26)],id.vars="pseudo_ID")
colnames(data_abs_lm22)<- c("pseudo_ID","cell","value")
data_abs_lm22 <- merge(data_abs_lm22,response_cat,by="pseudo_ID")
colnames(data_abs_lm22)<- c("subject","cell","value","ICI_type","pred_days","response","CIBERSORT_pass")

#Filter only high-quality baseline samples
data_filter_abs <- data_abs_lm22 %>% dplyr::filter(CIBERSORT_pass=="yes" & response != "NT" & pred_days == 0)
abs_lm22_cibersort <- ggplot(data_filter_abs,aes(reorder(subject,value,sum),value))+
  geom_bar(position="stack",stat="identity",aes(fill=cell))+
  scale_fill_manual(values=c(panel_22colors,response_colors))+
  geom_tile(aes(x=subject,y=-0.1,fill=response),width=0.9,height=0.1,color="white")+
  theme_classic()+
  xlab("")+
  ylab("CIBERSORT absolute abundance")
pdf("abs_lm22_cibersort_filtered_baseline.pdf",width=8,height=4)
plot(abs_lm22_cibersort)
dev.off()

#Cell abundance adjusted for ICI-type
cell_diff <- data_abs_lm22 %>% dplyr::filter(response!="NT" & CIBERSORT_pass=="yes" & pred_days==0)
cell_diff$ICI_main <- ifelse(cell_diff$ICI_type=="ipi_mono","ipi_nivo",
                             ifelse(cell_diff$ICI_type=="combination ICI","ipi_nivo","anti_PD_1"))
cell_diff$ICI_main <- factor(cell_diff$ICI_main,levels=c("anti_PD_1","ipi_nivo"))
cell_diff$response <- factor(cell_diff$response,levels=c("NT","R","NR"))
colors1 <- c("R"="#5ec962","NR"="#440154")
colors2 <- c("anti_PD_1"="#F98E09","ipi_nivo"="#BC3754")

df_adjusted <- data.frame(matrix(ncol=4))
colnames(df_adjusted) <- c("cell_type","coef_response","ci_lower","ci_upper")
r <- 1
for(c in unique(cell_diff$cell)){
  data <- cell_diff %>% dplyr::filter(cell==c)
  data$response <- factor(data$response,levels=c("R","NR"))
  model <- lm(value~response+ICI_main,data=data)
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
pdf("Adjusted_plot_CIBERSORT_baseline.pdf",width=5,height=6)
plot(adjusted_plot)
dev.off()

#Create gene set and perform GSVA
gene_sets <- read_excel("gene_sets.xlsx")

gene.sets <- list(bMIS_UC=gene_sets$bMIS_UC,
                  Th17=gene_sets$Th17,
                  Th2=gene_sets$Th2,
                  Th1=gene_sets$Th1,
                  Bmem_PC=gene_sets$Bmem_PC,
                  Th17_core=gene_sets$Th17_core,
                  Th17_selective=gene_sets$Th17_selective,
                  Th17.1_selective=gene_sets$Th17.1_selective)
gene.sets <- lapply(gene.sets,function(x) x[!is.na(x)])

#Log-normalize TPM data and merge again with metadata
merged_log <- merged
numeric_cols <- sapply(merged_log,is.numeric)
merged_log[,numeric_cols]<-lapply(merged_log[,numeric_cols]+0.001,log)
TPM_matrix <- merged_log

#Merge with metadata
setwd("~")
setwd(".../")
metadata <- read_excel("metadata.xlsx")
colnames_mapping <- setNames(metadata$pseudo_ID, metadata$GS_ID)
colnames(TPM_matrix) <- colnames_mapping[colnames(TPM_matrix)]
colnames(TPM_matrix)[1:2] <- c("Ensembl","Gene")
TPM_matrix_input <- TPM_matrix[,-1]
TPM_matrix_input <- TPM_matrix_input %>% dplyr::filter(Gene!="-")
TPM_matrix_input_new <- TPM_matrix_input %>% group_by(Gene) %>% summarise(across(everything(), sum))
TPM_matrix_input_new <- as.data.frame(TPM_matrix_input_new)
rownames(TPM_matrix_input_new) <- TPM_matrix_input_new$Gene
TPM_matrix_input_new <- TPM_matrix_input_new[,-1]

ColitisGSVA_param <- gsvaParam(as.matrix(TPM_matrix_input_new),
                               geneSets = gene.sets,
                               kcdf="Gaussian")
ColitisGSVA <- gsva(ColitisGSVA_param)

results <- as.data.frame(t(ColitisGSVA))
results$pseudo_ID <- rownames(results)
results <- merge(results,metadata,by=c("pseudo_ID"))
results$ICI_main <- ifelse(results$ICI_type=="combination ICI","ipi_nivo",
                           ifelse(results$ICI_type=="ipi_mono","ipi_nivo","anti_PD_1"))

#Figures
setwd("~")
setwd(".../")

results_2 <- results %>% dplyr::filter(response_pred!="NT" & n_days_pred==0)
results_2$response_pred <- factor(results_2$response_pred,levels=c("R","NR"))
colors1 <- c("R"="#5ec962","NR"="#440154") #"#FDE725"
results_2$ICI_main <- factor(results_2$ICI_main,levels=c("anti_PD_1","ipi_nivo"))
colors2 <- c("anti_PD_1"="#F98E09","ipi_nivo"="#BC3754")

#bMIS_UC score by response
plot1 <- ggplot(results_2,aes(factor(response_pred),bMIS_UC,fill=response_pred))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(pch=21,position=position_jitter(height=0))+
  scale_fill_manual(labels=c("Responder","Non-responder"),values=colors1)+
  geom_signif(annotation=paste0(expression(italic(P) == 0.93)),parse=T,y_position=0.7,xmin=1,xmax=2,tip_length=c(0,0))+
  xlab("")+
  ylab("GSVA score bMIS_UC")+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        axis.text.x=element_blank())
pdf("GSVA_bMIS_response_baseline.pdf",width=4,height=5)
plot(plot1)
dev.off()
wilcox.test(results_2$bMIS_UC~results_2$response_pred)

#Analyses in H&E-only cohort
setwd("~")
setwd(".../neutrophil_validation/results/")
list <- list.files(getwd())
file <- read.csv(list[1])

all_cols <- unique(unlist(lapply(list, function(f) names(read.csv(f, nrows = 1)))))
df <- as.data.frame(matrix(ncol=length(all_cols)+1,nrow=0))
colnames(df) <- c(all_cols,"sample")

rbind_fill <- function(df1, df2) {
  all_cols <- union(names(df1), names(df2))
  
  df1 <- df1 %>%
    dplyr::select(all_of(all_cols)) %>%
    dplyr::mutate(across(setdiff(all_cols, names(df1)), ~ NA))
  
  df2 <- df2 %>%
    dplyr::select(all_of(all_cols)) %>%
    dplyr::mutate(across(setdiff(all_cols, names(df2)), ~ NA))
  
  return(dplyr::bind_rows(df1, df2))
}

for(i in 1:length(list)) {
  file <- read.csv(list[i])
  file$sample <- sub(".csv","",list[i])
  missing_cols <- setdiff(names(df), names(file))
  file[missing_cols] <- NA
  df <- rbind(df,file)
  message(paste0("Done with sample ",i))
}
#Drop geomtry column
df <- df[,-1]

setwd("~")
setwd(".../")

n_tiles <- df %>%
  group_by(sample) %>%
  dplyr::summarize(count = n(), total_cells = sum(total_cell_cnt, na.rm = TRUE))

df$total_cnt_immune <- df$total_cnt_eosinophil+
  df$total_cnt_lymphocyte+
  df$total_cnt_neutrophil+
  df$total_cnt_plasma.cell
df$fraction_immune <- df$total_cnt_immune/df$total_cell_cnt

#Loading clinical metadata
clinical <- read_excel("complete_clinical.xlsx")
clinical$Lymphocytes <- clinical$`B cells naive`+
  clinical$`B cells memory`+
  clinical$`T cells CD8`+
  clinical$`T cells CD4 naive`+
  clinical$`T cells CD4 memory resting`+
  clinical$`T cells CD4 memory activated`+
  clinical$`T cells follicular helper`+
  clinical$`T cells regulatory (Tregs)`+
  clinical$`T cells gamma delta`
clinical <- clinical[,c("Tnr","Cohort","Response to steroids","Steroids duration at endoscopy (days)","Pre-steroids biopsy","pseudo_ID","Neutrophils","Lymphocytes","Plasma cells","Eosinophils")]
colnames(clinical)[1] <- "sample"
clinical$Response <- clinical$`Response to steroids`

merged <- merge(n_tiles,clinical,by="sample")

count <- df %>%
  dplyr::select(total_cell_cnt,total_cnt_neutrophil,total_cnt_lymphocyte,total_cnt_plasma.cell,sample,total_cnt_eosinophil)

df_fraction <- count %>%
  group_by(sample) %>%
  dplyr::summarize(
    total_neutrophils = sum(total_cnt_neutrophil, na.rm = TRUE),
    total_cells = sum(total_cell_cnt, na.rm = TRUE),
    total_lymphocytes = sum(total_cnt_lymphocyte, na.rm = TRUE),
    total_plasmacells = sum(total_cnt_plasma.cell, na.rm = TRUE),
    total_eosinophils = sum(total_cnt_eosinophil, na.rm = TRUE)
  ) %>%
  mutate(fraction_neutrophils = total_neutrophils / (total_lymphocytes+total_plasmacells+total_eosinophils+total_neutrophils),
         fraction_lymphocytes = total_lymphocytes / (total_lymphocytes+total_plasmacells+total_eosinophils+total_neutrophils),
         fraction_plasmacells = total_plasmacells / (total_lymphocytes+total_plasmacells+total_eosinophils+total_neutrophils),
         fraction_eosinophils = total_eosinophils / (total_lymphocytes+total_plasmacells+total_eosinophils+total_neutrophils),
         fraction_granulocytes = (total_eosinophils + total_neutrophils) / (total_lymphocytes+total_plasmacells+total_eosinophils+total_neutrophils),
         total_cell_cnt = total_cells,
         total_leukocyte_fraction = (total_neutrophils+total_lymphocytes+total_plasmacells+total_eosinophils)/total_cells
  ) %>%
  dplyr::select(sample,total_cell_cnt,total_leukocyte_fraction,fraction_neutrophils,fraction_lymphocytes,fraction_plasmacells,fraction_eosinophils,fraction_granulocytes,
                total_neutrophils,total_lymphocytes,total_plasmacells,total_eosinophils)

merged <- merge(merged,df_fraction,by="sample",all=T)
merged <- merged %>% drop_na(`Pre-steroids biopsy`)
treated <- merged %>% dplyr::filter(`Pre-steroids biopsy`!="no")
treated$Response <- treated$`Response to steroids`

treated$Response <- factor(treated$Response,levels=c("responder","non-responder"))

colors <- c("responder"="#5ec962","non-responder"="#440154") #"#FDE725"

treated_baseline <- treated %>% dplyr::filter(`Steroids duration at endoscopy (days)`==0)
plot <- ggplot(treated_baseline,aes(Response,100*total_neutrophils/total_cell_cnt,fill=Response))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(shape=21,position=position_jitter(height=0))+
  ylab("% neutrophils of all tissue cells")+
  xlab("")+
  labs(fill=NULL)+
  scale_fill_manual(values=colors)+
  geom_signif(annotation=paste0(expression(italic(P) == 0.56)),parse=T,y_position=3,xmin=1,xmax=2,tip_length=c(0,0))+
  theme_classic()
pdf("neutrophils_baseline.pdf",width=4,height=5)
plot(plot)
dev.off()
treated_baseline$fraction_neutros_total <- treated_baseline$total_neutrophils / treated_baseline$total_cell_cnt
wilcox_test(treated_baseline$fraction_neutros_total~treated_baseline$Response)

plot <- ggplot(treated_baseline,aes(Response,100*total_lymphocytes/total_cell_cnt,fill=Response))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(shape=21,position=position_jitter(height=0))+
  ylab("% lymphocytes of all tissue cells")+
  xlab("")+
  labs(fill=NULL)+
  scale_fill_manual(values=colors)+
  geom_signif(annotation=paste0(expression(italic(P) == 0.23)),parse=T,y_position=50,xmin=1,xmax=2,tip_length=c(0,0))+
  theme_classic()
pdf("lymphocytes_baseline.pdf",width=4,height=5)
plot(plot)
dev.off()
treated_baseline$fraction_lymphos_total <- treated_baseline$total_lymphocytes / treated_baseline$total_cell_cnt
wilcox.test(treated_baseline$fraction_lymphos_total~treated_baseline$Response)

plot <- ggplot(treated_baseline,aes(Response,100*total_leukocyte_fraction,fill=Response))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(shape=21,position=position_jitter(height=0))+
  ylab("% leukocytes of all tissue cells")+
  xlab("")+
  labs(fill=NULL)+
  scale_fill_manual(values=colors)+
  geom_signif(annotation=paste0(expression(italic(P) == 0.11)),parse=T,y_position=63,xmin=1,xmax=2,tip_length=c(0,0))+
  theme_classic()
pdf("leukocytes_baseline.pdf",width=4,height=5)
plot(plot)
dev.off()
wilcox.test(treated_baseline$total_leukocyte_fraction~treated_baseline$Response)

###################################
########### FIGURE 2 ##############
###################################
results_2 <- results %>% dplyr::filter(response_pred!="NT")
results_2$response_pred <- ifelse(results_2$response_pred=="R","Responder","Non-responder")
colors <- c("Responder"="#5ec962","Non-responder"="#440154") #"#FDE725"
results_3 <- results_2
results_3$pred_bin <- ifelse(results_3$n_days_pred==0,"0",
                             ifelse(results_3$n_days_pred==1,"1",
                                    ifelse(results_3$n_days_pred>7,">7","2-7")))
results_3$pred_bin <- factor(results_3$pred_bin,levels=c("0","1","2-7",">7"))
results_3$response_pred <- factor(results_3$response_pred,levels=c("Responder","Non-responder"))
plot <- ggplot(results_3,aes(pred_bin,bMIS_UC,fill=response_pred))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(shape=21,position=position_jitterdodge())+
  scale_fill_manual(values=colors)+
  theme_classic()+
  ylab("GSVA score bMIS_UC")+
  xlab("Days on steroids upon endoscopy")+
  labs(fill=NULL)
pdf("bMIS_UC_longitudinal.pdf",height=4,width=5)
plot(plot)
dev.off()

#H&E analyses. Time capped at 7 days!
treated_cont <- treated
treated_cont$`Steroids duration at endoscopy (days)` <- ifelse(treated_cont$`Steroids duration at endoscopy (days)`>7,7,
                                                               treated_cont$`Steroids duration at endoscopy (days)`)
treated_cont$shape <- ifelse(treated_cont$`Steroids duration at endoscopy (days)`==7,"capped","original")
treated_cont$Response <- ifelse(treated_cont$Response=="responder","Responder","Non-responder")
treated_cont$Response <- factor(treated_cont$Response,levels=c("Responder","Non-responder"))
plot <- ggplot(treated_cont,aes(`Steroids duration at endoscopy (days)`,100*total_neutrophils/total_cell_cnt,color=Response))+
  geom_smooth(method="loess",se=T,alpha=0.3,aes(fill=Response))+
  geom_point(aes(shape=shape),position=position_jitterdodge())+
  scale_shape_manual(values=c(17,21))+
  ylab("% neutrophils of all tissue cells")+
  xlab("Days on steroids upon endoscopy")+
  geom_segment(aes(x=-0.2,xend=0.8,y=4.5,yend=4.5),color="#5ec962")+
  annotate("text",x=0.3,y=4.8,label=expression(italic(P)~"=0.14"),parse=T,color="#5ec962")+
  geom_segment(aes(x=0.2,xend=1.2,y=5,yend=5),color="#440154")+
  annotate("text",x=0.7,y=5.3,label=expression(italic(P)~"=0.96"),parse=T,color="#440154")+
  labs(color=NULL)+
  labs(fill=NULL)+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  #geom_signif(annotation=paste0(expression(italic(P) == 0.83)),parse=T,y_position=3,xmin=1,xmax=2,tip_length=c(0,0))+
  theme_classic()
pdf("neutrophils_longitudinal_continuous.pdf",width=5,height=4)
plot(plot)
dev.off()

table(treated_cont$`Response to steroids`)

plot <- ggplot(treated_cont,aes(`Steroids duration at endoscopy (days)`,100*total_lymphocytes/total_cell_cnt,color=Response))+
  geom_smooth(method="loess",se=T,alpha=0.3,aes(fill=Response))+
  geom_point(aes(shape=shape),position=position_jitterdodge())+
  scale_shape_manual(values=c(17,21))+
  ylab("% lymphocytes of all tissue cells")+
  xlab("Days on steroids upon endoscopy")+
  geom_segment(aes(x=-0.2,xend=0.8,y=50,yend=50),color="#5ec962")+
  annotate("text",x=0.3,y=52,label=expression(italic(P)~"=0.0054"),parse=T,color="#5ec962")+
  geom_segment(aes(x=0.2,xend=1.2,y=54,yend=54),color="#440154")+
  annotate("text",x=0.7,y=56,label=expression(italic(P)~"=0.56"),parse=T,color="#440154")+
  labs(color=NULL)+
  labs(fill=NULL)+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  #geom_signif(annotation=paste0(expression(italic(P) == 0.83)),parse=T,y_position=3,xmin=1,xmax=2,tip_length=c(0,0))+
  theme_classic()
pdf("lymphocytes_longitudinal_continuous.pdf",width=5,height=4)
plot(plot)
dev.off()

test <- treated_cont %>% dplyr::filter(`Steroids duration at endoscopy (days)`<2 & Response=="Responder")
test$total_cell_cnt <- as.numeric(test$total_cell_cnt)
test$total_neutrophils <- as.numeric(test$total_neutrophils)
test$total_lymphocytes <- as.numeric(test$total_lymphocytes)
test$neutro_frac <- test$total_neutrophils/test$total_cell_cnt
test$lympho_frac <- test$total_lymphocytes/test$total_cell_cn
wilcox_test(test$neutro_frac~as.factor(test$`Steroids duration at endoscopy (days)`))
wilcox_test(test$lympho_frac~as.factor(test$`Steroids duration at endoscopy (days)`))

test <- treated_cont %>% dplyr::filter(`Steroids duration at endoscopy (days)`<2 & Response=="Non-responder")
test$total_cell_cnt <- as.numeric(test$total_cell_cnt)
test$total_neutrophils <- as.numeric(test$total_neutrophils)
test$total_lymphocytes <- as.numeric(test$total_lymphocytes)
test$neutro_frac <- test$total_neutrophils/test$total_cell_cnt
test$lympho_frac <- test$total_lymphocytes/test$total_cell_cn
wilcox_test(test$neutro_frac~as.factor(test$`Steroids duration at endoscopy (days)`))
wilcox_test(test$lympho_frac~as.factor(test$`Steroids duration at endoscopy (days)`))

####################################################
### Figure 3 and associated supplementary figure ###
####################################################

setwd("~")
setwd(".../")

#Filter only high-quality samples
data_all_abs <- data_abs_lm22 %>% dplyr::filter(CIBERSORT_pass=="yes" & response!="NT")
abs_lm22_cibersort_all <- ggplot(data_all_abs,aes(reorder(subject,value,sum),value))+
  geom_bar(position="stack",stat="identity",aes(fill=cell))+
  scale_fill_manual(values=c(panel_22colors,response_colors))+
  geom_tile(aes(x=subject,y=-0.1,fill=response),width=0.9,height=0.1,color="white")+
  theme_classic()+
  xlab("")+
  ylab("CIBERSORT absolute abundance")
pdf("abs_lm22_cibersort_all.pdf",width=12,height=4)
plot(abs_lm22_cibersort_all)
dev.off()

#bMIS_UC score by response including untreated
results_2 <- results %>% dplyr::filter(response_pred!="NT")
results_2$response_pred <- factor(results_2$response_pred,levels=c("R","NR"))
colors1 <- c("#5ec962","#9edea0","#440154","#8E6698") 
results_2$ICI_main <- factor(results_2$ICI_main,levels=c("anti_PD_1","ipi_nivo"))
colors2 <- c("anti_PD_1"="#F98E09","ipi_nivo"="#BC3754")
results_2$baseline <- ifelse(results_2$n_days_pred==0,"steroid-naive","treated")
results_2$combined <- paste0(results_2$response_pred,"_",results_2$baseline)

results_2$combined <- factor(results_2$combined,levels=c("R_steroid-naive",
                                                         "R_treated",
                                                         "NR_steroid-naive",
                                                         "NR_treated"))
plot1 <- ggplot(results_2,aes(combined,bMIS_UC,fill=combined))+
  geom_boxplot(outlier.shape=NA,alpha=0.7)+
  geom_point(pch=21,position=position_jitter(height=0))+
  scale_fill_manual(labels=c("Responder baseline",
                             "Responder treated",
                             "Non-responder baseline",
                             "Non-responder treated"),values=colors1)+
  geom_signif(annotation=paste0(expression(italic(P) == 0.096)),parse=T,y_position=0.83,xmin=1.5,xmax=3.5,tip_length=c(0.1,0.1))+
  geom_signif(annotation="",y_position=0.7,xmin=1,xmax=2,tip_length=c(0,0))+
  geom_signif(annotation="",y_position=0.7,xmin=3,xmax=4,tip_length=c(0,0))+
  xlab("")+
  ylab("GSVA score bMIS_UC")+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        axis.text.x = element_blank())
pdf("GSVA_bMIS_response_all_treated.pdf",width=6,height=5)
plot(plot1)
dev.off()
stat <- results_2
stat$group <- ifelse(stat$response_pred=="NR","NR","R")
wilcox.test(stat$bMIS_UC~stat$group)

colors3 <- c("#5ec962","#440154") 
plot2 <- ggplot(results_2,aes(mayo_score,bMIS_UC))+
  geom_smooth(method="lm",color="grey",se=T)+
  geom_point(aes(shape=baseline,fill=response_pred),position=position_jitter(height=0,width=0.1))+
  stat_cor(method="pearson",label.x=1,label.y=-0.2)+
  scale_fill_manual(labels=c("Responder",
                             "Non-responder"),
                    values=c(R="#5ec962",NR="#440154"))+
  scale_shape_manual(labels=c("Steroid-naive",
                             "Treated"),
                     values=c(21,24))+
  xlab("Mayo endoscopic score")+
  ylab("GSVA score bMIS_UC")+
  theme_classic()+
  theme(legend.title=element_blank())+
  guides(fill=guide_legend(override.aes=list(shape=21)))
pdf("GSVA_bMIS_mayo.pdf",width=5,height=4)
plot(plot2)
dev.off()

results_2_long <- melt(results_2,id.vars=c("pseudo_ID","response_pred","ICI_main"),
                       measure.vars=c("Th1","Th2","Th17","Bmem_PC"),
                       variable.name="Th",
                       value.name="score")

results_2_long$response_pred <- factor(results_2_long$response_pred,levels=c("NT","R","NR"))
results_2_long$Th <- factor(results_2_long$Th,levels=c("Th1","Th2","Th17","Bmem_PC"))
colors1 <- c("#5ec962","#440154")
plot3 <- ggplot(results_2_long,aes(Th,score,fill=response_pred))+
  geom_boxplot(outlier.shape=NA,alpha=0.7)+
  geom_point(pch=21,position=position_jitterdodge())+
  scale_fill_manual(labels=c("Responder","Non-responder"),values=colors1)+
  geom_signif(annotation=paste0(expression(italic(P) == 0.013)),parse=T,y_position=0.9,xmin=0.78,xmax=1.22,tip_length=c(0,0))+
  geom_signif(annotation=paste0(expression(italic(P) == 0.031)),parse=T,y_position=0.9,xmin=1.78,xmax=2.22,tip_length=c(0,0))+
  geom_signif(annotation=paste0(expression(italic(P) == 0.022)),parse=T,y_position=0.9,xmin=2.78,xmax=3.22,tip_length=c(0,0))+
  geom_signif(annotation=paste0(expression(italic(P) == 0.24)),parse=T,y_position=0.9,xmin=3.78,xmax=4.22,tip_length=c(0,0))+
  xlab("")+
  ylab("GSVA score immune response type")+
  theme_classic()+
  scale_x_discrete(labels=c("Th1"="Type 1","Th2"="Type 2","Th17"="Type 17","Bmem_PC"="Memory B/plasma cells"))+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())
pdf("GSVA_Th_response.pdf",width=6,height=5)
plot(plot3)
dev.off()
results_2_stat <- results_2_long %>% dplyr::filter(Th=="Bmem_PC")
wilcox.test(results_2_stat$score~results_2_stat$response_pred)

#Correlation plot CIBERSORT results (and bMIS_UC score)
results2 <- results %>% dplyr::filter(response_pred!="NT")
colnames(data_abs_lm22)[1] <- c("pseudo_ID")
corr_data <- merge(data_abs_lm22,results2[,c(1:2)],by="pseudo_ID")
corr_data <- corr_data %>% pivot_wider(names_from=cell,values_from = value)
corr_data$total_cells <- rowSums(corr_data[,c(6:ncol(corr_data))])
colors <- c("R"="#5ec962","NR"="#440154")
corr_plot_table <- as.matrix(corr_data[,c(6:(ncol(corr_data)-1))])
corr_plot_data <- cor(corr_plot_table,method="spearman")
corr_plot_data[is.na(corr_plot_data)] <- 0
testRes <- cor.mtest(corr_plot_table,method="spearman")
testRes$p[is.na(testRes$p)] <- 1
testRes$p.adj <- matrix(p.adjust(as.vector(as.matrix(testRes$p)),method="BH"),ncol=23)
rownames(testRes$p.adj) <- rownames(testRes$p)
colnames(testRes$p.adj) <- colnames(testRes$p)

setwd("~")
setwd(".../")

p_data <- testRes$p.adj
hclust_order <- corrMatOrder(corr_plot_data,order="hclust")
corr_plot_data <- corr_plot_data[hclust_order,hclust_order]
p_data <- p_data[hclust_order,hclust_order]

pdf("corr_plot_CIBERSORT.pdf",width=6,height=6)
corrplot(corr_plot_data,
         p.mat=p_data,
         method="color",
         tl.col="black",
         type="lower",
         diag = FALSE,
         sig.level=c(0.001,0.01,0.05),
         pch.cex=0.9,
         insig='label_sig',
         pch.col="black",
         col=colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(100))
dev.off()

results_2 <- results %>% dplyr::filter(response_pred!="NT")
results_2$response_pred <- ifelse(results_2$response_pred=="R","Responder","Non-responder")
colors <- c("Responder"="#5ec962","Non-responder"="#440154")
results_3 <- results_2
results_3$pred_bin <- ifelse(results_3$n_days_pred==0,"0",
                             ifelse(results_3$n_days_pred==1,"1",
                                    ifelse(results_3$n_days_pred>7,">7","2-7")))
results_3$pred_bin <- factor(results_3$pred_bin,levels=c("0","1","2-7",">7"))
results_3$response_pred <- factor(results_3$response_pred,levels=c("Responder","Non-responder"))
plot <- ggplot(results_3,aes(pred_bin,Th17.1_selective,fill=response_pred))+
  geom_boxplot(outlier.shape=NA,alpha=0.6)+
  geom_point(shape=21,position=position_jitterdodge())+
  scale_fill_manual(labels=c("Responder","Non-responder"),values=colors)+
  theme_classic()+
  ylab("GSVA score Th17.1 signature")+
  xlab("Days on steroids upon endoscopy")+
  labs(fill=NULL)
pdf("Th17_1_selective_longitudinal.pdf",height=4,width=5)
plot(plot)
dev.off()

#######################################################################
## ALL scRNAseq re-analyses have been performed in a separate script ##
#######################################################################





