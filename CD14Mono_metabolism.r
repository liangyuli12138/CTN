#!/usr/bin/Rscript 
# 2021-3-25
# metabolism analysis CTN T cell 
#=============================================================================================================================
################################### Metabolic#########################
library(Seurat)
library(scater)
library(stringr)
library("Rtsne")
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(scales)
library(ggplot2)
library(dplyr)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)

source('/public/workspace/lily/software/SingleCellMetabolic/utils.R')
source('/public/workspace/lily/software/SingleCellMetabolic/runGSEA_preRank.R')


# load a as_matrix function 

as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
    
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}


# read data
# all cell is tumor
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CD14Mono.RDS")
tmp_data=dat
all_data=as_matrix(tmp_data[['RNA']]@data) # make sure your data is correct size which can use as.matrix
cell_type=paste0("C",as.vector(tmp_data$seurat_clusters))  # cell type which means cluster or cell type 
tumor=unname(tmp_data$seurat_clusters) # this do not need to change ,this means all cells need to analysis


#=======================================================================================================================
col_data <- data.frame(tumor=tumor,cellType=as.character(cell_type),row.names=colnames(all_data))
pathways <- gmtPathways("/public/workspace/lily/software/SingleCellMetabolic/Data/KEGG_metabolism.gmt")
metabolics <- unique(as.vector(unname(unlist(pathways))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(all_data)),row.names = rownames(all_data))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE
#set SingleCellExperiment object which include expr matrix ,cluster info ,and gene info
sce <- SingleCellExperiment(
  assays = all_data,
  colData = col_data,
  rowData = row_data
)

selected_tumor_sce <- sce #the example code use "selected_tumor_sce" as the name
selected_tumor_metabolic_sce <- sce[rowData(sce)$metabolic,] # dims :1506 37637; 1506 metabolic genes


###################################### scRNA_pathway_activity ##########################################
pathway_file <- "/public/workspace/lily/software/SingleCellMetabolic/Data/KEGG_metabolism.gmt"
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
all_cell_types <- as.vector(selected_tumor_metabolic_sce$cellType)
cell_types <- unique(all_cell_types)

gene_pathway_number <- num_of_pathways(pathway_file,rownames(selected_tumor_metabolic_sce)[rowData(selected_tumor_metabolic_sce)$metabolic])
set.seed(123)
normalization_method <- "Deconvolution"

##Calculate the pathway activities
#mean ratio of genes in each pathway for each cell type
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
###calculate the pvalues using shuffle method
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = (list(pathway_names, cell_types)))

norm_tpm <- all_data


for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 5) next

  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))

  #remove genes which are zeros in any celltype to avoid extreme ratio value
  keep <- colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType>0.001)]

  if(length(keep)<3) next
  
  #using the loweset value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))

  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachCellType,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachCellType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 3) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  mean_exp_pathway <- apply(ratio_exp_eachCellType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle[p, ] <-  mean_exp_pathway[cell_types]
  mean_expression_noshuffle[p, ] <-  mean_exp_pathway[cell_types]
    
  ##shuffle 5000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(cell_types,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_cell_types_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachCellType_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:5000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_cell_types_list <- lapply(times,function(x) sample(all_cell_types)) 
  names(shuffle_cell_types_list) <- times
  mean_exp_eachCellType_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachCellType_list <- lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(cell_types),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- cell_types
  for(c in cell_types){
    if(is.na(mean_expression_shuffle[p,c])) next
    if(mean_expression_shuffle[p,c]>1){
      pval <- sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) / 5000 
    }else if(mean_expression_shuffle[p,c]<1){
      pval <- sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) / 5000
    }
    if(pval>0.01) mean_expression_shuffle[p, c] <- NA  ### NA is  blank in heatmap
    pvalues_mat[p,c] <- pval
  }
}
all_NA <- rowAlls(is.na(mean_expression_shuffle))
mean_expression_shuffle <- mean_expression_shuffle[!all_NA,]
#heatmap
dat <- mean_expression_shuffle

sort_row <- c()
sort_column <- c()

for(i in colnames(dat)){
  select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
  tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}
sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column <- names(sort_column)
dat[is.na(dat)] <- 1

#================================================================================
write.table(mean_expression_noshuffle,file="/public/workspace/lily/CTN/version_3_20/metabolism/CD14mono/KEGGpathway_activity_noshuffle.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(mean_expression_shuffle,file="/public/workspace/lily/CTN/version_3_20/metabolism/CD14mono/KEGGpathway_activity_shuffle.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(pvalues_mat,file="/public/workspace/lily/CTN/version_3_20/metabolism/CD14mono/KEGGpathway_activity_shuffle_pvalue.txt",row.names=T,col.names=T,quote=F,sep="\t")








