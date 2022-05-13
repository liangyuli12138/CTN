
# 2021-3-18
# delete one sample 
# and use double let 

# 1. delete sample
# CT1 sample is not OK
#=======================================================================================================================

filter <- function(dat){
    dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
    dat = NormalizeData(object = dat)
    dat <- FindVariableFeatures(object = dat)
    all.genes <- rownames(x = dat)
    dat <- ScaleData(object = dat, features = all.genes)
    # Run PCA
    dat <- RunPCA(object = dat, features = VariableFeatures(object = dat),npcs = 50)
    tmp_res <- 2
    if(ncol(dat)>10000){
        tmp_res <- 2
    }else{
        tmp_res <- 1
    }

    dat <- FindNeighbors(dat, dims = 1:20)
    dat <- FindClusters(dat, resolution = tmp_res)
# the tutorial says should use the same dims in Findcluster function
    dat <- RunTSNE(object = dat,dims = 1:20)
    return(dat)
}

# 0.1 PNAS samples 
library(data.table)
library(Seurat)
tmp.count <- as.data.frame(fread("/public/workspace/lily/CTN/PNAS_data/01.UMI.txt",sep="\t"))
tmp.ann <- as.data.frame(fread("/public/workspace/lily/CTN/PNAS_data/03.Cell.Barcodes.txt",sep="\t",header=F))

# gene name ann 
ann <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/hg19_ensemble.txt",sep="\t",header=T)
tmp.data <- merge(ann[,1:2],tmp.count,by.x="Gene.stable.ID",by.y="V1")
tmp.data$Gene.stable.ID <- NULL
tmp.res <- aggregate(.~Gene.name,data=tmp.data,FUN=max)
rownames(tmp.res) <- tmp.res$Gene.name
tmp.res$Gene.name <- NULL

colnames(tmp.ann) <- c("cells","sample","group")
rownames(tmp.ann) <- tmp.ann$cells
tmp.ann <- tmp.ann[colnames(tmp.res),]

# make a Seurat Object 

tmp.dat <- CreateSeuratObject(
counts=tmp.res,
project = "PNAS",
assay = "RNA",
names.field = 1,
names.delim = "_",
meta.data = tmp.ann
)

tmp.dat.f <- subset(tmp.dat,cells=which(tmp.dat$sample%in%c("CT2","CT3","CT4","CT5","SC1","SC2","SC3","SC4","SC5","SC6","SC7")))

prepare <- function(tmp_dat){

    tmp_dat[["percent.mt"]] <- PercentageFeatureSet(object = tmp_dat, pattern = "^MT-")
    tmp_dat = subset(x=tmp_dat,subset=nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
# seurat object
	tmp_dat = NormalizeData(object = tmp_dat)
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	# PCA
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	# clustering
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:10)
	# select proper resolution
	tmp_dat <- FindClusters(object = tmp_dat,resolution=0.8)
	# T-SNE
	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10,check_duplicates = FALSE)
	tmp_dat <- RunUMAP(tmp_dat,dims=1:10)
	return(tmp_dat)
}

 rs.f <- prepare(tmp.dat.f)
saveRDS(rs.f,file="/public/workspace/lily/CTN/version_3_20/PNAS.RDS")




# ourself data 
library(Seurat) 
filelist <- c("/public/workspace/lily/CTN/filtered_feature_bc_matrix_XWS3",
    "/public/workspace/lily/CTN/filtered_feature_bc_matrix_XWS1",
    "/public/workspace/lily/CTN/filtered_feature_bc_matrix_XWS2")
samplelist <- c("XWS3","XWS1","XWS2")
for(i in 1:length(filelist)){
    filepath <- filelist[i]
    tmp <- Read10X(data.dir = filepath)
    # dir.create(paste0("/public/workspace/lily/CTN/0_prepare/",samplelist[i]))
    # respath <- paste0("/public/workspace/lily/CTN/0_prepare/",samplelist[i],"/")
    # do some filter 
    name <- samplelist[i] # get file name
    dat <- CreateSeuratObject(counts = tmp,  project = name,min.cells = 3, min.features = 200)
    dat <- prepare(dat)
    saveRDS(dat,file=paste0("/public/workspace/lily/CTN/version_3_20/",samplelist[i],".RDS"))
    
}

#=====================================================================================================================
# integration 
#
library(Seurat)
XWS1 <- readRDS("/public/workspace/lily/CTN/version_3_20/XWS1.RDS")
XWS1$group <- "CTN"
XWS1$sample <- "XWS1"
XWS2 <- readRDS("/public/workspace/lily/CTN/version_3_20/XWS2.RDS")
XWS2$group <- "CTN"
XWS2$sample <- "XWS2"
XWS3 <- readRDS("/public/workspace/lily/CTN/version_3_20/XWS3.RDS")
XWS3$group <- "CTN"
XWS3$sample <- "XWS3"
PNAS <- readRDS("/public/workspace/lily/CTN/version_3_20/PNAS.RDS")



integration.anchors <- FindIntegrationAnchors(object.list = c(XWS1,XWS2,XWS3,PNAS))
inte <- IntegrateData(anchorset = integration.anchors)
#FindVariableFeatures
inte <- FindVariableFeatures(inte)
##Scaling the integrateda
all.genes <- rownames(inte)
inte <- ScaleData(inte, features = all.genes)
#PCA
inte <- RunPCA(inte)
#cluster
inte <- FindNeighbors(inte)
inte <- FindClusters(inte)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)
inte <- RunUMAP(inte,dims=1:10)
saveRDS(inte,file="/public/workspace/lily/CTN/version_3_20/CTN_inte.RDS")






##############################################################################################################################################
# 2021-3-22
# define cell type 
#=============================================================================================================================================
# some result plot should in rs_plot file 
# data file in data file 
#=============================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")
DefaultAssay(dat) <- "RNA"


# 2021-3-22
# run Doublets Finder to delete some cells that is unclassify 
#==============================================================================================================================================
# run_DoubletFinder <- function(seuratObj) {
#  set.seed(12345)
#  library(DoubletFinder)
#  sweep.res.list <- DoubletFinder::paramSweep_v3(seuratObj, PCs = 1:20, sct = FALSE)
#  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
#  bcmvn <- DoubletFinder::find.pK(sweep.stats)
#  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
#  homotypic.prop <- DoubletFinder::modelHomotypic(seuratObj$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
#  nExp_poi <- round(0.075*length(rownames(seuratObj@meta.data)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
#  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#  dat <- DoubletFinder::doubletFinder_v3(seuratObj, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#  return(dat)
# }

# dat.f <- run_DoubletFinder(dat)




pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/FeaturePlot_umap.pdf")
DimPlot(dat,label=T)
FeaturePlot(dat,features=c("CD3D","CD3E","CD4","CD8A"),label=T) # T cell 
FeaturePlot(dat,features=c("MS4A1","CD79A","CD19","IGHG1"),label=T) # B cell  Plasma
FeaturePlot(dat,features=c("CTSS","FCN1","S100A8","S100A9","LYZ","VCAN"),label=T) # Monocyte
FeaturePlot(dat,features=c("CLEC10A","CD1C","CLEC4C","LAMP3"),label=T) # DC 
FeaturePlot(dat,features=c("NCAM1","NKG7","GNLY","KLRD1"),label=T) # NK 
dev.off()


pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/FeaturePlot_tsne.pdf")
DimPlot(dat,label=T,reduction="tsne")
FeaturePlot(dat,features=c("CD3D","CD3E","CD4","CD8A"),label=T,reduction="tsne") # T cell 
FeaturePlot(dat,features=c("MS4A1","CD79A","CD19","IGHG1"),label=T,reduction="tsne") # B cell  Plasma
FeaturePlot(dat,features=c("CD14", "LYZ","FCGR3A", "MS4A7","CD14"),label=T,reduction="tsne") # Monocyte
FeaturePlot(dat,features=c("CLEC10A","CD1C","CLEC4C","LAMP3"),label=T,reduction="tsne") # DC 
FeaturePlot(dat,features=c("NCAM1","NKG7","GNLY","KLRD1"),label=T,reduction="tsne") # NK 
FeaturePlot(dat,features=c("TOP2A","MKI67"),label=T,reduction="tsne") # Prolifering Cells
FeaturePlot(dat,features=c("GP9","PF4"),label=T,reduction="tsne") # Platelet
dev.off()

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/VlnPlot_tsne.pdf",width=12)

VlnPlot(dat,features=c("CD3D","CD3E","CD4","CD8A"),pt.size=0) # T cell 
VlnPlot(dat,features=c("MS4A1","CD79A","CD19","IGHG1"),pt.size=0) # B cell  Plasma
VlnPlot(dat,features=c("CD14", "LYZ","FCGR3A", "MS4A7","CD14"),pt.size=0) # Monocyte
VlnPlot(dat,features=c("CLEC10A","CD1C","CLEC4C","LAMP3"),pt.size=0) # DC 
VlnPlot(dat,features=c("NCAM1","NKG7","GNLY","KLRD1"),pt.size=0) # NK 
VlnPlot(dat,features=c("TOP2A","MKI67"),pt.size=0) # Prolifering Cells
VlnPlot(dat,features=c("GP9","PF4"),pt.size=0) # Platelet
dev.off()





#==================================================================================================
# define celltype and plot result 
#==================================================================================================
dat$celltype <- "unclassify"
dat$celltype[which(dat$seurat_clusters%in%c(0,2,6,7,8,9,13,14,19,21))] <- "T cell"
dat$celltype[which(dat$seurat_clusters%in%c(11,12,22))] <- "B cell"
dat$celltype[which(dat$seurat_clusters%in%c(1,3))] <- "NK cell"
dat$celltype[which(dat$seurat_clusters%in%c(4,5,14,16,23))] <- "CD14 Mono"
dat$celltype[which(dat$seurat_clusters%in%c(10))] <- "CD16 Mono"
dat$celltype[which(dat$seurat_clusters%in%c(15,17))] <- "Platelet"
dat$celltype[which(dat$seurat_clusters%in%c(18))] <- "DC"
dat$celltype[which(dat$seurat_clusters%in%c(24))] <- "pDC"
dat$celltype[which(dat$seurat_clusters%in%c(20))] <- "Proliferating cell"


tmp <- table(dat$sample,dat$celltype)
tmp.res <- apply(tmp,1,function(x){x/sum(x)})
pdf("tmp.pdf")
for(i in 1:9){
    boxplot(tmp.res[i, 1:4], tmp.res[i, 12:14],tmp.res[i, 5:11],main=rownames(tmp.res)[i],names=c("Elder","CTN","SCT"))
}
dev.off()

apply(tmp.res,1,function(x){
    c(wilcox.test(x[1:4],x[12:14])$p.value,wilcox.test(x[1:4],x[5:11])$p.value)
})

apply(tmp.res,1,function(x){
    c(median(x[1:4]),median(x[12:14]),median(x[5:11]))
})









###################################################################################################################################################
# 2021-3-23
# 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")
DefaultAssay(dat) <- "RNA"

#=================================================================================================================================================
# my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
#          '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
#          '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
#          '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
#          '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
#          '#968175')

# library(Seurat)
# library(ggplot2)
# modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
#        p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
#                xlab("") + ylab(feature) + ggtitle("") +
#                theme(legend.position = "none",
#                axis.text.x = element_blank(),
#                axis.text.y = element_blank(),
#                axis.ticks.x = element_blank(),
#                axis.ticks.y = element_line(),
#                axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
#                plot.margin = plot.margin )
#        return(p)
# }

# ## main function
# StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
#        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
#             plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
#             theme(axis.text.x=element_text(), axis.ticks.x = element_line())
#        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
#        return(p)
# }

# #
# dat@active.ident <- factor(dat$celltype)
# markergene <- c("CD3D","CD3E","MS4A1","CD79A","CD14","MS4A7",
#     "CLEC10A","CLEC4C","SERPINF1","TCF4","NKG7","KLRD1",
#     "TOP2A","MKI67","GP9","PF4")

# # Fig1
# pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig1_violinplot.pdf",width=10,height=10)
# StackedVlnPlot(dat, markergene[1:8], pt.size=0, cols=my36colors)
# StackedVlnPlot(dat, markergene[9:16], pt.size=0, cols=my36colors)
# dev.off()















#=====================================================================================================================================
# 2021-3-24
# try to use doubletFinder
######################################################################################################################################

run_DoubletFinder <- function(seuratObj) {
 set.seed(12345)
 library(DoubletFinder)
 sweep.res.list <- DoubletFinder::paramSweep_v3(seuratObj, PCs = 1:20, sct = FALSE)
 sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
 bcmvn <- DoubletFinder::find.pK(sweep.stats)
 mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
 homotypic.prop <- DoubletFinder::modelHomotypic(seuratObj$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
 nExp_poi <- round(0.075*length(rownames(seuratObj@meta.data)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
 nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
 dat <- DoubletFinder::doubletFinder_v3(seuratObj, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
 return(dat)
}


library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/PNAS.RDS")
dat.f <- run_DoubletFinder(dat)

library(Seurat)
XWS1 <- readRDS("/public/workspace/lily/CTN/version_3_20/data/XWS1.RDS")
XWS1$group <- "CTN"
XWS1$sample <- "XWS1"
XWS1.f <- run_DoubletFinder(XWS1)
XWS2 <- readRDS("/public/workspace/lily/CTN/version_3_20/data/XWS2.RDS")
XWS2$group <- "CTN"
XWS2$sample <- "XWS2"
XWS2.f <- run_DoubletFinder(XWS2)
XWS3 <- readRDS("/public/workspace/lily/CTN/version_3_20/data/XWS3.RDS")
XWS3$group <- "CTN"
XWS3$sample <- "XWS3"
XWS3.f <- run_DoubletFinder(XWS3)


# > table(XWS1.f$DF.classifications_0.25_0.24_773)

# Doublet Singlet
#     773   10601

# > table(XWS2.f$DF.classifications_0.25_0.005_917)

# Doublet Singlet
#     917   12955
# > table(XWS3.f$DF.classifications_0.25_0.005_920)

# Doublet Singlet
#     920   12506

# table(dat.f$DF.classifications_0.25_0.005_3734)

# Doublet Singlet
#    3734   50910



#=========================================================================
# save data 
tmp.1 <- subset(XWS1.f,cells=which(XWS1.f$DF.classifications_0.25_0.24_773=="Singlet"))
saveRDS(tmp.1,file="/public/workspace/lily/CTN/version_3_20/tmp/XWS1.RDS")

tmp.2 <- subset(XWS2.f,cells=which(XWS2.f$DF.classifications_0.25_0.005_917=="Singlet"))
saveRDS(tmp.2,file="/public/workspace/lily/CTN/version_3_20/tmp/XWS2.RDS")

tmp.3 <- subset(XWS3.f,cells=which(XWS3.f$DF.classifications_0.25_0.005_920=="Singlet"))
saveRDS(tmp.3,file="/public/workspace/lily/CTN/version_3_20/tmp/XWS3.RDS")



####################################################################################################################################
library(Seurat)

# XWS3$group <- "CTN"
# XWS3$sample <- "XWS3"
# FindDoublelets in each sample 
PNAS <- readRDS("/public/workspace/lily/CTN/version_3_20/data/PNAS.RDS")

run_DoubletFinder <- function(seuratObj) {
 set.seed(12345)
 library(DoubletFinder)
 sweep.res.list <- DoubletFinder::paramSweep_v3(seuratObj, PCs = 1:20, sct = FALSE)
 sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
 bcmvn <- DoubletFinder::find.pK(sweep.stats)
 mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
 homotypic.prop <- DoubletFinder::modelHomotypic(seuratObj$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
 nExp_poi <- round(0.075*length(rownames(seuratObj@meta.data)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
 nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
 dat <- DoubletFinder::doubletFinder_v3(seuratObj, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
 return(dat)
}

inte.list <- list() 
samplelist <- unique(PNAS$sample)
for(i in 1:length(samplelist)){
	tmp <- subset(PNAS,cells=which(PNAS$sample==samplelist[i]))
	DefaultAssay(tmp) <- "RNA"
    tmp <- run_DoubletFinder(tmp)
    tmp.f <- subset(tmp,cells=which(tmp@meta.data[,ncol(tmp@meta.data)]=="Singlet"))
	inte.list[i] <- tmp.f
}


XWS1 <- readRDS("/public/workspace/lily/CTN/version_3_20/tmp/XWS1.RDS")
# XWS1$group <- "CTN"
# XWS1$sample <- "XWS1"
XWS2 <- readRDS("/public/workspace/lily/CTN/version_3_20/tmp/XWS2.RDS")
# XWS2$group <- "CTN"
# XWS2$sample <- "XWS2"
XWS3 <- readRDS("/public/workspace/lily/CTN/version_3_20/tmp/XWS3.RDS")

inte.list[12] <- XWS1
inte.list[13] <- XWS2
inte.list[14] <- XWS3

integration.anchors <- FindIntegrationAnchors(object.list = inte.list)
inte <- IntegrateData(anchorset = integration.anchors)
#FindVariableFeatures
inte <- FindVariableFeatures(inte)
##Scaling the integrateda
all.genes <- rownames(inte)
inte <- ScaleData(inte, features = all.genes)
#PCA
inte <- RunPCA(inte)
#cluster
inte <- FindNeighbors(inte)
inte <- FindClusters(inte)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)
# inte <- RunUMAP(inte,dims=1:10)
saveRDS(inte,file="/public/workspace/lily/CTN/version_3_20/tmp/CTN_inte.RDS")



# 2021-4-13 
# use resolution 1 to classify T cell and NK cell 
dat <- inte
DefaultAssay(dat) <- "RNA"

pdf("/public/workspace/lily/CTN/version_3_20/tmp_FeaturePlot_tsne.pdf")
DimPlot(dat,label=T,reduction="tsne")
FeaturePlot(dat,features=c("CD3D","CD3E","KLRB1","KLRD1"),label=T,reduction="tsne") # T cell 
FeaturePlot(dat,features=c("MS4A1","CD79A","CD19","IGHG1"),label=T,reduction="tsne") # B cell  Plasma
FeaturePlot(dat,features=c("CD14", "LYZ","FCGR3A", "MS4A7","CD14"),label=T,reduction="tsne") # Monocyte
FeaturePlot(dat,features=c("CLEC10A","CD1C","CLEC4C","LAMP3"),label=T,reduction="tsne") # DC 
FeaturePlot(dat,features=c("NCAM1","NKG7","GNLY","KLRD1"),label=T,reduction="tsne") # NK 
FeaturePlot(dat,features=c("TOP2A","MKI67"),label=T,reduction="tsne") # Prolifering Cells
FeaturePlot(dat,features=c("GP9","PF4"),label=T,reduction="tsne") # Platelet
dev.off()

pdf("/public/workspace/lily/CTN/version_3_20/tmp_VlnPlot_tsne.pdf",width=12)

VlnPlot(dat,features=c("CD3D","CD3E","KLRB1","KLRD1"),pt.size=0) # T cell 
VlnPlot(dat,features=c("MS4A1","CD79A","CD19","IGHG1"),pt.size=0) # B cell  Plasma
VlnPlot(dat,features=c("CD14", "LYZ","FCGR3A", "MS4A7","CD14"),pt.size=0) # Monocyte
VlnPlot(dat,features=c("CLEC10A","CD1C","CLEC4C","LAMP3"),pt.size=0) # DC 
VlnPlot(dat,features=c("NCAM1","KLRB1","GNLY","KLRD1"),pt.size=0) # NK 
VlnPlot(dat,features=c("TOP2A","MKI67"),pt.size=0) # Prolifering Cells
VlnPlot(dat,features=c("GP9","PF4"),pt.size=0) # Platelet
dev.off()


# 2021-4-13
# re classify cell types  [use seurat_cluster as resolution is 1]
dat$celltype <- "unclassify"
dat$celltype[which(dat$seurat_clusters%in%c(0,2,3,7,10,13,14,21))] <- "T cell"
dat$celltype[which(dat$seurat_clusters%in%c(11,15,23))] <- "B cell"
dat$celltype[which(dat$seurat_clusters%in%c(1,4,6,18))] <- "NK cell"
dat$celltype[which(dat$seurat_clusters%in%c(5,8,9,12,24))] <- "CD14 Mono"
dat$celltype[which(dat$seurat_clusters%in%c(16))] <- "CD16 Mono"
dat$celltype[which(dat$seurat_clusters%in%c(17,19))] <- "Platelet"
# dat$celltype[which(dat$seurat_clusters%in%c(18))] <- "DC"
dat$celltype[which(dat$seurat_clusters%in%c(22))] <- "pDC"
dat$celltype[which(dat$seurat_clusters%in%c(25))] <- "Proliferating cell"

# use single R to make sure 

# library(SingleR)

# # ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_HumanPrimaryCellAtlasData.RDS")
# # ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_NovershternHematopoieticData.RDS")

# ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_BluePrint.RDS")
# res.cluster<-SingleR(test=as_matrix(dat@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=dat$seurat_clusters,method="cluster")


# ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_MonacoImmuneData.RDS")
# res.cluster1<-SingleR(test=as_matrix(dat@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=dat$seurat_clusters,method="cluster")


#==========================================================================================================================================
# check some result 
        #     B cell          CD14 Mono          CD16 Mono            NK cell
        #       5163              15911               2391              21019
        #        pDC           Platelet Proliferating cell             T cell
        #        154               2300                150              39945
        # unclassify
        #        182
# col <- c("#f47b7b","#91be3e","#96cbb3","#836eaa","#39a6dd","#207c88","#e5352b","#ef9020","#dde2e0","#d4d7da")

# DimPlot(dat,group.by="celltype",cols=col)

tmp <- table(dat$sample,dat$celltype)
tmp.res <- apply(tmp,1,function(x){x/sum(x)})
pdf("tmp.pdf")
for(i in 1:9){
    boxplot(tmp.res[i, 1:4], tmp.res[i, 12:14],tmp.res[i, 5:11],main=rownames(tmp.res)[i],names=c("Elder","CTN","SCT"))
}
dev.off()

# # apply(tmp.res,1,function(x){
# #     c(wilcox.test(x[1:4],x[12:14])$p.value,wilcox.test(x[1:4],x[5:11])$p.value)
# # })

# rs.f <- apply(tmp.res,1,function(x){
#     c(median(x[1:4]),median(x[12:14]),median(x[5:11]))
# })


# dat.f <- rs.f[,c(1,2,5,9)]
# rownames(dat.f) <- c("Elder","CTN","SCT")
























##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################

# 2021-3-25
# Finally lly use this integration version to analysis 
# Find doublets in each sample and integration with each sample 
# and change tmp file into doublet find residentional file names as "doublets_rs"
dat@meta.data <- dat@meta.data[,c(1:9,38,39)]
saveRDS(dat,file="/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")














