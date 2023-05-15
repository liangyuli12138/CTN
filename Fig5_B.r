


# Fig4
# B cell analysis 
# 2021-3-24
#===================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")
Bcell <- subset(dat,cells=which(dat$celltype=="B cell"))

# divide into each sample and re-integration
#===================================================================================================================

inte.list <- list()
samplelist <- unique(Bcell$sample)
for(i in 1:length(samplelist)){
    tmp <- subset(Bcell,cells=which(Bcell$sample==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list,dims=1:20,k.filter = 10,k.score = 10)
inte <- IntegrateData(anchorset = integration.anchors,dims = 1:20,k.weight = 50)
#FindVariableFeatures
inte <- FindVariableFeatures(inte)
##Scaling the integrateda
all.genes <- rownames(inte)
inte <- ScaleData(inte, features = all.genes)
#PCA
inte <- RunPCA(inte)
#cluster
inte <- FindNeighbors(inte)
inte <- FindClusters(inte,resolution=0.5) # for B cell should change into 0.5
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)

saveRDS(inte,file="/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")







library(SingleR)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")
# ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_HumanPrimaryCellAtlasData.RDS")
# ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_NovershternHematopoieticData.RDS")

ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_BluePrint.RDS")
res.cluster<-SingleR(test=as.matrix(dat@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=dat$seurat_clusters,method="cluster")


ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_MonacoImmuneData.RDS")
res.cluster1<-SingleR(test=as.matrix(dat@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=dat$seurat_clusters,method="cluster")



#=========================================================================================================================================
# check sample group in each cluster
##########################################################################################################################################
dat$celltype.refine <- "unclassify"
dat$celltype.refine[which(dat$seurat_clusters%in%c(0,1,5,8))] <- "Naive B cell"
dat$celltype.refine[which(dat$seurat_clusters%in%c(2,4))] <- "Switched memory B cell"
dat$celltype.refine[which(dat$seurat_clusters%in%c(3,6,9))] <- "Memory B cell"
dat$celltype.refine[which(dat$seurat_clusters%in%c(7))] <- "Plasma"

saveRDS(dat,file="/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")






#########################################################################################################################################
# 2021-3-25
# 1. run trajectory 
#========================================================================================================================================
library(Seurat)
library(monocle)

dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")
sub.dat <- subset(dat,cells=colnames(dat)[which(dat$celltype.refine%in%c("Naive B cell","Switched memory B cell","Memory B cell","Plasma"))])

# run trajectory 
DefaultAssay(sub.dat) <- "RNA"
tmp.dat <- Seurat::as.CellDataSet(sub.dat)
tmp.dat <- estimateSizeFactors(tmp.dat)
tmp.dat <- detectGenes(tmp.dat, min_expr = 1)
# fData(dat)$use_for_ordering <-fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
expressed_genes <- row.names(subset(fData(tmp.dat),num_cells_expressed >= 10))

# clustering_DEG_genes <- differentialGeneTest(tmp.dat[expressed_genes,],
#        fullModelFormulaStr = '~celltype.refine',
#        cores = 10)
# ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:3000]

genes <- VariableFeatures(dat)[1:1000]
ordering_genes <- genes
tmp.dat <- setOrderingFilter(tmp.dat, ordering_genes = ordering_genes)
tmp.dat <- reduceDimension(tmp.dat, method = 'DDRTree')
tmp.dat <- orderCells(tmp.dat)

saveRDS(tmp.dat,file="/public/workspace/lily/CTN/version_3_20/trajectory/Bcell_pure_trajectory.RDS")


#=====================================================================================================================================
# 2021-4-22
# calculate trajectory result 
######################################################################################################################################
library(Seurat)
library(monocle)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")
DefaultAssay(dat) <- "integrated"
tmp.dat <- readRDS("/public/workspace/lily/CTN/version_3_20/monocle/Bcell_pure_trajectory.RDS")
genes <- VariableFeatures(dat)[1:1500]
ordering_genes <- genes
tmp.dat <- setOrderingFilter(tmp.dat, ordering_genes = ordering_genes)
tmp.dat <- reduceDimension(tmp.dat, method = 'DDRTree')
tmp.dat <- orderCells(tmp.dat)

# pdf("./tmp.pdf")
# plot_cell_trajectory(tmp.dat,color_by="celltype.refine")+scale_colour_manual(values=c("#ffc845","#00c16e","#f85a40","#7d3f98"))
# plot_cell_trajectory(tmp.dat,color_by="group")
# plot_cell_trajectory(tmp.dat,color_by="Pseudotime")
# plot_cell_trajectory(tmp.dat,color_by="sample") +facet_wrap(~sample, nrow = 3)
# dev.off()

# check trajectory result 
pdf("/public/workspace/lily/CTN/version_3_20/trajectory/Bcell_trajectory.pdf")
plot_cell_trajectory(tmp.dat,color_by="celltype.refine")+scale_colour_manual(values=c("#ffc845","#00c16e","#f85a40","#7d3f98"))
plot_cell_trajectory(tmp.dat,color_by="group")+scale_colour_manual(values=c("#003468","#5ec6f2","#4b1702"))
plot_cell_trajectory(tmp.dat,color_by="sample") +facet_wrap(~sample, nrow = 3)
dev.off()

saveRDS(tmp.dat,file="/public/workspace/lily/CTN/version_3_20/monocle/Bcell_pure_trajectory.RDS")











###########################################################################################################################################
# plot result 
# 2021-4-23
#==========================================================================================================================================
library(Seurat)

# 2021-4-23
# celltype.new 
# lly think we could cbind switch mem and mem B cell into memory B cell.
###########################################################################################################################################
dat$celltype.new <- as.vector(dat$celltype.refine)
dat$celltype.new[which(dat$celltype.new=="Switched memory B cell")] <- "Memory B cell"
saveRDS(dat,file="/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")


#==========================================================================================================================================
# 1. FigA B overview 
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")
dat$celltype.new <- factor(dat$celltype.new,levels=c("Naive B cell","Memory B cell","Plasma","unclassify"))
# naive memory switch-memory plasma
cols <- c("#009f4d","#f48924","#f85a40","#caccd1")
DefaultAssay(dat) <- "RNA"
# marker gene plot 

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig4_overview.pdf",useDingbats=F)
DimPlot(dat,group.by="celltype.new",cols=cols,reduction="tsne")
FeaturePlot(dat,features=c("CD22","TCL1A","MS4A1"),order=T) # naive B cell 
FeaturePlot(dat,features=c("CD38", "TNFRSF17","IGHG1","CD27"),order=T) # plasma and memeory B cell (CD27)
dev.off()


# 2023-1-3
# update to umap
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")
DefaultAssay(dat) <- "integrated"
dat <- RunUMAP(dat,dims=1:10)
DefaultAssay(dat) <- "RNA"

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig4_overview_umap.pdf",useDingbats=F)
DimPlot(dat,group.by="celltype.new",cols=cols,reduction="umap",raster=T)
FeaturePlot(dat,features=c("CD22","TCL1A","MS4A1"),order=T,raster=T) # naive B cell 
FeaturePlot(dat,features=c("CD38", "TNFRSF17","IGHG1","CD27"),order=T,raster=T) # plasma and memeory B cell (CD27)
dev.off()

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/supplmentary/Fig5_Bcell_cluster_umap.pdf",useDingbats=F)
DimPlot(dat,reduction="umap",raster=T)

dev.off()


#==========================================================================================================================================
# 2. different subtype in different group
# tmp <- table(dat$sample,dat$celltype.new)
# tmp.res <- apply(tmp,1,function(x){x/sum(x)})
# rs.f <- apply(tmp.res,1,function(x){
#     c(median(x[1:4]),median(x[12:14]),median(x[5:11]))
# })

library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")

library(reshape2)
library(ggplot2)
library(ggalluvial)
library(ggplot2)
library(RColorBrewer)

tmp.dat <- apply(table(dat$celltype.new,dat$group),2,function(x){x/sum(x)})
tmp.res <- melt(tmp.dat)
colnames(tmp.res) <- c("celltype.new","group","percentage")

cols.t <- c("#009f4d","#f48924","#f85a40","#caccd1")
tmp.res$group <- factor(tmp.res$group,levels=c("Elder","CTN","SC"))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig4_celltype_percent.pdf",useDingbats=F,height=10)
ggplot(tmp.res, aes(x = group, y = percentage, fill = celltype.new,stratum = celltype.new,alluvium = celltype.new)) +
geom_stratum() +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols.t)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()







#===========================================================================================================================================
# 3. trajectory analysis 
library(Seurat)
library(monocle)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/monocle/Bcell_pure_trajectory.RDS")
# dat$celltype.new <- as.vector(dat$celltype.refine)
# dat$celltype.new[which(dat$celltype.new=="Switched memory B cell")] <- "Memory B cell"
# saveRDS(dat,file="/public/workspace/lily/CTN/version_3_20/monocle/Bcell_pure_trajectory.RDS")

# 2021-9-16
# add group info 
dat$groupinfo <- "Unknow"
dat$groupinfo[which(dat$seurat_clusters%in%c(0,1,5,8))] <- "Elder"
dat$groupinfo[which(dat$seurat_clusters%in%c(6,9))] <- "CTN"
dat$groupinfo[which(dat$seurat_clusters%in%c(2,3,4,7))] <- "CTN.SC"




dat$celltype.new <- factor(dat$celltype.new,levels=c("Naive B cell","Memory B cell","Plasma"))
pdf("/public/workspace/lily/CTN/version_3_20/monocle/Bcell_trajectory.pdf",useDingbats=F)
plot_cell_trajectory(dat,color_by="celltype.new")+scale_colour_manual(values=c("#009f4d","#f48924","#f85a40"))
plot_cell_trajectory(dat,color_by="Pseudotime")
plot_cell_trajectory(dat,color_by="groupinfo")
plot_cell_trajectory(dat,color_by="group")
plot_cell_trajectory(dat,color_by="sample") +facet_wrap(~sample, nrow = 3)

dev.off()










#=========================================================================================================================================
# 4. calculate different cluster three group percentage and do DEG analysis 
# 
library(Seurat)
library(reshape2)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")
tmp.res <- apply(apply(table(dat$group,dat$seurat_clusters),1,function(x){x/sum(x)}),1,function(x){x/sum(x)})
colnames(tmp.res) <- paste0("C",colnames(tmp.res))
tmp.res.f <- tmp.res[,-11]
plot.data <- melt(tmp.res.f)
colnames(plot.data) <- c("group","Cluster","percentage")

library(ggplot2)
cols <- c("#5091cd","#7ac143","#f9a541")
plot.data$Cluster <- factor(plot.data$Cluster,levels=c("C0","C1","C5","C8","C2","C3","C4","C7","C6","C9"))
plot.data$group <- factor(plot.data$group,levels=c("Elder","CTN","SC"))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig4_cluster_percent.pdf",width=10,height=8)
ggplot(plot.data,aes(x=Cluster,y=percentage,fill=group,group=group))+ geom_bar(stat="identity",position="dodge")+
    scale_fill_manual(values=cols)+theme_bw()+ labs(y="percentage")


ggplot(plot.data,aes(x=Cluster,y=percentage,fill=group,group=group))+ geom_bar(stat="identity",position="stack")+
    scale_fill_manual(values=cols)+theme_bw()+ labs(y="percentage")

dev.off()



# dat$celltype.refine <- "unclassify"
# dat$celltype.refine[which(dat$seurat_clusters%in%c(0,1,5,8))] <- "Naive B cell"
# dat$celltype.refine[which(dat$seurat_clusters%in%c(3,6,9,2,4))] <- "Memory B cell"
# dat$celltype.refine[which(dat$seurat_clusters%in%c(7))] <- "Plasma"







#============================================================================================================================================
# 5. DEG analysis
# 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")
gene <- FindMarkers(dat,assay="RNA",group.by="seurat_clusters",ident.1 =c(2,3,4,7),ident.2=c(0,1,5,8))
tmp.gene <- gene[which(gene$p_val_adj<0.05),]
tmp.gene$logP <- tmp.gene$p_val_adj
tmp.gene$logP[which(tmp.gene$logP==0)] <- min(tmp.gene$logP[-which(tmp.gene$logP==0)])
tmp.gene$logP <- -log(tmp.gene$logP,100) # use log 100
tmp.gene$class <- "grey"
tmp.gene$class[which(tmp.gene$avg_logFC> 0.5)] <- "up"
tmp.gene$class[which(tmp.gene$avg_logFC< (-0.5))] <- "down"
tmp.gene$class <- factor(tmp.gene$class,levels=c("up","down","grey"))
library(ggplot2)
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig4_DEG.pdf",useDingbats=F)
ggplot(tmp.gene,aes(x=avg_logFC,y=logP,group=class,color=class)) + geom_point(size=1.5) + 
    scale_colour_manual(values=c("#ec2c22","#2e9df7","#d7d7d8")) +  geom_vline(aes(xintercept= -0.5),color="#2e9df7") +
     geom_vline(aes(xintercept= 0.5),color="#ec2c22")
dev.off()

gene.up <- rownames(tmp.gene)[which(tmp.gene$class=="up")]
write.table(gene.up,file="/public/workspace/lily/CTN/version_3_20/rs_plot/Fig4_DEG_SCT.CTN.txt",quote=F,row.names=F,col.names=F)
#============================================================================================================================================
# KEGG analysis by enrichR 
# different result maybe because 2018 Go database and 2021 Go database

dat <- read.table("/public/workspace/lily/CTN/version_3_20/EnrichR/GO_Biological_Process_2021_Bcell_SCT.CNT_elder.txt",sep="\t",header=T)
library(ggplot2)
dat$logP <- -log10(dat$Adjusted.P.value)
sub.dat <- head(dat,10)
sub.dat <- sub.dat[order(sub.dat$logP,decreasing=T),]
sub.dat$Term <- factor(sub.dat$Term,levels=rev(sub.dat$Term))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/SC_CTN_Bcell_DEG_GO_enrich.pdf",useDingbats=F,width=10)
ggplot(sub.dat,aes(x=Term,y=Odds.Ratio,fill=logP))+geom_bar(stat="identity",position = position_dodge(0.5),width=0.5) + coord_flip()
dev.off()

















#=============================================================================================================================================
# 6. transcript factor 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")
tf.dat <- read.table("/public/workspace/lily/CTN/version_3_20/SCENIC/Bcell/all_Bcell_auc_mtx.tsv",sep=",",header=T)
rownames(tf.dat) <- tf.dat$Cell
tf.dat$Cell <- NULL
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))
all(rownames(tf.dat)==rownames(dat@meta.data))

# add a column information and plot pheatmap 
#==============================================================================================================================================
#tf.dat$celltype <- dat$celltype.new
tf.dat$cluster <- paste0("C",dat$seurat_clusters)
tf.dat$celltype <- NULL
tmp.res <- aggregate(.~cluster,data=tf.dat,FUN=mean)
rownames(tmp.res) <- tmp.res$cluster
tmp.res$cluster <- NULL
tmp.res.f <- tmp.res[c("C0","C1","C5","C8","C2","C3","C4","C7","C6","C9"),] # C10 is unclassify 
# tmp.res.f <- tmp.res.f[,-which(colSums(tmp.res.f)==0)]  # some TF is all zero in all sample

annotation_row=data.frame(row.names=c("C0","C1","C5","C8","C2","C3","C4","C7","C6","C9"),
    group=c(rep("Elder",4),rep("SC.CTN",4),rep("CTN",2)),
    type=c(rep("Naive B",4),rep("Memory B cell",3),"Plasma",rep("Memory B cell",2)))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig4_heatmap_TFs.pdf",width=15,height=8)
pheatmap::pheatmap(tmp.res.f,scale="column",cluster_rows=F,annotation_row=annotation_row,color=colorRampPalette(c('steelblue','white','red'))(50))
dev.off()


# target gene 
p <- pheatmap::pheatmap(tmp.res.f,scale="column",cluster_rows=F,annotation_row=annotation_row,color=colorRampPalette(c('steelblue','white','red'))(50))
# just use SC.CTN up TFs 
tf.dat <- data.frame(TF=colnames(tmp.res.f)[p$tree_col$order][1:15],group=c(rep("SC.CTN",15)))

# TRUUST database 
trrust.dat <- read.table('/public/workspace/lily/metastasis/data/TRRUST/trrust_rawdata.human.tsv',sep="\t") 
colnames(trrust.dat) <- c("source","traget","effect","id")

tmp.tf <- merge(trrust.dat,tf.dat,by.x="source",by.y="TF")

#===============================================================================================================
tmp.up <- tmp.tf[which(tmp.tf$effect=="Activation"),]
write.table(as.vector(tmp.up$traget),file="/public/workspace/lily/CTN/version_3_20/rs_plot/Bcell_TF.up.gene.txt",sep="\t",row.names=F,col.names=T,quote=F)

# do enrich analysis in Enrich R target gene 


dat <- read.table("/public/workspace/lily/CTN/version_3_20/EnrichR/GO_Biological_Process_2021_Bcell_Up_TF_target_gene.txt",sep="\t",header=T)
library(ggplot2)
dat$logP <- -log10(dat$Adjusted.P.value)
sub.dat <- head(dat,10)
sub.dat <- sub.dat[order(sub.dat$logP,decreasing=T),]
sub.dat$Term <- factor(sub.dat$Term,levels=rev(sub.dat$Term))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/SC_CTN_Bcell_upTF_target_gene_GO_enrich.pdf",useDingbats=F,width=10)
ggplot(sub.dat,aes(x=Term,y=Odds.Ratio,fill=logP))+geom_bar(stat="identity",position = position_dodge(0.5),width=0.5) + coord_flip()
dev.off()





# 2021-9-7
# check which TFs in SCT.TF have highest 
tmp <- head(sort(table(tmp.up$traget),decreasing=T),20)

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Bcell_TFs_target_gene.pdf",useDingbats=F)
barplot(tmp,las=2,width=0.6)
dev.off()
















#=======================================================================================================================================
# 7. metabolism analysis
########################################################################################################################################
tmp <- read.table("/public/workspace/lily/CTN/version_3_20/metabolism/Bcell/KEGGpathway_activity_shuffle_OV.txt",sep="\t",header=T)
tmp.a <- tmp
tmp.a[is.na(tmp.a)] <- 1

tmp.f <- tmp.a[,-grep("C10",colnames(tmp.a))]
tmp.f <- tmp.f[,c("C0","C1","C5","C8","C2","C3","C4","C7","C6","C9")]
tmp.res <- t(apply(tmp.f,1,function(x){
    c(mean(x[1:4]),mean(x[c(5:8)]),mean(x[c(9:10)]))
}))
colnames(tmp.res) <- c("Elder","SC.CTN","CTN")


pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig4_Bcell_metabolism.pdf",height=12,width=8)
pheatmap::pheatmap(tmp.res,scale="row",color=colorRampPalette(c('steelblue','white','red'))(50))
dev.off()



# plot by another way by circle 
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig4_Bcell_metabolism_circle.pdf",useDingbats=F)
plot_Circ_heatmap(t(tmp.res))
dev.off()






























