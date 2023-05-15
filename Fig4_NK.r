

# Fig3 
# NK cell need to analysis
# 2021-3-24
#============================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")
NKcell <- subset(dat,cells=which(dat$celltype=="NK cell"))

# divide into each sample and re-integration
#============================================================================================================

inte.list <- list()
samplelist <- unique(NKcell$sample)
for(i in 1:length(samplelist)){
    tmp <- subset(NKcell,cells=which(NKcell$sample==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
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
inte <- FindClusters(inte,resolution=1)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)

saveRDS(inte,file="/public/workspace/lily/CTN/version_3_20/data/NKcell.RDS")











#############################################################################################################################################
# NK cell purefiy
# to study with pure NK cells
#============================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NKcell.RDS")
data <- dat[['RNA']]@data
sub.dat <- subset(dat,cells=colnames(dat)[-which(data["CD3E",]>0&data["CD3D",]>0)])
#==================================================================================
# new recluster 

inte.list <- list()
samplelist <- unique(sub.dat$sample)
for(i in 1:length(samplelist)){
    tmp <- subset(sub.dat,cells=which(sub.dat$sample==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
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
inte <- FindClusters(inte,resolution=1)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)

saveRDS(inte,file="/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")



# 2023-1-3
# plot umap result 
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
 DefaultAssay(dat) <- "integrated"
 dat <- RunUMAP(dat,dims=1:10)

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/NK_overview_umap.pdf",useDingbats=F)
DimPlot(dat,raster=T,label=T,label.size=8)
dev.off()

# plot CD56 signature score
library(ggplot2)
mod <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NKCD56_mod.RDS")
dat$CD56bri <- mod[,3]
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/NK_CD56score_umap.pdf",useDingbats=F)
FeaturePlot(dat,feature="CD56bri",raster=T)+
    scale_colour_gradientn(colours =  c("#3299cc","#3299cc","white","white","#E41A1C","#E41A1C"),values = c(0,0.5,0.5,0.55,1.0))
dev.off()


#####################################################################################################################################################
# 2021-4-13
# plot NK result 
#====================================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
tmp <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")
tmp.res <- as.numeric(table(dat$sample))/as.numeric(table(tmp$sample))
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/NK_percentage.pdf",useDingbats=F)
barplot(c(median(tmp.res[1:4]),median(tmp.res[12:14]),median(tmp.res[5:11])),names=c("Elder","CTN","SCT"),ylim=c(0,0.3),main="cell percentage of NK")
dev.off()




# define cell type 
# pending ....
#====================================================================================================================================================
# library(Seurat)
# dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
# # 2021-4-24 use zhangzemin ways to identify cluster type 
# gene <- FindAllMarkers(dat,assay="RNA",only.pos=T)
# # add some filter 
# gene <- gene[order(gene$avg_logFC,decreasing=T),]
# gene.f <- gene[which(gene$p_val_adj<0.05&gene$avg_logFC>0),]










# check which cluster are enrich 
library(Seurat)
library(reshape2)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
tmp.res <- apply(apply(table(dat$group,dat$seurat_clusters),1,function(x){x/sum(x)}),1,function(x){x/sum(x)})
colnames(tmp.res) <- paste0("C",colnames(tmp.res))
# tmp.res.f <- tmp.res[,-c(13,14)]
plot.data <- melt(tmp.res.f)
colnames(plot.data) <- c("group","Cluster","percentage")

library(ggplot2)
cols <- c("#5091cd","#7ac143","#f9a541")
plot.data$Cluster <- factor(plot.data$Cluster,levels=c("C0","C1","C5","C8","C2","C3","C4","C7","C6","C9"))
plot.data$group <- factor(plot.data$group,levels=c("Elder","CTN","SC"))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig3_cluster_percent.pdf",width=10,height=8)
ggplot(plot.data,aes(x=Cluster,y=percentage,fill=group,group=group))+ geom_bar(stat="identity",position="dodge")+
    scale_fill_manual(values=cols)+theme_bw()+ labs(y="percentage")


ggplot(plot.data,aes(x=Cluster,y=percentage,fill=group,group=group))+ geom_bar(stat="identity",position="stack")+
    scale_fill_manual(values=cols)+theme_bw()+ labs(y="percentage")

dev.off()






# use add module socre to identify CD56 dim and CD56 bright
# 2021-5-1
#=================================================================================================
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6269138/ 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
DefaultAssay(dat) <- "RNA"
cd56dim <- c("FGFBP2","GZMB","GZMA","SPON2","S100A4","CST7","FCGR3A",
"IGFBP7","GZMH","CFL1","ITGB2","CYBA","LGALS1",
"CD247","LAIR2","PRF1","DSTN","KIR2DL3","PTGDS",
"GTF3C1","EIF3B","KLRD1","TPM3","CTSD","LCP1",
"CD164","RBM38","SAR1A","ENC1","ZEB2","TLE1","EMP3","CD99")

cd56bri <- c("GZMK","CD44","PPP1R14B","CXCR3","RPL36A","SCML1",
"COTL1","NCF1","XCL1","HLA-DRB1","CD83","LTB","ZFP36L2","VIM",
"CD82","SPOCK2","CMC1","CNN2","PDE4B","IFRD1","RGCC","OCIAD2")

# dat <- AddModuleScore(dat,features=cd56dim,name="CD56dim")
# dat <- AddModuleScore(dat,features=cd56bri,name="CD56bri")

# use ssGSEA calculate 

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(cd56bri,"CD56bri",out="/public/workspace/lily/MOD_file/CD56bri.mod") 
mod.generate(cd56dim,"CD56dim",out="/public/workspace/lily/MOD_file/CD56dim.mod") 

mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("CD56bri","CD56dim"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

saveRDS(mod,file="/public/workspace/lily/CTN/version_3_20/data/NKCD56_mod.RDS")





#============================================================================
# 1. NCAM1 expression 
AverageExpression(dat,feature=c("NCAM1"),assay="RNA") -> tmp
dat.p <- data.frame(NCAM1 =as.numeric(as.vector(tmp$RNA)),Cluster=paste0("C",0:13))
dat.p <- dat.p[order(dat.p$NCAM1),]
dat.p$Cluster <- factor(dat.p$Cluster,levels=dat.p[order(dat.p$NCAM1,decreasing=T),"Cluster"])
library(ggplot2)
ggplot(dat.p,aes(x=Cluster,y=NCAM1,fill=NCAM1)) + geom_bar(stat="identity",position="dodge")+theme_bw()


# 2. classify 
dat$celltype.new <- "CD56dim"
dat$celltype.new[which(dat$seurat_clusters%in%c(11))] <- "CD56bri"
saveRDS(dat,file="/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")

DimPlot(dat,group.by="refine.type")

# 3. boxplot 
sub.dat <- subset(dat,cells=which(dat$refine.type%in%c("CD56dim","CD56bri")))
VlnPlot(sub.dat,feature="CD56dim1",group.by="refine.type",pt.size=0)
VlnPlot(sub.dat,feature="CD56bri1",group.by="refine.type",pt.size=0)



# 4. figure 
library(reshape2)
tmp <- melt(apply(apply(table(sub.dat$refine.type,sub.dat$group),2,function(x){x/sum(x)}),1,function(x){x/sum(x)}))
colnames(tmp) <- c("group","Type","value")
cols <- c("#5091cd","#7ac143","#f9a541")
tmp$group <- factor(tmp$group,levels=c("Elder","CTN","SC"))
ggplot(tmp,aes(x=Type,y=value,fill=group))+geom_bar(stat="identity",position="dodge")+theme_bw()+ scale_fill_manual(values=cols)





########################################################################################################################################
# 2021-5-22
# plot result 
# 0. DimPlot 
# 1. Boxplot show signature score
# 2. VlnPlot show FCGR3A and CD56? or show other markers
# 3. percentage difference
#=======================================================================================================================================
# 0. DimPlot
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig3_DimPlot.pdf",useDingbats=F)
DimPlot(dat)
DimPlot(dat,group.by="celltype.new",cols=c("#f9886c","#41afa5"))
dev.off()


# 1. barplot show signature score 
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
mod <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NKCD56_mod.RDS")
mod$Cluster <- paste0("C",dat$seurat_clusters)
tmp <- aggregate(.~Cluster,data=mod[,c(3:4,7)],FUN=median)
dat$CD56bri <- mod[,3]
tmp <- tmp[order(tmp$CD56bri_norm,decreasing=T),]
tmp$Cluster <- factor(tmp$Cluster,levels=tmp$Cluster)

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig3_score.pdf",useDingbats=F)
library(ggplot2)
ggplot(tmp,aes(x=Cluster,y=CD56bri_norm,fill=CD56bri_norm)) + geom_bar(stat="identity",position="dodge") + theme_bw() + 
    scale_fill_gradientn(colors=c("#3299cc","white","#ff4000"))

tmp <- tmp[order(tmp$CD56dim_norm,decreasing=T),]
tmp$Cluster <- factor(tmp$Cluster,levels=tmp$Cluster)
ggplot(tmp,aes(x=Cluster,y=CD56dim_norm,fill=CD56dim_norm)) + geom_bar(stat="identity",position="dodge",width=0.5) + theme_bw() +
    scale_fill_gradientn(colors=c("#3299cc","white","#ff4000"))

FeaturePlot(dat,feature="CD56bri")+
    scale_colour_gradientn(colours =  c("#3299cc","#3299cc","white","white","#E41A1C","#E41A1C"),values = c(0,0.5,0.5,0.55,1.0))
dev.off()




# 2. boxplot show NCAM1 and CD16 expression in dim and bri group
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
dat@active.ident <- factor(dat$celltype.new)
AverageExpression(dat,assays="RNA",features=c("NCAM1","FCGR3A"))

tmp.data <- data.frame(NCAM1=c(0.8227965,0.5743634),FCGR3A=c(2.9567476,8.5898070),group=c("CD56bri","CD56dim"))
library(ggplot2)
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig3_expression.pdf",useDingbats=F)
ggplot(tmp.data,aes(x=group,y=NCAM1)) + geom_bar(stat="identity",position="dodge") + theme_bw()
ggplot(tmp.data,aes(x=group,y=FCGR3A)) + geom_bar(stat="identity",position="dodge") + theme_bw()
dev.off()







# 3. percentage difference in different group 
#====================================================================================================================
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
tmp.res <- apply(apply(table(dat$group,dat$celltype.new),1,function(x){x/sum(x)}),1,function(x){x/sum(x)})

library(reshape2)
library(ggplot2)
res <- melt(tmp.res)
colnames(res) <- c("group","type","percentage")
cols <- c("#5091cd","#7ac143","#f9a541")
res$group <- factor(res$group,levels=c("Elder","CTN","SC"))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig3_percentage.pdf",useDingbats=F)
ggplot(res,aes(x=type,y=percentage,fill=group))+geom_bar(stat="identity",position="dodge")+theme_bw() +
scale_fill_manual(values=cols)
dev.off()





# 4. define CTN SCT enriched cluster 
library(reshape2)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
tmp.res <- apply(apply(table(dat$group,dat$seurat_clusters),1,function(x){x/sum(x)}),1,function(x){x/sum(x)})
colnames(tmp.res) <- paste0("C",colnames(tmp.res))
# tmp.res.f <- tmp.res[,-c(13,14)]
plot.data <- melt(tmp.res)
colnames(plot.data) <- c("group","Cluster","percentage")

library(ggplot2)
cols <- c("#5091cd","#7ac143","#f9a541")
plot.data$Cluster <- factor(plot.data$Cluster,levels=c("C0","C5","C8","C10","C2","C3","C11","C13","C1","C4","C6","C7","C9","C12"))
plot.data$group <- factor(plot.data$group,levels=c("Elder","CTN","SC"))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig3_cluster_percent.pdf",width=10,height=8)
ggplot(plot.data,aes(x=Cluster,y=percentage,fill=group,group=group))+ geom_bar(stat="identity",position="dodge")+
    scale_fill_manual(values=cols)+theme_bw()+ labs(y="percentage")


ggplot(plot.data,aes(x=Cluster,y=percentage,fill=group,group=group))+ geom_bar(stat="identity",position="stack")+
    scale_fill_manual(values=cols)+theme_bw()+ labs(y="percentage")

dev.off()



# 5. Vaccon plot show different 
# vaccon plot result is not verey good
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
gene <- FindMarkers(dat,assay="RNA",group.by="seurat_clusters",ident.1 =c(0,5,8,10),ident.2 =c(1,4,6,7),logfc.threshold=0.1)
# tmp.gene <- gene[which(gene$p_val_adj<0.05),]
# tmp.gene$logP <- tmp.gene$p_val_adj
# tmp.gene$logP[which(tmp.gene$logP==0)] <- min(tmp.gene$logP[-which(tmp.gene$logP==0)])
# tmp.gene$logP <- -log(tmp.gene$logP,100) # use log 100
# tmp.gene$class <- "grey"
# tmp.gene$class[which(tmp.gene$avg_logFC> 0.1)] <- "up"
# tmp.gene$class[which(tmp.gene$avg_logFC< (-0.1))] <- "down"
# tmp.gene$class <- factor(tmp.gene$class,levels=c("up","down","grey"))
# library(ggplot2)
# #pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig4_DEG.pdf",useDingbats=F)
# ggplot(tmp.gene,aes(x=avg_logFC,y=logP,group=class,color=class)) + geom_point(size=1.5) + 
#     scale_colour_manual(values=c("#ec2c22","#2e9df7","#d7d7d8")) +  geom_vline(aes(xintercept= -0.2),color="#2e9df7") +
#      geom_vline(aes(xintercept= 0.2),color="#ec2c22")
# dev.off()

rownames(gene)[which(gene$avg_logFC>0.1&gene$p_val_adj<0.05)] -> tmp
write.table(tmp,file="/public/workspace/lily/CTN/version_3_20/rs_plot/Fig3_DEG.txt",sep="\t",col.names=F,quote=F,row.names=F)


# 2021-9-7 
# get gene and then do enrichment in enrichr 
#=====================================================================================================================================================
dat <- read.table("/public/workspace/lily/CTN/version_3_20/EnrichR/GO_Biological_Process_NK_SCT.CNT_Elder.txt",sep="\t",header=T)

library(ggplot2)
dat$logP <- -log10(dat$Adjusted.P.value)
sub.dat <- head(dat,10)
sub.dat <- sub.dat[order(sub.dat$logP,decreasing=T),]
sub.dat$Term <- factor(sub.dat$Term,levels=rev(sub.dat$Term))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/SC_CTN_NKcell_DEG_GO_enrich.pdf",useDingbats=F,width=10)
ggplot(sub.dat,aes(x=Term,y=Odds.Ratio,fill=logP))+geom_bar(stat="identity",position = position_dodge(0.5),width=0.5) + coord_flip()
dev.off()





# 6. TFs and Metabolism 
# pySCENIC in pipeline 
# metabolism 
#==========================================================================================================================
# TFs 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
tf.dat <- read.table("/public/workspace/lily/CTN/version_3_20/SCENIC/NK/step3.auc_mtx.csv",sep=",",header=T)
rownames(tf.dat) <- tf.dat$Cell
tf.dat$Cell <- NULL
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))
all(rownames(tf.dat)==rownames(dat@meta.data))

# add a column information and plot pheatmap 
#==============================================================================================================================================
# tf.dat$celltype <- dat$celltype.new
tf.dat$cluster <- paste0("C",dat$seurat_clusters)
tf.dat$celltype <- NULL
tmp.res <- aggregate(.~cluster,data=tf.dat,FUN=mean)
rownames(tmp.res) <- tmp.res$cluster
tmp.res$cluster <- NULL
tmp.res.f <- tmp.res[c("C0","C5","C8","C10","C2","C3","C11","C13","C1","C4","C6","C7","C9","C12"),] # C10 is unclassify 
# tmp.res.f <- tmp.res.f[,-which(colSums(tmp.res.f)==0)]  # some TF is all zero in all sample

annotation_row=data.frame(row.names=c("C0","C5","C8","C10","C2","C3","C11","C13","C1","C4","C6","C7","C9","C12"),
    group=c(rep("SC.CTN",4),rep("Uniform",4),rep("Elder",4),rep("CTN",2)))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig3_heatmap_TFs.pdf",width=15,height=8)
pheatmap::pheatmap(tmp.res.f,scale="column",cluster_rows=F,annotation_row=annotation_row,color=colorRampPalette(c('steelblue','white','red'))(50))
dev.off()

p <- pheatmap::pheatmap(tmp.res.f,scale="column",cluster_rows=F,annotation_row=annotation_row,color=colorRampPalette(c('steelblue','white','red'))(50))
# just use SC.CTN up TFs 
tf.dat <- data.frame(TF=colnames(tmp.res.f)[p$tree_col$order][89:116],group=c(rep("SC.CTN",28)))

#write.table(tf.dat$TF,file="./tmp.txt",row.names=F,col.names=F,quote=F)
trrust.dat <- read.table('/public/workspace/lily/metastasis/data/TRRUST/trrust_rawdata.human.tsv',sep="\t") 
colnames(trrust.dat) <- c("source","traget","effect","id")

tmp.tf <- merge(trrust.dat,tf.dat,by.x="source",by.y="TF")

#===============================================================================================================
tmp.up <- tmp.tf[which(tmp.tf$effect=="Activation"),]
write.table(as.vector(tmp.up$traget),file="/public/workspace/lily/CTN/version_3_20/rs_plot/NK_TF.up.gene.txt",sep="\t",row.names=F,col.names=F,quote=F)

##########################################################################
# metabolism 
tmp <- read.table("/public/workspace/lily/CTN/version_3_20/metabolism/NK/KEGGpathway_activity_shuffle.txt",sep="\t",header=T)
tmp.a <- tmp
tmp.a[is.na(tmp.a)] <- 1


tmp.f <- tmp.a[,c("C0","C5","C8","C10","C2","C3","C11","C13","C1","C4","C6","C7","C9","C12")]
tmp.res <- t(apply(tmp.f,1,function(x){
    c(mean(x[1:4]),mean(x[c(5:8)]),mean(x[c(9:12)]),mean(x[c(13,14)]))
}))
colnames(tmp.res) <- c("SC.CTN","Uniform","Elder","CTN")


# 2021-9-7 
# another way to plot result 


pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig2_NKcell_metabolism_circle.pdf",useDingbats=F)
plot_Circ_heatmap(t(tmp.res))
dev.off()




pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig3_NK_metabolism.pdf",height=12,width=8)
pheatmap::pheatmap(tmp.res,scale="row",color=colorRampPalette(c('steelblue','white','red'))(50))
dev.off()













