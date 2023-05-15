
#!/usr/bin/Rscript

# 2021-4-23 
# analysis all monocyte cell in PBMC 
# decide to use CD14 mono or all mono 
# check data and find CD16 is always low than 5% ,not good . so analysis just use CD14 monocyte 
# CD14 Mono  CD16 Mono  
# 0.1359392 0.03029750 
# 0.1861506 0.02686222 
# 0.1735324 0.02453782 
#############################################################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")
sub.dat <- subset(dat,cells=which(dat$celltype=="CD14 Mono"))
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

saveRDS(inte,file="/public/workspace/lily/CTN/version_3_20/data/CD14Mono.RDS")




#==============================================================================================================================================
# 2021-5-21
# it semms no subtype for CD14 Monocyte 
# https://www.frontiersin.org/articles/10.3389/fimmu.2019.02035/full said major subtype 

# 1. define cluster 
# seem not every interseting
library(reshape2)
library(ggplot2)

dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CD14Mono.RDS")
tmp.dat <- apply(apply(table(dat$seurat_clusters,dat$group),2,function(x){x/sum(x)}),1,function(x){x/sum(x)})
tmp.res <- melt(tmp.dat)
colnames(tmp.res) <- c("group","Cluster","percentage")
tmp.res$Cluster <- paste0("C",tmp.res$Cluster)
cols <- c("#5091cd","#7ac143","#f9a541")
tmp.res$group <- factor(tmp.res$group,levels=c("Elder","CTN","SC"))
tmp.res$Cluster <- factor(tmp.res$Cluster,levels=c("C12","C14","C7","C9","C0","C1","C2","C3","C4","C5","C11","C6","C8","C13","C10","C15"))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig5_celltype_percent.pdf",useDingbats=F,height=10)
ggplot(tmp.res,aes(x=Cluster,y=percentage,fill=group,group=group))+ geom_bar(stat="identity",position="stack")+
    scale_fill_manual(values=cols)+theme_bw()+ labs(y="percentage")
dev.off()



pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig5_TSNE_plot.pdf",useDingbats=F)
DimPlot(dat)
dev.off()

# 2023-1-3
# set umap 
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CD14Mono.RDS")
dat <- RunUMAP(dat,dims=1:10)

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig6_Mono_umap_plot.pdf",useDingbats=F)
DimPlot(dat,label=T,label.size=8,raster=T)
dev.off()



# 2021-11-1
# set new cell group 
#=====================================================================================================================================================
library(reshape2)
library(ggplot2)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CD14Mono.RDS")
dat$type.new <- "Unknow"
dat$type.new[which(dat$seurat_clusters%in%c(0,1,2,3,4,5,11))] <- "Uniform"
dat$type.new[which(dat$seurat_clusters%in%c(6,8,13))] <- "CTN.SC"
dat$type.new[which(dat$seurat_clusters%in%c(12,14))] <- "Elder"
dat$type.new[which(dat$seurat_clusters%in%c(10,15))] <- "CTN"
dat$type.new[which(dat$seurat_clusters%in%c(7,9))] <- "Elder.SC"

saveRDS(dat,file="/public/workspace/lily/CTN/version_3_20/data/CD14Mono.RDS")


#================================================================================================================================================
# now find DEGS 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CD14Mono.RDS")
DefaultAssay(dat) <- "RNA"
tmp.gene <- FindMarkers(dat,ident.1="CTN.SC",ident.2="Elder",assay="RNA",group.by="type.new")
tmp.gene <- tmp.gene[which(tmp.gene$p_val_adj<0.05),]

#===============================================================================================================================================
tmp.gene$logP <- tmp.gene$p_val_adj
tmp.gene$logP[which(tmp.gene$logP==0)] <- min(tmp.gene$logP[-which(tmp.gene$logP==0)])
tmp.gene$logP <- -log(tmp.gene$logP,100) # use log 100
tmp.gene$class <- "grey"
tmp.gene$class[which(tmp.gene$avg_logFC> 0.5)] <- "up"
tmp.gene$class[which(tmp.gene$avg_logFC< (-0.5))] <- "down"
tmp.gene$class <- factor(tmp.gene$class,levels=c("up","down","grey"))

# add some gene names 
tmp.gene$geneName <- "NA"
tmp.gene$geneName[which(tmp.gene$avg_logFC>0.7)] <- rownames(tmp.gene[which(tmp.gene$avg_logFC>0.7),]) # top5
tmp.gene$geneName[which(tmp.gene$avg_logFC< (-1.9))] <- rownames(tmp.gene[which(tmp.gene$avg_logFC<(-1.9)),]) # bottom5
tmp.gene$geneName[which(tmp.gene$geneName=="NA")] <- NA


library(ggplot2)
library(ggrepel)

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig5_DEG_vaconao.pdf",useDingbats=F)
ggplot(tmp.gene,aes(x=avg_logFC,y=logP,group=class,color=class)) + geom_point(size=3) + 
    scale_colour_manual(values=c("#ec2c22","#2e9df7","#d7d7d8")) +  geom_vline(aes(xintercept= -0.5),color="#2e9df7") +
     geom_vline(aes(xintercept= 0.5),color="#ec2c22") + xlim(c(-3,3)) + geom_text(aes(label = geneName), size = 3)+ labs(title="CTN.SC vs Elder")

dev.off()




# get gene name to do enrichment analysis 
gene.up <- rownames(tmp.gene[which(tmp.gene$avg_logFC>0),])
write.table(gene.up,file="~/tmp/gene.up.txt",quote=F,row.names=F,col.names=F,sep="\t")
gene.dn <- rownames(tmp.gene[which(tmp.gene$avg_logFC<0),])
write.table(gene.dn,file="~/tmp/gene.dn.txt",quote=F,row.names=F,col.names=F,sep="\t")


# then do enrichment analsis in EnrichR 
# now plot enrichment analysis result
#====================================================================================================================================================
dat <- read.table("/public/workspace/lily/CTN/version_3_20/EnrichR/GO_Biological_Process_2021_CD14mono_SC_CTN_vs_Elder_up.txt",sep="\t",header=T)

library(ggplot2)
dat$logP <- -log10(dat$Adjusted.P.value)
sub.dat <- head(dat,10)
sub.dat <- sub.dat[order(sub.dat$logP,decreasing=T),]
sub.dat$Term <- factor(sub.dat$Term,levels=rev(sub.dat$Term))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/SC_CTN_CD14mono_up_DEG_GO_enrich.pdf",useDingbats=F,width=10)
ggplot(sub.dat,aes(x=Term,y=Odds.Ratio,fill=logP))+
geom_bar(stat="identity",position = position_dodge(0.5),width=0.5) + coord_flip() 
# + scale_fill_gradientn(colours=c("grey","red"))
dev.off()




#====================================================================================================================================================
dat <- read.table("/public/workspace/lily/CTN/version_3_20/EnrichR/GO_Biological_Process_2021_CD14mono_SC_CTN_vs_Elder_dn.txt",sep="\t",header=T)

library(ggplot2)
dat$logP <- -log10(dat$Adjusted.P.value)
sub.dat <- head(dat,10)
sub.dat <- sub.dat[order(sub.dat$logP,decreasing=T),]
sub.dat$Term <- factor(sub.dat$Term,levels=rev(sub.dat$Term))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/SC_CTN_CD14mono_dn_DEG_GO_enrich.pdf",useDingbats=F,width=10)
ggplot(sub.dat,aes(x=Term,y=Odds.Ratio,fill=logP))+geom_bar(stat="identity",position = position_dodge(0.5),width=0.5) + coord_flip()
dev.off()







#=============================================================================================================================================
# 2021-11-2
# analysis about metabolism 
#=============================================================================================================================================


tmp <- read.table("/public/workspace/lily/CTN/version_3_20/metabolism/CD14mono/KEGGpathway_activity_shuffle.txt",sep="\t",header=T)
tmp.a <- tmp
tmp.a[is.na(tmp.a)] <- 1

tmp.f <- tmp.a[,c("C12","C14","C7","C9","C0","C1","C2","C3","C4","C5","C11","C6","C8","C13","C10","C15")]
tmp.res <- t(apply(tmp.f,1,function(x){
    c(mean(x[1:2]),mean(x[3:4]),mean(x[c(5:11)]),mean(x[c(12:14)]),mean(x[c(15:16)]))  
}))
colnames(tmp.res) <- c("Elder","Elder.SC","Uniform","CTN.SC","CTN")


pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig5_CD14mono_metabolism.pdf",height=12,width=8)
pheatmap::pheatmap(tmp.res,scale="row",color=colorRampPalette(c('steelblue','white','red'))(50))
dev.off()


# plot by another way by circle 
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig5_CD14mono_metabolism_circle.pdf",useDingbats=F)
plot_Circ_heatmap(t(tmp.res))
dev.off()














#=============================================================================================================================================
# 2021-11-2
# analysis about TFs 
#=============================================================================================================================================
library(Seurat)
library(reshape2)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CD14Mono.RDS")
sub.dat <- dat
#===================================================================================================================================================
# regultons 
reg <- read.delim2("/public/workspace/lily/CTN/version_3_20/SCENIC/CD14mono/step3.auc_mtx.csv",sep=",",header=T)
rownames(reg) <- reg$Cell
reg$Cell <- NULL
gsub("\\.\\.\\.","",colnames(reg)) -> colnames(reg)
reg.m <- apply(reg,2,function(x){as.numeric(as.vector(x))})
reg.m <- data.frame(reg.m)
rownames(reg.m) <- rownames(reg)
reg.f <- reg.m[which(rownames(reg.m)%in%colnames(sub.dat)),] # subdat cells

reg.f$Cluster <- paste0("C",sub.dat$seurat_clusters)
tmp.res <- aggregate(.~Cluster,data=reg.f,FUN=mean)
rownames(tmp.res) <- tmp.res$Cluster
tmp.res$Cluster <- NULL

tmp.res.f <- tmp.res[c("C12","C14","C7","C9","C0","C1","C2","C3","C4","C5","C11","C6","C8","C13","C10","C15"),]

annotation_row=data.frame(row.names=c("C12","C14","C7","C9","C0","C1","C2","C3","C4","C5","C11","C6","C8","C13","C10","C15"),
    group=c(rep("Elder",2),rep("Elder.SC",2),rep("Uniform",7),rep("CTN.SC",3),rep("CTN",2))
)



pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig5_heatmap_TFs.pdf",width=30,height=8,useDingbats=F)
pheatmap::pheatmap(tmp.res.f,scale="column",cluster_rows=F,
    annotation_row=annotation_row,color=colorRampPalette(c('steelblue','white','red'))(50),cellwidth=10,cellheight=10)
dev.off()


# circle
# library(circlize)
# col_fun <- colorRamp2(c(-4, 0, 4), c("steelblue", "white", "red"))
# circos.heatmap(t(tmp.res.f), col = col_fun,track.height = 0.5,dend.side = "inside",rownames.side = "outside")
# circos.clear()


pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig5_heatmap_TFs_circlize.pdf",width=30,height=8,useDingbats=F)
plot_Circ_heatmap(tmp.res.f)
dev.off()












































