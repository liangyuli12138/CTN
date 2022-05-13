

# Fig2 
# 2021-3-24
# plot T cell analysis 
#============================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")
Tcell <- subset(dat,cells=which(dat$celltype=="T cell"))

# divide into each sample and re-integration
#============================================================================================================

inte.list <- list()
samplelist <- unique(Tcell$sample)
for(i in 1:length(samplelist)){
    tmp <- subset(Tcell,cells=which(Tcell$sample==samplelist[i]))
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

saveRDS(inte,file="/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")


#=================================================================================================================================
# use Single R function to dentify samples 
# 2021-3-24
##################################################################################################################################


library(SingleR)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
# ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_HumanPrimaryCellAtlasData.RDS")
# ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_NovershternHematopoieticData.RDS")

ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_BluePrint.RDS")
res.cluster<-SingleR(test=as.matrix(dat@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=dat$seurat_clusters,method="cluster")


# ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_MonacoImmuneData.RDS")
# res.cluster1<-SingleR(test=as.matrix(dat@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=dat$seurat_clusters,method="cluster")


# ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_HumanPrimaryCellAtlasData.RDS")
# res.cluster2<-SingleR(test=as.matrix(dat@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=dat$seurat_clusters,method="cluster")

#####################################################################################################################################
# 2021-3-26
# use SingleR and marker gene to make sure 
# not every difference with SingleR result 
# use Blue print as the reference 
#####################################################################################################################################
# make sure
dat$celltype.refine <- "unclassify"
dat$celltype.refine[which(dat$seurat_clusters%in%c(0,1,2,4,8,10,11,13))] <- "CD8+ Tem"
dat$celltype.refine[which(dat$seurat_clusters%in%c(3,6))] <- "Naive T"
dat$celltype.refine[which(dat$seurat_clusters%in%c(7,12))] <- "CD8+ Tcm"
dat$celltype.refine[which(dat$seurat_clusters%in%c(5))] <- "CD4+ Tem"
dat$celltype.refine[which(dat$seurat_clusters%in%c(9))] <- "Tregs"
dat$celltype.refine[which(dat$seurat_clusters%in%c(14))] <- "Unclassify"

saveRDS(dat,file="/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
# # try to use Find marker to find each cluster marker gene 
# gene.tcell <- FindAllMarkers(dat,assay="RNA")

# tmp <- table(dat$sample,dat$seurat_clusters)
# tmp.res <- apply(tmp,1,function(x){x/sum(x)})


#                  0          1         2         3         4         5
#   CTN   0.54588542 0.48265664 0.5198954 0.1651989 0.3039248 0.2230200
#   Elder 0.01611098 0.08100185 0.1102143 0.5916688 0.2575138 0.5185489
#   SC    0.43800360 0.43634151 0.3698902 0.2431324 0.4385614 0.2584312

#                  6         7         8         9        10        11        12
#   CTN   0.06631604 0.2177113 0.4133861 0.3293772 0.2940330 0.2394032 0.6073908
#   Elder 0.76812796 0.4138678 0.2058944 0.4276933 0.3363261 0.4030696 0.2740023
#   SC    0.16555601 0.3684209 0.3807195 0.2429294 0.3696408 0.3575272 0.1186069

#                13          14
#   CTN   0.3695887 0.974114851
#   Elder 0.3075970 0.021383770
#   SC    0.3228143 0.004501379














##################################################################################################################################################
# T cell TSNE Plot result 
# 2021-3-26
#=================================================================================================================================================
# CD4+ Tem   CD8+ Tcm   CD8+ Tem    Naive T      Tregs Unclassify
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
dat <- RunUMAP(dat,dims=1:10)
cols.t <- c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3")
DimPlot(dat,group.by="celltype.refine",cols=cols.t)




#=================================================================================================================================================
# plot Figure 2 
# result plot 
##################################################################################################################################################

# 1. overview plot 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
dat <- RunUMAP(dat,dims=1:10)
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig2_overview.pdf",width=10,height=10,useDingbats=F)
cols.t <- c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3")
DimPlot(dat,group.by="celltype.refine",cols=cols.t)
dev.off()


# 2. stack bar plot 
library(reshape2)
tmp.dat <- apply(table(dat$celltype.refine,dat$group),2,function(x){x/sum(x)})
tmp.res <- melt(tmp.dat)
colnames(tmp.res) <- c("celltype.refine","group","percentage")

# plot result 
library(ggalluvial)
library(ggplot2)
library(RColorBrewer)

# plot river plot 
cols.t <- c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3")
tmp.res$group <- factor(tmp.res$group,levels=c("Elder","CTN","SC"))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig2_celltype_percent.pdf",useDingbats=F,height=10)
ggplot(tmp.res, aes(x = group, y = percentage, fill = celltype.refine,stratum = celltype.refine,alluvium = celltype.refine)) +
geom_stratum() +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols.t)+
theme_bw() + theme(panel.grid.major = element_blank())

dev.off()





# 2. add some bar plot to show 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")

#========================================================================================================
library(reshape2)
tmp.res <- apply(apply(table(dat$group,dat$seurat_clusters),1,function(x){x/sum(x)}),1,function(x){x/sum(x)})
colnames(tmp.res) <- paste0("C",colnames(tmp.res))
tmp.res.f <- tmp.res[,-15]
plot.data <- melt(tmp.res.f)
colnames(plot.data) <- c("group","Cluster","percentage")

library(ggplot2)
cols <- c("#7ac143","#5091cd","#f9a541")
plot.data$Cluster <- factor(plot.data$Cluster,levels=c("C0","C1","C2","C3","C5","C6","C7","C8","C9","C4","C10","C11","C13","C12"))


pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig2_cluster_percent.pdf",width=10,height=8)
ggplot(plot.data,aes(x=Cluster,y=percentage,fill=group,group=group))+ geom_bar(stat="identity",position="dodge")+
    scale_fill_manual(values=cols)+theme_bw()+ labs(y="percentage")

ggplot(plot.data,aes(x=Cluster,y=percentage,fill=group,group=group))+ geom_bar(stat="identity",position="stack")+
    scale_fill_manual(values=cols)+theme_bw()+ labs(y="percentage")

dev.off()








# 3. trajectory 
##########################################################################################################
library(Seurat)
library(monocle)


#=========================================================================================================
# 2021-3-27
# use integration to find order gene 
# 2021-9-1 
# delet Treg 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
sub.dat <- subset(dat,cells=which(dat$celltype.refine%in%c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T")))

inte.list <- list() 
samplelist <- unique(sub.dat$sample)
for(i in 1:length(samplelist)){
	tmp <- subset(sub.dat,cells=which(sub.dat$sample==samplelist[i]))
	DefaultAssay(tmp) <- "RNA"
	inte.list[i] <- tmp
}

# try to another way to integration 
#===========================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
sub.dat <- subset(dat,cells=which(dat$celltype.refine!="Unclassify"))

inte.list <- list() 
samplelist <- unique(sub.dat$celltype.refine)
for(i in 1:length(samplelist)){
	tmp <- subset(sub.dat,cells=which(sub.dat$celltype.refine==samplelist[i]))
	DefaultAssay(tmp) <- "RNA"
	inte.list[i] <- tmp
}









integration.anchors <- FindIntegrationAnchors(object.list = inte.list)
inte <- IntegrateData(anchorset = integration.anchors)
#FindVariableFeatures
inte <- FindVariableFeatures(inte)
############################################################################################################
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

genes <- VariableFeatures(inte)[1:500]
ordering_genes <- genes
tmp.dat <- setOrderingFilter(tmp.dat, ordering_genes = ordering_genes)
tmp.dat <- reduceDimension(tmp.dat, method = 'DDRTree')
tmp.dat <- orderCells(tmp.dat)

saveRDS(tmp.dat,file="/public/workspace/lily/CTN/version_3_20/trajectory/Tcell_pure_trajectory.RDS")


# check result 
library(monocle)
library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/CTN/version_3_20/trajectory/Tcell_pure_trajectory.RDS")

plot_cell_trajectory(tmp.dat,color_by="celltype.refine")+scale_colour_manual(values=c("#ffc845","#00c16e","#f85a40","#7d3f98","red"))
plot_cell_trajectory(tmp.dat,color_by="group")
plot_cell_trajectory(tmp.dat,color_by="sample") +facet_wrap(~sample, nrow = 3)








#======================================================================================================================================================
# no re - inte 

dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
sub.dat <- subset(dat,cells=colnames(dat)[which(dat$celltype.refine%in%c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"))])

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

DefaultAssay(dat) <- "integrated"
genes <- VariableFeatures(dat)[1:100]
ordering_genes <- genes
tmp.dat <- setOrderingFilter(tmp.dat, ordering_genes = ordering_genes)
tmp.dat <- reduceDimension(tmp.dat, method = 'DDRTree')
tmp.dat <- orderCells(tmp.dat)
saveRDS(tmp.dat,file="/public/workspace/lily/CTN/version_3_20/trajectory/Tcell_pure_trajectory.RDS")

pdf("/public/workspace/lily/CTN/version_3_20/monocle/Tcell_pure_trajectory.pdf",useDingbats=F)

plot_cell_trajectory(tmp.dat,color_by="celltype.refine")+scale_colour_manual(values=c("#ffc845","#00c16e","#f85a40","#7d3f98","red"))
plot_cell_trajectory(tmp.dat,color_by="group")
plot_cell_trajectory(tmp.dat,color_by="sample") +facet_wrap(~sample, nrow = 3)
dev.off()






# 4. show different Genes in CD8+ Tem
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
sub.dat <- subset(dat,cells=which(dat$celltype.refine=="CD8+ Tem"))
DefaultAssay(sub.dat) <- "RNA"
sub.dat@active.ident <- factor(sub.dat$group)

# should use CTN and SC to find different genes with 
# gene <- 

CTN.gene <- FindMarkers(sub.dat,ident.1="CTN",ident.2="Elder",assay="RNA",group.by="group")
SC.gene <- FindMarkers(sub.dat,ident.1="SC",ident.2="Elder",assay="RNA",group.by="group")


#===================================================================================================================================================
# plot result 
# use avg log FC and pct.1-pct.2 to calculate result
####################################################################################################################################################
CTN.gene.f <- CTN.gene[which(CTN.gene$p_val_adj<0.05),]
SC.gene.f <- SC.gene[which(SC.gene$p_val_adj<0.05),]

# calculate pct difference 
CTN.gene.f$pct.diff <- CTN.gene.f$pct.1 - CTN.gene.f$pct.2
SC.gene.f$pct.diff <- SC.gene.f$pct.1 - SC.gene.f$pct.2

# add grou information 
CTN.gene.f$group <- "CTN"
SC.gene.f$group <- "SC"

# add gene names 
CTN.gene.f$gene <- rownames(CTN.gene.f)
SC.gene.f$gene <- rownames(SC.gene.f)

#===================================================================================================================================================
tmp.1 <- CTN.gene.f[,c("avg_logFC","pct.diff","group","gene","p_val_adj")]
tmp.2 <- SC.gene.f[,c("avg_logFC","pct.diff","group","gene","p_val_adj")]

tmp.1$logP <- tmp.1$p_val_adj
tmp.1$logP[which(tmp.1$logP==0)] <- min(tmp.1$logP[-which(tmp.1$logP==0)])
tmp.1$logP <- -log(tmp.1$logP,100) # use log 100

tmp.2$logP <- tmp.2$p_val_adj
tmp.2$logP[which(tmp.2$logP==0)] <- min(tmp.2$logP[-which(tmp.2$logP==0)])
tmp.2$logP <- -log(tmp.2$logP,100) # use log 100

tmp.res <- rbind(tmp.1,tmp.2)
rownames(tmp.res) <- NULL
tmp.res$type <- "down"
tmp.res$type[which(tmp.res$avg_logFC>0)] <- "up"


tmp.res$class <- paste0(tmp.res$group,"_",tmp.res$type)
tmp.res$class[which(tmp.res$gene%in%names(which(table(tmp.res$gene)==2))&tmp.res$type=="up")] <- "co-up"
tmp.res$class[which(tmp.res$gene%in%names(which(table(tmp.res$gene)==2))&tmp.res$type=="dn")] <- "co-dn"


write.table(tmp.res,file="/public/workspace/lily/CTN/version_3_20/rs_plot/SC_CTN_DEG.txt",sep="\t",col.names=T,quote=F,row.names=F)

#===================================================================================================================================================
# not suit to plot vacano plot

library(ggplot2)
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig2_vacano_plot.pdf",useDingbats=F)
ggplot(tmp.res,aes(x=avg_logFC,y=logP,group=class,color=class)) + geom_point() + 
    scale_colour_manual(values=c("#ec2c22","#2e9df7","#edd812","#109dc0","#f38020"))
dev.off()
##############################################################################################
# use co expression gene to do KEGG analysis 
write.table(unique(tmp.res$gene[which(tmp.res$group=="Co_up")]),file="./tmp.txt",row.names=F,col.names=F,quote=F)

# Panther 2016 GO 
# have some interesting pathway 








# 5. TFs 
#===================================================================================================================================================
library(Seurat)
library(reshape2)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
sub.dat <- subset(dat,cells=which(dat$celltype.refine!="Unclassify"))
#===================================================================================================================================================
# regultons 
reg <- read.delim2("/public/workspace/lily/CTN/version_3_20/SCENIC/Tcell/all_Tcell_auc_mtx.tsv",sep="\t",header=T)
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

tmp.res.f <- tmp.res[c("C0","C1","C2","C3","C5","C6","C7","C8","C9","C4","C10","C11","C13","C12"),]

annotation_row=data.frame(row.names=c("C0","C1","C2","C3","C5","C6","C7","C8","C9","C4","C10","C11","C13","C12"),
    group=c(rep("SC.CTN",3),rep("Elder",3),rep("Uniform",7),rep("CTN",1)),
    type=c(rep("CD8+ Tem",3),"naive T","CD4+ Tem","naive T","CD8+ Tcm","CD8+ Tem","Tregs","CD8+ Tem","CD8+ Tem","CD8+ Tem","CD8+ Tem","CD8+ Tcm"))



pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig2_heatmap_TFs.pdf",width=15,height=8,useDingbats=F)
pheatmap::pheatmap(tmp.res.f,scale="column",cluster_rows=F,annotation_row=annotation_row,color=colorRampPalette(c('steelblue','white','red'))(50))
dev.off()


# use TFs to plot result 
#===================================================================================================================================================
p <- pheatmap::pheatmap(tmp.res.f,scale="column",cluster_rows=F,annotation_row=annotation_row,color=colorRampPalette(c('steelblue','white','red'))(50))

tf.dat <- data.frame(TF=colnames(tmp.res.f)[p$tree_col$order],group=c(rep("Elder",13),rep("CTN.SC",33),rep("ungroup",6),rep("CTN",26)))

# TRUUST database 
trrust.dat <- read.table('/public/workspace/lily/metastasis/data/TRRUST/trrust_rawdata.human.tsv',sep="\t") 
colnames(trrust.dat) <- c("source","traget","effect","id")

tmp.tf <- merge(trrust.dat,tf.dat,by.x="source",by.y="TF")

#===============================================================================================================
tmp.up <- tmp.tf[which(tmp.tf$effect=="Activation"),]

tmp.dn <- tmp.tf[which(tmp.tf$effect=="Repression"),]

write.table(tmp.up,file="/public/workspace/lily/CTN/version_3_20/rs_plot/TF.up.txt",sep="\t",row.names=F,col.names=T,quote=F)
write.table(tmp.dn,file="/public/workspace/lily/CTN/version_3_20/rs_plot/TF.dn.txt",sep="\t",row.names=F,col.names=T,quote=F)
# plot result 
#################################################################################################################
tmp.tf.f <- tmp.up[which(tmp.up$group=="CTN.SC"),]
links.up <- tmp.tf.f[,c("source","traget","effect")]
links.up$effect <- 1
node.up <- unique(c(as.vector(tmp.tf.f$source),as.vector(tmp.tf.f$traget)))
node.up <- as.data.frame(node.up)






#=====================================================================================================================================================
# 2021-9-2
# analysis TFs target genes 
tmp.tf.f$traget <- as.vector(tmp.tf.f$traget)
sort(table(tmp.tf.f$traget)) # found IFNG gene 

tmp <- head(sort(table(tmp.tf.f$traget),decreasing=T),17)

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/TFs_target_gene.pdf",useDingbats=F)
barplot(tmp,las=2)
dev.off()




#====================================================================================================================================================
# 2021-9-2
gene <- unique(as.vector(tmp.tf.f$traget))
write.table(gene,file="/public/workspace/lily/CTN/version_3_20/rs_plot/SC_CTN_TF_act_target_gene.txt",row.names=F,col.names=F,quote=F)

# run in enrichR and plot 
dat <- read.table("/public/workspace/lily/CTN/version_3_20/EnrichR/GO_Biological_Process_Tcell_SC_CTN_TF_act_target_gene.txt",header=T,sep="\t")

library(ggplot2)
dat$logP <- -log10(dat$Adjusted.P.value)
sub.dat <- head(dat,10)
sub.dat <- sub.dat[order(sub.dat$logP,decreasing=T),]
sub.dat$Term <- factor(sub.dat$Term,levels=rev(sub.dat$Term))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/SC_CTN_TF_act_target_gene_GO_enrich.pdf",useDingbats=F)
ggplot(sub.dat,aes(x=Term,y=Odds.Ratio,fill=logP))+geom_bar(stat="identity",position = position_dodge(0.5),width=0.5) + coord_flip()
dev.off()








# up 
library(ggalluvial)
library(ggplot2)
library(RColorBrewer)

links <- links.up
node <- node.up
tmp <- to_lodes_form(links,axes = 1:2,id = "Cohort") 
tmp$stratum <- factor(tmp$stratum,levels=unique(tmp$stratum))
# plot 
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175')
cols <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75')
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig2_SC.CTN.up_TF.pdf",useDingbats=F,height=10)
ggplot(tmp,aes(x =factor(x,level = c("source","traget")),y=effect,stratum = stratum, alluvium = Cohort,fill = stratum, label =stratum)) +
geom_flow( width = 1/3) + # flow
geom_stratum( width = 1/3,linetype=0,size=0.5,alpha =0.5,color = "black")+ 
geom_text(stat ="stratum" , size =3) + #添加名字
scale_x_discrete(limits = c() )+ #去掉横坐标轴
theme_bw() +
theme(legend.position="none") +
scale_fill_manual(values = c(my36colors[1:19],rep("grey",179)))
dev.off()





#================================================================================================================================================
# 2021-6-21
# analysis about T cell TFs 
#================================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
sub.dat <- subset(dat,cells=which(dat$celltype.refine!="Unclassify"))
#===================================================================================================================================================
# regultons 
reg <- read.delim2("/public/workspace/lily/CTN/version_3_20/SCENIC/Tcell/all_Tcell_auc_mtx.tsv",sep="\t",header=T)
rownames(reg) <- reg$Cell
reg$Cell <- NULL
gsub("\\.\\.\\.","",colnames(reg)) -> colnames(reg)
reg.m <- apply(reg,2,function(x){as.numeric(as.vector(x))})
reg.m <- data.frame(reg.m)
rownames(reg.m) <- rownames(reg)
reg.f <- reg.m[which(rownames(reg.m)%in%colnames(sub.dat)),] # subdat cells

# add group info 
sub.dat$group.fine <- "Unknow"
sub.dat$group.fine[which(sub.dat$seurat_clusters%in%c(3,5,6))] <- "Elder"
sub.dat$group.fine[which(sub.dat$seurat_clusters%in%c(0,1,2))] <- "SCT.CTN"


reg.f$Cluster <- sub.dat$group.fine
tmp.res <- aggregate(.~Cluster,data=reg.f,FUN=mean)
rownames(tmp.res) <- tmp.res$Cluster
tmp.res$Cluster <- NULL 
tmp.res <- t(tmp.res)
tmp.res <- data.frame(tmp.res)
tmp.res$df <- tmp.res$SCT.CTN - tmp.res$Elder
tmp.res <- tmp.res[order(tmp.res$df,decreasing=T),]
####################################################################################################################################################
tmp.res$df <- tmp.res$SCT.CTN -tmp.res$Elder




# 5. metabolism 
#===================================================================================================================================================

tmp <- read.table("/public/workspace/lily/CTN/version_3_20/metabolism/Tcell/KEGGpathway_activity_shuffle.txt",sep="\t",header=T)
tmp.a <- tmp
tmp.a[is.na(tmp.a)] <- 1

tmp.f <- tmp.a[,-grep("C14",colnames(tmp.a))]
tmp.f <- tmp.f[,c("C0","C1","C2","C4","C8","C3","C5","C6","C7","C9","C10","C11","C13","C12")]
tmp.res <- t(apply(tmp.f,1,function(x){
    c(mean(x[1:5]),mean(x[c(6:10)]),mean(x[c(11:13)]),mean(x[c(14)]))  
}))
colnames(tmp.res) <- c("SC.CTN","Elder","Uniform","CTN")
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig2_Tcell_metabolism.pdf",height=12,width=8)
pheatmap::pheatmap(tmp.res,scale="row",color=colorRampPalette(c('steelblue','white','red'))(50))
dev.off()


# plot by another way by circle 
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig2_Tcell_metabolism_circle.pdf",useDingbats=F)
plot_Circ_heatmap(t(tmp.res))
dev.off()














# 2021-9-3 
# just plot TCA and OXPH by boxplot 
tmp <- read.table("/public/workspace/lily/CTN/version_3_20/metabolism/Tcell/KEGGpathway_activity_shuffle.txt",sep="\t",header=T)
tmp.a <- tmp
tmp.a[is.na(tmp.a)] <- 1

tmp.f <- tmp.a[,-grep("C14",colnames(tmp.a))]
tmp.f <- tmp.f[,c("C0","C1","C2","C4","C8","C3","C5","C6","C7","C9","C10","C11","C13","C12")]

boxplot(as.numeric(tmp.f[2,1:5]),as.numeric(tmp.f[2,c(6:10)]),as.numeric(tmp.f[2,c(11:13)]),as.numeric(tmp.f[2,c(14)]))

boxplot(as.numeric(tmp.f[15,1:5]),as.numeric(tmp.f[15,c(6:10)]),as.numeric(tmp.f[15,c(11:13)]),as.numeric(tmp.f[15,c(14)]))


























#===================================================================================================================================================
# add some result 
# 2021-9-1
# 1. add plot for 
#===================================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
sub.dat <- subset(dat,cells=which(dat$celltype.refine!="Unclassify"))
DefaultAssay(sub.dat) <- "RNA"
sub.dat@active.ident <- factor(sub.dat$group)

# calculate average expression 
tmp <- AverageExpression(sub.dat,assays="RNA")$RNA
tmp.f <- tmp[which(tmp$CTN>0&tmp$Elder>0&tmp$SC>0),]
tmp.f$log.fc.CTN <- log10(tmp.f$CTN/tmp.f$Elder)
tmp.f$log.fc.SCT <- log10(tmp.f$SC/tmp.f$Elder)

plot(tmp.f$log.fc.CTN,tmp.f$log.fc.SCT)

gene <- rownames(tmp.f[which(tmp.f$log.fc.CTN>1&tmp.f$log.fc.SCT>1),])
write.table(gene,file="~/tmp/tmp.gene.txt",row.names=F,col.names=F,quote=F)

rownames(tmp.f[which(tmp.f$log.fc.CTN<(-1)&tmp.f$log.fc.SCT<(-1)),])




# use ssGSEA to calculate score 
# sub.dat do not have Unclassify t cell 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.rs <- mod.analyze2(as.matrix(sub.dat[['RNA']]@data),c("Tcell_act"),"/public/workspace/lily/MOD_file/",permN=0)
mod.rs <- as.data.frame(mod.rs)
saveRDS(mod.rs,file="/public/workspace/lily/CTN/version_3_20/data/Tcell_tcell_act.RDS")

# plot result 
mod.rs$group <- sub.dat$group
mod.rs$group <- factor(mod.rs$group,levels=c("Elder","SC","CTN"))
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig2_Tcell_act.pdf",useDingbats=F)
boxplot(Tcell_act_norm~group,data=mod.rs,mian="MSGDB T cell act signature",outline=F,ylim=c(0,1))
dev.off()




#============================================================================================================================================
# use half violin plot to show result
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig2_GZMB_RidgePlot.pdf",useDingbats=F)
RidgePlot(sub.dat,assay="RNA",features=c("GZMA","GZMB"))
dev.off()




# 2022-3-19 
# test use cytotoxic gene in PNAS 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(c("GZMA","GZMB","GZMH","PRF1"),"Tcell_cyto",out="/public/workspace/lily/CTN/version_3_20/Signature/Tcell_cyto.mod") 

library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
sub.dat <- subset(dat,cells=which(dat$celltype.refine!="Unclassify"))

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.rs <- mod.analyze2(as.matrix(sub.dat[['RNA']]@data),c("Tcell_cyto"),"/public/workspace/lily/CTN/version_3_20/Signature/",permN=0)
mod.rs <- as.data.frame(mod.rs)
saveRDS(mod.rs,file="/public/workspace/lily/CTN/version_3_20/data/Tcell_cytotoxic_pnas_mod.RDS")



















