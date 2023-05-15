
# Figure1 plot 

# 0. maybe a work flow picture 








# 1. overview big figure 
#===============================================================================================
# 2021-3-23
################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")

# B cell CD14 Mono CD16 Mono DC NK pDC Platelet Proliferating T cell T cell
col <- c("#f47b7b","#91be3e","#96cbb3","#836eaa","#39a6dd","#207c88","#e5352b","#ef9020","#dde2e0","#d4d7da")

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig1_overview.pdf",width=8,useDingbats=F)
DimPlot(dat,group.by="celltype",cols=col,raster=F,reduction="tsne")
dev.off()


# 2022-12-29
# plot by umap 
dat <- RunUMAP(dat,dims=1:10)

col <- c("#f47b7b","#91be3e","#96cbb3","#836eaa","#39a6dd","#207c88","#e5352b","#ef9020","#dde2e0","#d4d7da")

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig1_overview_umap.pdf",width=8,useDingbats=F)
DimPlot(dat,group.by="celltype",cols=col,raster=T,reduction="umap")
dev.off()


# 2023-1-4
# supplementary plot 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")
dat <- RunUMAP(dat,dims=1:10)

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/supplmentary/CTN_inte_FeaturePlot_umap.pdf",useDingbats=F)
DimPlot(dat,raster=T,reduction="umap")
DimPlot(dat,group.by="sample",raster=T,reduction="umap")
DimPlot(dat,group.by="group",cols=c("#003468","#5ec6f2","#4b1702"),raster=T,reduction="umap")

FeaturePlot(dat,features=c("CD3D","CD3E","MS4A1","CD79A"),raster=T,reduction="umap") # T cell 
FeaturePlot(dat,features=c("KLRD1","NKG7","CD14","MS4A7"),raster=T,reduction="umap") # B cell  Plasma
FeaturePlot(dat,features=c("TOP2A","MKI67","GP9","PF4"),raster=T,reduction="umap") # Prolifering Cells
dev.off()


# 2. cell type in different group cells
#===============================================================================================
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig1_celltype_group.pdf",width=8,useDingbats=F)
DimPlot(dat,cells=which(dat$group=="Elder"),reduction="tsne",group.by="celltype",cols=col,raster=F)+labs(title="Elder") # control
DimPlot(dat,cells=which(dat$group=="CTN"),reduction="tsne",group.by="celltype",cols=col,raster=F)+labs(title="Centenarian") # cetenraian
DimPlot(dat,cells=which(dat$group=="SC"),reduction="tsne",group.by="celltype",cols=col,raster=F)+labs(title="Supercentenarian") # supercentenarain
dev.off()


pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig1_celltype_group_umap.pdf",width=8,useDingbats=F)
DimPlot(dat,cells=which(dat$group=="Elder"),reduction="umap",group.by="celltype",cols=col,raster=T)+labs(title="Elder") # control
DimPlot(dat,cells=which(dat$group=="CTN"),reduction="umap",group.by="celltype",cols=col,raster=T)+labs(title="Centenarian") # cetenraian
DimPlot(dat,cells=which(dat$group=="SC"),reduction="umap",group.by="celltype",cols=col,raster=T)+labs(title="Supercentenarian") # supercentenarain
dev.off()





# 3. vlnplot of marker genes 
#==============================================================================================
# my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
#          '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
#          '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
#          '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
#          '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
#          '#968175')
col <- c("#f47b7b","#91be3e","#96cbb3","#836eaa","#39a6dd","#207c88","#e5352b","#ef9020","#dde2e0","#d4d7da")

library(Seurat)
library(ggplot2)
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +
               theme(legend.position = "none",
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_line(),
               axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
               plot.margin = plot.margin )
       return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
       plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
            plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
            theme(axis.text.x=element_text(), axis.ticks.x = element_line())
       p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
       return(p)
}

#
dat@active.ident <- factor(dat$celltype)
markergene <- c("CD3D","CD3E","MS4A1","CD79A","CD14","MS4A7","CLEC4C","SERPINF1","NKG7","KLRD1",
    "TOP2A","MKI67","GP9","PF4")

# Fig1
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig1_violinplot.pdf",width=10,height=10)
StackedVlnPlot(dat, markergene[1:7], pt.size=0, cols=col)
StackedVlnPlot(dat, markergene[8:14], pt.size=0, cols=col)
dev.off()




# barplot show median per group 
#========================================================================================================
tmp <- table(dat$sample,dat$celltype)
tmp.res <- apply(tmp,1,function(x){x/sum(x)})
# pdf("tmp.pdf")
# for(i in 1:9){
#     boxplot(tmp.res[i, 1:4], tmp.res[i, 12:14],tmp.res[i, 5:11],main=rownames(tmp.res)[i],names=c("Elder","CTN","SCT"))
# }
# dev.off()

# apply(tmp.res,1,function(x){
#     c(wilcox.test(x[1:4],x[12:14])$p.value,wilcox.test(x[1:4],x[5:11])$p.value)
# })

rs.f <- apply(tmp.res,1,function(x){
    c(median(x[1:4]),median(x[12:14]),median(x[5:11]))
})


dat.f <- rs.f[,c(1,2,4,8)]
rownames(dat.f) <- c("Elder","CTN","SC")

# plot result 
library(reshape2)
library(ggplot2)
dat.f <- melt(dat.f)
colnames(dat.f) <- c("group","celltype","value")
cols <- c("#7ac143","#5091cd","#f9a541")
dat.f$celltype <- factor(dat.f$celltype,levels=c("T cell","NK cell","CD14 Mono","B cell"))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig1_barplot.pdf")
ggplot(dat.f,aes(x=celltype,y=value,fill=group,group=group))+ geom_bar(stat="identity",position="dodge")+
    scale_fill_manual(values=cols)+theme_bw()+ labs(y="percentage")
dev.off()



# 2022-12-29
# change group color 
# young, Elder,CTN,SC
# "#d4c78c","#5ec6f2","#003468","#4b1702"
# and subsample to calculate 

library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")
dat <- RunUMAP(dat,dims=1:10)

Tsub <- subset(dat,cells=which(dat$celltype=="T cell"))

dat1 <- c()
dat2 <- c()
dat3 <- c()

i=1
while(i <= 10){
    set.seed(123*i)
    cellsElder <- sample(which(Tsub$group=="Elder"),5000)
    # print(head(cellsElder))
    TsubElder <- subset(dat,cells=cellsElder)
    set.seed(1234*i)
    cellsCTN <- sample(which(Tsub$group=="CTN"),5000)
    TsubCTN <- subset(dat,cells=cellsCTN)
    set.seed(12345*i)
    cellsSC <- sample(which(Tsub$group=="SC"),5000)
    TsubSC <- subset(dat,cells=cellsSC)

    dat1 <- c(dat1,length(which(TsubElder$seurat_clusters%in%c(3,8)))/5000)
    dat2 <- c(dat2,length(which(TsubCTN$seurat_clusters%in%c(3,8)))/5000)
    dat3 <- c(dat3,length(which(TsubSC$seurat_clusters%in%c(3,8)))/5000)
    i=i+1
}

pdat <- data.frame(percent=c(dat1,dat2,dat3),group=c(rep(c("Elder","CTN","SC"),each=10)))
pdat$group <- factor(pdat$group,levels=c("Elder","CTN","SC"))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig1_sample_Tcellsubset.pdf",useDingbats=F)
# boxplot(percent~group,data=pdat,names=c("Elder","CTN","SC"),main="cell percentage in TG1",ylims=c(0.05,0.25),outline=F)
barplot(c(mean(dat1),mean(dat2),mean(dat3)),names=c("Elder","CTN","SC"),main="cell percentage in TG1",ylim=c(0,1))
beeswarm::beeswarm(percent~group,data=pdat,col = 4, pch = 16,add = TRUE)
dev.off()


pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig1_sample_Tcellsubset_long.pdf",useDingbats=F)
# boxplot(percent~group,data=pdat,names=c("Elder","CTN","SC"),main="cell percentage in TG1",ylims=c(0.05,0.25),outline=F)
barplot(c(mean(1-dat1),mean(1-dat2),mean(1-dat3)),names=c("Elder","CTN","SC"),main="cell percentage in TG1",ylim=c(0,1))
beeswarm::beeswarm((1-percent)~group,data=pdat,col = 4, pch = 16,add = TRUE)
dev.off()










