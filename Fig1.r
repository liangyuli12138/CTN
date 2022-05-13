
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






# 2. cell type in different group cells
#===============================================================================================
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig1_celltype_group.pdf",width=8,useDingbats=F)
DimPlot(dat,cells=which(dat$group=="Elder"),reduction="tsne",group.by="celltype",cols=col,raster=F)+labs(title="Elder") # control
DimPlot(dat,cells=which(dat$group=="CTN"),reduction="tsne",group.by="celltype",cols=col,raster=F)+labs(title="Centenarian") # cetenraian
DimPlot(dat,cells=which(dat$group=="SC"),reduction="tsne",group.by="celltype",cols=col,raster=F)+labs(title="Supercentenarian") # supercentenarain
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
























