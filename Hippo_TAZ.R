# !/usr/bin/Rscript
# 2021-10-6
# this program is used to analysis T cell (CD4)
#====================================================================================================================================================
# 
#  1. calculate Hippo signaling pathway 

library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")

res.list <- list()
type <- unique(dat$group)  # except unclassify group 
for(i in 1:length(type)){
	sub.dat <- subset(dat,cells=which(dat$group==type[i]))
	source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
	mod <- mod.analyze2(as.matrix(sub.dat[['RNA']]@data),c("KEGG_hippo","Reactome_hippo"),"/public/workspace/lily/CTN/",permN=0)
	mod <- as.data.frame(mod)
	mod$type <- sub.dat$celltype.refine
	mod$cluster <- sub.dat$seurat_clusters
	res.list[[i]] <- mod
	names(res.list)[i] <- type[i]
}




saveRDS(res.list,file="/public/workspace/lily/CTN/version_3_20/data/Tcell_Hippo.RDS")


#==================================================================================================================================================
# plot result 
library(ggplot2)
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Tcell_Hippo_3group.pdf",useDingbats=F)
barplot(aggregate(KEGG_hippo_norm~type,data=res.list[[2]],FUN=median)[,2][-6],main="KEGG_hippo pathway (Elder)",ylab="Signature score",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,0.5))
barplot(aggregate(KEGG_hippo_norm~type,data=res.list[[3]],FUN=median)[,2][-6],main="KEGG_hippo pathway (CTN)",ylab="Signature score",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,0.5))
barplot(aggregate(KEGG_hippo_norm~type,data=res.list[[1]],FUN=median)[,2][-6],main="KEGG_hippo pathway (SC)",ylab="Signature score",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,0.5))
dev.off()

normal <- function(x){
	(x-min(x))/max(x)
}

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Tcell_Hippo_3group.normalized.pdf",useDingbats=F)
barplot(scale(aggregate(KEGG_hippo_norm~type,data=res.list[[2]],FUN=median)[,2][-6],center=F)[,1],main="KEGG_hippo pathway (Elder)",ylab="Signature score",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,1))
barplot(scale(aggregate(KEGG_hippo_norm~type,data=res.list[[3]],FUN=median)[,2][-6],center=F)[,1],main="KEGG_hippo pathway (CTN)",ylab="Signature score",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,1))
barplot(scale(aggregate(KEGG_hippo_norm~type,data=res.list[[1]],FUN=median)[,2][-6],center=F)[,1],main="KEGG_hippo pathway (SC)",ylab="Signature score",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,1))
dev.off()


res.dat <- data.frame(Elder=scale(aggregate(KEGG_hippo_norm~type,data=res.list[[2]],FUN=median)[,2][-6],center=F)[,1],
	CTN=scale(aggregate(KEGG_hippo_norm~type,data=res.list[[3]],FUN=median)[,2][-6],center=F)[,1],
	SC=scale(aggregate(KEGG_hippo_norm~type,data=res.list[[1]],FUN=median)[,2][-6],center=F)[,1])

rownames(res.dat) <- c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs")

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Tcell_Hippo_3group.normalized.heatmap.pdf",useDingbats=F)
pheatmap::pheatmap(res.dat,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),display_numbers=T,number_color="black",
	fontsize_number=12,main="KEGG Hippo signaling")
dev.off()








# 2. calculate YAP and TAFAZZIN gene 
library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
DefaultAssay(tmp.dat) <- "RNA"
dat <- subset(tmp.dat,cells=colnames(tmp.dat)[-which(tmp.dat$celltype.refine=="Unclassify")])
dat$type.group <- paste0(dat$group,"_",dat$celltype.refine)


percent_feature <- function(dat,genelist,group){
    res.list <- c()
    for(i in 1:length(genelist)){
        dat$tmp_gene <- ifelse(dat[["RNA"]]@data[genelist[i],]>0,"Y","N")
        if(all(dat$tmp_gene=="N")){
            tmpp <- c(0,0,0,0)
            names(tmpp) <- names(table(dat@meta.data[,group]))
           res.list[[i]] <- tmpp
        }else{
           res.list[[i]] <- apply(table(dat$tmp_gene,dat@meta.data[,group]),2,function(x){x/sum(x)})[2,] 
        }
        
        names(res.list)[i] <- genelist[i]
    }
    return(res.list)
}

res <- percent_feature(dat,genelist=c("TAZ"),group="type.group")




# plot result 
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Tcell_TAZ_3group.normalized.pdf",useDingbats=F)
barplot(scale(res$TAZ[6:10],center=F)[,1],main="TAZ expression percentage (Elder)",ylab="Scaled percentage",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,1.5))

barplot(scale(res$TAZ[1:5],center=F)[,1],main="TAZ expression percentage (CTN)",ylab="Scaled percentage",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,1.5))

barplot(scale(res$TAZ[11:15],center=F)[,1],main="TAZ expression percentage (SC)",ylab="Scaled percentage",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,1.5))
dev.off()



# another ways to plot 
res.dat <- data.frame(Elder=scale(res$TAZ[6:10],center=F)[,1],
	CTN=scale(res$TAZ[1:5],center=F)[,1],
	SC =scale(res$TAZ[11:15],center=F)[,1])
rownames(res.dat)=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs")

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Tcell_TAZ_3group.normalized.heatmap.pdf",useDingbats=F)
pheatmap::pheatmap(res.dat,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),display_numbers=T,number_color="black",
	fontsize_number=12,main="TAZ expressed percentage")
dev.off()






















#=========================================================================================================================================
# 2021-11-4
# calculate IL-19 and TGF-Beta signaling pathway 
# all IL-19 expression is 0


library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
DefaultAssay(tmp.dat) <- "RNA"
dat <- subset(tmp.dat,cells=colnames(tmp.dat)[-which(tmp.dat$celltype.refine=="Unclassify")])
dat$type.group <- paste0(dat$group,"_",dat$celltype.refine)


percent_feature <- function(dat,genelist,group){
    res.list <- c()
    for(i in 1:length(genelist)){
        dat$tmp_gene <- ifelse(dat[["RNA"]]@data[genelist[i],]>0,"Y","N")
        if(all(dat$tmp_gene=="N")){
            tmpp <- rep(0,length(table(dat@meta.data[,group])))
            names(tmpp) <- names(table(dat@meta.data[,group]))
           res.list[[i]] <- tmpp
        }else{
           res.list[[i]] <- apply(table(dat$tmp_gene,dat@meta.data[,group]),2,function(x){x/sum(x)})[2,] 
        }
        
        names(res.list)[i] <- genelist[i]
    }
    return(res.list)
}

res <- percent_feature(dat,genelist=c("IL19"),group="type.group")


#==================================================================================================================================================
tmp.dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
DefaultAssay(tmp.dat) <- "RNA"
dat <- subset(tmp.dat,cells=colnames(tmp.dat)[-which(tmp.dat$celltype.refine=="Unclassify")])

res.list <- list()
type <- unique(dat$group)  # except unclassify group 
for(i in 1:length(type)){
	sub.dat <- subset(dat,cells=which(dat$group==type[i]))
	source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
	mod <- mod.analyze2(as.matrix(sub.dat[['RNA']]@data),c("HALLMARK_TGF_BETA_SIGNALING"),"/public/workspace/lily/MOD_file/HALLMARK/",permN=0)
	mod <- as.data.frame(mod)
	mod$type <- sub.dat$celltype.refine
	mod$cluster <- sub.dat$seurat_clusters
	res.list[[i]] <- mod
	names(res.list)[i] <- type[i]
}

saveRDS(res.list,file="/public/workspace/lily/CTN/version_3_20/data/Tcell_TGFB.RDS")





normal <- function(x){
	(x-min(x))/max(x)
}

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Tcell_Hippo_3group.normalized.pdf",useDingbats=F)
barplot(scale(aggregate(KEGG_hippo_norm~type,data=res.list[[2]],FUN=median)[,2][-6],center=F)[,1],main="KEGG_hippo pathway (Elder)",ylab="Signature score",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,1))
barplot(scale(aggregate(KEGG_hippo_norm~type,data=res.list[[3]],FUN=median)[,2][-6],center=F)[,1],main="KEGG_hippo pathway (CTN)",ylab="Signature score",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,1))
barplot(scale(aggregate(KEGG_hippo_norm~type,data=res.list[[1]],FUN=median)[,2][-6],center=F)[,1],main="KEGG_hippo pathway (SC)",ylab="Signature score",xlab="T cells",
	col=c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3"),names=c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs"),ylim=c(0,1))
dev.off()


res.dat <- data.frame(Elder=scale(aggregate(HALLMARK_TGF_BETA_SIGNALING_norm~type,data=res.list[[2]],FUN=median)[,2][-6],center=F)[,1],
	CTN=scale(aggregate(HALLMARK_TGF_BETA_SIGNALING_norm~type,data=res.list[[3]],FUN=median)[,2][-6],center=F)[,1],
	SC=scale(aggregate(HALLMARK_TGF_BETA_SIGNALING_norm~type,data=res.list[[1]],FUN=median)[,2][-6],center=F)[,1])

rownames(res.dat) <- c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs")

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Tcell_Hippo_3group.normalized.heatmap.pdf",useDingbats=F)
pheatmap::pheatmap(res.dat,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),display_numbers=T,number_color="black",
	fontsize_number=12,main="TGF-beta signaling")
dev.off()


























# 2021-11-26 
# analysis about GZMK+ PD1+ T cells
#==================================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
DefaultAssay(dat) <- "RNA"
dat <- RunUMAP(dat,dims=1:10)

FeaturePlot(dat,features=c("PDCD1","GZMK"))
DimPlot(dat,label=T,label.size=4)




# 2021-12-9
# analysis about Tem about GZMK GZMA GZMB and so on 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
DefaultAssay(dat) <- "RNA"
dat <- RunUMAP(dat,dims=1:10)

sub.dat <- subset(dat,cells=which(dat$celltype.refine=="CD8+ Tem"))
sub.dat@active.ident <- factor(sub.dat$group)

RidgePlot(sub.dat,assay="RNA",features=c("GZMA","GZMB"))


FeaturePlot(dat,features=c("PDCD1","GZMK"))
DimPlot(dat,label=T,label.size=4)














#=====================================================================================================================================================
# 2022-1-4 
# analysis about SASP, TERT and some signature in MSGDB
#
#=====================================================================================================================================================
# make mod 
library(qusage)
tmp.sig <- read.gmt("/public/workspace/lily/CTN/version_3_20/Signature.gmt")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
for(i in 1:length(tmp.sig)){
	mod.generate(tmp.sig[[i]],names(tmp.sig)[[i]],out=paste0("/public/workspace/lily/CTN/version_3_20/Signature/",names(tmp.sig)[[i]],".mod")) # make a mod file 
}



# calculate result 
# for scRNA 
#====================================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
modfile <- gsub("\\.mod$","",dir("/public/workspace/lily/CTN/version_3_20/Signature/"))
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),modfile,"/public/workspace/lily/CTN/version_3_20/Signature/",permN=0)
mod <- as.data.frame(mod)
saveRDS(mod,file="/public/workspace/lily/CTN/version_3_20/data/Tcel_multiSig_mod.RDS")


mod$group <- dat$group
mod$celltype.refine <- dat$celltype.refine
mod$sample <- dat$sample




#====================================================================================================================================================
# plot result 
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/scRNA_Signature_sample.pdf",useDingbats=F)
tmp.dat <- mod[which(mod$celltype.refine=="CD8+ Tem"),]
for(i in 9:16){
	tmp.res <- aggregate(tmp.dat[,i]~sample,data=tmp.dat,FUN=median)
	colnames(tmp.res)[2] <- colnames(tmp.dat)[i]
	tmp.res$group <- c(rep("Elder",4),rep("SC",7),rep("CTN",3))
	tmp.res$group <- factor(tmp.res$group,levels=c("Elder","CTN","SC"))
	boxplot(tmp.res[,2]~group,data=tmp.res,main=gsub("_norm$","",colnames(tmp.res)[2]),outline=F,ylim=c(0,0.5))
	legend("topright",legend=(paste0("p=",kruskal.test(tmp.res[,2]~group,data=tmp.res)$p.value)))
	beeswarm::beeswarm(tmp.res[,2]~group,data=tmp.res,col = 4, pch = 16,add = TRUE)
}
dev.off()



# plot result by group 
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/scRNA_Signature_group.pdf",useDingbats=F)
tmp.dat <- mod[which(mod$celltype.refine=="CD8+ Tem"),]
for(i in 9:16){
	boxplot(tmp.dat[,i]~group,data=tmp.dat,main=gsub("_norm$","",colnames(tmp.dat)[i]),outline=F,ylim=c(0,1))
	legend("topright",legend=(paste0("p=",kruskal.test(tmp.dat[,i]~group,data=tmp.dat)$p.value)))
	# beeswarm::beeswarm(tmp.dat[,i]~group,data=tmp.dat,col = 4, pch = 16,add = TRUE)
}
dev.off()



dat$longevity_score <- mod$BIOCARTA_LONGEVITY_PATHWAY_norm
dat$reactiom_Hippo_score <- mod$REACTOME_SIGNALING_BY_HIPPO_norm












