

# use python to load dta
conda activate scanpy
python

import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read("/public/workspace/lily/CTN/version_3_20/Response/GSE158055/COVID19_ALL.h5ad")

adata_sub = adata[adata.obs['Outcome'].isin(['control'])]
adata_sub.write_h5ad("/public/workspace/lily/CTN/version_3_20/Response/GSE158055/All_Control.h5ad")
adata_sub.write_csvs("/public/workspace/lily/CTN/version_3_20/Response/GSE158055/All_Control_meta.csv", skip_data=False)

# bytlib load r/4.1.2
library(SeuratDisk)
library(Seurat)
Convert('/public/workspace/lily/CTN/version_3_20/Response/GSE158055/All_Control.h5ad', "h5seurat",
        overwrite = TRUE,assay = "RNA")
data <- LoadH5Seurat("/public/workspace/lily/CTN/version_3_20/Response/GSE158055/All_Control.h5seurat")
metadata <- read.csv("~/CTN/version_3_20/Response/GSE158055/All_Control_meta/obs.csv")
rownames(metadata) <- metadata$X
metadata <- as.data.frame(metadata)
all(rownames(metadata)==colnames(data))
data@meta.data <- metadata
saveRDS(data,file="~/CTN/version_3_20/Response/GSE158055/All_Control.RDS")
#===========================================================================================================================================
# GSE158055
dat <- readRDS("~/CTN/version_3_20/Response/GSE158055/All_Control.RDS")
subdat <- subset(dat,cells=which(dat$Sample.type%in%c("fresh PBMC","frozen PBMC"))) # do not use FACS
subdat$CellType <- "unknow"
subdat$CellType[which(subdat$majorType%in%c("B","Plasma"))] <- "Bcell"
subdat$CellType[which(subdat$majorType%in%c("CD4","CD8"))] <- "Tcell"
subdat$CellType[which(subdat$majorType%in%c("Mono"))] <- "CD14+ Monocytes"
subdat$CellType[which(subdat$celltype%in%c("Mono_c5-CD16"))] <- "CD16+ Monocytes"
subdat$CellType[which(subdat$majorType%in%c("DC"))] <- "DC"
subdat$CellType[which(subdat$majorType%in%c("NK"))] <- "NK"


tmp.dat <- data.frame(reshape2::melt(apply(table(subdat$CellType,subdat$sampleID),2,function(x){x/sum(x)})))
colnames(tmp.dat) <- c("Celltype","Sample","Percent")
tmpsex <- data.frame(reshape2::melt(table(subdat$sampleID,subdat$Sex)))
tmpsex <- tmpsex[which(tmpsex$value>0),]
tmp.dat <- merge(tmp.dat,tmpsex[,c(1,2)],by.x="Sample",by.y="Var1")
colnames(tmp.dat)[4] <- "Sex"
tmp.dat$Sex <- as.character(unlist(tmp.dat$Sex))
tmp.dat$Sex[which(tmp.dat$Sex%in%c("F"))] <- "Female"
tmp.dat$Sex[which(tmp.dat$Sex%in%c("M"))] <- "Male"

saveRDS(tmp.dat,file="~/CTN/version_3_20/Response/GSE158055/GSE158055_analysis_data.RDS")





# GSE190992
# analysis 
library(Seurat)
dat <- readRDS("~/CTN/version_3_20/Response/GSE190992/GSE190992_AIFI-scRNA-PBMC-FinalData.RDS")
dat$CellType <- "unknow"
dat$CellType[which(dat$seurat_pbmc_type%in%c("B cell progenitor","pre-B cell"))] <- "Bcell"
dat$CellType[which(dat$seurat_pbmc_type%in%c("CD14+ Monocytes"))] <- "CD14+ Monocytes"
dat$CellType[which(dat$seurat_pbmc_type%in%c("CD16+ Monocytes"))] <- "CD16+ Monocytes"
dat$CellType[which(dat$seurat_pbmc_type%in%c("CD4 Memory","CD4 Naive","CD8 effector","CD8 Naive","Double negative T cell"))] <- "Tcell"
dat$CellType[which(dat$seurat_pbmc_type%in%c("Dendritic cell"))] <- "DC"
dat$CellType[which(dat$seurat_pbmc_type%in%c("pDC"))] <- "pDC"
dat$CellType[which(dat$seurat_pbmc_type%in%c("NK cell"))] <- "NK"
dat$CellType[which(dat$seurat_pbmc_type%in%c("Platelets"))] <- "Platelets"

tmp.dat <- data.frame(reshape2::melt(apply(table(dat$CellType,dat$Sample),2,function(x){x/sum(x)})))
colnames(tmp.dat) <- c("Celltype","Sample","Percent")
tmp.dat$ID <- sapply(strsplit(as.character(tmp.dat$Sample),"W"),function(x){x[[1]]})
tmp.dat$Sex <- "Male"
tmp.dat$Sex[which(tmp.dat$ID%in%c("PB5206"))] <- "Female"

saveRDS(tmp.dat,file="~/CTN/version_3_20/Response/GSE190992/GSE190992_analysis_data.RDS")


#======================================================================================================================
# plot result 
library(ggplot2)
library(ggpubr)
dat1 <- readRDS("~/CTN/version_3_20/Response/GSE158055/GSE158055_analysis_data.RDS")
dat1$Group <- "China"
dat2 <- readRDS("~/CTN/version_3_20/Response/GSE190992/GSE190992_analysis_data.RDS")
dat2$Group <- "America"

subdat1 <- dat1[which(dat1$Celltype%in%c("Bcell","Tcell","CD14+ Monocytes","CD16+ Monocytes","NK","DC")),]
subdat2 <- dat2[which(dat2$Celltype%in%c("Bcell","Tcell","CD14+ Monocytes","CD16+ Monocytes","NK","DC")),]
subdat2$ID <- NULL

dat <- rbind(subdat1,subdat2)

pdf("~/CTN/version_3_20/Response/locality_diff_barplot.pdf",useDingbats=F)
ggplot(dat,aes(x=Celltype,y=Percent,fill=Group))+
	geom_bar(position=position_dodge(0.8),stat="summary",width=0.5,size=1)+
	stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black",width = 0.2,position = position_dodge(0.8))+
	geom_point(data = dat, aes(x=Celltype,y = Percent,fill=Group),position = position_dodge(0.8),size = 3, shape = 21)+
	scale_fill_manual(values = c('#76b852','#fd9f3e'))+ theme_bw() +
	theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
dev.off()
wilcox.test(dat$Percent[which(dat$Celltype=="Tcell"&dat$Group=="China")],dat$Percent[which(dat$Celltype=="Tcell"&dat$Group=="America")])
wilcox.test(dat$Percent[which(dat$Celltype=="Bcell"&dat$Group=="China")],dat$Percent[which(dat$Celltype=="Bcell"&dat$Group=="America")])
wilcox.test(dat$Percent[which(dat$Celltype=="CD14+ Monocytes"&dat$Group=="China")],dat$Percent[which(dat$Celltype=="CD14+ Monocytes"&dat$Group=="America")])
wilcox.test(dat$Percent[which(dat$Celltype=="CD16+ Monocytes"&dat$Group=="China")],dat$Percent[which(dat$Celltype=="CD16+ Monocytes"&dat$Group=="America")])
wilcox.test(dat$Percent[which(dat$Celltype=="NK"&dat$Group=="China")],dat$Percent[which(dat$Celltype=="NK"&dat$Group=="America")])
wilcox.test(dat$Percent[which(dat$Celltype=="DC"&dat$Group=="China")],dat$Percent[which(dat$Celltype=="DC"&dat$Group=="America")])



# plot for sex
pdf("~/CTN/version_3_20/Response/sex_diff_barplot_GSE158055.pdf",useDingbats=F)

ggplot(subdat1,aes(x=Celltype,y=Percent,fill=Sex))+
	geom_bar(position=position_dodge(0.8),stat="summary",width=0.5,size=1)+
	stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black",width = 0.2,position = position_dodge(0.8))+
	geom_point(data = subdat1, aes(x=Celltype,y = Percent,fill=Sex),position = position_dodge(0.8),size = 3, shape = 21)+
	scale_fill_manual(values = c('#e18283','#5494cc'))+ theme_bw() +
	theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
dev.off()
wilcox.test(subdat1$Percent[which(subdat1$Celltype=="Tcell"&subdat1$Sex=="Male")],subdat1$Percent[which(subdat1$Celltype=="Tcell"&subdat1$Sex=="Female")])
wilcox.test(subdat1$Percent[which(subdat1$Celltype=="Bcell"&subdat1$Sex=="Male")],subdat1$Percent[which(subdat1$Celltype=="Bcell"&subdat1$Sex=="Female")])
wilcox.test(subdat1$Percent[which(subdat1$Celltype=="CD14+ Monocytes"&subdat1$Sex=="Male")],subdat1$Percent[which(subdat1$Celltype=="CD14+ Monocytes"&subdat1$Sex=="Female")])
wilcox.test(subdat1$Percent[which(subdat1$Celltype=="CD16+ Monocytes"&subdat1$Sex=="Male")],subdat1$Percent[which(subdat1$Celltype=="CD16+ Monocytes"&subdat1$Sex=="Female")])
wilcox.test(subdat1$Percent[which(subdat1$Celltype=="NK"&subdat1$Sex=="Male")],subdat1$Percent[which(subdat1$Celltype=="NK"&subdat1$Sex=="Female")])
wilcox.test(subdat1$Percent[which(subdat1$Celltype=="DC"&subdat1$Sex=="Male")],subdat1$Percent[which(subdat1$Celltype=="DC"&subdat1$Sex=="Female")])

pdf("~/CTN/version_3_20/Response/sex_diff_barplot_GSE190992.pdf",useDingbats=F)

ggplot(subdat2,aes(x=Celltype,y=Percent,fill=Sex))+
	geom_bar(position=position_dodge(0.8),stat="summary",width=0.5,size=1)+
	stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black",width = 0.2,position = position_dodge(0.8))+
	geom_point(data = subdat2, aes(x=Celltype,y = Percent,fill=Sex),position = position_dodge(0.8),size = 3, shape = 21)+
	scale_fill_manual(values = c('#e18283','#5494cc'))+ theme_bw() +
	theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
dev.off()
wilcox.test(subdat2$Percent[which(subdat2$Celltype=="Tcell"&subdat2$Sex=="Male")],subdat2$Percent[which(subdat2$Celltype=="Tcell"&subdat2$Sex=="Female")])
wilcox.test(subdat2$Percent[which(subdat2$Celltype=="Bcell"&subdat2$Sex=="Male")],subdat2$Percent[which(subdat2$Celltype=="Bcell"&subdat2$Sex=="Female")])
wilcox.test(subdat2$Percent[which(subdat2$Celltype=="CD14+ Monocytes"&subdat2$Sex=="Male")],subdat2$Percent[which(subdat2$Celltype=="CD14+ Monocytes"&subdat2$Sex=="Female")])
wilcox.test(subdat2$Percent[which(subdat2$Celltype=="CD16+ Monocytes"&subdat2$Sex=="Male")],subdat2$Percent[which(subdat2$Celltype=="CD16+ Monocytes"&subdat2$Sex=="Female")])
wilcox.test(subdat2$Percent[which(subdat2$Celltype=="NK"&subdat2$Sex=="Male")],subdat2$Percent[which(subdat2$Celltype=="NK"&subdat2$Sex=="Female")])
wilcox.test(subdat2$Percent[which(subdat2$Celltype=="DC"&subdat2$Sex=="Male")],subdat2$Percent[which(subdat2$Celltype=="DC"&subdat2$Sex=="Female")])


pdf("~/CTN/version_3_20/Response/sex_diff_barplot.pdf",useDingbats=F)
ggplot(dat,aes(x=Celltype,y=Percent,fill=Sex))+
	geom_bar(position=position_dodge(0.8),stat="summary",width=0.5,size=1)+
	stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black",width = 0.2,position = position_dodge(0.8))+
	geom_point(data = dat, aes(x=Celltype,y = Percent,fill=Sex),position = position_dodge(0.8),size = 3, shape = 21)+
	scale_fill_manual(values = c('#e18283','#5494cc'))+ theme_bw() +
	theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
dev.off()

wilcox.test(dat$Percent[which(dat$Celltype=="Tcell"&dat$Sex=="Male")],dat$Percent[which(dat$Celltype=="Tcell"&dat$Sex=="Female")])
wilcox.test(dat$Percent[which(dat$Celltype=="Bcell"&dat$Sex=="Male")],dat$Percent[which(dat$Celltype=="Bcell"&dat$Sex=="Female")])
wilcox.test(dat$Percent[which(dat$Celltype=="CD14+ Monocytes"&dat$Sex=="Male")],dat$Percent[which(dat$Celltype=="CD14+ Monocytes"&dat$Sex=="Female")])
wilcox.test(dat$Percent[which(dat$Celltype=="CD16+ Monocytes"&dat$Sex=="Male")],dat$Percent[which(dat$Celltype=="CD16+ Monocytes"&dat$Sex=="Female")])
wilcox.test(dat$Percent[which(dat$Celltype=="NK"&dat$Sex=="Male")],dat$Percent[which(dat$Celltype=="NK"&dat$Sex=="Female")])
wilcox.test(dat$Percent[which(dat$Celltype=="DC"&dat$Sex=="Male")],dat$Percent[which(dat$Celltype=="DC"&dat$Sex=="Female")])


















#===============================================================================================================================================
# 2023-2-13
# anlysis GSE158055 as young 
# bytlib load r/4.1.2
library(Seurat)
dat <- readRDS("~/CTN/version_3_20/Response/GSE158055/All_Control.RDS")
subdat <- subset(dat,cells=which(dat$Sample.type%in%c("fresh PBMC","frozen PBMC") & dat$Age < 50)) # do not use FACS
subdat$CellType <- "unknow"
subdat$CellType[which(subdat$majorType%in%c("B","Plasma"))] <- "Bcell"
subdat$CellType[which(subdat$majorType%in%c("CD4","CD8"))] <- "Tcell"
subdat$CellType[which(subdat$majorType%in%c("Mono"))] <- "CD14+ Monocytes"
subdat$CellType[which(subdat$celltype%in%c("Mono_c5-CD16"))] <- "CD16+ Monocytes"
subdat$CellType[which(subdat$majorType%in%c("DC"))] <- "DC"
subdat$CellType[which(subdat$majorType%in%c("NK"))] <- "NK"

# centenarians data
CTN <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")
tmp <- table(CTN$sample,CTN$celltype)
tmp.res <- apply(tmp,1,function(x){x/sum(x)})

################
# get percentage 
young <- apply(table(subdat$sampleID,subdat$CellType),1,function(x){x/sum(x)})[c(1,2,5,6),]
old <- tmp.res[c(1,2,4,8),]

# cbind data
pdat <- data.frame(reshape2::melt(cbind(young,old)),stringsAsFactors=F)
colnames(pdat) <- c("Celltype","Sample","Percent")
pdat$Group <- "Young"
pdat$Group[grep("^XWS",pdat$Sample)] <- "CTN"
pdat$Group[grep("^SC",pdat$Sample)] <- "SC"
pdat$Group[grep("^CT",pdat$Sample)] <- "Elder"

pdat$Group <- factor(pdat$Group,levels=c("Young","Elder","CTN","SC"))


# 1. cell type proportion
library(ggplot2)
pdf("/public/workspace/lily/CTN/version_3_20/Response/Young_CTN_cell_proportion.pdf",useDingbats=F)
ggplot(pdat,aes(x=Celltype,y=Percent,fill=Group))+
	geom_bar(fun="median",position=position_dodge(0.8),stat="summary",width=0.5,size=1)+
	stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black",width = 0.2,position = position_dodge(0.8))+
	geom_point(data = pdat, aes(x=Celltype,y = Percent,fill=Group),position = position_dodge(0.8),size = 3, shape = 21)+
	scale_fill_manual(values = c("#d4c78c","#5ec6f2","#003468","#4b1702"))+ theme_bw() +
	theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
dev.off()






# 2. T cell gene marker and subset proportion
dat <- readRDS("~/CTN/version_3_20/Response/GSE158055/All_Control.RDS")
subdat <- subset(dat,cells=which(dat$Sample.type%in%c("fresh PBMC","frozen PBMC") & dat$Age < 50)) # do not use FACS
subdat$CellType <- "unknow"
subdat$CellType[which(subdat$majorType%in%c("B","Plasma"))] <- "Bcell"
subdat$CellType[which(subdat$majorType%in%c("CD4","CD8"))] <- "Tcell"
subdat$CellType[which(subdat$majorType%in%c("Mono"))] <- "CD14+ Monocytes"
subdat$CellType[which(subdat$celltype%in%c("Mono_c5-CD16"))] <- "CD16+ Monocytes"
subdat$CellType[which(subdat$majorType%in%c("DC"))] <- "DC"
subdat$CellType[which(subdat$majorType%in%c("NK"))] <- "NK"
Tdat <- subset(subdat,cells=which(subdat$CellType=="Tcell"))

percent_feature <- function(dat,genelist,group){
    res.list <- c()
    for(i in 1:length(genelist)){
        dat$tmp_gene <- ifelse(dat[["RNA"]]@data[genelist[i],]>0,"Y","N")
        if(all(dat$tmp_gene=="N")){
           res.list[[i]] <- rep(0,length=length(table(dat@meta.data[,group])))
        }else{
           res.list[[i]] <- apply(table(dat$tmp_gene,dat@meta.data[,group]),2,function(x){x/sum(x)})[2,] 
        }
        
        names(res.list)[i] <- genelist[i]
    }
    return(res.list)
}


percent_feature(Tdat,c("GZMA","GZMB","GZMK","GZMH"),group="sampleID")

# CTN data
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
percent_feature(dat,c("GZMA","GZMB","GZMK","GZMH"),group="sample")

Youngdat <- matrix(unlist(percent_feature(Tdat,c("GZMA","GZMB","GZMK","GZMH"),group="sampleID")),ncol=4)
CTNdat <- matrix(unlist(percent_feature(dat,c("GZMA","GZMB","GZMK","GZMH"),group="sample")),ncol=4)

pdat <- data.frame(rbind(Youngdat,CTNdat),stringsAsFactors=F)
colnames(pdat) <- c("GZMA","GZMB","GZMK","GZMH")
pdat$Sample <- c(names(table(Tdat$sampleID)),names(table(dat$sample)))

pdat <- reshape2::melt(pdat)
colnames(pdat) <- c("Sample","Gene","Proportion")
pdat$Group <- "Young"
pdat$Group[grep("^XWS",pdat$Sample)] <- "CTN"
pdat$Group[grep("^SC",pdat$Sample)] <- "SC"
pdat$Group[grep("^CT",pdat$Sample)] <- "Elder"

pdat$Group <- factor(pdat$Group,levels=c("Young","Elder","CTN","SC"))

library(ggplot2)
pdf("/public/workspace/lily/CTN/version_3_20/Response/Young_CTN_Tcell_gene_proportion.pdf",useDingbats=F)
ggplot(pdat,aes(x=Gene,y=Proportion,fill=Group))+
	geom_bar(fun="median",position=position_dodge(0.8),stat="summary",width=0.5,size=1)+
	stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black",width = 0.2,position = position_dodge(0.8))+
	geom_point(data = pdat, aes(x=Gene,y = Proportion,fill=Group),position = position_dodge(0.8),size = 3, shape = 21)+
	scale_fill_manual(values = c("#d4c78c","#5ec6f2","#003468","#4b1702"))+ theme_bw() +
	theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
dev.off()








# 3. T cell subset proportion 
# use Average expression of GZMK to identify 
# use Average GZMK > 1
library(Seurat)

dat <- readRDS("~/CTN/version_3_20/Response/GSE158055/All_Control.RDS")
subdat <- subset(dat,cells=which(dat$Sample.type%in%c("fresh PBMC","frozen PBMC") & dat$Age < 50)) # do not use FACS
subdat$CellType <- "unknow"
subdat$CellType[which(subdat$majorType%in%c("B","Plasma"))] <- "Bcell"
subdat$CellType[which(subdat$majorType%in%c("CD4","CD8"))] <- "Tcell"
subdat$CellType[which(subdat$majorType%in%c("Mono"))] <- "CD14+ Monocytes"
subdat$CellType[which(subdat$celltype%in%c("Mono_c5-CD16"))] <- "CD16+ Monocytes"
subdat$CellType[which(subdat$majorType%in%c("DC"))] <- "DC"
subdat$CellType[which(subdat$majorType%in%c("NK"))] <- "NK"
Tdat <- subset(subdat,cells=which(subdat$CellType=="Tcell"))

# new celltype refine
# FeaturePlot(Tdat,features=c("GZMK","CD8A","NKG7","GZMH","CD4","CD27","IL7R"),label=T) 
# ZNF683 zeming zhang recognized as CD8Tm , but show NKG7 and GZMH 
Tdat$celltype.refine <- "unknow"
Tdat$celltype.refine[which(Tdat$celltype%in%c("T_CD8_c02-GPR183","T_CD8_c03-GZMK",
	"T_CD8_c06-TNF","T_CD8_c09-SLC4A10","T_CD8_c10-MKI67-GZMK","T_CD8_c05-ZNF683"))] <- "CD8Tem"

YCD8TEM <- apply(table(Tdat$celltype.refine,Tdat$sampleID),2,function(x){x/sum(x)})[1,]
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
CD8TEM <- apply(table(dat$sample,dat$celltype.refine),1,function(x){x/sum(x)})[5,]



pdf("/public/workspace/lily/CTN/version_3_20/Response/Young_CTN_CD8Tem_proportion.pdf",useDingbats=F)
boxplot(YCD8TEM,CD8TEM[1:4],CD8TEM[12:14],CD8TEM[5:11])
stripchart(list(YCD8TEM,CD8TEM[1:4],CD8TEM[12:14],CD8TEM[5:11]),vertical = T,method = "jitter",cex = 0.8,
           pch = 20,col = "red",add=T)
dev.off()













# genes <- c('ABL1','ADA','ADAM17','ADAM8','AGER','AIF1','AIRE','AKT1','ANXA1','AP3B1','AP3D1','APBB1IP','ARG1','ARG2','ATG5','ATP7A','AZI2','B2M',
# 	'BAD','BATF','BAX','BCL10','BCL11B','BCL2','BCL3','BCL6','BMP4','BTLA','BTN2A2','BTN3A1','CAMK4','CARD11','CASP3','CASP8','CAV1','CBFB',
# 	'CBLB','CCDC88B','CCL19','CCL2','CCL21','CCL5','CCND3','CCR2','CCR6','CCR7','CCR9','CD151','CD160','CD1C','CD1D','CD2','CD209','CD24',
# 	'CD27','CD274','CD276','CD28','CD300A','CD3D','CD3E','CD3G','CD4','CD40LG','CD44','CD46','CD47','CD5','CD55','CD6','CD7','CD70','CD74',
# 	'CD80','CD81','CD83','CD86','CD8A','CD8B','CDC42','CDH26','CEACAM1','CEBPB','CGAS','CHD7','CLC','CLEC4A','CLEC4G','CLEC7A','CLECL1',
# 	'CLPTM1','CORO1A','CR1','CRTAM','CSK','CTLA4','CTPS1','CTSL','CXADR','CYP26B1','CYRIB','DDOST','DLG1','DLG5','DLL4','DNAJA3','DOCK2',
# 	'DOCK8','DPP4','DROSHA','DTX1','DUSP10','DUSP3','EBI3','EFNB1','EFNB3','EGR1','EGR3','EIF2AK4','ELF4','ENTPD7','EOMES','EPO','ERBB2',
# 	'F2RL1','FADD','FANCA','FANCD2','FCER1G','FCGR2B','FGL1','FGL2','FKBP1A','FKBP1B','FLOT2','FOXJ1','FOXN1','FOXP1','FOXP3','FUT7','FYN',
# 	'FZD5','FZD7','GATA3','GBA','GLI2','GLI3','GLMN','GNRH1','GPAM','GPNMB','GPR183','GPR89A','GPR89B','GRAP2','GRB2','GSN','HAVCR2','HES1',
# 	'HFE','HHLA2','HLA-A','HLA-DMB','HLA-DOA','HLA-DPA1','HLA-DPB1','HLA-DRA','HLA-DRB1','HLA-DRB3','HLA-E','HLA-G','HLX','HMGB1','HSPD1',
# 	'HSPH1','ICAM1','ICOS','ICOSLG','IDO1','IFNA1','IFNA10','IFNA13','IFNA14','IFNA16','IFNA17','IFNA2','IFNA21','IFNA4','IFNA5','IFNA6','IFNA7',
# 	'IFNA8','IFNB1','IFNE','IFNG','IFNK','IFNL1','IFNW1','IGF1','IGF2','IGFBP2','IHH','IL10','IL12A','IL12B','IL12RB1','IL15','IL18','IL18R1',
# 	'IL1A','IL1B','IL1RL2','IL2','IL20RB','IL21','IL23A','IL23R','IL27','IL27RA','IL2RA','IL36B','IL4','IL4R','IL6','IL6R','IL6ST','IL7','IL7R',
# 	'INS','IRF1','IRF4','ITGAL','ITK','ITPKB','JAG2','JAK3','JAML','JMJD6','KAT2A','KIF13B','KIT','KLRC1','KLRC4-KLRK1','KLRK1','LAG3','LAPTM5',
# 	'LAT','LAX1','LCK','LCP1','LEF1','LEP','LEPR','LFNG','LGALS1','LGALS3','LGALS7B','LGALS9','LGALS9B','LGALS9C','LIG4','LILRB1','LILRB2',
# 	'LILRB4','LMBR1L','LMO1','LOXL3','LRRC32','LY9','LYN','MAD1L1','MAFB','MALT1','MAP3K8','MAPK8IP1','MARCHF7','MDK','METTL3','MICA','MICB',
# 	'MIR181C','MIR21','MIR27A','MIR30B','MR1','MSN','MTOR','MYB','MYH9','NCAPH2','NCK1','NCK2','NCKAP1L','NCSTN','NDFIP1','NFATC2','NFKBID',
# 	'NFKBIZ','NHEJ1','NKAP','NKX2-3','NLRC3','NLRP3','NOD2','NRARP','PAG1','PAK1','PAK2','PAK3','PATZ1','PAWR','PAX1','PCK1','PDCD1','PDCD1LG2',
# 	'PDPK1','PELI1','PIK3CA','PIK3CD','PIK3CG','PIK3R1','PIK3R6','PKNOX1','PLA2G2D','PLA2G2F','PNP','PPP3CA','PPP3CB','PRDM1','PRDX2','PRELID1',
# 	'PREX1','PRKAR1A','PRKCQ','PRKCZ','PRKDC','PRNP','PRR7','PSEN1','PSMB10','PSMB11','PTGER4','PTPN11','PTPN2','PTPN22','PTPN6','PTPRC',
# 	'PYCARD','RAB27A','RAB29','RABL3','RAC1','RAC2','RAG1','RAG2','RARA','RASAL3','RASGRP1','RC3H1','RC3H2','RELB','RHOH','RIPK2','RIPK3',
# 	'RIPOR2','RORA','RORC','RPL22','RPS3','RPS6','RSAD2','RUNX1','RUNX3','SART1','SASH3','SCGB1A1','SCRIB','SDC4','SELENOK','SEMA4A','SFTPD',
# 	'SH3RF1','SHH','SIRPA','SIRPB1','SIRPG','SIT1','SLA2','SLAMF6','SLC11A1','SLC46A2','SLC7A1','SMAD3','SMAD7','SOCS1','SOCS5','SOCS6','SOD1',
# 	'SOS1','SOS2','SOX12','SOX13','SOX4','SP3','SPINK5','SPN','SPTA1','SRC','SRF','STAT3','STAT5B','STK11','STOML2','SYK','TARM1','TBX21','TCF7',
# 	'TCIRG1','TESPA1','TFRC','TGFBR2','THEMIS','THY1','TIGIT','TMEM131L','TMEM98','TMIGD2','TNFAIP8L2','TNFRSF13C','TNFRSF14','TNFRSF18',
# 	'TNFRSF1B','TNFRSF21','TNFRSF4','TNFSF11','TNFSF13B','TNFSF14','TNFSF18','TNFSF4','TNFSF8','TNFSF9','TOX','TP53','TRAF6','TREML2','TREX1',
# 	'TSC1','TWSG1','VAV1','VCAM1','VNN1','VSIG4','VSIR','VTCN1','WAS','WDFY4','WNT1','WNT4','XBP1','XCL1','YES1','ZAP70','ZBTB1','ZBTB16',
# 	'ZBTB7B','ZC3H12A','ZC3H8','ZFP36L1','ZFP36L2','ZFPM1','ZMIZ1','ZNF683','ZP3','ZP4')

subdatT <- AddModuleScore(subdatT,features=list(genes),name="Tcell_act")














#################################################################################################################################################
# 2023-2-14
# Double lets Finder result check 
# data/XWS，PNAS are samples before Doublets 
# CTN data accroding to prepare.R line 336
# > table(XWS1.f$DF.classifications_0.25_0.24_773)

# Doublet Singlet
#     773   10601

# > table(XWS2.f$DF.classifications_0.25_0.005_917)

# Doublet Singlet
#     917   12955
# > table(XWS3.f$DF.classifications_0.25_0.005_920)

# Doublet Singlet
#     920   12506
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS") # after Double and inte
tmp <- readRDS("/public/workspace/lily/CTN/version_3_20/data/PNAS.RDS") # before Double

pdat <- data.frame(singlets=table(dat$sample),raw=c(table(tmp$sample),11374,13872,13426),stringsAsFactors=F)
colnames(pdat)<- c("Sample","Singlets","Raw")
pdat$Doublets <- pdat$Raw -pdat$Singlets

pdat <- reshape2::melt(pdat[,c(1,2,4)])
colnames(pdat)<- c("Sample","Type","Num")
pdat$Sample <- factor(pdat$Sample,levels=c("CT2","CT3","CT4","CT5","XWS1","XWS2","XWS3","SC1","SC2","SC3","SC4","SC5","SC6","SC7"))
library(ggplot2)

pdf("/public/workspace/lily/CTN/version_3_20/Response/CTN_Doublets.pdf",useDingbats=F)

ggplot(pdat,aes(x=Sample,y=Num,fill=Type))+geom_bar(stat = "identity",position = "fill")+
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
dev.off()






# calculate inflammation signature 

library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS") # after Double and inte
gene <- c("IL1B","IL6","IFNB1","CCL2","CCL5","CXCL1","CXCL8","TGFB1","MMP1","GDF15","PTGES2") # https://www.nature.com/articles/s41577-021-00646-4#Sec4
dat <- AddModuleScore(dat,assay="RNA",features=list(gene),name="inflammation")

pdf("/public/workspace/lily/CTN/version_3_20/Response/CTN_inflammation.pdf",useDingbats=F)
boxplot(inflammation1~group,data=dat@meta.data,outline=F,names=c("CTN","Elder","SC"))
dev.off()



# anti-inflammation
# not use not ok
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS") # after Double and inte

gene <- c("IL1RA","IL10","IL4","IL13","MERTK","FCER1A")
dat <- AddModuleScore(dat,assay="RNA",features=list(gene),name="anti_inflammation")

boxplot(anti_inflammation1~group,data=dat@meta.data,outline=F,names=c("CTN","Elder","SC"))














#################################################################################################################################################
# 2023-2-14
library(Seurat)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
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


dat$FOXP3_act <- as.numeric(as.character(scale(reg.m$FOXP3,center=F,scale=T)))
dat$LEF1_act <- as.numeric(as.character(scale(reg.m$LEF1,center=F,scale=T)))
dat$TCF7_act <- as.numeric(as.character(scale(reg.m$TCF7,center=F,scale=T)))
# scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


pdf("/public/workspace/lily/CTN/version_3_20/Response/CTN_TF_act_gene.pdf",useDingbats=F,height=10,width=10)
DimPlot(dat,cells.highlight=colnames(dat)[which(dat$seurat_clusters==12)],raster=T)
FeaturePlot(dat,features=c("FOXP3_act","LEF1_act","TCF7_act"),raster=T) & scale_colour_gradientn(colours=c("grey","lightgrey","orange","red"))
FeaturePlot(dat,features=c("FOXP3","LEF1","TCF7","IFNG"),raster=T)
FeaturePlot(dat,features=c("CD4","CD8B","CD8A","IFNG"),raster=T)

dev.off()




# 2023-02-16
# classify CD4T and CD8T
dat$tmp <- dat$celltype.refine
dat$tmp[grep("CD4+",dat$tmp)] <- "CD4 T"
dat$tmp[grep("CD8+",dat$tmp)] <- "CD8 T"
dat$tmp[grep("Tregs",dat$tmp)] <- "CD4 T"
dat$tmp[grep("γδT",dat$tmp)] <- "Unclassify"
pdf("/public/workspace/lily/CTN/version_3_20/Response/CTN_T_CD4_CD8.pdf",useDingbats=F)
DimPlot(dat,group.by="tmp",cols=c("#b5c327","#109dc0","grey"),raster=T)
dev.off()

