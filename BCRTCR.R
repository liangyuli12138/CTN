

# 2021-9-18
# analysis TCR and BCR 
#============================================================================================================================================
library(Seurat)

# dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN_inte.RDS")
Tdat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
Bdat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")

# Sample1 
# test
# tmp.1 <- subset(dat,cells=which(dat$orig.ident=="XWS1"))
# tmp.1$cells <- sapply(strsplit(colnames(tmp.1),"_"),function(x){x[[1]]})
# tmp <- read.table("/public/workspace/lily/CTN/version_3_20/TRUST/XWS1/XWS1_report.tsv",header=F)
# tmp$cells <- sapply(strsplit(as.vector(tmp$V9),"_"),function(x){x[[1]]})
# tmp.t <- tmp.f[grep("^T",tmp$V8),]

# tmp.1$group <- "others"
# tmp.1$group[which(tmp.1$cells %in% tmp.t$cells)] <- "T "

# T cell
#============================================================================================================================================
S1t <- subset(Tdat,cells=which(Tdat$orig.ident=="XWS1"))
S1t$cells <- sapply(strsplit(colnames(S1t),"_"),function(x){x[[1]]})
tmp <- read.table("/public/workspace/lily/CTN/version_3_20/TRUST/XWS1/XWS1_report.tsv",header=F)
tmp$cells <- sapply(strsplit(as.vector(tmp$V9),"_"),function(x){x[[1]]})


tmp.f <- tmp[which(tmp$V1>1),] # expand  cells
tmp.f.t <- tmp.f[grep("^T",tmp.f$V8),] # expanded T cells 
# this is calculate 
S1t@meta.data[which(S1t$cells %in% tmp.f.t$cells),"celltype.refine"] 

S1.t.res <- table(S1t@meta.data[which(S1t$cells %in% tmp.f.t$cells),"celltype.refine"] )


# B cells 
#============================================================================================================================================
S1b <- subset(Bdat,cells=which(Bdat$orig.ident=="XWS1")) # all
S1b$cells <- sapply(strsplit(colnames(S1b),"_"),function(x){x[[1]]})
tmp.f.b <- tmp.f[grep("^I",tmp.f$V8),] # expanded T cells 
# this is calculate 
S1b@meta.data[which(S1b$cells %in% tmp.f.b$cells),"celltype.new"] 

S1.b.res <- table(S1b@meta.data[which(S1b$cells %in% tmp.f.b$cells),"celltype.new"])





# sample2 
#############################################################################################################################################
S2t <- subset(Tdat,cells=which(Tdat$orig.ident=="XWS2"))
S2t$cells <- sapply(strsplit(colnames(S2t),"_"),function(x){x[[1]]})
tmp <- read.table("/public/workspace/lily/CTN/version_3_20/TRUST/XWS2/XWS2_report.tsv",header=F)
tmp$cells <- sapply(strsplit(as.vector(tmp$V9),"_"),function(x){x[[1]]})


tmp.f <- tmp[which(tmp$V1>1),] # expand  cells
tmp.f.t <- tmp.f[grep("^T",tmp.f$V8),] # expanded T cells 
# this is calculate 
S2t@meta.data[which(S2t$cells %in% tmp.f.t$cells),"celltype.refine"] 

S2.t.res <- table(S2t@meta.data[which(S2t$cells %in% tmp.f.t$cells),"celltype.refine"])


# B cells 
#============================================================================================================================================
S2b <- subset(Bdat,cells=which(Bdat$orig.ident=="XWS2")) # all
S2b$cells <- sapply(strsplit(colnames(S2b),"_"),function(x){x[[1]]})
tmp.f.b <- tmp.f[grep("^I",tmp.f$V8),] # expanded T cells 
# this is calculate 
S2b@meta.data[which(S2b$cells %in% tmp.f.b$cells),"celltype.new"] 

S2.b.res <- table(S2b@meta.data[which(S2b$cells %in% tmp.f.b$cells),"celltype.new"])




# sample3 
#############################################################################################################################################
S3t <- subset(Tdat,cells=which(Tdat$orig.ident=="XWS3"))
S3t$cells <- sapply(strsplit(colnames(S3t),"_"),function(x){x[[1]]})
tmp <- read.table("/public/workspace/lily/CTN/version_3_20/TRUST/XWS3/XWS3_report.tsv",header=F)
tmp$cells <- sapply(strsplit(as.vector(tmp$V9),"_"),function(x){x[[1]]})


tmp.f <- tmp[which(tmp$V1>1),] # expand  cells
tmp.f.t <- tmp.f[grep("^T",tmp.f$V8),] # expanded T cells 
# this is calculate 
S3t@meta.data[which(S3t$cells %in% tmp.f.t$cells),"celltype.refine"] 

S3.t.res <- table(S3t@meta.data[which(S3t$cells %in% tmp.f.t$cells),"celltype.refine"])



# B cells 
#============================================================================================================================================
S3b <- subset(Bdat,cells=which(Bdat$orig.ident=="XWS3")) # all
S3b$cells <- sapply(strsplit(colnames(S3b),"_"),function(x){x[[1]]})
tmp.f.b <- tmp.f[grep("^I",tmp.f$V8),] # expanded T cells 
# this is calculate 
S3b@meta.data[which(S3b$cells %in% tmp.f.b$cells),"celltype.new"] 


S3.b.res <- table(S3b@meta.data[which(S3b$cells %in% tmp.f.b$cells),"celltype.new"])







#============================================================================================================================================
# plot result 

# T cell 
S1.t.res/sum(S1.t.res)
S2.t.res/sum(S2.t.res)
S3.t.res/sum(S3.t.res)
T.res <- data.frame(S1=c(0.019607843,0.049019608,0.921568627,0,0.009803922,0),
		S2= c(0,0,0.94666667,0,0.05333333,0),
		S3= c(0.024,0.040,0.904,0.008,0.008,0.016)
	)
rownames(T.res) <- c("CD4+ Tem","CD8+ Tcm","CD8+ Tem","Naive T","Tregs","Unclassify")
T.res$type <- rownames(T.res)
library(ggplot2)
T.dat <- reshape2::melt(T.res,id.var="type")

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/TRUSR4_TCR.pdf",useDingbats=F)
cols.t <- c("#4a8594","#faae40","#598c14","#b5c327","#3a4958","#dbe0e3")
ggplot(T.dat,aes(x=variable,y=value,fill=type))+ geom_bar(stat="identity",position="stack") +
	scale_fill_manual(values=cols.t)+theme_bw()

dev.off()


# 2023-1-3
# new cell type 
# T cell 
S1.t.res/sum(S1.t.res)
S2.t.res/sum(S2.t.res)
S3.t.res/sum(S3.t.res)
T.res <- data.frame(S1=c(0.147058824,0,0.019607843,0.715686275,0.049019608,0.009803922,0,0.058823529),
		S2= c(0.18666667,0,0,0.74666667,0,0.05333333,0,0.01333333),
		S3= c(0.224,0.008,0.024,0.584,0.040,0.008,0.016,0.096)
	)
rownames(T.res) <- c("CD4+ CTL","CD4+ Tcm","CD4+ Tm","CD8+ Tem","CD8+ Tn","Tregs","Unclassify","γδT")
T.res$type <- rownames(T.res)
library(ggplot2)
T.dat <- reshape2::melt(T.res,id.var="type")
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/TRUSR4_TCR_new.pdf",useDingbats=F)
cols.t <- c("#00b2a9","#0066a1","#a626aa","#00ad45","#70b29c","#3a4958","#bfbfbf","#aea400") # 2022-12-31

ggplot(T.dat,aes(x=variable,y=value,fill=type))+ geom_bar(stat="identity",position="stack") +
	scale_fill_manual(values=cols.t)+theme_bw()

dev.off()





S1.b.res/sum(S1.b.res)
S2.b.res/sum(S2.b.res)
S3.b.res/sum(S3.b.res)

B.res <- data.frame(S1=c(0,0,1,0),
		S2= c(0.13157895,0.10526316,0.73684211,0.02631579),
		S3= c( 0.44,0.16,0.36,0.04)
	)
rownames(B.res) <- c("Memory B cell","Naive B cell","Plasma","unclassify")
B.res$type <- rownames(B.res)
library(ggplot2)
B.dat <- reshape2::melt(B.res,id.var="type")

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/TRUSR4_BCR.pdf",useDingbats=F)
B.dat$type <- factor(B.dat$type,levels=c("Naive B cell","Memory B cell","Plasma","unclassify"))
# naive memory switch-memory plasma
cols <- c("#009f4d","#f48924","#f85a40","#caccd1")
ggplot(B.dat,aes(x=variable,y=value,fill=type))+ geom_bar(stat="identity",position="stack") +
	scale_fill_manual(values=cols)+theme_bw()

dev.off()





























