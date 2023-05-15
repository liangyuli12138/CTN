
# 2021-11-25 
# analysis bluk RNA-seq 
#===============================================================================================================================================
# 0. prepare data 
# dat <- data.table::fread("/public/workspace/lily/CTN/version_3_20/data/genes_fpkm_expression.txt",header=T,sep="\t")
# dat.res <- dat[,c(6,14:28)]
# colnames(dat.res) <- sapply(strsplit(colnames(dat.res),"\\."),function(x){x[[length(x)]]})
# tmp.res <- aggregate(.~gene_name,data=dat.res,FUN=median)
# tmp.res <- data.frame(tmp.res)
# rownames(tmp.res) <- tmp.res$gene_name
# tmp.res$gene_name <- NULL

# saveRDS(tmp.res,file="/public/workspace/lily/CTN/version_3_20/data/Mice_Bulk_RNAseq.RDS")



# transate human hippo gene into mouse 
tmp <- read.table("/public/workspace/lily/REF/Mouse_Human_gene.txt",sep="\t",header=T)
gene.m <- as.vector(unique(tmp$Gene.name[which(tmp$Human.gene.name%in%gene)]))
mod.generate(gene.m,"KEGG_hippo.m",out="/public/workspace/lily/MOD_file/KEGG_hippo.m.mod") # make a mod file 
mod.generate(gene.m,"Reactome_hippo.m",out="/public/workspace/lily/MOD_file/Reactome_hippo.m.mod") # make a mod file 
mod.generate(gene.m,"Tcell_act.m",out="/public/workspace/lily/MOD_file/Tcell_act.m.mod") # make a mod file from MSGDB GO 



# 1. calculate some reuslt 
#==============================================================================================================================================
dat <-readRDS("/public/workspace/lily/CTN/version_3_20/data/Mice_Bulk_RNAseq.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("Reactome_hippo.m","KEGG_hippo.m","Tcell_act.m"),"/public/workspace/lily/MOD_file/",permN=0)

mod <- mod.analyze2(as.matrix(dat[,c(1:5,11:15)]),c("Reactome_hippo.m","KEGG_hippo.m","Tcell_act.m"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$group <- sapply(strsplit(rownames(mod),"_"),function(x){x[[1]]})



boxplot(KEGG_hippo.m_norm~group,data=mod,FUN=median,main=" Long vs. Age (KEGG hippo)",ylab="KEGG Hippo signature score")
legend("topleft",legend=paste0("P=",round(wilcox.test(KEGG_hippo.m_norm~group,data=mod)$p.value,2)),bty="n")

boxplot(Reactome_hippo.m_norm~group,data=mod,FUN=median,main=" Long vs. Age (Reactome hippo)",ylab="Reactome Hippo signature score")
legend("topleft",legend=paste0("P=",round(wilcox.test(Reactome_hippo.m_norm~group,data=mod)$p.value,2)),bty="n")

boxplot(Tcell_act.m_norm~group,data=mod,FUN=median,main=" Long vs. Age (Tact. sig.)",ylab="T cells activate signature score")
legend("topleft",legend=paste0("P=",round(wilcox.test(Tcell_act.m_norm~group,data=mod)$p.value,2)),bty="n")


apply(mod[,4:6],2,function(x){c(median(x[1:5]),median(x[6:10]),median(x[11:15]))})
apply(dat["Gzma",],1,function(x){c(median(x[1:5]),median(x[6:10]),median(x[11:15]))})





boxplot(mod[6:10,6],mod[1:5,6],mod[11:15,6],names=c("infant","aged","long"),main="T cell activate signature score")

boxplot(as.numeric(dat["Il10",6:10]),as.numeric(dat["Il10",1:5]),as.numeric(dat["Il10",11:15]),names=c("infant","aged","long"),main="Il10 expression")

boxplot(as.numeric(dat["Ifng",6:10]),as.numeric(dat["Ifng",1:5]),as.numeric(dat["Ifng",11:15]),names=c("infant","aged","long"),main="Ifng expression")

boxplot(as.numeric(dat["Gzma",6:10]),as.numeric(dat["Gzma",1:5]),as.numeric(dat["Gzma",11:15]),names=c("infant","aged","long"),main="Gzma expression")

boxplot(as.numeric(dat["Gzmb",6:10]),as.numeric(dat["Gzmb",1:5]),as.numeric(dat["Gzmb",11:15]),names=c("infant","aged","long"),main="Gzmb expression")

boxplot(as.numeric(dat["Gzmk",6:10]),as.numeric(dat["Gzmk",1:5]),as.numeric(dat["Gzmk",11:15]),names=c("infant","aged","long"),main="Gzmk expression")






# 2. DEG for aged and long 
tmp.res <- data.frame(t(apply(dat,1,function(x){
 	c(
 		mean(x[1:5]),
 		mean(x[11:15]),
 		mean(x[11:15])/mean(x[1:5]),
 		wilcox.test(x[1:5],x[11:15])$p.value
 	)
 })))
colnames(tmp.res) <- c("Exp_aged","Exp_log","FC","pvalue")
tmp.res$log2FC <- log2(tmp.res$FC)

tmp.res <- tmp.res[which(tmp.res$Exp_aged>0&tmp.res$Exp_log>0),]

res <- tmp.res[which(tmp.res$FC>2&tmp.res$pvalue<0.05),]

# no Gm gene (Pesudogene) and rik (lncRNA)
tmp.gene <- grep("^Gm|Rik$",rownames(res),value=T,invert=T)

write.table(tmp.gene,file="~/tmp.txt",row.names=F,col.names=F,quote=F)




# plot volcano
tmp.res$change <- "unkonw"
tmp.res$change[which(tmp.res$FC<0.5&tmp.res$pvalue<0.05)] <- "Dn"
tmp.res$change[which(tmp.res$FC>2&tmp.res$pvalue<0.05)] <- "Up"

tmp.res$log2P <- -log10(tmp.res$pvalue)

tmp.res$label <- NA
tmp.res$label[which(rownames(tmp.res)%in%c("Gzma","Gzmb"))] <- rownames(tmp.res)[which(rownames(tmp.res)%in%c("Gzma","Gzmb"))]

tmp.res <- tmp.res[which(tmp.res$pvalue<0.05),]
library(ggplot2)
library(ggrepel)

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Mice_Bulk_RNAseq_aged_long_DEG.pdf",useDingbats=F)
ggplot(data=tmp.res,aes(x=log2FC,y=log2P,color=change))+
	geom_point(alpha=0.4, size=3.5) + scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
	geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
	geom_hline(yintercept = 5,lty=4,col="black",lwd=0.8) +
	# 坐标轴
	labs(x="Average log(fold change)",
	   y="-log2 (p-value)")+
	theme_bw()+
	# 图例
	theme(plot.title = element_text(hjust = 0.5), 
	    legend.position="right", 
	    legend.title = element_blank()) +
	geom_label(aes(label=label))

dev.off()











# 2021-12-22 
# plot a heamtap 

dat <-readRDS("/public/workspace/lily/CTN/version_3_20/data/Mice_Bulk_RNAseq.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("Reactome_hippo.m","KEGG_hippo.m","Tcell_act.m"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- data.frame(mod)

mod$Ifng <- as.numeric(dat["Ifng",])
mod$Il10 <- as.numeric(dat["Il10",])
mod$Cdkn1a <- as.numeric(dat["Cdkn1a",])
mod$Gzma <- as.numeric(dat["Gzma",])
mod$Gzmk <- as.numeric(dat["Gzmk",])
mod$Gzmb <- as.numeric(dat["Gzmb",])
# get plot data 
plot.dat <- mod[,c("Tcell_act.m_norm","Ifng","Il10","Cdkn1a","Gzmb","Gzma")]

# heatmap is not OK
# pheatmap::pheatmap(plot.dat[-c(6:10),],scale="column",cluster_rows=F,color=colorRampPalette(c('steelblue','white','red'))(100))


# use boxplot 
plot.dat$group <- sapply(strsplit(rownames(plot.dat),"_"),function(x){x[[1]]})
tmp.res <- reshape2::melt(plot.dat,id="group")
colnames(tmp.res) <- c("group","type","value")
tmp.res$group <- factor(tmp.res$group,levels=c("infan","aged","long"))

library(ggplot2)
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Mice_Bulk_RNAseq_factor.pdf",useDingbats=F)
ggplot(tmp.res,aes(x=group,y=value))+geom_boxplot(aes(group=group,color=group)) + geom_point()+facet_wrap( ~ type, ncol=2,scales="free")
dev.off()

















#====================================================================================================================================================
# analysis some signature 
# 2022-1-4
# need to trans into mouse 
# /public/workspace/lily/CTN/version_3_20/Signature_MM/
#====================================================================================================================================================
library(qusage)
tmp.sig <- read.gmt("/public/workspace/lily/CTN/version_3_20/Signature.gmt")
trans.gene <- read.table("~/REF/Mouse_Human_gene.txt",sep="\t",header=T)
tmp.sig <- lapply(tmp.sig,function(x){as.vector(unique(trans.gene$Gene.name[which(trans.gene$Human.gene.name%in%x)]))})

# write out gmt 
conn=file("/public/workspace/lily/CTN/version_3_20/Signature_MM.gmt",'w')
for(i in 1:length(tmp.sig)){
	cat(names(tmp.sig)[i],names(tmp.sig)[i],tmp.sig[[i]],file=conn,sep='\t')	
	cat("\n",file=conn)
}


source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
for(i in 1:length(tmp.sig)){
	mod.generate(tmp.sig[[i]],names(tmp.sig)[[i]],out=paste0("/public/workspace/lily/CTN/version_3_20/Signature_MM/",names(tmp.sig)[[i]],".mod")) # make a mod file 
}


#====================================================================================================================================================
dat <-readRDS("/public/workspace/lily/CTN/version_3_20/data/Mice_Bulk_RNAseq.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
modfile <- gsub("\\.mod$","",dir("/public/workspace/lily/CTN/version_3_20/Signature_MM/"))
mod <- mod.analyze2(as.matrix(dat),modfile,"/public/workspace/lily/CTN/version_3_20/Signature_MM/",permN=0)
mod <- as.data.frame(mod)
mod$group <- sapply(strsplit(rownames(mod),"_"),function(x){x[[1]]})



# plot result 
library(ggplot2)
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Mice_Bulk_RNAseq_Signature.pdf",useDingbats=F)

for(i in 9:16){
	boxplot(mod[,i]~group,data=mod,main=gsub("_norm$","",colnames(mod)[i]),ylim=c(0,1))
	legend("topright",legend=(paste0("p=",kruskal.test(mod[,i]~group,data=mod)$p.value)))
	beeswarm::beeswarm(mod[,i]~group,data=mod,col = 4, pch = 16,add = TRUE)
}


dev.off()











#============================================================================================================================================
# do GSEA for some pathway 
# 2022-1-5
library(clusterProfiler)
library(enrichplot)
dat <-readRDS("/public/workspace/lily/CTN/version_3_20/data/Mice_Bulk_RNAseq.RDS")
# DEG for aged and long 
tmp.res <- data.frame(t(apply(dat,1,function(x){
 	c(
 		mean(x[1:5]),
 		mean(x[11:15]),
 		mean(x[11:15])/mean(x[1:5]),
 		wilcox.test(x[1:5],x[11:15])$p.value
 	)
 })))
colnames(tmp.res) <- c("Exp_aged","Exp_log","FC","pvalue")
tmp.res$log2FC <- log2(tmp.res$FC)
tmp.res <- tmp.res[which(tmp.res$Exp_aged>0&tmp.res$Exp_log>0),]


#==============================================================================================================================================
tmp.res <- tmp.res[order(tmp.res$FC,decreasing=T),]

genelist <- as.numeric(tmp.res$FC)
names(genelist)<- rownames(tmp.res)

tmp.sig <- read.gmt("/public/workspace/lily/CTN/version_3_20/Signature_MM.gmt")
gsea <- GSEA(genelist,TERM2GENE = tmp.sig,pvalueCutoff=1) #GSEA分析

gseaplot2(gsea,2)








#=============================================================================================================================================
# do CD8+Tem analysis for CTN vs Elder
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
sub.dat <- subset(dat,cells=which(dat$celltype.refine=="CD8+ Tem"))
DefaultAssay(sub.dat) <- "RNA"
sub.dat@active.ident <- factor(sub.dat$group)

CTN.gene <- FindMarkers(sub.dat,ident.1="CTN",ident.2="Elder",assay="RNA",group.by="group",logfc.threshold = 0,min.pct = 0)
SC.gene <- FindMarkers(sub.dat,ident.1="SC",ident.2="Elder",assay="RNA",group.by="group")

CTN.gene.f <- CTN.gene[which(CTN.gene$p_val_adj<0.05),]
SC.gene.f <- SC.gene[which(SC.gene$p_val_adj<0.05),]

# for CTN.gene

genelist <- as.numeric(CTN.gene.f$avg_logFC)
names(genelist)<- rownames(CTN.gene.f)
genelist <- sort(genelist,decreasing= T)

tmp.sig <- read.gmt("/public/workspace/lily/CTN/version_3_20/Signature.gmt")
gsea <- clusterProfiler::GSEA(genelist,TERM2GENE = tmp.sig,pvalueCutoff=1) #GSEA分析

gseaplot2(gsea,)




























#=============================================================================================================================================
# 2022-1-5
# calculate metabolism pathways
# trans data into Human genes

dat <-readRDS("/public/workspace/lily/CTN/version_3_20/data/Mice_Bulk_RNAseq.RDS")
trans.gene <- read.table("~/REF/Mouse_Human_gene.txt",sep="\t",header=T)
dat$gene <- rownames(dat)
dat.t <- merge(dat,trans.gene,by.x="gene",by.y="Gene.name")
dat.t$gene <- NULL
dat.res <- aggregate(.~Human.gene.name,data=dat.t,FUN=median)
dat.res <- dat.res[-1,]
rownames(dat.res) <- dat.res$Human.gene.name
dat.res$Human.gene.name <- NULL

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
filename <- gsub("\\.mod","",dir("/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism"))
mod.rs <- mod.analyze2(as.matrix(dat.res),filename,"/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/",permN=0)
mod.rs <- as.data.frame(mod.rs)

res1 <- data.frame(t(apply(mod.rs[,86:170],2,function(x){
	c(
		median(x[1:5]),
		median(x[6:10]),
		median(x[11:15])
	)})))

colnames(res1) <- c("aged","infant","long")
res1$diff <- res1$long - res1$aged
res1 <- res1[order(res1$diff,decreasing=T),]

#==============================================================================================================================================
# get scRNA metabolism data 
tmp <- read.table("/public/workspace/lily/CTN/version_3_20/metabolism/Tcell/KEGGpathway_activity_shuffle.txt",sep="\t",header=T)
tmp.a <- tmp
tmp.a[is.na(tmp.a)] <- 1

tmp.f <- tmp.a[,-grep("C14",colnames(tmp.a))]
tmp.f <- tmp.f[,c("C0","C1","C2","C4","C8","C3","C5","C6","C7","C9","C10","C11","C13","C12")]
tmp.res <- t(apply(tmp.f,1,function(x){
    c(mean(x[1:5]),mean(x[c(6:10)]),mean(x[c(11:13)]),mean(x[c(14)]))  
}))
colnames(tmp.res) <- c("SC.CTN","Elder","Uniform","CTN")
tmp.res <- data.frame(tmp.res)
tmp.res$diff <- tmp.res$SC.CTN - tmp.res$Elder
tmp.res <- tmp.res[order(tmp.res$diff,decreasing=T),]


# then got top 10 diff pathway 
intersect(gsub(" norm","",gsub("_"," ",head(rownames(res1),10))),head(rownames(tmp.res),10))
# and wirte result 
write.table(head(res1,10),file="/public/workspace/lily/CTN/version_3_20/rs_plot/Mice_Bulk_RNAseq_metabolism_top10.txt",quote=F,sep="\t")
write.table(head(tmp.res,10),file="/public/workspace/lily/CTN/version_3_20/rs_plot/Tcell_scRNA_metabolism_top10.txt",quote=F,sep="\t")





#=============================================================================================================================================
# now plot result for Mice data 
tmp.dat <- mod.rs[,c("Fatty_acid_degradation_norm","Oxidative_phosphorylation_norm")]
tmp.dat$group <- sapply(strsplit(rownames(tmp.dat),"_"),function(x){x[[1]]})

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Mice_Bulk_RNAseq_metabolism.pdf",useDingbats=F)

	boxplot(Fatty_acid_degradation_norm~group,data=tmp.dat,main="Fatty acid degradation",ylim=c(0,1))
	legend("topright",legend=(paste0("p=",kruskal.test(Fatty_acid_degradation_norm~group,data=tmp.dat)$p.value)))
	beeswarm::beeswarm(Fatty_acid_degradation_norm~group,data=tmp.dat,col = 4, pch = 16,add = TRUE)

	boxplot(Oxidative_phosphorylation_norm~group,data=tmp.dat,main="Oxidative phosphorylation",ylim=c(0,1))
	legend("topright",legend=(paste0("p=",kruskal.test(Oxidative_phosphorylation_norm~group,data=tmp.dat)$p.value)))
	beeswarm::beeswarm(Oxidative_phosphorylation_norm~group,data=tmp.dat,col = 4, pch = 16,add = TRUE)

dev.off()





























