
# 2022-1-4
# this program is used to analysis 
#===============================================================================================================================================
library(readxl)

tmp1 <- read_xlsx("/public/workspace/lily/CTN/version_3_20/Cytometry/Panel1_data.xlsx")

data1 <- data.frame(
	Treg = apply(tmp1,1,function(x){as.numeric(x[10])/as.numeric(x[5])}), # Tregs
	CD8Tem = apply(tmp1,1,function(x){as.numeric(x[17])/as.numeric(x[5])}), # CD8+ Tem 
	CD4Tem = apply(tmp1,1,function(x){as.numeric(x[13])/as.numeric(x[5])}), # CD4+ Tem
	CD8Tcm = apply(tmp1,1,function(x){as.numeric(x[16])/as.numeric(x[5])}), # CD8+ Tcm  
	naiveT = apply(tmp1,1,function(x){(as.numeric(x[14])+as.numeric(x[19]))/as.numeric(x[5])}), # naive T cells
	GZMBCD8T = apply(tmp1,1,function(x){as.numeric(x[18])/as.numeric(x[17])}), # GZMB+ CD8T cells percentage
	group=c(rep("CTN",5),rep("Elder",5),rep("Young",6))
	)


tmp.res <- reshape2::melt(aggregate(.~group,data=data1,FUN=median))
colnames(tmp.res)[3] <- "median"
tmp.res$sd <- as.numeric(reshape2::melt(aggregate(.~group,data=data1,FUN=sd))[,3])
# plot result 
library(ggplot2)
library(ggpubr)
tmp.res$group <- factor(tmp.res$group,levels=c("Young","Elder","CTN"))

pdf("/public/workspace/lily/CTN/version_3_20/Cytometry/Panel1_res.pdf",useDingbats=F)
ggplot(data=tmp.res,aes(x=group,y=median,group=group)) + geom_bar(stat = "identity",aes(fill=group)) +
	geom_errorbar(aes(ymin = median, ymax = median + sd), width = 0.2, position = position_dodge(0.9)) +
	theme_classic() +
	facet_wrap(.~variable,scale="free")+ ylab(label="% T cells") +
	stat_compare_means(comparisons=list(c("Elder","CTN")))
dev.off()



pdf("/public/workspace/lily/CTN/version_3_20/Cytometry/Panel1_res.pdf",useDingbats=F)
ggplot(data=tmp.res,aes(x=group,y=median,group=group)) + geom_bar(stat = "identity",aes(fill=group)) +
	geom_errorbar(aes(ymin = median, ymax = median + sd), width = 0.2, position = position_dodge(0.9)) +
	theme_classic() +
	facet_wrap(.~variable,scale="free")+ ylab(label="% T cells") +
	stat_compare_means(comparisons=list(c("Elder","CTN")))
dev.off()















#===========================================================================================================================================
library(readxl)

tmp2 <- read_xlsx("/public/workspace/lily/CTN/version_3_20/Cytometry/Panel2_data.xlsx")

data2 <- data.frame(
	Tcell = apply(tmp2,1,function(x){as.numeric(x[16])/as.numeric(x[3])}), # T cells
	Bcell = apply(tmp2,1,function(x){as.numeric(x[9])/as.numeric(x[3])}), # B cells
	MemBcell = apply(tmp2,1,function(x){as.numeric(x[10])/as.numeric(x[9])}), # Memery B cells 
	CD14Mono = apply(tmp2,1,function(x){as.numeric(x[4])/as.numeric(x[3])}), # CD14+ monocyte
	NK = apply(tmp2,1,function(x){as.numeric(x[12])/as.numeric(x[3])}), # NK cells
	non_MemB = apply(tmp2,1,function(x){(as.numeric(x[9])-as.numeric(x[10]))/as.numeric(x[9])}), # non-Memory B cells
	group=c(rep("CTN",5),rep("Elder",5),rep("Young",6))
	)


tmp.res <- reshape2::melt(aggregate(.~group,data=data2,FUN=median))
colnames(tmp.res)[3] <- "median"
tmp.res$sd <- as.numeric(reshape2::melt(aggregate(.~group,data=data2,FUN=sd))[,3])
# plot result 
library(ggplot2)
library(ggpubr)
tmp.res$group <- factor(tmp.res$group,levels=c("Young","Elder","CTN"))

pdf("/public/workspace/lily/CTN/version_3_20/Cytometry/Panel2_res.pdf",useDingbats=F)
ggplot(data=tmp.res,aes(x=group,y=median,group=group)) + geom_bar(stat = "identity",aes(fill=group)) +
	geom_errorbar(aes(ymin = median, ymax = median + sd), width = 0.2, position = position_dodge(0.9)) +
	theme_classic() +
	facet_wrap(.~variable,scale="free")+ ylab(label="% all cells") +
	stat_compare_means(comparisons=list(c("Elder","CTN")))
dev.off()

























