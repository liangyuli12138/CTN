
# 2021-5-22
# this program is used to run Cellchat 
#=======================================================================================================
library(Seurat)
library(CellChat)

B <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Bcell.RDS")
T <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
T$celltype.new <- T$celltype.refine
NK <- readRDS("/public/workspace/lily/CTN/version_3_20/data/NK.pure.RDS")
Mono <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CD14Mono.RDS")

tmp <- merge(x=B,y=c(T,NK,Mono))

saveRDS(tmp,file="/public/workspace/lily/CTN/version_3_20/data/merge.celltype.RDS")

##########################################################################################################
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(CellChat)
options(future.globals.maxSize= 891289600) 

dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/merge.celltype.RDS")
# create a Cellchat object
data.input <- GetAssayData(dat, assay = "RNA", slot = "data") # normalized data matrix
meta <- data.frame(dat@meta.data) # create a dataframe of the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
# change ident and add meta information
#cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

CellChatDB <- CellChatDB.human
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB) 
# use all DB 
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

# 
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel

#####################################################################################################
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (should use)
cellchat <- projectData(cellchat, PPI.human)
# Calculate cell cell communication 
####################################################################################################
cellchat <- computeCommunProb(cellchat) # 1W cells 5min
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

saveRDS(cellchat,file="/public/workspace/lily/CTN/version_3_20/data/merge.celltype.cellchat.RDS")



######################################################################################################
# do some prepare 
# 2021-5-31 
# re-run analysis 
cellchat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/merge.celltype.cellchat.RDS")
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/merge.celltype.RDS")
colnames(cellchat@meta) <- colnames(dat@meta.data)
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
saveRDS(cellchat,file="/public/workspace/lily/CTN/version_3_20/data/merge.celltype.cellchat.RDS")

#=====================================================================================================
# plot result
cellchat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/merge.celltype.cellchat.RDS")
groupSize <- as.numeric(table(cellchat@idents))


pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/Fig6_TME.pdf",useDingbats=F)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",arrow.size=1,arrow.width=1)

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

tmp.dat <- subsetCommunication(cellchat)
tmp.dat$source.target <- paste0(tmp.dat$source," -> ",tmp.dat$target)
df <- tmp.dat[which(tmp.dat$prob>0.01&tmp.dat$pval<0.01),]
ggplot(df, aes(x = source.target, y = interaction_name_2,
	color = prob, size = pval)) + geom_point(pch = 16) +
	theme_linedraw() + theme(panel.grid.major = element_blank()) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1,
	vjust = 0.5), axis.title.x = element_blank(),
	axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
dev.off()









#=================================================================================================================================================
# 2021-6-21
# analysis cell cell communication 
# maybe should classify into 3 group to analysis 
#=================================================================================================================================================
library(Seurat)
library(CellChat)
library(ggplot2)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/merge.celltype.RDS")

sub1 <- subset(dat,cells=which(dat$group=="CTN"))
sub2 <- subset(dat,cells=which(dat$group=="Elder"))
sub3 <- subset(dat,cells=which(dat$group=="SC"))

data <- list(CTN=sub1,Elder=sub2,SCT=sub3)

# analysis 
#=================================================================================================================================================
for(i in 1:length(data)){
	library(patchwork)
	options(stringsAsFactors = FALSE)
	library(Seurat)
	library(CellChat)
	options(future.globals.maxSize= 891289600) 

# run cellchat
dat <- data[[i]]
# create a Cellchat object
data.input <- GetAssayData(dat, assay = "RNA", slot = "data") # normalized data matrix
meta <- data.frame(dat@meta.data) # create a dataframe of the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
# change ident and add meta information
#cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

CellChatDB <- CellChatDB.human
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB) 
# use all DB 
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

# 
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 8) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (should use)
cellchat <- projectData(cellchat, PPI.human)
# Calculate cell cell communication 
####################################################################################################
cellchat <- computeCommunProb(cellchat) # 1W cells 5min
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

saveRDS(cellchat,file=paste0("/public/workspace/lily/CTN/version_3_20/data/",names(data)[i],".cellchat.RDS"))


}



# check result 
#=====================================================================================================================================================
library(Seurat)
library(CellChat)
CTN <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN.cellchat.RDS")
groupSize <- table(CTN@idents)
# 1. plot cell cell communication numbers
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/CellChat/CTN.cellchat.pdf",useDingbats=F)
netVisual_circle(CTN@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.off()

# Elder 
Elder <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Elder.cellchat.RDS")
groupSize <- table(CTN@idents)
# 1. plot cell cell communication numbers
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/CellChat/Elder.cellchat.pdf",useDingbats=F)
netVisual_circle(Elder@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.off()


# SCT 
SCT <- readRDS("/public/workspace/lily/CTN/version_3_20/data/SCT.cellchat.RDS")
groupSize <- table(SCT@idents)
# 1. plot cell cell communication numbers
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/CellChat/SCT.cellchat.pdf",useDingbats=F)
netVisual_circle(SCT@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.off()






#============================================================================================================================================
# 2021-9-17
# add some figures for cell chat 
library(Seurat)
library(CellChat)

# 1.0 plot for merge data  
cellchat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/merge.celltype.cellchat.RDS")
groupSize <- as.numeric(table(cellchat@idents))

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/CellChat/merge.data.plot.pdf",useDingbats=F)
netVisual_heatmap(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = c("IL1"), width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = c("IL6"), width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = c("IGF"), width = 8, height = 2.5, font.size = 10)
dev.off()




#==========================================================================================================================================
# 2.0 add to a list 
CTN <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN.cellchat.RDS")
Elder <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Elder.cellchat.RDS")
SCT <- readRDS("/public/workspace/lily/CTN/version_3_20/data/SCT.cellchat.RDS")

object.list <- list(CTN = CTN, Elder = Elder,SCT=SCT)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat



# 2.1 interaction weight 
pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/CellChat/interaction_weight_3group.pdf",useDingbats=F)
compareInteractions(cellchat, show.legend = F,measure="weight") # show higher in SCT and CTN
dev.off()


# 2.2  not ok 
pdf("~/tmp.1.pdf")
netVisual_heatmap(Elder,measure="weight")
netVisual_heatmap(CTN,measure="weight")
netVisual_heatmap(SCT,measure="weight")
dev.off()


# library(ComplexHeatmap)
# cellchat <- netAnalysis_computeCentrality(cellchat)
# i = 1
# # combining all the identified signaling pathways from different datasets 
# pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
# ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
# ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
# ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 6)

# draw(ht1 + ht2 +ht3, ht_gap = unit(0.5, "cm"))








pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/CellChat/signaling_Supplementary.pdf",useDingbats=F)
CTN <- netAnalysis_computeCentrality(CTN)
netAnalysis_signalingRole_network(CTN, signaling = c("IL1"), width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(CTN, signaling = c("IL6"), width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(CTN, signaling = c("IGF"), width = 8, height = 2.5, font.size = 10)

Elder <- netAnalysis_computeCentrality(Elder)
netAnalysis_signalingRole_network(Elder, signaling = c("IL1"), width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(Elder, signaling = c("IL6"), width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(Elder, signaling = c("IGF"), width = 8, height = 2.5, font.size = 10)

SCT <- netAnalysis_computeCentrality(SCT)
netAnalysis_signalingRole_network(SCT, signaling = c("IL1"), width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(SCT, signaling = c("IL6"), width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(SCT, signaling = c("IGF"), width = 8, height = 2.5, font.size = 10)
dev.off()


























#========================================================================================================================================
# analysis about T cell 
# 2021-6-21
library(Seurat)
library(CellChat)

CTN <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN.cellchat.RDS")
Elder <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Elder.cellchat.RDS")
SCT <- readRDS("/public/workspace/lily/CTN/version_3_20/data/SCT.cellchat.RDS")


tmp <- netVisual_bubble(CTN, sources.use = 2, targets.use = 4, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "CTN"
tmp1 <- tmp.res

# Elder 
tmp <- netVisual_bubble(Elder, sources.use = 2, targets.use = 4, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "Elder"
tmp2 <- tmp.res

# SCT 
tmp <- netVisual_bubble(SCT, sources.use = 2, targets.use = 4, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "SCT"
tmp3 <- tmp.res


tmp.f <- merge(tmp1,tmp2,all.x=T,all.y=T,by="interaction_name_2")
colnames(tmp.f)[2:5] <- c("prob.CTN","group.CTN","prob.Elder","group.Elder")
tmp.rs <- merge(tmp.f,tmp3,all.x=T,all.y=T,by="interaction_name_2")

# just use prob info 
res <- tmp.rs[,c(1,2,4,6)]
colnames(res)[4] <- "prob.SCT"

# check result 
tmp.1 <- res[which(res$prob.SCT>0&res$prob.CTN>0&is.na(res$prob.Elder)),]
tmp.2 <- res[which(res$prob.Elder>0&is.na(res$prob.CTN)&is.na(res$prob.SCT)),]


rownames(tmp.1) <- tmp.1$interaction_name_2
tmp.1$interaction_name_2 <- NULL
tmp.1 <- tmp.1[order(tmp.1$prob.SCT,decreasing=T),]


rownames(tmp.2) <- tmp.2$interaction_name_2
tmp.2$interaction_name_2 <- NULL
tmp.2 <- tmp.2[order(tmp.2$prob.Elder,decreasing=T),]


pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/CellChat/Mono.to.Tcell.pdf",useDingbats=F)
pheatmap::pheatmap(tmp.1,cluster_row=F,cluster_col=F)
pheatmap::pheatmap(tmp.2,cluster_row=F,cluster_col=F)
dev.off()


#==========================================================================================================================================
# Venn plot show different pathway 

library(Seurat)
library(CellChat)

CTN <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN.cellchat.RDS")
Elder <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Elder.cellchat.RDS")
SCT <- readRDS("/public/workspace/lily/CTN/version_3_20/data/SCT.cellchat.RDS")

tmp.sct <- subsetCommunication(SCT,slot.name="netP")
tmp.sct <- tmp.sct[order(tmp.sct$prob,decreasing=T),]
sct <- unique(head(tmp.sct,20)$pathway_name)

tmp.ctn <- subsetCommunication(CTN,slot.name="netP")
tmp.ctn <- tmp.ctn[order(tmp.ctn$prob,decreasing=T),]
ctn <- unique(head(tmp.ctn,20)$pathway_name)

tmp.eld <- subsetCommunication(Elder,slot.name="netP")
tmp.eld <- tmp.eld[order(tmp.eld$prob,decreasing=T),]
eld <- unique(head(tmp.eld,20)$pathway_name)

dat <- list(SCT=sct,CTN=ctn,Elder=eld)

pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/CellChat/Top.pathway.pdf",useDingbats=F)
gplots::venn(dat)
dev.off()







sct <- unique(tmp.sct$pathway_name)
ctn <- unique(tmp.ctn$pathway_name)
eld <- unique(tmp.eld$pathway_name)

tmp <- intersect(head(ctn,40),head(sct,40))

tmp[-which(tmp%in%head(eld,40))]












#=========================================================================================================================================================
# analysis about Monocyte to B cell and NK cells
# maybe not suitable 

library(Seurat)
library(CellChat)

CTN <- readRDS("/public/workspace/lily/CTN/version_3_20/data/CTN.cellchat.RDS")
Elder <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Elder.cellchat.RDS")
SCT <- readRDS("/public/workspace/lily/CTN/version_3_20/data/SCT.cellchat.RDS")


tmp <- netVisual_bubble(CTN, sources.use = 2, targets.use = 3, remove.isolate = T,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "CTN"
tmp1 <- tmp.res

# Elder 
tmp <- netVisual_bubble(Elder, sources.use = 2, targets.use = 3, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "Elder"
tmp2 <- tmp.res

# SCT 
tmp <- netVisual_bubble(SCT, sources.use = 2, targets.use = 3, remove.isolate = T,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "SCT"
tmp3 <- tmp.res


tmp.f <- merge(tmp1,tmp2,all.x=T,all.y=T,by="interaction_name_2")
colnames(tmp.f)[2:5] <- c("prob.CTN","group.CTN","prob.Elder","group.Elder")
tmp.rs <- merge(tmp.f,tmp3,all.x=T,all.y=T,by="interaction_name_2")

# just use prob info 
res <- tmp.rs[,c(1,2,4,6)]
colnames(res)[4] <- "prob.SCT"

# check result 
tmp.1 <- res[which(res$prob.SCT>0&res$prob.CTN>0&is.na(res$prob.Elder)),]
tmp.2 <- res[which(res$prob.Elder>0&is.na(res$prob.CTN)&is.na(res$prob.SCT)),]


rownames(tmp.1) <- tmp.1$interaction_name_2
tmp.1$interaction_name_2 <- NULL
tmp.1 <- tmp.1[order(tmp.1$prob.SCT,decreasing=T),]


rownames(tmp.2) <- tmp.2$interaction_name_2
tmp.2$interaction_name_2 <- NULL
tmp.2 <- tmp.2[order(tmp.2$prob.Elder,decreasing=T),]


pdf("/public/workspace/lily/CTN/version_3_20/rs_plot/CellChat/Mono.to.NKcell.pdf",useDingbats=F)
pheatmap::pheatmap(tmp.1,cluster_row=F,cluster_col=F)
pheatmap::pheatmap(tmp.2,cluster_row=F,cluster_col=F)
dev.off()













tmp1 <- netVisual_bubble(CTN, sources.use = 2, targets.use = c(1,3,4), remove.isolate = T,
	return.data=T,signaling=c("MHC-II","MHC-I","CCL","CLEC","GALECTIN","CD22","CD45","CD99"))$communication

tmp2 <- netVisual_bubble(Elder, sources.use = 2, targets.use = c(1,3,4), remove.isolate = T,
	return.data=T,signaling=c("MHC-II","MHC-I","CCL","CLEC","GALECTIN","CD22","CD45","CD99"))$communication

tmp3 <- netVisual_bubble(SCT, sources.use = 2, targets.use = c(1,3,4), remove.isolate = T,
	return.data=T,signaling=c("MHC-II","MHC-I","CCL","CLEC","GALECTIN","CD22","CD45","CD99"))$communication

 # cast 
rs1 <- reshape2::dcast(tmp1,source.target~interaction_name_2,value.var="prob")
rownames(rs1) <- rs1$source.target
rs1$source.target <- NULL
rs1 <- t(rs1)
colnames(rs1) <- paste0(colnames(rs1),".CTN")

rs2 <- reshape2::dcast(tmp2,source.target~interaction_name_2,value.var="prob")
rownames(rs2) <- rs2$source.target
rs2$source.target <- NULL
rs2 <- t(rs2)
colnames(rs2) <- paste0(colnames(rs2),".Elder")

rs3 <- reshape2::dcast(tmp3,source.target~interaction_name_2,value.var="prob")
rownames(rs3) <- rs3$source.target
rs3$source.target <- NULL
rs3 <- t(rs3)
colnames(rs3) <- paste0(colnames(rs3),".SCT")


a <- merge(rs1,rs2,all.x=T,all.y=T,by="row.names")
rownames(a) <- a$Row.names
a$Row.names <- NULL
b <- merge(a,rs3,all.x=T,all.y=T,by="row.names")
rownames(b) <- b$Row.names
b$Row.names <- NULL
pheatmap::pheatmap(b,cluster_rows=F,cluster_cols=F)


























































