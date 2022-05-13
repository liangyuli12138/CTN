#!/usr/bin/Rscript
# 2021-2-19
# TFs analysis in T cells 
#==================================================================================================


# first you have to generate a csv file 
# run in R 
library(Seurat)
dat <- readRDS("/public/workspace/lily/CTN/version_3_20/data/Tcell.RDS")
# sub.data=subset(dat,cells=which(dat$refine.group.new%in%c("CT.SC.dominant","CTN.dominant","NC.dominant")))
# if the result is all cells with TFs ,so do not subset cells 
write.csv(t(as.matrix(dat[['RNA']]@data)),file="/public/workspace/lily/CTN/version_3_20/trajectory/Tcell/all_Tcell.csv")

#===============================================================================================
# run in shell use pyscenic 
#===============================================================================================
#!/bin/sh


bytlib load python-3.6.6
/public/workspace/lily/Lung2Brain/inte7/Pyscenic/Scenic/arboreto_with_multiprocessing.py \
    -o /public/workspace/lily/CTN/version_3_20/trajectory/Tcell/all_Tcell.adjacencies.csv \
    --method grnboost2 \
    --seed 12345 \
    --num_workers 5 \
    /public/workspace/lily/CTN/version_3_20/trajectory/Tcell/all_Tcell.csv \
    /public/workspace/zhumy/ref/SCENIC/hs_hgnc_curated_tfs.txt 

echo "step1 done"
sleep 10

pyscenic ctx \
/public/workspace/lily/CTN/version_3_20/trajectory/Tcell/all_Tcell.adjacencies.csv \
/public/workspace/zhumy/ref/SCENIC/hg19-500bp-upstream-7species.mc9nr.feather \
/public/workspace/zhumy/ref/SCENIC/hg19-tss-centered-10kb-7species.mc9nr.feather \
--annotations_fname /public/workspace/lily/Lung2Brain/inte7/Pyscenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname /public/workspace/lily/CTN/version_3_20/trajectory/Tcell/all_Tcell.csv \
--output /public/workspace/lily/CTN/version_3_20/trajectory/Tcell/all_Tcell_regulons.tsv \
--num_workers 6 \
--mode "custom_multiprocessing"

#=====================================================================================================
echo "step2 done"
sleep 10

pyscenic aucell \
/public/workspace/lily/CTN/version_3_20/trajectory/Tcell/all_Tcell.csv /public/workspace/lily/CTN/version_3_20/trajectory/Tcell/all_Tcell_regulons.tsv \
--output /public/workspace/lily/CTN/version_3_20/trajectory/Tcell/all_Tcell_auc_mtx.tsv \
--num_workers 5 --seed 12345


# pyscenic aucell \
#     $expression_mtx_fname $ctx_output \
#     --output $aucell_output \
#     --num_workers 5 --seed 12345










































