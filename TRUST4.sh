
# 2021-9-8
# calculate TCR by 10X Genomic RNA-seq
#================================================================================================================================
/public/workspace/lily/software/TRUST4/run-trust4 \
-b /public/workspace/lily/CTN/Raw_data/XWS3/SI-GA-E1/outs/possorted_genome_bam.bam \
-f /public/workspace/lily/software/TRUST4/hg38_bcrtcr.fa \
--ref /public/workspace/lily/software/TRUST4/human_IMGT+C.fa \
-o XWS3 \
--od /public/workspace/lily/CTN/TRUST4/XWS3 \
-t 6 \
--barcode CB











































