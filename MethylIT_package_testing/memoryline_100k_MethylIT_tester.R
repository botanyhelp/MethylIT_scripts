
# .libPaths(c("/usr/lib64/R/library","/usr/share/R/library","~/R/x86_64-redhat-linux-gnu-library/3.6"))
library(BiocManager)
library(MethylIT)
library(MethylIT.utils)

pryr::mem_used()

date()

githubURL <- "https://github.com/botanyhelp/MethylIT_scripts/raw/master/MethylIT_package_testing/RData/memoryLine_100k_samples_GR_CG.RData"
load(url(githubURL))

names(memoryLine_100k_samples_GR_CG)
#  [1] "M_2_1" "M_2_2" "M_2_3" "M_2_4" "M_2_5" "W_2_1" "W_2_2" "W_2_3" "W_2_4" "W_2_5" "M_3_1" "M_3_2" "M_3_3" "M_3_4" "M_3_5" "W_3_1" "W_3_2" "W_3_3" "W_3_4" "W_3_5" "M_4_1" "M_4_2" "M_4_3" "M_4_4" "M_4_5" "W_4_1" "W_4_2" "W_4_3" "W_4_4"
# [30] "W_4_5" "M_5_1" "M_5_2" "M_5_3" "M_5_4" "M_5_5" "W_5_1" "W_5_2" "W_5_3" "W_5_4" "W_5_5" "M_6_1" "M_6_2" "M_6_3" "M_6_4" "M_6_5" "W_6_1" "W_6_2" "W_6_3" "W_6_4" "W_6_5"

ml_samples_W31W32W33_Ref_CG_sum <- poolFromGRlist(list(memoryLine_100k_samples_GR_CG$W_3_1, memoryLine_100k_samples_GR_CG$W_3_2, memoryLine_100k_samples_GR_CG$W_3_3), stat = "sum", num.cores = 12L)
ml_samples_W34W35M31M32M33M34M35_LR <- memoryLine_100k_samples_GR_CG[c("W_3_4","W_3_5","M_3_1","M_3_2","M_3_3","M_3_4","M_3_5")]
names(ml_samples_W34W35M31M32M33M34M35_LR)
# [1] "W_3_4" "W_3_5" "M_3_1" "M_3_2" "M_3_3" "M_3_4" "M_3_5"

ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR = estimateDivergence(ref = ml_samples_W31W32W33_Ref_CG_sum, 
                                                     indiv = ml_samples_W34W35M31M32M33M34M35_LR, 
                                                     Bayesian = TRUE, 
                                                     min.coverage = 4, 
                                                     high.coverage = 300, 
                                                     percentile = 0.999, 
                                                     num.cores = 64L, tasks = 20L, verbose = TRUE )

nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull2P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="Weibull2P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull3P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="Weibull3P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma2P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="Gamma2P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma3P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="Gamma3P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma3P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="GGamma3P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma4P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="GGamma4P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_ECDF = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="ECDF", num.cores = 64L, verbose = TRUE)


PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull2P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull2P, dist.name="Weibull2P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
#NO: PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull3P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull3P, dist.name="Weibull3P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
# Error in pweibull(q - m[3], shape = m[1], scale = m[2], lower.tail = FALSE) : 
#   Non-numeric argument to mathematical function
PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma2P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma2P, dist.name="Gamma2P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma3P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma3P, dist.name="Gamma3P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma3P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma3P, dist.name="GGamma3P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma4P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma4P, dist.name="GGamma4P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_ECDF = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_ECDF, dist.name="ECDF", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)

cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR = estimateCutPoint(LR = PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, control.names = c( "W_3_4", "W_3_5"),treatment.names = c("M_3_1","M_3_2","M_3_3","M_3_4","M_3_5"),div.col = 9, verbose = FALSE)
cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR$cutpoint
# [1] 2.514052

DIMPs_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR <- selectDIMP(PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, div.col = 9, cutpoint=cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR$cutpoint)

AG_gff3 = import("ftp://ftp.ensemblgenomes.org/pub/plants/current/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.44.gff3.gz")
# AG_gff3 = import("/data/TAIR10_gff3/Arabidopsis_thaliana.TAIR10.38.gff3.gz")
dim(as.data.frame(AG_gff3[AG_gff3$type=="gene"]))
# [1] 27655    25

genes <-AG_gff3[AG_gff3$type=="gene"]

Genes2kb = GeneUpDownStream(genes, upstream = 2000, downstream = 2000)

for(i in names(ml_samples_W34W35M31M32M33M34M35_LR)){assign(paste0("DIMPs_ml_",i,"_gene"),getDIMPatGenes(GR = DIMPs_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR[[i]],GENES=Genes2kb))}

head(DIMPs_ml_M_3_1_gene)
# GRanges object with 6 ranges and 2 metadata columns:
#       seqnames      ranges strand |    GeneID     DIMPs
#          <Rle>   <IRanges>  <Rle> |  <factor> <integer>
#   [1]        1   1631-7899      + | AT1G01010         2
#   [2]        1  4788-11130      - | AT1G01020         2
#   [3]        1  9649-15714      - | AT1G01030         1
#   [4]        1 21121-33227      + | AT1G01040        11
#   [5]        1 29170-35171      - | AT1G01050         7
#   [6]        1 31365-39871      - | AT1G01060         2

Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35 = uniqueGRanges(list(DIMPs_ml_M_3_1_gene[,2],DIMPs_ml_M_3_2_gene[,2],DIMPs_ml_M_3_3_gene[,2],DIMPs_ml_M_3_4_gene[,2],DIMPs_ml_M_3_5_gene[,2],DIMPs_ml_W_3_4_gene[,2],DIMPs_ml_W_3_5_gene[,2]), type = "equal", verbose = TRUE, chromosomes = c("1", "2", "3", "4", "5"), ignore.strand = TRUE )

colnames( mcols(Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35)) <- c("M_3_1","M_3_2","M_3_3","M_3_4","M_3_5","W_3_4","W_3_5")

# Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35_ID = subsetByOverlaps(Genes2kb, Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35, type = "equal", ignore.strand = FALSE)
# 
# dmps = data.frame( mcols( Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35 ) )
# dmps = apply( dmps, 2, as.numeric )
# rownames(dmps) <- Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35_ID$gene_id

condition = data.frame(condition = factor(c("TT","TT","TT","TT","TT","CT","CT"), levels = c("CT", "TT")))
rownames(condition)
# [1] "1" "2" "3" "4" "5" "6" "7"
names(mcols(Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35))
# [1] "M_3_1" "M_3_2" "M_3_3" "M_3_4" "M_3_5" "W_3_4" "W_3_5"

rownames(condition) <- names(mcols(Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35))

DIMR_ml_samples_W34W35M31M32M33M34M35 <- glmDataSet(GR = Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35, colData = condition)

DMGs_ml_samples_W34W35M31M32M33M34M35_countTest2 = countTest2(DIMR_ml_samples_W34W35M31M32M33M34M35, num.cores = 24L, countFilter = T)
dim(mcols(DMGs_ml_samples_W34W35M31M32M33M34M35_countTest2))
# [1] 29 14
# [1] 31 14    ## 20190730 1657pm found 2 more DMGs because it used the latest version of MethylIT:
