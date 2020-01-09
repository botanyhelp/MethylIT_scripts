# soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200109.R


.libPaths()

library(BiocManager)
library(MethylIT)
library(MethylIT.utils)

rm(list=ls())
length(ls()) 
# [1] 0

# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F4_P37_CG_4_29_2019.RData")
# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F6_P37_CG_4_29_2019.RData")

# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F4_P37_CHH_4_29_2019.RData")
# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F6_P37_CHH_4_29_2019.RData")


load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F4_P37_CHH_4_29_2019.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F6_P37_CHH_4_29_2019.RData")

names(F4_P37_CHH)
# [1] "G13A_CHH" "G14A_CHH" "G15A_CHH"

names(F6_P37_CHH)
# [1] "G22A_CHH" "G23A_CHH" "G24A_CHH"

ls()
# [1] "F4_P37_CHH" "F6_P37_CHH"

seqnames(F4_P37_CHH$G14A_CHH)
# factor-Rle of length 254667767 with 20 runs
#   Lengths: 15619475 13855537  9164311 10652937 11910229 13405530 ... 11297345 13530428 11878478 12636861 13507325
#   Values :        1       10       11       12       13       14 ...        5        6        7        8        9
# Levels(20): 1 10 11 12 13 14 15 16 17 18 19 2 20 3 4 5 6 7 8 9

G22A_CHH_grl <- split(F6_P37_CHH$G22A_CHH,seqnames(F4_P37_CHH$G14A_CHH))
G23A_CHH_grl <- split(F6_P37_CHH$G23A_CHH,seqnames(F4_P37_CHH$G14A_CHH))
G24A_CHH_grl <- split(F6_P37_CHH$G24A_CHH,seqnames(F4_P37_CHH$G14A_CHH))

G13A_CHH_grl <- split(F4_P37_CHH$G13A_CHH,seqnames(F4_P37_CHH$G13A_CHH))
G14A_CHH_grl <- split(F4_P37_CHH$G14A_CHH,seqnames(F4_P37_CHH$G13A_CHH))
G15A_CHH_grl <- split(F4_P37_CHH$G15A_CHH,seqnames(F4_P37_CHH$G13A_CHH))

names(G22A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

names(G14A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

## 20190109 saved:
#
# save(F4_P37_CHH,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/F4_P37_CHH.RData")
# save(F6_P37_CHH,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/F6_P37_CHH.RData")
# save(G13A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G13A_CHH_grl.RData")
# save(G14A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G14A_CHH_grl.RData")
# save(G15A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G15A_CHH_grl.RData")
# save(G22A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G22A_CHH_grl.RData")
# save(G23A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G23A_CHH_grl.RData")
# save(G24A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G24A_CHH_grl.RData")


date()
# [1] "Thu Jan  9 10:44:19 2020"

# 32min runtime:
for(item in names(G14A_CHH_grl)) {
	    assign(
		           paste0("soy_WT_G222324_Ref_CHH_sum_chr", item),
			           poolFromGRlist(
						              list(G22A_CHH_grl[[item]], G23A_CHH_grl[[item]], G24A_CHH_grl[[item]]),
							                  stat = "sum",
									              num.cores = 12L
							              ))}

date()
# [1] "Thu Jan  9 11:16:09 2020"

ls()
#  [1] "F4_P37_CHH"                       "F6_P37_CHH"                       "G13A_CHH_grl"                    
#  [4] "G14A_CHH_grl"                     "G15A_CHH_grl"                     "G22A_CHH_grl"                    
#  [7] "G23A_CHH_grl"                     "G24A_CHH_grl"                     "item"                            
# [10] "soy_WT_G222324_Ref_CHH_sum_chr1"  "soy_WT_G222324_Ref_CHH_sum_chr10" "soy_WT_G222324_Ref_CHH_sum_chr11"
# [13] "soy_WT_G222324_Ref_CHH_sum_chr12" "soy_WT_G222324_Ref_CHH_sum_chr13" "soy_WT_G222324_Ref_CHH_sum_chr14"
# [16] "soy_WT_G222324_Ref_CHH_sum_chr15" "soy_WT_G222324_Ref_CHH_sum_chr16" "soy_WT_G222324_Ref_CHH_sum_chr17"
# [19] "soy_WT_G222324_Ref_CHH_sum_chr18" "soy_WT_G222324_Ref_CHH_sum_chr19" "soy_WT_G222324_Ref_CHH_sum_chr2" 
# [22] "soy_WT_G222324_Ref_CHH_sum_chr20" "soy_WT_G222324_Ref_CHH_sum_chr3"  "soy_WT_G222324_Ref_CHH_sum_chr4" 
# [25] "soy_WT_G222324_Ref_CHH_sum_chr5"  "soy_WT_G222324_Ref_CHH_sum_chr6"  "soy_WT_G222324_Ref_CHH_sum_chr7" 
# [28] "soy_WT_G222324_Ref_CHH_sum_chr8"  "soy_WT_G222324_Ref_CHH_sum_chr9" 

length(ls(pattern="soy_WT_G222324_Ref_CHH_sum_chr"))
# [1] 20

## 20190109 saved:
#
# save(soy_WT_G222324_Ref_CHH_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr1.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr10.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr11.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr12.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr13.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr14.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr15.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr16.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr17.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr18.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr19.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr2.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr20.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr3.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr4.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr5.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr6.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr7.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr8.RData")
# save(soy_WT_G222324_Ref_CHH_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr9.RData")

date()

## 20190109: this code was NOT run..instead, the 4 code blocks shown farther below 
#  ..was run concurrently in 4 R consoles simultaneously.  That is why code is 
#  ..shown farther down that uses load() to get those objects available here.
#  ..this code works but was not run and will take a long time if it does run. 
#
for(item in names(G14A_CHH_grl)) {
	    assign(
		           paste0(
				              "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr",
					                  item
					              ),
		           estimateDivergence(
					                  ref = eval(str2expression(paste0("soy_WT_G222324_Ref_CHH_sum_chr",item))),
							              indiv = list(G22A_CHH_grl[[item]], G23A_CHH_grl[[item]], G24A_CHH_grl[[item]],G13A_CHH_grl[[item]], G14A_CHH_grl[[item]], G15A_CHH_grl[[item]]),
								                  Bayesian = TRUE,
										              min.coverage = c(12,4),
											                  min.meth = 3,
													              high.coverage = 300,
														                  percentile = 0.999,
																              num.cores = 64L,
																	                  tasks = 20L,
																			              verbose = FALSE
							          )
			       )
}

################################################################################
################################################################################
# START multi-R console to cut runtime by 75%:
################################################################################
################################################################################

# Here we see 4 blocks of code that process 5 of the 20 soybean chromoaomes. 
# These code blocks can be run simultaneously (lets hope!) and cut the time by 4x.



# ################################################################################
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G13A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G14A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G15A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G22A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G23A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G24A_CHH_grl.RData")
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr1.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr10.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr11.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr12.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr13.RData")
# 
# for(item in c("1","10","11","12","13")) {
#     assign(
#         paste0(
#             "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr",
#             item
#         ),
#         estimateDivergence(
#             ref = eval(str2expression(paste0("soy_WT_G222324_Ref_CHH_sum_chr",item))),
#             indiv = list(G22A_CHH_grl[[item]], G23A_CHH_grl[[item]], G24A_CHH_grl[[item]],G13A_CHH_grl[[item]], G14A_CHH_grl[[item]], G15A_CHH_grl[[item]]),
#             Bayesian = TRUE,
#             min.coverage = c(12,4),
#             min.meth = 3,
#             high.coverage = 300,
#             percentile = 0.999,
#             num.cores = 64L,
#             tasks = 20L,
#             verbose = FALSE
#         )
#     )
# }
# 
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13.RData")
# ################################################################################
# 
# 
# 
# 
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G13A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G14A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G15A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G22A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G23A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G24A_CHH_grl.RData")
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr14.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr15.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr16.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr17.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr18.RData")
# 
# for(item in c("14","15","16","17","18")) {
#     assign(
#         paste0(
#             "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr",
#             item
#         ),
#         estimateDivergence(
#             ref = eval(str2expression(paste0("soy_WT_G222324_Ref_CHH_sum_chr",item))),
#             indiv = list(G22A_CHH_grl[[item]], G23A_CHH_grl[[item]], G24A_CHH_grl[[item]],G13A_CHH_grl[[item]], G14A_CHH_grl[[item]], G15A_CHH_grl[[item]]),
#             Bayesian = TRUE,
#             min.coverage = c(12,4),
#             min.meth = 3,
#             high.coverage = 300,
#             percentile = 0.999,
#             num.cores = 64L,
#             tasks = 20L,
#             verbose = FALSE
#         )
#     )
# }
# 
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18.RData")
# ################################################################################
# 
# 
# 
# 
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G13A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G14A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G15A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G22A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G23A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G24A_CHH_grl.RData")
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr19.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr2.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr20.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr3.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr4.RData")
# 
# for(item in c("19","2","20","3","4")) {
#     assign(
#         paste0(
#             "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr",
#             item
#         ),
#         estimateDivergence(
#             ref = eval(str2expression(paste0("soy_WT_G222324_Ref_CHH_sum_chr",item))),
#             indiv = list(G22A_CHH_grl[[item]], G23A_CHH_grl[[item]], G24A_CHH_grl[[item]],G13A_CHH_grl[[item]], G14A_CHH_grl[[item]], G15A_CHH_grl[[item]]),
#             Bayesian = TRUE,
#             min.coverage = c(12,4),
#             min.meth = 3,
#             high.coverage = 300,
#             percentile = 0.999,
#             num.cores = 64L,
#             tasks = 20L,
#             verbose = FALSE
#         )
#     )
# }
# 
# 
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4.RData")
# ################################################################################
# 
# 
# 
# 
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G13A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G14A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G15A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G22A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G23A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G24A_CHH_grl.RData")
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr5.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr6.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr7.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr8.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr9.RData")
# 
# for(item in c("5","6","7","8","9")) {
#     assign(
#         paste0(
#             "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr",
#             item
#         ),
#         estimateDivergence(
#             ref = eval(str2expression(paste0("soy_WT_G222324_Ref_CHH_sum_chr",item))),
#             indiv = list(G22A_CHH_grl[[item]], G23A_CHH_grl[[item]], G24A_CHH_grl[[item]],G13A_CHH_grl[[item]], G14A_CHH_grl[[item]], G15A_CHH_grl[[item]]),
#             Bayesian = TRUE,
#             min.coverage = c(12,4),
#             min.meth = 3,
#             high.coverage = 300,
#             percentile = 0.999,
#             num.cores = 64L,
#             tasks = 20L,
#             verbose = FALSE
#         )
#     )
# }
# 
# 
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9.RData")

################################################################################
################################################################################
################################################################################
# END multi-R console to cut runtime by 75%:
################################################################################
################################################################################
################################################################################

date()
# [1] "Thu Jan  9 16:40:03 2020"

# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F6_R10_CHH_4_29_2019.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9.RData")


length(ls(pattern="HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr"))
# [1] 20

names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1)
# NULL

names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")

names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1)
# [1] "G22A_CHH" "G23A_CHH" "G24A_CHH" "G13A_CHH" "G14A_CHH" "G15A_CHH"


for(item in names(G14A_CHH_grl)) {
	    assign(paste0("critical.val_G222324refctrl_vs_G131415_sum_chr", item),
		       do.call(rbind, lapply(eval(str2expression(paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr",item))), function(x) {
					         hd.95 = quantile(x$hdiv, 0.95)
						     baytv.95 = quantile(abs(x$bay.TV), 0.95)
						     tv.95 = quantile(abs(x$TV), 0.95)
						         return(c(tv = tv.95, baytv = baytv.95, hd = hd.95, num.sites.hd95 = sum(x$hdiv > hd.95), num.sites.tv95 = sum(x$TV > tv.95), num.sites.baytv95 = sum(x$bay.TV > baytv.95)))}))
		       )
}
    
for(item in names(G14A_CHH_grl)) {
	 print(eval(str2expression(paste0("critical.val_G222324refctrl_vs_G131415_sum_chr",item))))
}
## 20190109: 
#
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1702128 0.08572858 0.3633274          77216          28090             13260
# G23A_CHH 0.1775873 0.08798573 0.3667813          71467          40638             45736
# G24A_CHH 0.1958042 0.10661883 0.3897832          48326          35543             42428
# G13A_CHH 0.2500000 0.13016440 0.6369731          58507          18344             10921
# G14A_CHH 0.1875862 0.11691683 1.0376375          93064          53812             59550
# G15A_CHH 0.2500000 0.13731769 0.7895883          64005          39620             48166
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1666667 0.08508452 0.3842621          63988          22807             10078
# G23A_CHH 0.1700000 0.08571900 0.3877397          61553          34482             37290
# G24A_CHH 0.1851852 0.10479647 0.4567538          48466          36212             41919
# G13A_CHH 0.2352941 0.12836617 0.7039856          55242          17529              8897
# G14A_CHH 0.1858696 0.11640307 1.0663130          74443          42940             47361
# G15A_CHH 0.2361111 0.13346783 0.8398077          58101          37511             42743
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1739130 0.08764946 0.3711059          34124          12882              6221
# G23A_CHH 0.1818182 0.08980278 0.3733259          31833          17771             19778
# G24A_CHH 0.1929825 0.10691250 0.4177441          24051          17657             20845
# G13A_CHH 0.2541273 0.13422005 0.6508937          25608           8599              4765
# G14A_CHH 0.1906308 0.11977002 1.0606248          41539          23959             26410
# G15A_CHH 0.2500000 0.13883708 0.7905970          28403          17884             20786
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1685079 0.08695618 0.3826133          46358          17440              7762
# G23A_CHH 0.1750000 0.08762936 0.3896761          44672          24681             25863
# G24A_CHH 0.1875000 0.10570427 0.4633457          36424          26544             30680
# G13A_CHH 0.2400000 0.13065409 0.6983777          39852          12893              6587
# G14A_CHH 0.1904762 0.11845125 1.0621110          53855          30561             33824
# G15A_CHH 0.2393939 0.13458041 0.8399924          42755          26683             30272
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1764706 0.09035498 0.3958731          38686          15180              7069
# G23A_CHH 0.1818182 0.09043757 0.3917350          36779          20187             21230
# G24A_CHH 0.1923077 0.10710433 0.4377925          28233          20397             23845
# G13A_CHH 0.2500000 0.13586566 0.7227406          32575          10428              5459
# G14A_CHH 0.1971417 0.12353670 1.0984441          45466          26740             29541
# G15A_CHH 0.2500000 0.13861246 0.8431039          35280          21030             24806
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1720779 0.08395616 0.3454061          61812          22767             12129
# G23A_CHH 0.1789773 0.08735663 0.3520917          56173          31877             37705
# G24A_CHH 0.2000000 0.10679407 0.3540858          33841          25053             30595
# G13A_CHH 0.2525253 0.12822551 0.5901472          43189          14878              9028
# G14A_CHH 0.1884615 0.11603940 1.0129287          75229          43375             48434
# G15A_CHH 0.2536484 0.13811607 0.7465016          46650          30427             36231
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1739130 0.09535415 0.4239546          53578          19737              6870
# G23A_CHH 0.1785714 0.09455904 0.4248049          52699          29001             25635
# G24A_CHH 0.1818182 0.11032908 0.5711890          50984          36751             34050
# G13A_CHH 0.2375000 0.14079851 0.8100267          51743          16278              6326
# G14A_CHH 0.2000000 0.12538484 1.1030601          61731          33590             37681
# G15A_CHH 0.2389937 0.13931037 0.9025650          53677          32556             33952
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1703297 0.08572355 0.3646013          42293          15678              7619
# G23A_CHH 0.1769042 0.08780623 0.3676937          39315          22361             24714
# G24A_CHH 0.1876209 0.10512335 0.4342014          31776          23537             27358
# G13A_CHH 0.2417582 0.12916599 0.6736894          35713          11500              6315
# G14A_CHH 0.1875000 0.11697468 1.0638778          51850          29499             33242
# G15A_CHH 0.2467532 0.13578253 0.8001238          36480          22670             26644
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1666667 0.08642830 0.4038205          49946          17889              7668
# G23A_CHH 0.1666667 0.08570133 0.4130120          49907          27987             28446
# G24A_CHH 0.1717172 0.10057944 0.5840247          50023          37150             40732
# G13A_CHH 0.2222222 0.12809469 0.8012522          50736          15434              6789
# G14A_CHH 0.1875000 0.11827447 1.0854928          55951          32307             35302
# G15A_CHH 0.2268908 0.12975165 0.9067941          51582          31780             34975
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1733333 0.09003565 0.3941853          68518          25262             11096
# G23A_CHH 0.1764706 0.09072788 0.3989271          66258          37021             36121
# G24A_CHH 0.1863636 0.10768586 0.4880410          57100          41314             43846
# G13A_CHH 0.2400000 0.13468682 0.7445786          62492          19913              9279
# G14A_CHH 0.1944444 0.12154212 1.0807811          80398          46003             49988
# G15A_CHH 0.2411765 0.13663610 0.8503343          65122          39998             44217
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1666667 0.08430907 0.3711623          65546          23841             10652
# G23A_CHH 0.1724138 0.08517116 0.3754506          62142          34598             38743
# G24A_CHH 0.1868132 0.10456967 0.4461310          49057          36429             43177
# G13A_CHH 0.2381234 0.12726001 0.6665924          54186          17603              9532
# G14A_CHH 0.1857143 0.11523623 1.0484332          77605          44728             49497
# G15A_CHH 0.2380952 0.13303457 0.8200026          58593          37115             43219
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1677019 0.08579754 0.3769762          53872          20337              9488
# G23A_CHH 0.1750000 0.08685048 0.3825409          51495          28691             31438
# G24A_CHH 0.1875000 0.10503918 0.4371538          39390          29092             34103
# G13A_CHH 0.2400000 0.12916016 0.6879230          45764          14895              7844
# G14A_CHH 0.1875000 0.11707088 1.0608481          63965          36469             40713
# G15A_CHH 0.2380952 0.13330817 0.8275686          49764          31302             35866
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1657143 0.08365454 0.3808387          65523          23371             10341
# G23A_CHH 0.1666667 0.08381870 0.3879392          63885          36143             39290
# G24A_CHH 0.1833333 0.10354016 0.4594700          49996          37277             44020
# G13A_CHH 0.2307692 0.12509868 0.7048490          58636          17843              9068
# G14A_CHH 0.1833333 0.11417040 1.0590002          75835          42880             47182
# G15A_CHH 0.2272727 0.12842397 0.8605776          64709          40647             46710
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1764706 0.08800224 0.3563486          46820          18145              9197
# G23A_CHH 0.1861472 0.09126118 0.3597739          42312          23559             26610
# G24A_CHH 0.2032258 0.10881540 0.3617690          27213          19722             23538
# G13A_CHH 0.2631579 0.13409670 0.6151753          33347          11867              6869
# G14A_CHH 0.1964286 0.12060142 1.0332247          57065          32725             36214
# G15A_CHH 0.2631579 0.14255461 0.7518475          35879          23023             27041
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1785714 0.09278200 0.3745215          55307          20738              9003
# G23A_CHH 0.1875000 0.09524377 0.3716191          50532          28401             29717
# G24A_CHH 0.2083333 0.11388867 0.3810290          33023          23339             26224
# G13A_CHH 0.2666667 0.14058430 0.6411083          39612          13413              6958
# G14A_CHH 0.2000000 0.12361135 1.0457716          67879          37556             42215
# G15A_CHH 0.2666667 0.14702685 0.7718042          42722          27194             31163
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1666667 0.08426926 0.3827010          51402          18390              8454
# G23A_CHH 0.1681818 0.08497054 0.3849733          49409          27933             30029
# G24A_CHH 0.1837838 0.10386990 0.4651248          39885          29680             34441
# G13A_CHH 0.2323232 0.12654806 0.6989990          45020          14203              7216
# G14A_CHH 0.1851852 0.11576469 1.0627181          59591          33887             37588
# G15A_CHH 0.2307692 0.13040271 0.8493567          49013          30887             35291
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1666667 0.08534436 0.3836610          58186          21621              9877
# G23A_CHH 0.1724942 0.08593526 0.3868555          55986          31212             33725
# G24A_CHH 0.1833333 0.10367979 0.4726601          47142          35096             40692
# G13A_CHH 0.2352941 0.12842016 0.7075758          51325          16581              8349
# G14A_CHH 0.1875000 0.11748964 1.0706983          67466          38736             43074
# G15A_CHH 0.2352941 0.13174956 0.8427399          54865          34217             39234
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1724138 0.08721289 0.3744809          50774          19448              9179
# G23A_CHH 0.1785714 0.08862225 0.3757187          47613          26433             28997
# G24A_CHH 0.1950000 0.10699484 0.4114233          34238          25099             29772
# G13A_CHH 0.2500000 0.13196276 0.6610193          40056          12837              7231
# G14A_CHH 0.1907692 0.11961368 1.0627798          60325          35019             38039
# G15A_CHH 0.2500000 0.13668481 0.8028904          44078          25923             31692
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1800000 0.09307513 0.3825358          39217          15202              6935
# G23A_CHH 0.1889401 0.09484867 0.3758664          35635          19873             20390
# G24A_CHH 0.2028986 0.11228535 0.3964220          24794          17904             19938
# G13A_CHH 0.2631579 0.14058003 0.6702746          29421           9896              5237
# G14A_CHH 0.2000000 0.12570978 1.0813262          47392          27535             30895
# G15A_CHH 0.2653061 0.14560256 0.7768702          30957          19292             22117
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1696429 0.08522688 0.3667863          59499          21887             10396
# G23A_CHH 0.1764706 0.08721077 0.3717428          55824          31290             35157
# G24A_CHH 0.1904762 0.10574493 0.4359485          43923          32565             38890
# G13A_CHH 0.2490842 0.12944654 0.6540627          47118          15259              8674
# G14A_CHH 0.1875000 0.11643146 1.0538785          72009          40441             45475
# G15A_CHH 0.2428571 0.13466023 0.8047050          52522          33016             38646

date()
# [1] "Thu Jan  9 16:45:15 2020"


names(G14A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

names(G14A_CHH_grl[2:20])
# [1] "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

date()


## 20190109: this code was NOT run..instead, the 4 code blocks shown farther below 
#  ..was run concurrently in 4 R consoles simultaneously.  That is why code is 
#  ..shown farther down that uses load() to get those objects available here.
#  ..this code works but was not run and will take a long time if it does run. 
#
# for(item in names(G14A_CHH_grl[2:20])) {
for(item in names(G14A_CHH_grl)) {
	    assign(
		           paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item),
			           gofReport(
					                 eval(str2expression(
									                     paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr", item)
											                 )),
					                 model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
							             column = 9,
								                 absolute = FALSE,
										             output = c("best.model"),
											                 num.cores = 64L,
													             verbose = FALSE
							         )
			       )
    print(eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
}

################################################################################
################################################################################
################################################################################
# START multi-R console to cut runtime by 75%:
################################################################################
################################################################################
################################################################################

# Here we see 4 blocks of code that process 5 of the 20 soybean chromoaomes. 
# These code blocks can be run simultaneously (lets hope!) and cut the time by 4x.



#
# library(BiocManager)
# library(MethylIT)
# library(MethylIT.utils)
# 
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13.RData")
# 
# for(item in c("1","10","11","12","13")) {
#   assign(
#     paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item),
#     gofReport(
#       eval(str2expression(
#         paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr", item)
#       )),
#       model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
#       column = 9,
#       absolute = FALSE,
#       output = c("best.model"),
#       num.cores = 64L,
#       verbose = FALSE
#     )
#   )
#   print(eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
# }
# 
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
# 
# 
# library(BiocManager)
# library(MethylIT)
# library(MethylIT.utils)
# 
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18.RData")
# 
# for(item in c("14","15","16","17","18")) {
#   assign(
#     paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item),
#     gofReport(
#       eval(str2expression(
#         paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr", item)
#       )),
#       model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
#       column = 9,
#       absolute = FALSE,
#       output = c("best.model"),
#       num.cores = 64L,
#       verbose = FALSE
#     )
#   )
#   print(eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
# }
# 
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
# 
# 
# library(BiocManager)
# library(MethylIT)
# library(MethylIT.utils)
# 
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4.RData")
# 
# for(item in c("19","2","20","3","4")) {
#   assign(
#     paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item),
#     gofReport(
#       eval(str2expression(
#         paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr", item)
#       )),
#       model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
#       column = 9,
#       absolute = FALSE,
#       output = c("best.model"),
#       num.cores = 64L,
#       verbose = FALSE
#     )
#   )
#   print(eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
# }
# 
# 
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
# 
# 
# library(BiocManager)
# library(MethylIT)
# library(MethylIT.utils)
# 
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9.RData")
# 
# for(item in c("5","6","7","8","9")) {
#   assign(
#     paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item),
#     gofReport(
#       eval(str2expression(
#         paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr", item)
#       )),
#       model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
#       column = 9,
#       absolute = FALSE,
#       output = c("best.model"),
#       num.cores = 64L,
#       verbose = FALSE
#     )
#   )
#   print(eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
# }
# 
# 
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")
# 
################################################################################
################################################################################
################################################################################
# END multi-R console to cut runtime by 75%:
################################################################################
################################################################################
################################################################################

load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")


date()
# [1] "Thu Jan  9 16:46:20 2020"

ls()
#  [1] "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1"                        "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10"                      
#  [3] "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11"                       "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12"                      
#  [5] "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13"                       "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14"                      
#  [7] "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15"                       "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16"                      
#  [9] "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17"                       "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18"                      
# [11] "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19"                       "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2"                       
# [13] "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20"                       "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3"                       
# [15] "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4"                        "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5"                       
# [17] "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6"                        "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7"                       
# [19] "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8"                        "bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9"                       
# [21] "critical.val_G222324refctrl_vs_G131415_sum_chr1"                        "critical.val_G222324refctrl_vs_G131415_sum_chr10"                      
# [23] "critical.val_G222324refctrl_vs_G131415_sum_chr11"                       "critical.val_G222324refctrl_vs_G131415_sum_chr12"                      
# [25] "critical.val_G222324refctrl_vs_G131415_sum_chr13"                       "critical.val_G222324refctrl_vs_G131415_sum_chr14"                      
# [27] "critical.val_G222324refctrl_vs_G131415_sum_chr15"                       "critical.val_G222324refctrl_vs_G131415_sum_chr16"                      
# [29] "critical.val_G222324refctrl_vs_G131415_sum_chr17"                       "critical.val_G222324refctrl_vs_G131415_sum_chr18"                      
# [31] "critical.val_G222324refctrl_vs_G131415_sum_chr19"                       "critical.val_G222324refctrl_vs_G131415_sum_chr2"                       
# [33] "critical.val_G222324refctrl_vs_G131415_sum_chr20"                       "critical.val_G222324refctrl_vs_G131415_sum_chr3"                       
# [35] "critical.val_G222324refctrl_vs_G131415_sum_chr4"                        "critical.val_G222324refctrl_vs_G131415_sum_chr5"                       
# [37] "critical.val_G222324refctrl_vs_G131415_sum_chr6"                        "critical.val_G222324refctrl_vs_G131415_sum_chr7"                       
# [39] "critical.val_G222324refctrl_vs_G131415_sum_chr8"                        "critical.val_G222324refctrl_vs_G131415_sum_chr9"                       
# [41] "F4_P37_CHH"                                                             "F6_P37_CHH"                                                            
# [43] "G13A_CHH_grl"                                                           "G14A_CHH_grl"                                                          
# [45] "G15A_CHH_grl"                                                           "G22A_CHH_grl"                                                          
# [47] "G23A_CHH_grl"                                                           "G24A_CHH_grl"                                                          
# [49] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1"  "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10"
# [51] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12"
# [53] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14"
# [55] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16"
# [57] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18"
# [59] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2" 
# [61] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3" 
# [63] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4"  "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5" 
# [65] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6"  "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7" 
# [67] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8"  "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9" 
# [69] "item"                                                                   "soy_WT_G222324_Ref_CHH_sum_chr1"                                       
# [71] "soy_WT_G222324_Ref_CHH_sum_chr10"                                       "soy_WT_G222324_Ref_CHH_sum_chr11"                                      
# [73] "soy_WT_G222324_Ref_CHH_sum_chr12"                                       "soy_WT_G222324_Ref_CHH_sum_chr13"                                      
# [75] "soy_WT_G222324_Ref_CHH_sum_chr14"                                       "soy_WT_G222324_Ref_CHH_sum_chr15"                                      
# [77] "soy_WT_G222324_Ref_CHH_sum_chr16"                                       "soy_WT_G222324_Ref_CHH_sum_chr17"                                      
# [79] "soy_WT_G222324_Ref_CHH_sum_chr18"                                       "soy_WT_G222324_Ref_CHH_sum_chr19"                                      
# [81] "soy_WT_G222324_Ref_CHH_sum_chr2"                                        "soy_WT_G222324_Ref_CHH_sum_chr20"                                      
# [83] "soy_WT_G222324_Ref_CHH_sum_chr3"                                        "soy_WT_G222324_Ref_CHH_sum_chr4"                                       
# [85] "soy_WT_G222324_Ref_CHH_sum_chr5"                                        "soy_WT_G222324_Ref_CHH_sum_chr6"                                       
# [87] "soy_WT_G222324_Ref_CHH_sum_chr7"                                        "soy_WT_G222324_Ref_CHH_sum_chr8"                                       
# [89] "soy_WT_G222324_Ref_CHH_sum_chr9"                                        "soybean_gff3"      

## 20200109: these were saved in separate R session:
#
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")

## 20190109 saved:
#
# save(critical.val_G222324refctrl_vs_G131415_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr1.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr10.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr11.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr12.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr13.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr14.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr15.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr16.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr17.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr18.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr19.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr2.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr20.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr3.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr4.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr5.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr6.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr7.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr8.RData")
# save(critical.val_G222324refctrl_vs_G131415_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr9.RData")

## 20200109: these were saved in separate R session:
#
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9.RData")

















































date()

covr <-
	    lapply(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1, function(x) {
		           cov <- x$c2 + x$t2
			           return(cov)
			       })

do.call(rbind, lapply(covr, function(x) {
		          q60 <- quantile(x, 0.6)
			      q9999 <- quantile(x, 0.9999)
			      idx1 <- which(x >= q60)
			          idx2 <- which(x <= 500)
			          q95 <- quantile(x, 0.95)
				      idx <- intersect(idx1, idx2)
				      return(c(
					               round(summary(x)),
						               q60,
							               quantile(x, c(0.95, 0.99, 0.999, 0.9999)),
								               num.siteGreater_8 = sum(x >= 8),
									               q60_to_500 = sum((x >= q60) & (x <= 500)),
										               num.siteGreater_500 = sum(x > 500)
						           ))}))

## 20190109: we see:
#   - num.siteGreater_500 are all 0, confirming that our "high.coverage=300" filter in estimateDivergence worked
#   - 1st quartile values of about 10, showing that our data has good coverage, >10, for 75% of sites
#   - our data is extremely similar for all 6 samples in terms of coverage
# See also soybean_G222324refctrl_vs_G1415_sepChroms_20200102.R for more discussion.
#          Min. 1st Qu. Median Mean 3rd Qu. Max. 60% 95% 99% 99.9% 99.99% num.siteGreater_8 q60_to_500 num.siteGreater_500
# G22A_CHH    0      12     20   21      29   67  23  41  48    54     59           1378919     645309                   0
# G23A_CHH    0      11     17   18      24   53  19  32  37    42     46           1256964     625871                   0
# G24A_CHH    0      10     13   13      16   37  14  21  24    28     32            850429     413332                   0
# G13A_CHH    0       9     13   14      18   55  15  25  31    37     42            957931     509097                   0
# G14A_CHH    0      21     31   33      43  122  36  61  73    86     97           1814093     749869                   0
# G15A_CHH    0      10     15   15      20   60  17  27  32    39     43           1101701     514599                   0

length(mcols(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1$G13A_CHH)$hdiv)
# [1] 1170149

critical.val_G222324refctrl_vs_G131415_sum_chr1["G13A_CHH",][4]
# num.sites.hd95 
#          58507 

58507/1170149
# [1] 0.04999962

critical.val_G222324refctrl_vs_G131415_sum_chr1
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1702128 0.08572858 0.3633274          77216          28090             13260
# G23A_CHH 0.1775873 0.08798573 0.3667813          71467          40638             45736
# G24A_CHH 0.1958042 0.10661883 0.3897832          48326          35543             42428
# G13A_CHH 0.2500000 0.13016440 0.6369731          58507          18344             10921
# G14A_CHH 0.1875862 0.11691683 1.0376375          93064          53812             59550
# G15A_CHH 0.2500000 0.13731769 0.7895883          64005          39620             48166

critical.val_G222324refctrl_vs_G131415_sum_chr10
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1666667 0.08508452 0.3842621          63988          22807             10078
# G23A_CHH 0.1700000 0.08571900 0.3877397          61553          34482             37290
# G24A_CHH 0.1851852 0.10479647 0.4567538          48466          36212             41919
# G13A_CHH 0.2352941 0.12836617 0.7039856          55242          17529              8897
# G14A_CHH 0.1858696 0.11640307 1.0663130          74443          42940             47361
# G15A_CHH 0.2361111 0.13346783 0.8398077          58101          37511             42743

critical.val_G222324refctrl_vs_G131415_sum_chr13
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHH 0.1764706 0.09035498 0.3958731          38686          15180              7069
# G23A_CHH 0.1818182 0.09043757 0.3917350          36779          20187             21230
# G24A_CHH 0.1923077 0.10710433 0.4377925          28233          20397             23845
# G13A_CHH 0.2500000 0.13586566 0.7227406          32575          10428              5459
# G14A_CHH 0.1971417 0.12353670 1.0984441          45466          26740             29541
# G15A_CHH 0.2500000 0.13861246 0.8431039          35280          21030             24806


date()
# [1] "Thu Jan  9 16:49:35 2020"

## 20200109:
#
# We call getPotentialDIMP with tv.cut values shown here:
#   - tv.col= 7L,
#   - tv.cut=0.2,
# ..normally, per a best practice, we would use a calculated value like this for chr1:
#
max(critical.val_G222324refctrl_vs_G131415_sum_chr1[,"tv.95%"])
# [1] 0.25

#
# ..but that will cause too many DIMPs to be filtered out.  Therefore, the best 
# ..practice may be to use a hardcoded 0.2 in cases where tv.95% values are 
# ..lower than, say 0.25. 
#
# UPDATE 20200109, we will use 0.2 anyway, even though maybe 0.25 is better
#

for(item in names(G14A_CHH_grl)) {
	    assign(
		           paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr", item),
			           getPotentialDIMP(
						                eval(str2expression(
										                    paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr", item)
												                )),
						                nlms = eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$nlms"))),
								            div.col = 9L,
									                tv.col= 7L,
											            tv.cut=0.2,
												                dist.name = eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel")))
								        )
			       )
    #print(eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
}

## 20190109 saved:
#
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr10.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr11.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr12.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr13.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr14.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr15.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr16.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr17.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr18.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr19.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr2.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr20.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr3.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr4.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr5.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr6.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr7.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr8.RData")
save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr9.RData")

head(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1$G14A_CHH,n=2)
# GRanges object with 2 ranges and 10 metadata columns:
#       seqnames    ranges strand |        c1        t1        c2        t2                p1                 p2                TV             bay.TV             hdiv
#          <Rle> <IRanges>  <Rle> | <numeric> <numeric> <numeric> <numeric>         <numeric>          <numeric>         <numeric>          <numeric>        <numeric>
#   [1]        1       419      - |        11        26         2        24 0.247992098311154 0.0956613892055564 -0.22037422037422 -0.152330709105597 1.34566796501388
#   [2]        1      1241      + |         4        28        10        18 0.126116237506308  0.295236553029761 0.232142857142857  0.169120315523454 1.37421937340063
#                      wprob
#                  <numeric>
#   [1] 0.000566233109272829
#   [2] 0.000505528132568407
#   -------
#   seqinfo: 20 sequences from an unspecified genome; no seqlengths

date()
# [1] "Thu Jan  9 16:57:22 2020"

names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1)
# [1] "G22A_CHH" "G23A_CHH" "G24A_CHH" "G13A_CHH" "G14A_CHH" "G15A_CHH"

ls(pattern="^PS")
#  [1] "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1"  "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr10" "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr11"
#  [4] "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr12" "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr13" "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr14"
#  [7] "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr15" "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr16" "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr17"
# [10] "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr18" "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr19" "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr2" 
# [13] "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr20" "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr3"  "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr4" 
# [16] "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr5"  "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr6"  "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr7" 
# [19] "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr8"  "PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr9" 

mcols(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1$G14A_CHH)
# DataFrame with 52403 rows and 10 columns
#              c1        t1        c2        t2                 p1                 p2                 TV             bay.TV              hdiv                wprob
#       <numeric> <numeric> <numeric> <numeric>          <numeric>          <numeric>          <numeric>          <numeric>         <numeric>            <numeric>
# 1            11        26         2        24  0.247992098311154 0.0956613892055564  -0.22037422037422 -0.152330709105597  1.34566796501388 0.000566233109272829
# 2             4        28        10        18  0.126116237506308  0.295236553029761  0.232142857142857  0.169120315523454  1.37421937340063 0.000505528132568407
# 3            10        24         2        29  0.242794864056033 0.0843073253680471 -0.229601518026565 -0.158487538687986   1.6209580070877 0.000194168874521359
# 4            12        59        14        21  0.161838613161413  0.337150511778247  0.230985915492958  0.175311898616834  2.02000753842289 4.44510535793226e-05
# 5            14        13         1        13  0.377658987690913  0.101549064887709 -0.447089947089947 -0.276109922803205  2.20382939337605  2.3106223118918e-05
# ...         ...       ...       ...       ...                ...                ...                ...                ...               ...                  ...
# 52399         6        19         0        12  0.197718428151061 0.0670906754323647              -0.24 -0.130627752718696 0.682671705744405  0.00984072845790378
# 52400         1        34         5        16 0.0589008385415779  0.203930487209133   0.20952380952381  0.145029648667555  1.35733260592977 0.000540558529992592
# 52401         1        19         3         7 0.0839568303710847  0.215443770709155               0.25   0.13148694033807 0.512650806543467   0.0227142972896739
# 52402         1        21         3         7 0.0794504758000514  0.215443770709155  0.254545454545455  0.135993294909103 0.575371499929727   0.0165583552745657
# 52403         1        20         3         6 0.0816415162325848   0.22614825452633  0.285714285714286  0.144506738293745 0.580432540272651   0.0161480276196027

treatment_names_CHH <- names(F4_P37_CHH)

control_names_CHH <- names(F6_P37_CHH)

control_names_CHH
# [1] "G22A_CHH" "G23A_CHH" "G24A_CHH"

treatment_names_CHH
# [1] "G13A_CHH" "G14A_CHH" "G15A_CHH"

names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr2) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr3) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr4) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr5) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr6) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr7) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr8) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr9) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr10) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr11) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr12) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr13) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr14) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr15) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr16) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr17) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr18) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr19) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr20) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")


















date()
# [1] "Thu Jan  9 16:58:18 2020"

#
# cutpoint with machine learning:
#
for(item in names(G14A_CHH_grl)) {
	    assign(
		           paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item),
			           estimateCutPoint(
						                eval(str2expression(
										                    paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
												                )),
						                simple = FALSE,
								            control.names = control_names_CHH,
									                treatment.names = treatment_names_CHH,
											            column = c(hdiv = TRUE, bay.TV = TRUE, wprob = TRUE, pos = TRUE),
												                div.col = 9, clas.perf = TRUE,
														            classifier1 = "pca.qda", n.pc = 4, 
															                classifier2 = "pca.lda",
																	            center = TRUE, scale = TRUE,
																		                verbose = TRUE)
			           )

}
## 20190109: here we see that:
#  - 40 calls are made because we have:
#  - 20 samples
#  - 2 classifiers
#
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos

cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr1$cutpoint
# 0.3651608
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr2$cutpoint
# 0.378848
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr3$cutpoint
# 0.3507506
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr4$cutpoint
# 0.3741857
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr5$cutpoint
# 0.3831403
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr6$cutpoint
# 0.3851075
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr7$cutpoint
# 0.3768692
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr8$cutpoint
# 0.3749555
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr9$cutpoint
# 0.3686761
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr10$cutpoint
# 0.3850044
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr11$cutpoint
# 0.3652215
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr12$cutpoint
# 0.3738825
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr13$cutpoint
# 0.3883132
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr14$cutpoint
# 0.3191474
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr15$cutpoint
# 0.4222549
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr16$cutpoint
# 0.3585283
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr17$cutpoint
# 0.3954204
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr18$cutpoint
# 0.3948675
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr19$cutpoint
# 0.3726884
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr20$cutpoint
# 0.3830666



for(item in names(G14A_CHH_grl)) {
	    print(paste0("chr",item))
    print(eval(str2expression(
			              paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$testSetPerformance$table")
				          )))
        print(eval(str2expression(
				          paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$testSetPerformance$overall")
					      )))
}
# [1] "chr1"
#           Reference
# Prediction    CT    TT
#         CT     0     3
#         TT  5501 35538
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   0.8658934750  -0.0001461334   0.8625586139   0.8691764792   0.8659665708   0.5209241158   0.0000000000 
# [1] "chr10"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  4376 27556
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8629588      0.0000000      0.8591382      0.8667132      0.8629588      0.5040313      0.0000000 
# [1] "chr11"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  2758 16785
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8588753      0.0000000      0.8539152      0.8637284      0.8588753      0.5050787      0.0000000 
# [1] "chr12"
#           Reference
# Prediction    CT    TT
#         CT     0     1
#         TT  3539 21550
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   8.589079e-01  -7.969687e-05   8.545390e-01   8.631936e-01   8.589478e-01   5.117171e-01   0.000000e+00 
# [1] "chr13"
#           Reference
# Prediction    CT    TT
#         CT     0     1
#         TT  3058 19553
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   8.647178e-01  -8.842752e-05   8.601909e-01   8.691508e-01   8.647621e-01   5.125771e-01   0.000000e+00 
# [1] "chr14"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  4652 28642
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8602751      0.0000000      0.8565046      0.8639827      0.8602751      0.5039103      0.0000000 
# [1] "chr15"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  4080 25165
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8604890      0.0000000      0.8564651      0.8644411      0.8604890      0.5041754      0.0000000 
# [1] "chr16"
#           Reference
# Prediction    CT    TT
#         CT     0     1
#         TT  3264 19927
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   8.592187e-01  -8.621765e-05   8.546760e-01   8.636712e-01   8.592618e-01   5.121990e-01   0.000000e+00 
# [1] "chr17"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  3658 21488
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8545295      0.0000000      0.8501109      0.8588661      0.8545295      0.5044108      0.0000000 
# [1] "chr18"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  5112 31809
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8615422      0.0000000      0.8579773      0.8650502      0.8615422      0.5037301      0.0000000 
# [1] "chr19"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  4625 28348
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8597337      0.0000000      0.8559386      0.8634653      0.8597337      0.5039218      0.0000000 
# [1] "chr2"
#           Reference
# Prediction    CT    TT
#         CT     0     1
#         TT  3918 24497
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   8.620847e-01  -7.036988e-05   8.580210e-01   8.660743e-01   8.621199e-01   5.111229e-01   0.000000e+00 
# [1] "chr20"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  4406 26889
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8592107      0.0000000      0.8553080      0.8630467      0.8592107      0.5040182      0.0000000 
# [1] "chr3"
#           Reference
# Prediction    CT    TT
#         CT     0     4
#         TT  3801 23982
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   0.8630654623  -0.0002876845   0.8589674162   0.8670874324   0.8632094145   0.5321441465   0.0000000000 
# [1] "chr4"
#           Reference
# Prediction    CT    TT
#         CT     0     1
#         TT  4542 28777
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   0.8636554622  -0.0000600144   0.8599242009   0.8673232084   0.8636854742   0.5103248619   0.0000000000 
# [1] "chr5"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  3541 21749
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8599842      0.0000000      0.8556464      0.8642391      0.8599842      0.5044820      0.0000000 
# [1] "chr6"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  4214 25332
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8573749      0.0000000      0.8533356      0.8613439      0.8573749      0.5041091      0.0000000 
# [1] "chr7"
#           Reference
# Prediction    CT    TT
#         CT     0     2
#         TT  3826 23920
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8620441     -0.0001441      0.8579306      0.8660817      0.8621162      0.5181963      0.0000000 
# [1] "chr8"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  3303 20919
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.8636364      0.0000000      0.8592504      0.8679349      0.8636364      0.5046400      0.0000000 
# [1] "chr9"
#           Reference
# Prediction    CT    TT
#         CT     0     1
#         TT  4480 27179
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   8.584649e-01  -6.316109e-05   8.545767e-01   8.622873e-01   8.584965e-01   5.104163e-01   0.000000e+00 


# 
# cutpoint simple
#
for(item in names(G14A_CHH_grl)) {
	    assign(
		           paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item),
			           estimateCutPoint(
						                eval(str2expression(
										                    paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
												                )),
						                simple = TRUE,
								            control.names = control_names_CHH,
									                treatment.names = treatment_names_CHH,
											            column = c(hdiv = TRUE, bay.TV = TRUE, wprob = TRUE, pos = TRUE),
												                div.col = 9,
														            clas.perf = TRUE,
															                classifier1 = "qda",
																	            prop = 0.6,
																		                verbose = TRUE)
			           )

}
## 20190105: here we see that:
#  - 20 calls are made because we have:
#  - 20 samples
#  - 1 classifiers
#
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos

cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr1$cutpoint
# 0.7381446
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr2$cutpoint
# 0.7648554
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr3$cutpoint
# 0.7105012
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr4$cutpoint
# 0.720265
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr5$cutpoint
# 0.7593088
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr6$cutpoint
# 0.7932721
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr7$cutpoint
# 0.7421807
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr8$cutpoint
# 0.7555899
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr9$cutpoint
# 0.732585
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr10$cutpoint
# 0.7911854
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr11$cutpoint
# 0.7811807
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr12$cutpoint
# 0.7528698
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr13$cutpoint
# 0.7574951
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr14$cutpoint
# 0.6855358
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr15$cutpoint
# 0.8216133
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr16$cutpoint
# 0.7595319
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr17$cutpoint
# 0.785654
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr18$cutpoint
# 0.7876163
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr19$cutpoint
# 0.7356416
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr20$cutpoint
# 0.7852044


for(item in names(G14A_CHH_grl)) {
	    print(paste0("chr",item))
    print(eval(str2expression(
			              paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$testSetPerformance$table")
				          )))
        print(eval(str2expression(
				          paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$testSetPerformance$overall")
					      )))
}
# [1] "chr1"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  1412 20163
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.345539e-01   0.000000e+00   9.311717e-01   9.378184e-01   9.345539e-01   5.070813e-01  1.410131e-308 
# [1] "chr10"
#           Reference
# Prediction    CT    TT
#         CT   467  1809
#         TT   685 13804
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   8.512377e-01   1.994122e-01   8.457603e-01   8.565929e-01   9.312854e-01   1.000000e+00  5.563545e-112 
# [1] "chr11"
#           Reference
# Prediction   CT   TT
#         CT    0    0
#         TT  649 8959
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.324521e-01   0.000000e+00   9.272495e-01   9.373907e-01   9.324521e-01   5.104442e-01  1.002811e-142 
# [1] "chr12"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  1003 12284
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.245127e-01   0.000000e+00   9.198895e-01   9.289489e-01   9.245127e-01   5.084033e-01  1.088141e-219 
# [1] "chr13"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT   903 11428
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.267699e-01   0.000000e+00   9.220298e-01   9.313075e-01   9.267699e-01   5.088558e-01  5.945552e-198 
# [1] "chr14"
#           Reference
# Prediction    CT    TT
#         CT   266     5
#         TT   897 16789
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.497689e-01   3.552062e-01   9.464728e-01   9.529184e-01   9.352342e-01   1.217013e-16  2.043958e-193 
# [1] "chr15"
#           Reference
# Prediction    CT    TT
#         CT   504  2044
#         TT   722 12217
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   8.213986e-01   1.793708e-01   8.152738e-01   8.274024e-01   9.208368e-01   1.000000e+00  3.201201e-139 
# [1] "chr16"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT   796 11102
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.330980e-01   0.000000e+00   9.284596e-01   9.375231e-01   9.330980e-01   5.094309e-01  1.087182e-174 
# [1] "chr17"
#           Reference
# Prediction   CT   TT
#         CT  886 4622
#         TT  228 8128
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.6501731      0.1545872      0.6421676      0.6581156      0.9196480      1.0000000      0.0000000 
# [1] "chr18"
#           Reference
# Prediction    CT    TT
#         CT   135   408
#         TT  1251 17614
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.145198e-01   1.039429e-01   9.104981e-01   9.184168e-01   9.285862e-01   1.000000e+00   6.150784e-95 
# [1] "chr19"
#           Reference
# Prediction    CT    TT
#         CT   508  2026
#         TT   832 14843
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   8.430446e-01   1.836725e-01   8.376790e-01   8.483005e-01   9.264100e-01   1.000000e+00  2.604931e-110 
# [1] "chr2"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  1074 13850
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.280354e-01   0.000000e+00   9.237717e-01   9.321313e-01   9.280354e-01   5.081202e-01  4.021553e-235 
# [1] "chr20"
#           Reference
# Prediction    CT    TT
#         CT   494  1799
#         TT   700 13585
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   8.492581e-01   2.083509e-01   8.437202e-01   8.546731e-01   9.279768e-01   1.000000e+00  6.309408e-107 
# [1] "chr3"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT   921 13556
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.363818e-01   0.000000e+00   9.322830e-01   9.403042e-01   9.363818e-01   5.087672e-01  7.265417e-202 
# [1] "chr4"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  1265 16794
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.299518e-01   0.000000e+00   9.261320e-01   9.336325e-01   9.299518e-01   5.074820e-01  1.240961e-276 
# [1] "chr5"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  1018 12646
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.254977e-01   0.000000e+00   9.209670e-01   9.298461e-01   9.254977e-01   5.083410e-01  5.973880e-223 
# [1] "chr6"
#           Reference
# Prediction    CT    TT
#         CT   478  1858
#         TT   632 12270
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   8.365927e-01   1.982411e-01   8.306253e-01   8.424314e-01   9.271558e-01   1.000000e+00  4.416622e-133 
# [1] "chr7"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  1037 13778
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.300034e-01   0.000000e+00   9.257761e-01   9.340608e-01   9.300034e-01   5.082635e-01  4.430445e-227 
# [1] "chr8"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT   856 11671
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.316676e-01   0.000000e+00   9.271081e-01   9.360253e-01   9.316676e-01   5.090947e-01  9.810843e-188 
# [1] "chr9"
#           Reference
# Prediction    CT    TT
#         CT     0     0
#         TT  1271 15711
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.251560e-01   0.000000e+00   9.210957e-01   9.290701e-01   9.251560e-01   5.074650e-01  6.163791e-278 


for(item in names(G14A_CHH_grl)) {
	    assign(
		           paste0("DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr", item),
			           selectDIMP(
					                  eval(str2expression(
									                      paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
											                  )),
					                  div.col = 9,
							              cutpoint = eval(str2expression(
												                     paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item, '$cutpoint')
														                 )),
					                  )
			           )

}

ls(pattern="DIMPsYI")
#  [1] "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr1"  "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr10" "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr11"
#  [4] "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr12" "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr13" "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr14"
#  [7] "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr15" "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr16" "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr17"
# [10] "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr18" "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr19" "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr2" 
# [13] "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr20" "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr3"  "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr4" 
# [16] "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr5"  "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr6"  "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr7" 
# [19] "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr8"  "DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr9" 

unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr1, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    3546     1534     2418    25677    30986    34859 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr2, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2419     1075     2195    19523    20692    23703 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr3, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2251     1099     1324    15981    21785    21643 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr4, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2781     1519     1868    19670    26698    26632 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr5, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2147     1010     2367    18045    18952    22250 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr6, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2350     1053     2513    20235    21234    23954 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr7, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2393     1099     1947    18178    21414    23012 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr8, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1929     1003     1305    14278    18881    17498 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr9, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2697     1314     2815    20907    24106    27179 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr10, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2614     1165     2485    22175    23044    26757 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr11, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1466      652     1182    11370    14049    14343 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr12, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2065     1012     2266    17037    18404    20703 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr13, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1968      896     1731    15469    17454    18254 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr14, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2976     1458     1731    20303    26976    28481 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr15, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2263     1137     3060    20025    20896    23122 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr16, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1792      814     1722    15195    17101    18638 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr17, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2183      994     2996    19569    18213    21690 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr18, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2697     1370     3289    25355    27075    29863 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr19, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    3019     1306     2980    23051    25296    29431 
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr20, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2660     1127     2611    22659    22800    27134 

date()
# [1] "Thu Jan  9 17:19:16 2020"

DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- data.frame("G22A_CHH"=0, "G23A_CHH"=0, "G24A_CHH"=0, "G13A_CHH"=0, "G14A_CHH"=0, "G15A_CHH"=0)
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr1, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr2, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr3, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr4, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr5, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr6, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr7, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr8, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr9, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr10, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr11, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr12, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr13, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr14, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr15, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr16, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr17, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr18, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr19, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr20, length)))

colSums(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS)
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    48216    22637    44805   384702   436056   479146

for(item in names(G14A_CHH_grl)) {
	    assign(
		           paste0("DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr", item),
			           selectDIMP(
					                  eval(str2expression(
									                      paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
											                  )),
					                  div.col = 9,
							              cutpoint = eval(str2expression(
												                     paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item, '$cutpoint')
														                 )),
					                  )
			           )

}

ls(pattern="DIMPsML")
#  [1] "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr1"  "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr10" "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr11"
#  [4] "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr12" "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr13" "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr14"
#  [7] "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr15" "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr16" "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr17"
# [10] "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr18" "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr19" "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr2" 
# [13] "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr20" "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr3"  "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr4" 
# [16] "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr5"  "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr6"  "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr7" 
# [19] "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr8"  "DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr9" 


unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr1, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    9774     8068    12663    50217    52403    63489 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr2, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6801     5546     9456    36125    35536    43512 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr3, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6495     5752     7737    32186    36336    40675 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr4, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    7591     6673     9243    37998    43622    48040 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr5, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    5890     4745     9231    32588    31813    39414 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr6, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6920     5766    10710    36947    37518    45504 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr7, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6561     5428     8802    34346    35524    41596 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr8, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    5653     4815     6840    28246    32108    33784 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr9, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    7490     6113    11324    39262    40134    48546 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr10, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    7394     6113    11257    39951    40271    50201 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr11, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    4691     4001     6354    23901    24999    28611 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr12, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    5843     4747     8823    31558    31484    37862 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr13, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    5387     4278     6904    28007    29044    32697 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr14, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    7919     7887     9531    39266    43428    50413 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr15, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6435     5450    10339    34994    37139    42728 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr16, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    5410     4673     8039    29421    29344    35222 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr17, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    5888     4923     9887    31899    31024    38830 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr18, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    8178     7088    12723    44912    47637    55604 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr19, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    7906     6136    11962    40496    41575    52049 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr20, length))
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    7423     5809    11396    39992    39416    50179


date()
# [1] "Thu Jan  9 17:21:01 2020"

DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- data.frame("G22A_CHH"=0, "G23A_CHH"=0, "G24A_CHH"=0, "G13A_CHH"=0, "G14A_CHH"=0, "G15A_CHH"=0)
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr1, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr2, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr3, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr4, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr5, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr6, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr7, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr8, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr9, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr10, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr11, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr12, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr13, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr14, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr15, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr16, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr17, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr18, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr19, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr20, length)))

colSums(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS)
# G22A_CHH G23A_CHH G24A_CHH G13A_CHH G14A_CHH G15A_CHH 
#   135649   114011   193221   712312   740355   878956 

pryr::mem_used()
# 47.9 GB

length(ls(pattern="best"))
# [1] 20
length(ls(pattern="critical"))
# [1] 20
length(ls(pattern="cutpointsML"))
# [1] 20
length(ls(pattern="cutpointsYI"))
# [1] 20
length(ls(pattern="DIMPsML"))
# [1] 21
length(ls(pattern="DIMPsYI"))
# [1] 21
length(ls(pattern="HD"))
# [1] 20
length(ls(pattern="PS"))
# [1] 60
length(ls(pattern="^PS_"))
# [1] 20

## 20190109:
# #
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")
save(control_names_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/control_names_CHH.RData")
save(treatment_names_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/treatment_names_CHH.RData")
save(covr, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/covr.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr1.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr10.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr11.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr12.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr13.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr14.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr15.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr16.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr17.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr18.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr19.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr2.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr20.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr3.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr4.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr5.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr6.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr7.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr8.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr9.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS.RData")
save(F4_P37_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/F4_P37_CHH.RData")
save(F6_P37_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/F6_P37_CHH.RData")


# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9.RData")
# save(item, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/item.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr10.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr11.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr12.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr13.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr14.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr15.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr16.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr17.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr18.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr19.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr2.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr20.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr3.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr4.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr5.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr6.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr7.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr8.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr9.RData")



date()







# soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200109.R



.libPaths()

library(BiocManager)
library(MethylIT)
library(MethylIT.utils)

rm(list=ls())
length(ls()) 
# [1] 0

# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F4_P37_CG_4_29_2019.RData")
# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F6_P37_CG_4_29_2019.RData")

# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F4_P37_CHH_4_29_2019.RData")
# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F6_P37_CHH_4_29_2019.RData")

# TODO 20190109

load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F4_P37_CHH_4_29_2019.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F6_P37_CHH")

names(F4_P37_CHH)
# [1] "G13A_CHH" "G14A_CHH" "G15A_CHH"

names(F6_P37_CHH)
# [1] "G22A_CHH" "G23A_CHH" "G24A_CHH"

ls()
# [1] "F4_P37_CHH" "F6_P37_CHH"

seqnames(F4_P37_CHH$G14A_CHH)


G22A_CHH_grl <- split(F6_P37_CHH$G22A_CHH,seqnames(F4_P37_CHH$G14A_CHH))
G23A_CHH_grl <- split(F6_P37_CHH$G23A_CHH,seqnames(F4_P37_CHH$G14A_CHH))
G24A_CHH_grl <- split(F6_P37_CHH$G24A_CHH,seqnames(F4_P37_CHH$G14A_CHH))

G13A_CHH_grl <- split(F4_P37_CHH$G13A_CHH,seqnames(F4_P37_CHH$G13A_CHH))
G14A_CHH_grl <- split(F4_P37_CHH$G14A_CHH,seqnames(F4_P37_CHH$G13A_CHH))
G15A_CHH_grl <- split(F4_P37_CHH$G15A_CHH,seqnames(F4_P37_CHH$G13A_CHH))

names(G22A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

names(G14A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 


date()

for(item in names(G14A_CHH_grl)) {
    assign(
        paste0("soy_WT_G222324_Ref_CHH_sum_chr", item),
        poolFromGRlist(
            list(G22A_CHH_grl[[item]], G23A_CHH_grl[[item]], G24A_CHH_grl[[item]]),
            stat = "sum",
            num.cores = 12L
        ))}

date()

ls()


length(ls(pattern="soy_WT_G222324_Ref_CHH_sum_chr"))
# [1] 20

## 20190109:
#
save(F4_P37_CHH,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/F4_P37_CHH.RData")
save(F6_P37_CHH,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/F6_P37_CHH.RData")
save(G13A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G13A_CHH_grl.RData")
save(G14A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G14A_CHH_grl.RData")
save(G15A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G15A_CHH_grl.RData")
save(G22A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G22A_CHH_grl.RData")
save(G23A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G23A_CHH_grl.RData")
save(G24A_CHH_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/G24A_CHH_grl.RData")
save(item,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/item.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr1.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr10.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr11.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr12.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr13.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr14.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr15.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr16.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr17.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr18.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr19.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr2.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr20.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr3.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr4.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr5.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr6.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr7.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr8.RData")
save(soy_WT_G222324_Ref_CHH_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/soy_WT_G222324_Ref_CHH_sum_chr9.RData")

date()

## 20190107: 61 minute runtime:
#
for(item in names(G14A_CHH_grl)) {
    assign(
        paste0(
            "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr",
            item
        ),
        estimateDivergence(
            ref = eval(str2expression(paste0("soy_WT_G222324_Ref_CHH_sum_chr",item))),
            indiv = list(G22A_CHH_grl[[item]], G23A_CHH_grl[[item]], G24A_CHH_grl[[item]],G13A_CHH_grl[[item]], G14A_CHH_grl[[item]], G15A_CHH_grl[[item]]),
            Bayesian = TRUE,
            min.coverage = c(12,4),
            min.meth = 3,
            high.coverage = 300,
            percentile = 0.999,
            num.cores = 64L,
            tasks = 20L,
            verbose = FALSE
        )
    )
}

date()


length(ls(pattern="HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr"))
# [1] 20

names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1)
# NULL

names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")

names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1)
# [1] "G22A_CHH" "G23A_CHH" "G24A_CHH" "G13A_CHH" "G14A_CHH" "G15A_CHH"


for(item in names(G14A_CHH_grl)) {
    assign(paste0("critical.val_G222324refctrl_vs_G131415_sum_chr", item),
    do.call(rbind, lapply(eval(str2expression(paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr",item))), function(x) {
    hd.95 = quantile(x$hdiv, 0.95)
    baytv.95 = quantile(abs(x$bay.TV), 0.95)
    tv.95 = quantile(abs(x$TV), 0.95)
    return(c(tv = tv.95, baytv = baytv.95, hd = hd.95, num.sites.hd95 = sum(x$hdiv > hd.95), num.sites.tv95 = sum(x$TV > tv.95), num.sites.baytv95 = sum(x$bay.TV > baytv.95)))}))
)
}
    
for(item in names(G1A_CHH_grl)) {
 print(eval(str2expression(paste0("critical.val_G222324refctrl_vs_G131415_sum_chr",item))))
}
## 20190109: 
#

date()

## 20190107: almost 3 hours runtime:
#


names(G14A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

names(G14A_CHH_grl[2:20])
# [1] "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

date()
# [1] "Wed Jan  8 07:54:00 2020"

# for(item in names(G14A_CHH_grl[2:20])) {
for(item in names(G14A_CHH_grl)) {
    assign(
        paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item),
        gofReport(
            eval(str2expression(
                paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr", item)
            )),
            model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
            column = 9,
            absolute = FALSE,
            output = c("best.model"),
            num.cores = 64L,
            verbose = FALSE
        )
    )
    print(eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
}



date()



ls()

## 20200109: these were saved in separate R session:
#
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
# save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")

save(critical.val_G222324refctrl_vs_G131415_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr1.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr10.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr11.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr12.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr13.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr14.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr15.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr16.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr17.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr18.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr19.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr2.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr20.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr3.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr4.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr5.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr6.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr7.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr8.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr9.RData")

## 20200109: these were saved in separate R session:
#
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9.RData")

















































date()

covr <-
    lapply(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1, function(x) {
        cov <- x$c2 + x$t2
        return(cov)
    })

do.call(rbind, lapply(covr, function(x) {
    q60 <- quantile(x, 0.6)
    q9999 <- quantile(x, 0.9999)
    idx1 <- which(x >= q60)
    idx2 <- which(x <= 500)
    q95 <- quantile(x, 0.95)
    idx <- intersect(idx1, idx2)
    return(c(
        round(summary(x)),
        q60,
        quantile(x, c(0.95, 0.99, 0.999, 0.9999)),
        num.siteGreater_8 = sum(x >= 8),
        q60_to_500 = sum((x >= q60) & (x <= 500)),
        num.siteGreater_500 = sum(x > 500)
    ))}))

## 20190109: we see:
#   - num.siteGreater_500 are all 0, confirming that our "high.coverage=300" filter in estimateDivergence worked
#   - 1st quartile values of about 10, showing that our data has good coverage, >10, for 75% of sites
#   - our data is extremely similar for all 6 samples in terms of coverage
# See also soybean_G222324refctrl_vs_G1415_sepChroms_20200102.R for more discussion.

length(mcols(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1$G13A_CHH)$hdiv)
# [1] 1095092

critical.val_G222324refctrl_vs_G131415_sum_chr1["G13A_CHH",][4]
# num.sites.hd95 
#          54755 

54755/1095092
# [1] 0.05000037

critical.val_G222324refctrl_vs_G131415_sum_chr1

critical.val_G222324refctrl_vs_G131415_sum_chr10

critical.val_G222324refctrl_vs_G131415_sum_chr13


date()

## 20200107:
#
# We call getPotentialDIMP with tv.cut values shown here:
#   - tv.col= 7L,
#   - tv.cut=0.2,
# ..normally, per a best practice, we would use a calculated value like this for chr1:
#
max(critical.val_G222324refctrl_vs_G131415_sum_chr1[,"tv.95%"])

#
# ..but that will cause too many DIMPs to be filtered out.  Therefore, the best 
# ..practice may be to use a hardcoded 0.2 in cases where tv.95% values are 
# ..lower than, say 0.25. 
#
# UPDATE 20200109, we will use 0.2 anyway, even though maybe 0. is better
#

for(item in names(G14A_CHH_grl)) {
    assign(
        paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr", item),
        getPotentialDIMP(
            eval(str2expression(
                paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr", item)
            )),
            nlms = eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$nlms"))),
            div.col = 9L,
            tv.col= 7L,
            tv.cut=0.2,
            dist.name = eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel")))
        )
    )
    #print(eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
}

## 20190109:
#
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr10.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr11.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr12.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr13.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr14.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr15.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr16.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr17.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr18.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr19.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr2.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr20.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr3.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr4.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr5.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr6.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr7.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr8.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr9.RData")

head(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1$G14A_CHH,n=2)

date()
# [1] "Wed Jan  8 13:45:42 2020"

names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1)
# [1] "G22A_CHH" "G23A_CHH" "G24A_CHH" "G13A_CHH" "G14A_CHH" "G15A_CHH"

ls(pattern="^PS")


mcols(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1$G14A_CHH)

treatment_names_CHH <- names(F4_P37_CHH)

control_names_CHH <- names(F6_P37_CHH)

control_names_CHH
# [1] "G22A_CHH" "G23A_CHH" "G24A_CHH"

treatment_names_CHH
# [1] "G13A_CHH" "G14A_CHH" "G15A_CHH"

names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr2) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr3) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr4) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr5) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr6) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr7) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr8) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr9) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr10) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr11) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr12) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr13) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr14) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr15) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr16) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr17) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr18) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr19) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr20) <- c("G22A_CHH","G23A_CHH","G24A_CHH","G13A_CHH","G14A_CHH","G15A_CHH")


















date()

#
# cutpoint with machine learning:
#
for(item in names(G14A_CHH_grl)) {
    assign(
        paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item),
        estimateCutPoint(
            eval(str2expression(
                paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
            )),
            simple = FALSE,
            control.names = control_names_CHH,
            treatment.names = treatment_names_CHH,
            column = c(hdiv = TRUE, bay.TV = TRUE, wprob = TRUE, pos = TRUE),
            div.col = 9, clas.perf = TRUE,
            classifier1 = "pca.qda", n.pc = 4, 
            classifier2 = "pca.lda",
            center = TRUE, scale = TRUE,
            verbose = TRUE)
        )

}
## 20190105: here we see that:
#  - 40 calls are made because we have:
#  - 20 samples
#  - 2 classifiers
#
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos

date()

cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr1$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr2$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr3$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr4$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr5$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr6$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr7$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr8$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr9$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr10$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr11$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr12$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr13$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr14$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr15$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr16$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr17$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr18$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr19$cutpoint
cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr20$cutpoint



for(item in names(G14A_CHH_grl)) {
    print(paste0("chr",item))
    print(eval(str2expression(
        paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$testSetPerformance$table")
    )))
    print(eval(str2expression(
        paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$testSetPerformance$overall")
    )))
}


# 
# cutpoint simple
#
for(item in names(G14A_CHH_grl)) {
    assign(
        paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item),
        estimateCutPoint(
            eval(str2expression(
                paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
            )),
            simple = TRUE,
            control.names = control_names_CHH,
            treatment.names = treatment_names_CHH,
            column = c(hdiv = TRUE, bay.TV = TRUE, wprob = TRUE, pos = TRUE),
            div.col = 9,
            clas.perf = TRUE,
            classifier1 = "qda",
            prop = 0.6,
            verbose = TRUE)
        )

}
## 20190105: here we see that:
#  - 20 calls are made because we have:
#  - 20 samples
#  - 1 classifiers
#
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
# Model: treat ~ hdiv + bay.TV + logP + pos
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr1$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr2$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr3$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr4$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr5$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr6$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr7$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr8$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr9$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr10$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr11$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr12$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr13$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr14$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr15$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr16$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr17$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr18$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr19$cutpoint
cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr20$cutpoint


for(item in names(G14A_CHH_grl)) {
    print(paste0("chr",item))
    print(eval(str2expression(
        paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$testSetPerformance$table")
    )))
    print(eval(str2expression(
        paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item,"$testSetPerformance$overall")
    )))
}



for(item in names(G14A_CHH_grl)) {
    assign(
        paste0("DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr", item),
        selectDIMP(
            eval(str2expression(
                paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
            )),
            div.col = 9,
            cutpoint = eval(str2expression(
                paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item, '$cutpoint')
            )),
            )
        )

}

ls(pattern="DIMPsYI")

unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr1, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr2, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr3, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr4, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr5, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr6, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr7, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr8, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr9, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr10, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr11, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr12, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr13, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr14, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr15, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr16, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr17, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr18, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr19, length))
unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr20, length))

date()

DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- data.frame("G22A_CHH"=0, "G23A_CHH"=0, "G24A_CHH"=0, "G13A_CHH"=0, "G14A_CHH"=0, "G15A_CHH"=0)
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr1, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr2, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr3, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr4, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr5, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr6, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr7, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr8, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr9, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr10, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr11, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr12, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr13, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr14, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr15, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr16, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr17, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr18, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr19, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr20, length)))

colSums(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS)

for(item in names(G14A_CHH_grl)) {
    assign(
        paste0("DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr", item),
        selectDIMP(
            eval(str2expression(
                paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
            )),
            div.col = 9,
            cutpoint = eval(str2expression(
                paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr", item, '$cutpoint')
            )),
            )
        )

}

ls(pattern="DIMPsML")


unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr1, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr2, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr3, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr4, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr5, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr6, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr7, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr8, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr9, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr10, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr11, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr12, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr13, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr14, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr15, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr16, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr17, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr18, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr19, length))
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr20, length))



date()

DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- data.frame("G22A_CHH"=0, "G23A_CHH"=0, "G24A_CHH"=0, "G13A_CHH"=0, "G14A_CHH"=0, "G15A_CHH"=0)
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr1, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr2, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr3, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr4, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr5, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr6, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr7, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr8, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr9, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr10, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr11, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr12, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr13, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr14, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr15, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr16, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr17, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr18, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr19, length)))
DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr20, length)))

colSums(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS)

date()

pryr::mem_used()

length(ls(pattern="best"))
# [1] 20
length(ls(pattern="critical"))
# [1] 20
length(ls(pattern="cutpointsML"))
# [1] 20
length(ls(pattern="cutpointsYI"))
# [1] 20
length(ls(pattern="DIMPsML"))
# [1] 21
length(ls(pattern="DIMPsYI"))
# [1] 21
length(ls(pattern="HD"))
# [1] 20
length(ls(pattern="PS"))
# [1] 60
length(ls(pattern="^PS_"))
# [1] 20

## 20190109:
# #
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/bestFits_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")
save(control_names_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/control_names_CHH.RData")
save(treatment_names_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/treatment_names_CHH.RData")
save(covr, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/covr.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr1.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr10.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr11.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr12.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr13.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr14.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr15.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr16.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr17.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr18.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr19.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr2.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr20.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr3.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr4.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr5.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr6.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr7.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr8.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr9.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHH_sum_COUNTS.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr1.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr10.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr11.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr12.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr13.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr14.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr15.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr16.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr17.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr18.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr19.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr2.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr20.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr3.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr4.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr5.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr6.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr7.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr8.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_chr9.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHH_sum_COUNTS.RData")
save(F4_P37_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/F4_P37_CHH.RData")
save(F6_P37_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/F6_P37_CHH.RData")


# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr1.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr10.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr11.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr12.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr13.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr14.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr15.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr16.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr17.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr18.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr19.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr2.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr20.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr3.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr4.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr5.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr6.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr7.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr8.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHH_sum_ID_chr9.RData")
# save(item, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/item.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr1.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr10.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr11.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr12.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr13.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr14.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr15.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr16.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr17.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr18.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr19.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr2.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr20.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr3.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr4.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr5.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr6.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr7.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr8.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHH_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHH_sum_chr9.RData")



date()






