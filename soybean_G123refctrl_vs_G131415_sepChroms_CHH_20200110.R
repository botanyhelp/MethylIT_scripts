# soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110.R

.libPaths()

library(BiocManager)
library(MethylIT)
library(MethylIT.utils)

rm(list = ls())
length(ls())
# [1] 0

## 20200110: warning: these RData files are 14GB each:

# G123 CHH:
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F4_WT_CHH_4_29_2019.RData")
names(F4_WT_CHH)
# [1] "G1A_CHH" "G2A_CHH" "G3A_CHH"

# G131415 CHH:
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F4_P37_CHH_4_29_2019.RData")
names(F4_P37_CHH)
# [1] "G13A_CHH" "G14A_CHH" "G15A_CHH"

ls()
# [1] "F4_P37_CHH"   "F4_WT_CHH"    "G13A_CHH_grl" "G14A_CHH_grl" "G15A_CHH_grl" "item"        

seqnames(F4_P37_CHH$G14A_CHH)
# factor-Rle of length 254667767 with 20 runs
#   Lengths: 15619475 13855537  9164311 10652937 11910229 13405530 ... 14344213 11297345 13530428 11878478 12636861 13507325
#   Values :        1       10       11       12       13       14 ...        4        5        6        7        8        9
# Levels(20): 1 10 11 12 13 14 15 16 17 18 19 2 20 3 4 5 6 7 8 9

G1A_CHH_grl <-
    split(F4_WT_CHH$G1A_CHH, seqnames(F4_P37_CHH$G14A_CHH))
G2A_CHH_grl <-
    split(F4_WT_CHH$G2A_CHH, seqnames(F4_P37_CHH$G14A_CHH))
G3A_CHH_grl <-
    split(F4_WT_CHH$G3A_CHH, seqnames(F4_P37_CHH$G14A_CHH))

G13A_CHH_grl <-
    split(F4_P37_CHH$G13A_CHH, seqnames(F4_P37_CHH$G13A_CHH))
G14A_CHH_grl <-
    split(F4_P37_CHH$G14A_CHH, seqnames(F4_P37_CHH$G13A_CHH))
G15A_CHH_grl <-
    split(F4_P37_CHH$G15A_CHH, seqnames(F4_P37_CHH$G13A_CHH))

names(G1A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

names(G14A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

date()
# [1] "Fri Jan 10 06:55:20 2020"

## 20200110: 30 minutes runtime

for (item in names(G14A_CHH_grl)) {
    assign(
        paste0("soy_WT_G123_Ref_CHH_sum_chr", item),
        poolFromGRlist(
            list(G1A_CHH_grl[[item]], G2A_CHH_grl[[item]], G3A_CHH_grl[[item]]),
            stat = "sum",
            num.cores = 12L
        )
    )
}

date()

ls()
#  [1] "F4_P37_CHH"                    "F4_WT_CHH"                     "G13A_CHH_grl"                 
#  [4] "G14A_CHH_grl"                  "G15A_CHH_grl"                  "G1A_CHH_grl"                  
#  [7] "G2A_CHH_grl"                   "G3A_CHH_grl"                   "item"                         
# [10] "soy_WT_G123_Ref_CHH_sum_chr1"  "soy_WT_G123_Ref_CHH_sum_chr10" "soy_WT_G123_Ref_CHH_sum_chr11"
# [13] "soy_WT_G123_Ref_CHH_sum_chr12" "soy_WT_G123_Ref_CHH_sum_chr13" "soy_WT_G123_Ref_CHH_sum_chr14"
# [16] "soy_WT_G123_Ref_CHH_sum_chr15" "soy_WT_G123_Ref_CHH_sum_chr16" "soy_WT_G123_Ref_CHH_sum_chr17"
# [19] "soy_WT_G123_Ref_CHH_sum_chr18" "soy_WT_G123_Ref_CHH_sum_chr19" "soy_WT_G123_Ref_CHH_sum_chr2" 
# [22] "soy_WT_G123_Ref_CHH_sum_chr20" "soy_WT_G123_Ref_CHH_sum_chr3"  "soy_WT_G123_Ref_CHH_sum_chr4" 
# [25] "soy_WT_G123_Ref_CHH_sum_chr5"  "soy_WT_G123_Ref_CHH_sum_chr6"  "soy_WT_G123_Ref_CHH_sum_chr7" 
# [28] "soy_WT_G123_Ref_CHH_sum_chr8"  "soy_WT_G123_Ref_CHH_sum_chr9" 

length(ls(pattern = "soy_WT_G123_Ref_CHH_sum_chr"))
# [1] 20

date()
# [1] "Fri Jan 10 07:24:55 2020"

## 20200110: 2hrs runtime:

for (item in names(G14A_CHH_grl)) {
    assign(
        paste0(
            "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr",
            item
        ),
        estimateDivergence(
            ref = eval(str2expression(
                paste0("soy_WT_G123_Ref_CHH_sum_chr", item)
            )),
            indiv = list(
                G1A_CHH_grl[[item]],
                G2A_CHH_grl[[item]],
                G3A_CHH_grl[[item]],
                G13A_CHH_grl[[item]],
                G14A_CHH_grl[[item]],
                G15A_CHH_grl[[item]]
            ),
            Bayesian = TRUE,
            min.coverage = c(12, 4),
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
# [1] "Fri Jan 10 09:25:33 2020"


################################################################################
################################################################################
################################################################################
# estimateDivergence blocks to run in 4 separate R consoles concurrently:
################################################################################
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G13A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G14A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G15A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G22A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G23A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G24A_CHH_grl.RData")
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr14.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr15.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr16.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr17.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr18.RData")
# 
# for(item in c("14","15","16","17","18")) {
#     assign(
#         paste0(
#             "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr",
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
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr14.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr15.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr16.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr17.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr18.RData")
# 
# 
# 
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G13A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G14A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G15A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G22A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G23A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G24A_CHH_grl.RData")
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr19.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr2.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr20.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr3.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr4.RData")
# 
# for(item in c("19","2","20","3","4")) {
#     assign(
#         paste0(
#             "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr",
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
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr19.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr20.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr2.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr3.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr4.RData")
# 
# 
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G13A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G14A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G15A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G22A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G23A_CHH_grl.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/G24A_CHH_grl.RData")
# 
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr5.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr6.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr7.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr8.RData")
# load(file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/soy_WT_G222324_Ref_CHH_sum_chr9.RData")
# 
# for(item in c("5","6","7","8","9")) {
#     assign(
#         paste0(
#             "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr",
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
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr5.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr6.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr7.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr8.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr9.RData")
# 
# 
################################################################################
################################################################################
################################################################################


length(
    ls(pattern = "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr")
)
# [1] 20

names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr1)
# NULL

names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr1) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr2) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr3) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr4) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr5) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr6) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr7) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr8) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr9) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr10) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr11) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr12) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr13) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr14) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr15) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr16) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr17) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr18) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr19) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr20) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")

names(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr1)
# [1] "G1A_CHH"  "G2A_CHH"  "G3A_CHH"  "G13A_CHH" "G14A_CHH" "G15A_CHH"


save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr1,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr1.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr10,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr10.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr11,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr11.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr12,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr12.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr13,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr13.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr14,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr14.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr15,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr15.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr16,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr16.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr17,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr17.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr18,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr18.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr19,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr19.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr20,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr20.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr2,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr2.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr3,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr3.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr4,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr4.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr5,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr5.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr6,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr6.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr7,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr7.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr8,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr8.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr9,file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr9.RData")

date()


for (item in names(G14A_CHH_grl)) {
    assign(
        paste0("critical.val_G123refctrl_vs_G131415_sum_chr", item),
        do.call(rbind, lapply(eval(
            str2expression(
                paste0(
                    "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr",
                    item
                )
            )
        ), function(x) {
            hd.95 = quantile(x$hdiv, 0.95)
            baytv.95 = quantile(abs(x$bay.TV), 0.95)
            tv.95 = quantile(abs(x$TV), 0.95)
            return(
                c(
                    tv = tv.95,
                    baytv = baytv.95,
                    hd = hd.95,
                    num.sites.hd95 = sum(x$hdiv > hd.95),
                    num.sites.tv95 = sum(x$TV > tv.95),
                    num.sites.baytv95 = sum(x$bay.TV > baytv.95)
                )
            )
        }))
    )
}

for (item in names(G14A_CHH_grl)) {
    print(eval(str2expression(
        paste0("critical.val_G123refctrl_vs_G131415_sum_chr", item)
    )))
}
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1875000 0.08791114 0.3403946          65929          23965             12459
# G2A_CHH  0.1683064 0.08829886 0.3436536          56749          34456             40669
# G3A_CHH  0.1846154 0.10048581 0.3967412          51303          37345             45167
# G13A_CHH 0.2428571 0.12460086 0.6257777          57092          22099             17786
# G14A_CHH 0.1907692 0.12020499 1.1264722          87948          58504             69389
# G15A_CHH 0.2500000 0.14373325 0.9214695          62877          44486             54685
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1811594 0.08695748 0.3643958          56198          20516              9205
# G2A_CHH  0.1590909 0.08505201 0.3791735          53435          32310             36834
# G3A_CHH  0.1727273 0.09636966 0.4402549          50303          36513             43251
# G13A_CHH 0.2299950 0.12137768 0.6845995          53252          20194             15298
# G14A_CHH 0.1904762 0.12046218 1.1512402          68821          45425             54233
# G15A_CHH 0.2403101 0.13968271 0.9824801          56043          40253             48324
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1923077 0.09032063 0.3445514          28779          10776              5756
# G2A_CHH  0.1666667 0.08879422 0.3572528          26649          16292             18971
# G3A_CHH  0.1846154 0.10054783 0.4042878          24157          17295             20984
# G13A_CHH 0.2500000 0.12792552 0.6373847          25065           9530              7990
# G14A_CHH 0.1948052 0.12340258 1.1505485          39098          26116             30773
# G15A_CHH 0.2504496 0.14527800 0.9271446          27948          20037             23906
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1842105 0.08823508 0.3645239          40980          15191              7316
# G2A_CHH  0.1632653 0.08672823 0.3739692          37798          22743             25482
# G3A_CHH  0.1754386 0.09745910 0.4350219          36016          25582             30095
# G13A_CHH 0.2333333 0.12360289 0.6883084          38609          14989             11246
# G14A_CHH 0.1944444 0.12192475 1.1397496          49987          33069             38890
# G15A_CHH 0.2420263 0.13937863 0.9656145          41330          28734             34510
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1916996 0.09291280 0.3805988          33611          12753              5641
# G2A_CHH  0.1666667 0.08918706 0.3806444          30914          18731             20791
# G3A_CHH  0.1818182 0.09991615 0.4337257          28991          20417             24273
# G13A_CHH 0.2424242 0.12924463 0.7140562          31349          12187              9097
# G14A_CHH 0.2000000 0.12793898 1.2098107          42058          28603             33629
# G15A_CHH 0.2500000 0.14443869 0.9892908          33901          23325             28419
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1904762 0.08626285 0.3180217          50764          18785             11852
# G2A_CHH  0.1754386 0.08899423 0.3146377          40279          24022             30437
# G3A_CHH  0.1911422 0.10217352 0.3662187          35866          26328             32613
# G13A_CHH 0.2500000 0.12357281 0.5797997          41966          16180             14218
# G14A_CHH 0.1904762 0.11879476 1.0995879          71964          47554             57004
# G15A_CHH 0.2564103 0.14520602 0.8719483          45724          32971             40388
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1875000 0.09847774 0.4169291          48216          17794              5912
# G2A_CHH  0.1619484 0.09229520 0.4332579          48460          28607             24690
# G3A_CHH  0.1717647 0.10149848 0.5054394          48078          33728             32405
# G13A_CHH 0.2341880 0.13360363 0.7868020          48605          18269             10913
# G14A_CHH 0.2028825 0.12917853 1.1898432          57125          36824             42920
# G15A_CHH 0.2460317 0.14447564 1.0248450          50298          33733             38622
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1875000 0.08780829 0.3495608          37070          13334              6754
# G2A_CHH  0.1627119 0.08611169 0.3695709          35656          21654             24938
# G3A_CHH  0.1750000 0.09725933 0.4317359          34032          24701             29338
# G13A_CHH 0.2352941 0.12351180 0.6706789          34760          13509             10802
# G14A_CHH 0.1920737 0.12145005 1.1683362          48270          32353             38179
# G15A_CHH 0.2500000 0.14284085 0.9498092          35719          24014             30614
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1785714 0.08813235 0.3896892          44172          15732              6350
# G2A_CHH  0.1555556 0.08395398 0.4144580          44278          26633             29215
# G3A_CHH  0.1666667 0.09342373 0.4944482          44286          31340             36644
# G13A_CHH 0.2222222 0.12102990 0.7567466          45255          16715             11717
# G14A_CHH 0.1931818 0.12214432 1.1619366          51186          34316             40243
# G15A_CHH 0.2333333 0.13671346 1.0270212          46338          32382             38317
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1875000 0.09243325 0.3858615          61470          22436              9613
# G2A_CHH  0.1612903 0.08824082 0.4029735          60512          36063             36099
# G3A_CHH  0.1727273 0.09816959 0.4649386          58652          41191             44280
# G13A_CHH 0.2333333 0.12738287 0.7265106          59869          22689             15741
# G14A_CHH 0.2000000 0.12500069 1.1621578          74527          48267             57415
# G15A_CHH 0.2456140 0.14154388 0.9843770          62473          42910             50735
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1818182 0.08585982 0.3489073          57151          20968             10329
# G2A_CHH  0.1600000 0.08469924 0.3662877          54029          32920             38268
# G3A_CHH  0.1742424 0.09623080 0.4334537          51278          37313             44810
# G13A_CHH 0.2309524 0.12101497 0.6580579          52870          20336             15804
# G14A_CHH 0.1889405 0.11881498 1.1277427          71960          47935             56639
# G15A_CHH 0.2400000 0.13836819 0.9479937          57122          40649             49073
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1844926 0.08758194 0.3575419          47271          17303              8254
# G2A_CHH  0.1612903 0.08536939 0.3717785          44556          26984             31142
# G3A_CHH  0.1755218 0.09695985 0.4262397          41183          29586             35688
# G13A_CHH 0.2333333 0.12261291 0.6734806          44488          17240             13216
# G14A_CHH 0.1904762 0.12053081 1.1448139          59268          39715             46717
# G15A_CHH 0.2400000 0.13860475 0.9601910          48224          34255             41125
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1764706 0.08483111 0.3630203          58233          20850              9544
# G2A_CHH  0.1558442 0.08268465 0.3785351          55690          33710             39174
# G3A_CHH  0.1691057 0.09432223 0.4461685          52821          38799             46387
# G13A_CHH 0.2222222 0.11777257 0.6842494          56142          21327             16075
# G14A_CHH 0.1875000 0.11802063 1.1331402          69714          45875             54812
# G15A_CHH 0.2307692 0.13474905 0.9974581          60924          43471             52177
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1996727 0.09103323 0.3323589          38479          14268              8382
# G2A_CHH  0.1775899 0.09176763 0.3306823          32117          19383             22882
# G3A_CHH  0.1923077 0.10359744 0.3763480          28993          20716             25058
# G13A_CHH 0.2555192 0.12949908 0.6034294          32147          13385             10758
# G14A_CHH 0.1988304 0.12340165 1.1289412          54564          35825             42953
# G15A_CHH 0.2655172 0.14963917 0.8821794          34909          24946             30371
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.2000000 0.09662453 0.3516032          45301          16627              8051
# G2A_CHH  0.1818182 0.09631120 0.3442900          38014          22574             24547
# G3A_CHH  0.1965812 0.10807698 0.3899137          34189          24316             27316
# G13A_CHH 0.2608696 0.13542524 0.6266558          38000          15149             11158
# G14A_CHH 0.2000000 0.12626437 1.1397316          64888          42266             49865
# G15A_CHH 0.2698413 0.15446560 0.8992651          41392          29352             35160
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1781004 0.08565258 0.3619281          45187          16462              7650
# G2A_CHH  0.1586022 0.08437894 0.3801917          42789          26075             30006
# G3A_CHH  0.1727941 0.09552173 0.4349179          39744          28586             34149
# G13A_CHH 0.2271075 0.11987554 0.6884040          43242          16107             12415
# G14A_CHH 0.1894737 0.11960964 1.1487646          54801          36396             43057
# G15A_CHH 0.2352941 0.13678272 0.9862991          46632          32985             39662
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1818182 0.08698647 0.3653826          51433          19111              9101
# G2A_CHH  0.1593407 0.08463120 0.3802626          48812          29414             33560
# G3A_CHH  0.1724138 0.09520831 0.4406823          46518          33549             40272
# G13A_CHH 0.2307692 0.12144544 0.6937035          49350          18647             14180
# G14A_CHH 0.1923077 0.12149326 1.1538538          62188          41663             49301
# G15A_CHH 0.2380952 0.13747692 0.9772733          52482          37105             44698
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1894737 0.08941058 0.3534747          43492          16560              8318
# G2A_CHH  0.1666667 0.08828781 0.3587121          39030          23423             27604
# G3A_CHH  0.1833333 0.09966280 0.4074846          35753          25561             31053
# G13A_CHH 0.2424242 0.12585303 0.6525738          38963          15380             12043
# G14A_CHH 0.1956522 0.12319526 1.1509771          56306          37619             44244
# G15A_CHH 0.2500000 0.14253851 0.9353918          42909          28942             36472
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.2000000 0.09616889 0.3615868          33102          12073              5985
# G2A_CHH  0.1764706 0.09446074 0.3549505          28338          17079             18454
# G3A_CHH  0.1932773 0.10651560 0.4012141          25572          18062             20581
# G13A_CHH 0.2564103 0.13597230 0.6711545          28491          11313              8299
# G14A_CHH 0.2023786 0.12879263 1.1723971          44723          30098             35351
# G15A_CHH 0.2666667 0.15273188 0.9103824          30157          21262             25085
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1875000 0.08730268 0.3462276          51227          18656              9671
# G2A_CHH  0.1666667 0.08660777 0.3577828          47002          27428             33639
# G3A_CHH  0.1790659 0.09826972 0.4153036          43811          32198             38829
# G13A_CHH 0.2380952 0.12308839 0.6456746          45900          17967             14509
# G14A_CHH 0.1904762 0.12029184 1.1428297          66985          44460             52566
# G15A_CHH 0.2457002 0.14062165 0.9455568          51232          36170             44114

date()
# [1] "Fri Jan 10 10:23:58 2020"


names(G14A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 


## 20200110_1100am: these gofReport calls are being run simultaneously in four 
#  ..separate terminals.  This code works but will take very long and so we 
#  ..run concurrently to save time.  Therefore, we will use load() below, in 
#  ..this program to load() the RData objects save()-ed elsewhere.  The four 
#  ..code blocks that were run in another R console are pasted in the comments 
#  ..below. 

for (item in names(G14A_CHH_grl)) {
    assign(
        paste0("bestFits_G123refctrl_vs_G131415_CHH_sum_chr", item),
        gofReport(
            eval(str2expression(
                paste0(
                    "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr",
                    item
                )
            )),
            model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
            column = 9,
            absolute = FALSE,
            output = c("best.model"),
            num.cores = 64L,
            verbose = FALSE
        )
    )
    print(eval(str2expression(
        paste0(
            "bestFits_G123refctrl_vs_G131415_CHH_sum_chr",
            item,
            "$bestModel"
        )
    )))
}

## 20200110_1100am: these gofReport calls are being run simultaneously in four 
#  ..separate terminals.  This code works but will take very long and so we 
#  ..run concurrently to save time.  Therefore, we will use load() below, in 
#  ..this program to load() the RData objects save()-ed elsewhere.  The four 
#  ..code blocks that were run in another R console are pasted in the comments 
#  ..below. 
################################################################################
## gofReport blocks to run in 4 separate R consoles concurrently:
################################################################################
# 
# library(BiocManager)
# library(MethylIT)
# library(MethylIT.utils)
# 
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr1.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr10.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr11.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr12.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr13.RData")
# 
# for(item in c("1","10","11","12","13")) {
#   assign(
#     paste0("bestFits_G123refctrl_vs_G131415_CHH_sum_chr", item),
#     gofReport(
#       eval(str2expression(
#         paste0("HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr", item)
#       )),
#       model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
#       column = 9,
#       absolute = FALSE,
#       output = c("best.model"),
#       num.cores = 64L,
#       verbose = FALSE
#     )
#   )
#   print(eval(str2expression(paste0("bestFits_G123refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
# }
# 
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr1.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr10.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr11.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr12.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr13.RData")
# 
# 
# library(BiocManager)
# library(MethylIT)
# library(MethylIT.utils)
# 
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr14.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr15.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr16.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr17.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr18.RData")
# 
# for(item in c("14","15","16","17","18")) {
#   assign(
#     paste0("bestFits_G123refctrl_vs_G131415_CHH_sum_chr", item),
#     gofReport(
#       eval(str2expression(
#         paste0("HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr", item)
#       )),
#       model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
#       column = 9,
#       absolute = FALSE,
#       output = c("best.model"),
#       num.cores = 64L,
#       verbose = FALSE
#     )
#   )
#   print(eval(str2expression(paste0("bestFits_G123refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
# }
# 
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr14.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr15.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr16.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr17.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr18.RData")
# 
# 
# library(BiocManager)
# library(MethylIT)
# library(MethylIT.utils)
# 
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr19.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr20.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr2.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr3.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr4.RData")
# 
# for(item in c("19","2","20","3","4")) {
#   assign(
#     paste0("bestFits_G123refctrl_vs_G131415_CHH_sum_chr", item),
#     gofReport(
#       eval(str2expression(
#         paste0("HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr", item)
#       )),
#       model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
#       column = 9,
#       absolute = FALSE,
#       output = c("best.model"),
#       num.cores = 64L,
#       verbose = FALSE
#     )
#   )
#   print(eval(str2expression(paste0("bestFits_G123refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
# }
# 
# 
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr19.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr2.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr20.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr3.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr4.RData")
# 
# 
# library(BiocManager)
# library(MethylIT)
# library(MethylIT.utils)
# 
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr5.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr6.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr7.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr8.RData")
# load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr9.RData")
# 
# for(item in c("5","6","7","8","9")) {
#   assign(
#     paste0("bestFits_G123refctrl_vs_G131415_CHH_sum_chr", item),
#     gofReport(
#       eval(str2expression(
#         paste0("HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr", item)
#       )),
#       model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
#       column = 9,
#       absolute = FALSE,
#       output = c("best.model"),
#       num.cores = 64L,
#       verbose = FALSE
#     )
#   )
#   print(eval(str2expression(paste0("bestFits_G123refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
# }
# 
# 
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr5.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr6.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr7.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr8.RData")
# save(bestFits_G123refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr9.RData")
# 
################################################################################
################################################################################
################################################################################







load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr1.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr10.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr11.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr12.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr13.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr14.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr15.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr16.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr17.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr18.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr19.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr2.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr20.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr3.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr4.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr5.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr6.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr7.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr8.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/bestFits_G123refctrl_vs_G131415_CHH_sum_chr9.RData")


date()
# [1] "Fri Jan 10 15:01:56 2020"

ls()
#  [1] "bestFits_G123refctrl_vs_G131415_CHH_sum_chr1"                        "bestFits_G123refctrl_vs_G131415_CHH_sum_chr10"                      
#  [3] "bestFits_G123refctrl_vs_G131415_CHH_sum_chr11"                       "bestFits_G123refctrl_vs_G131415_CHH_sum_chr12"                      
#  [5] "bestFits_G123refctrl_vs_G131415_CHH_sum_chr13"                       "bestFits_G123refctrl_vs_G131415_CHH_sum_chr14"                      
#  [7] "bestFits_G123refctrl_vs_G131415_CHH_sum_chr15"                       "bestFits_G123refctrl_vs_G131415_CHH_sum_chr16"                      
#  [9] "bestFits_G123refctrl_vs_G131415_CHH_sum_chr17"                       "bestFits_G123refctrl_vs_G131415_CHH_sum_chr18"                      
# [11] "bestFits_G123refctrl_vs_G131415_CHH_sum_chr19"                       "bestFits_G123refctrl_vs_G131415_CHH_sum_chr2"                       
# [13] "bestFits_G123refctrl_vs_G131415_CHH_sum_chr20"                       "bestFits_G123refctrl_vs_G131415_CHH_sum_chr3"                       
# [15] "bestFits_G123refctrl_vs_G131415_CHH_sum_chr4"                        "bestFits_G123refctrl_vs_G131415_CHH_sum_chr5"                       
# [17] "bestFits_G123refctrl_vs_G131415_CHH_sum_chr6"                        "bestFits_G123refctrl_vs_G131415_CHH_sum_chr7"                       
# [19] "bestFits_G123refctrl_vs_G131415_CHH_sum_chr8"                        "bestFits_G123refctrl_vs_G131415_CHH_sum_chr9"                       
# [21] "critical.val_G123refctrl_vs_G131415_sum_chr1"                        "critical.val_G123refctrl_vs_G131415_sum_chr10"                      
# [23] "critical.val_G123refctrl_vs_G131415_sum_chr11"                       "critical.val_G123refctrl_vs_G131415_sum_chr12"                      
# [25] "critical.val_G123refctrl_vs_G131415_sum_chr13"                       "critical.val_G123refctrl_vs_G131415_sum_chr14"                      
# [27] "critical.val_G123refctrl_vs_G131415_sum_chr15"                       "critical.val_G123refctrl_vs_G131415_sum_chr16"                      
# [29] "critical.val_G123refctrl_vs_G131415_sum_chr17"                       "critical.val_G123refctrl_vs_G131415_sum_chr18"                      
# [31] "critical.val_G123refctrl_vs_G131415_sum_chr19"                       "critical.val_G123refctrl_vs_G131415_sum_chr2"                       
# [33] "critical.val_G123refctrl_vs_G131415_sum_chr20"                       "critical.val_G123refctrl_vs_G131415_sum_chr3"                       
# [35] "critical.val_G123refctrl_vs_G131415_sum_chr4"                        "critical.val_G123refctrl_vs_G131415_sum_chr5"                       
# [37] "critical.val_G123refctrl_vs_G131415_sum_chr6"                        "critical.val_G123refctrl_vs_G131415_sum_chr7"                       
# [39] "critical.val_G123refctrl_vs_G131415_sum_chr8"                        "critical.val_G123refctrl_vs_G131415_sum_chr9"                       
# [41] "F4_P37_CHH"                                                          "F4_WT_CHH"                                                          
# [43] "G13A_CHH_grl"                                                        "G14A_CHH_grl"                                                       
# [45] "G15A_CHH_grl"                                                        "G1A_CHH_grl"                                                        
# [47] "G2A_CHH_grl"                                                         "G3A_CHH_grl"                                                        
# [49] "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr1"  "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr10"
# [51] "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr11" "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr12"
# [53] "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr13" "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr14"
# [55] "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr15" "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr16"
# [57] "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr17" "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr18"
# [59] "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr19" "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr2" 
# [61] "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr20" "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr3" 
# [63] "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr4"  "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr5" 
# [65] "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr6"  "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr7" 
# [67] "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr8"  "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr9" 
# [69] "item"                                                                "soy_WT_G123_Ref_CHH_sum_chr1"                                       
# [71] "soy_WT_G123_Ref_CHH_sum_chr10"                                       "soy_WT_G123_Ref_CHH_sum_chr11"                                      
# [73] "soy_WT_G123_Ref_CHH_sum_chr12"                                       "soy_WT_G123_Ref_CHH_sum_chr13"                                      
# [75] "soy_WT_G123_Ref_CHH_sum_chr14"                                       "soy_WT_G123_Ref_CHH_sum_chr15"                                      
# [77] "soy_WT_G123_Ref_CHH_sum_chr16"                                       "soy_WT_G123_Ref_CHH_sum_chr17"                                      
# [79] "soy_WT_G123_Ref_CHH_sum_chr18"                                       "soy_WT_G123_Ref_CHH_sum_chr19"                                      
# [81] "soy_WT_G123_Ref_CHH_sum_chr2"                                        "soy_WT_G123_Ref_CHH_sum_chr20"                                      
# [83] "soy_WT_G123_Ref_CHH_sum_chr3"                                        "soy_WT_G123_Ref_CHH_sum_chr4"                                       
# [85] "soy_WT_G123_Ref_CHH_sum_chr5"                                        "soy_WT_G123_Ref_CHH_sum_chr6"                                       
# [87] "soy_WT_G123_Ref_CHH_sum_chr7"                                        "soy_WT_G123_Ref_CHH_sum_chr8"                                       
# [89] "soy_WT_G123_Ref_CHH_sum_chr9"                                       


covr <-
    lapply(HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr1, function(x) {
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
    return(
        c(
            round(summary(x)),
            q60,
            quantile(x, c(0.95, 0.99, 0.999, 0.9999)),
            num.siteGreater_8 = sum(x >= 8),
            q60_to_500 = sum((x >= q60) & (x <= 500)),
            num.siteGreater_500 = sum(x > 500)
        )
    )
}))
#          Min. 1st Qu. Median Mean 3rd Qu. Max. 60% 95% 99% 99.9% 99.99% num.siteGreater_8 q60_to_500 num.siteGreater_500
# G1A_CHH     0      10     17   17      24   52  19  33  38    44     48           1127306     564902                   0
# G2A_CHH     0      11     16   16      20   42  17  26  30    34     37           1041260     517733                   0
# G3A_CHH     0      11     14   14      18   38  16  23  27    30     33            926655     425097                   0
# G13A_CHH    0       9     14   14      19   50  16  26  31    37     42            950349     458078                   0
# G14A_CHH    0      21     31   33      44  124  36  63  77    92    103           1718486     714279                   0
# G15A_CHH    0      10     15   16      20   60  17  28  33    39     45           1095208     537077                   0


length(
    mcols(
        HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr1$G13A_CHH
    )$hdiv
)
# [1] 1141870

critical.val_G123refctrl_vs_G131415_sum_chr1["G13A_CHH", ][4]
# num.sites.hd95 
#          57092 

57092 / 1141870
# [1] 0.04999869

critical.val_G123refctrl_vs_G131415_sum_chr1
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1875000 0.08791114 0.3403946          65929          23965             12459
# G2A_CHH  0.1683064 0.08829886 0.3436536          56749          34456             40669
# G3A_CHH  0.1846154 0.10048581 0.3967412          51303          37345             45167
# G13A_CHH 0.2428571 0.12460086 0.6257777          57092          22099             17786
# G14A_CHH 0.1907692 0.12020499 1.1264722          87948          58504             69389
# G15A_CHH 0.2500000 0.14373325 0.9214695          62877          44486             54685

critical.val_G123refctrl_vs_G131415_sum_chr10
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1811594 0.08695748 0.3643958          56198          20516              9205
# G2A_CHH  0.1590909 0.08505201 0.3791735          53435          32310             36834
# G3A_CHH  0.1727273 0.09636966 0.4402549          50303          36513             43251
# G13A_CHH 0.2299950 0.12137768 0.6845995          53252          20194             15298
# G14A_CHH 0.1904762 0.12046218 1.1512402          68821          45425             54233
# G15A_CHH 0.2403101 0.13968271 0.9824801          56043          40253             48324

critical.val_G123refctrl_vs_G131415_sum_chr13
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G1A_CHH  0.1916996 0.09291280 0.3805988          33611          12753              5641
# G2A_CHH  0.1666667 0.08918706 0.3806444          30914          18731             20791
# G3A_CHH  0.1818182 0.09991615 0.4337257          28991          20417             24273
# G13A_CHH 0.2424242 0.12924463 0.7140562          31349          12187              9097
# G14A_CHH 0.2000000 0.12793898 1.2098107          42058          28603             33629
# G15A_CHH 0.2500000 0.14443869 0.9892908          33901          23325             28419

max(critical.val_G123refctrl_vs_G131415_sum_chr1[, "tv.95%"])
# [1] 0.25

for (item in names(G14A_CHH_grl)) {
    assign(
        paste0("PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr", item),
        getPotentialDIMP(
            eval(str2expression(
                paste0(
                    "HD_mc12_4_mm3_hc300_p999_WT_G123refctrl_vs_G131415_CHH_sum_ID_chr",
                    item
                )
            )),
            nlms = eval(str2expression(
                paste0(
                    "bestFits_G123refctrl_vs_G131415_CHH_sum_chr",
                    item,
                    "$nlms"
                )
            )),
            div.col = 9L,
            tv.col = 7L,
            tv.cut = 0.2,
            dist.name = eval(str2expression(
                paste0(
                    "bestFits_G123refctrl_vs_G131415_CHH_sum_chr",
                    item,
                    "$bestModel"
                )
            ))
        )
    )
    #print(eval(str2expression(paste0("bestFits_G123refctrl_vs_G131415_CHH_sum_chr", item,"$bestModel"))))
}

## 20200110_1508pm saved:
#
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr1.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr10.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr11.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr12.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr13.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr14.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr15.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr16.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr17.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr18.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr19.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr2.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr20.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr3.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr4.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr5.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr6.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr7.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr8.RData")
# save(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr9.RData")

head(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr1$G14A_CHH, n = 2)
# GRanges object with 2 ranges and 10 metadata columns:
#       seqnames    ranges strand |        c1        t1        c2        t2                 p1                p2                TV            bay.TV             hdiv
#          <Rle> <IRanges>  <Rle> | <numeric> <numeric> <numeric> <numeric>          <numeric>         <numeric>         <numeric>         <numeric>        <numeric>
#   [1]        1      5248      - |         6        54        14        21  0.105206384184913 0.336897833227396               0.3 0.231691449042484 3.75351745767507
#   [2]        1      5931      - |         1        35         4         8 0.0575925675921796 0.241042880945057 0.305555555555556 0.183450313352878  1.4027596083637
#                      wprob
#                  <numeric>
#   [1] 0.000268067649449869
#   [2]    0.028338455585167
#   -------
#   seqinfo: 20 sequences from an unspecified genome; no seqlengths

names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr1)
# [1] "G1A_CHH"  "G2A_CHH"  "G3A_CHH"  "G13A_CHH" "G14A_CHH" "G15A_CHH"

ls(pattern = "^PS")
#  [1] "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr1"  "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr10" "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr11"
#  [4] "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr12" "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr13" "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr14"
#  [7] "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr15" "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr16" "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr17"
# [10] "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr18" "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr19" "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr2" 
# [13] "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr20" "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr3"  "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr4" 
# [16] "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr5"  "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr6"  "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr7" 
# [19] "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr8"  "PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr9" 

mcols(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr1$G14A_CHH)
# DataFrame with 19356 rows and 10 columns
#              c1        t1        c2        t2                 p1                p2                TV            bay.TV             hdiv                wprob
#       <numeric> <numeric> <numeric> <numeric>          <numeric>         <numeric>         <numeric>         <numeric>        <numeric>            <numeric>
# 1             6        54        14        21  0.105206384184913 0.336897833227396               0.3 0.231691449042484 3.75351745767507 0.000268067649449869
# 2             1        35         4         8 0.0575925675921796 0.241042880945057 0.305555555555556 0.183450313352878  1.4027596083637    0.028338455585167
# 3             7        24         9         5  0.191431143526062 0.418306026237778  0.41705069124424 0.226874882711716 1.27489785653584   0.0368395103814135
# 4             2        38         6        15 0.0713404500396498 0.235828928291945 0.235714285714286 0.164488478252295 1.59679759353086    0.019098150975153
# 5            24        52        20        18  0.282597095147611 0.437718284901343 0.210526315789474 0.155121189753732 1.36365435384998   0.0306992879122935
# ...         ...       ...       ...       ...                ...               ...               ...               ...              ...                  ...
# 19352         2        45         6        11 0.0634114292505699 0.268942774393751 0.310387984981227 0.205531345143181 2.19628509719486  0.00574465888801232
# 19353         1        40        10        33 0.0525389648538767 0.214021242378413 0.208167895632445 0.161482277524537 2.66535171081031  0.00226984107222526
# 19354         9        31        19        23  0.196381328045228 0.386287644058239 0.227380952380952  0.18990631601301 1.87199195593272   0.0109727140649642
# 19355         6        30        13        21  0.153780278034479 0.322320009420001 0.215686274509804 0.168539731385521 1.44530940167387   0.0259806370601423
# 19356         0        30         5        16 0.0433598467664298 0.205047430093628 0.238095238095238 0.161687583327198 1.73211267888778    0.014532752363519

treatment_names_CHH <- names(F4_P37_CHH)

control_names_CHH <- names(F4_WT_CHH)

control_names_CHH
# [1] "G1A_CHH" "G2A_CHH" "G3A_CHH"

treatment_names_CHH
# [1] "G13A_CHH" "G14A_CHH" "G15A_CHH"

names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr1)
# [1] "G1A_CHH"  "G2A_CHH"  "G3A_CHH"  "G13A_CHH" "G14A_CHH" "G15A_CHH"

names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr1) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr2) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr3) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr4) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr5) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr6) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr7) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr8) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr9) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr10) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr11) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr12) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr13) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr14) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr15) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr16) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr17) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr18) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr19) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")
names(PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr20) <-
    c("G1A_CHH",
        "G2A_CHH",
        "G3A_CHH",
        "G13A_CHH",
        "G14A_CHH",
        "G15A_CHH")


for (item in names(G14A_CHH_grl)) {
    assign(
        paste0(
            "cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr",
            item
        ),
        estimateCutPoint(
            eval(str2expression(
                paste0("PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
            )),
            simple = FALSE,
            control.names = control_names_CHH,
            treatment.names = treatment_names_CHH,
            column = c(
                hdiv = TRUE,
                bay.TV = TRUE,
                wprob = TRUE,
                pos = TRUE
            ),
            div.col = 9,
            clas.perf = TRUE,
            classifier1 = "pca.qda",
            n.pc = 4,
            classifier2 = "pca.lda",
            center = TRUE,
            scale = TRUE,
            verbose = TRUE
        )
    )
    
}

cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr1$cutpoint
# 0.343572
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr2$cutpoint
# 0.4299256
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr3$cutpoint
# 0.3775363
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr4$cutpoint
# 0.369985
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr5$cutpoint
# 0.4434331
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr6$cutpoint
# 0.443965
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr7$cutpoint
# 0.4126668
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr8$cutpoint
# 0.370007
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr9$cutpoint
# 0.4198379
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr10$cutpoint
# 0.6523464
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr11$cutpoint
# 0.6006159
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr12$cutpoint
# 0.3909242
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr13$cutpoint
# 0.4435291
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr14$cutpoint
# 0.3689361
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr15$cutpoint
# 0.510349
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr16$cutpoint
# 0.374841
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr17$cutpoint
# 0.4978735
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr18$cutpoint
# 0.5218214
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr19$cutpoint
# 0.4337695
cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr20$cutpoint
# 0.4470871



for (item in names(G14A_CHH_grl)) {
    print(paste0("chr", item))
    print(eval(str2expression(
        paste0(
            "cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr",
            item,
            "$testSetPerformance$table"
        )
    )))
    print(eval(str2expression(
        paste0(
            "cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr",
            item,
            "$testSetPerformance$overall"
        )
    )))
}
# [1] "chr1"
#           Reference
# Prediction    CT    TT
#         CT  3513     0
#         TT  1145 17886
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.492104e-01   8.295954e-01   9.462636e-01   9.520407e-01   7.933818e-01   0.000000e+00  1.489538e-250 
# [1] "chr10"
#           Reference
# Prediction    CT    TT
#         CT  1425     4
#         TT    50 13560
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.964093e-01   9.794183e-01   9.953175e-01   9.973015e-01   9.019217e-01   0.000000e+00   9.141298e-10 
# [1] "chr11"
#           Reference
# Prediction   CT   TT
#         CT  939    1
#         TT   71 8426
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.923705e-01   9.588287e-01   9.904014e-01   9.940257e-01   8.929745e-01   0.000000e+00   4.232137e-16 
# [1] "chr12"
#           Reference
# Prediction   CT   TT
#         CT 1850    0
#         TT  492 9916
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.598629e-01   8.588269e-01   9.562366e-01   9.632688e-01   8.089411e-01   0.000000e+00  1.423541e-108 
# [1] "chr13"
#           Reference
# Prediction   CT   TT
#         CT 1636    0
#         TT  417 9066
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.624966e-01   8.648238e-01   9.587970e-01   9.659515e-01   8.153611e-01   0.000000e+00   2.987093e-92 
# [1] "chr14"
#           Reference
# Prediction    CT    TT
#         CT  2732     0
#         TT   828 14701
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.546575e-01   8.415853e-01   9.515393e-01   9.576298e-01   8.050490e-01   0.000000e+00  1.199616e-181 
# [1] "chr15"
#           Reference
# Prediction    CT    TT
#         CT  1302    14
#         TT  1120 11919
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.210031e-01   6.557298e-01   9.164700e-01   9.253648e-01   8.312783e-01  3.338896e-217  3.749811e-236 
# [1] "chr16"
#           Reference
# Prediction   CT   TT
#         CT 1784    0
#         TT  585 9792
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.518954e-01   8.308243e-01   9.479415e-01   9.556312e-01   8.051969e-01   0.000000e+00  8.339706e-129 
# [1] "chr17"
#           Reference
# Prediction   CT   TT
#         CT 1622    0
#         TT  337 9871
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.715131e-01   8.892834e-01   9.683545e-01   9.744367e-01   8.344041e-01   0.000000e+00   7.819137e-75 
# [1] "chr18"
#           Reference
# Prediction    CT    TT
#         CT  1709    10
#         TT   843 15014
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.514679e-01   7.738486e-01   9.481871e-01   9.545983e-01   8.548020e-01   0.000000e+00  1.690754e-178 
# [1] "chr19"
#           Reference
# Prediction    CT    TT
#         CT  2433     0
#         TT   547 14028
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.678387e-01   8.800549e-01   9.650765e-01   9.704393e-01   8.247883e-01   0.000000e+00  1.539234e-120 
# [1] "chr2"
#           Reference
# Prediction    CT    TT
#         CT  2070     0
#         TT   542 11765
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.623009e-01   8.620812e-01   9.590577e-01   9.653552e-01   8.183209e-01   0.000000e+00  1.883784e-119 
# [1] "chr20"
#           Reference
# Prediction    CT    TT
#         CT  2220     0
#         TT   554 13093
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.650848e-01   8.686507e-01   9.621105e-01   9.678870e-01   8.251717e-01   0.000000e+00  4.618677e-122 
# [1] "chr3"
#           Reference
# Prediction    CT    TT
#         CT  2071     0
#         TT   810 11872
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.450959e-01   8.044961e-01   9.412976e-01   9.487175e-01   8.047177e-01   0.000000e+00  9.827871e-178 
# [1] "chr4"
#           Reference
# Prediction    CT    TT
#         CT  2944     0
#         TT   812 14458
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.554189e-01   8.519816e-01   9.523207e-01   9.583707e-01   7.937850e-01   0.000000e+00  3.611022e-178 
# [1] "chr5"
#           Reference
# Prediction    CT    TT
#         CT  1412     0
#         TT   808 10367
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.358068e-01   7.421763e-01   9.313836e-01   9.400271e-01   8.236276e-01  6.105165e-301  2.674792e-177 
# [1] "chr6"
#           Reference
# Prediction    CT    TT
#         CT  2121     0
#         TT   493 12321
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.669903e-01   8.765198e-01   9.639990e-01   9.697978e-01   8.249749e-01   0.000000e+00  8.625467e-109 
# [1] "chr7"
#           Reference
# Prediction    CT    TT
#         CT  2156     0
#         TT   575 11732
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.602434e-01   8.588188e-01   9.569306e-01   9.633694e-01   8.111733e-01   0.000000e+00  1.248420e-126 
# [1] "chr8"
#           Reference
# Prediction    CT    TT
#         CT  1887     0
#         TT   733 10396
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.436847e-01   8.043940e-01   9.395856e-01   9.475839e-01   7.987093e-01   0.000000e+00  5.425746e-161 
# [1] "chr9"
#           Reference
# Prediction    CT    TT
#         CT  2317     0
#         TT   628 13463
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.617260e-01   8.582477e-01   9.586769e-01   9.646100e-01   8.205144e-01   0.000000e+00  3.702054e-138 


for (item in names(G14A_CHH_grl)) {
    assign(
        paste0(
            "cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr",
            item
        ),
        estimateCutPoint(
            eval(str2expression(
                paste0("PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
            )),
            simple = TRUE,
            control.names = control_names_CHH,
            treatment.names = treatment_names_CHH,
            column = c(
                hdiv = TRUE,
                bay.TV = TRUE,
                wprob = TRUE,
                pos = TRUE
            ),
            div.col = 9,
            clas.perf = TRUE,
            classifier1 = "qda",
            prop = 0.6,
            verbose = TRUE
        )
    )
    
}

cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr1$cutpoint
# 0.9031353
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr2$cutpoint
# 0.9350445
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr3$cutpoint
# 0.8590336
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr4$cutpoint
# 0.8780821
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr5$cutpoint
# 0.9704326
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr6$cutpoint
# 0.9549699
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr7$cutpoint
# 0.913291
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr8$cutpoint
# 0.8910695
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr9$cutpoint
# 0.9201756
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr10$cutpoint
# 0.9649035
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr11$cutpoint
# 0.9071639
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr12$cutpoint
# 0.9464278
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr13$cutpoint
# 0.9645466
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr14$cutpoint
# 0.8501074
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr15$cutpoint
# 0.9936481
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr16$cutpoint
# 0.9335336
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr17$cutpoint
# 0.999052
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr18$cutpoint
# 0.9576507
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr19$cutpoint
# 0.9282379
cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr20$cutpoint
# 0.9736832


for (item in names(G14A_CHH_grl)) {
    print(paste0("chr", item))
    print(eval(str2expression(
        paste0(
            "cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr",
            item,
            "$testSetPerformance$table"
        )
    )))
    print(eval(str2expression(
        paste0(
            "cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr",
            item,
            "$testSetPerformance$overall"
        )
    )))
}
# [1] "chr1"
#           Reference
# Prediction    CT    TT
#         CT   760     0
#         TT     0 15234
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9997694      1.0000000      0.9524822      0.0000000            NaN 
# [1] "chr10"
#           Reference
# Prediction    CT    TT
#         CT   580     0
#         TT     0 11592
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.996970e-01   1.000000e+00   9.523497e-01  8.130759e-259            NaN 
# [1] "chr11"
#           Reference
# Prediction   CT   TT
#         CT  384    0
#         TT    0 7109
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.995078e-01   1.000000e+00   9.487522e-01  6.394296e-172            NaN 
# [1] "chr12"
#           Reference
# Prediction   CT   TT
#         CT  441    0
#         TT    0 8513
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.995881e-01   1.000000e+00   9.507483e-01  3.970773e-197            NaN 
# [1] "chr13"
#           Reference
# Prediction   CT   TT
#         CT  421    0
#         TT    0 7799
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.995513e-01   1.000000e+00   9.487835e-01  2.058662e-188            NaN 
# [1] "chr14"
#           Reference
# Prediction    CT    TT
#         CT   610     0
#         TT     0 12362
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.997157e-01   1.000000e+00   9.529756e-01  4.457588e-272            NaN 
# [1] "chr15"
#           Reference
# Prediction    CT    TT
#         CT   611     0
#         TT     0 10576
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.996703e-01   1.000000e+00   9.453830e-01  1.331999e-273            NaN 
# [1] "chr16"
#           Reference
# Prediction   CT   TT
#         CT  433    0
#         TT    0 8402
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.995826e-01   1.000000e+00   9.509904e-01  1.534793e-193            NaN 
# [1] "chr17"
#           Reference
# Prediction   CT   TT
#         CT  476    0
#         TT    0 8628
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.995949e-01   1.000000e+00   9.477153e-01  4.736818e-213            NaN 
# [1] "chr18"
#           Reference
# Prediction    CT    TT
#         CT   696     0
#         TT     0 13104
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.997327e-01   1.000000e+00   9.495652e-01  6.953646e-311            NaN 
# [1] "chr19"
#           Reference
# Prediction    CT    TT
#         CT   625     0
#         TT     0 11989
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.997076e-01   1.000000e+00   9.504519e-01  4.080416e-279            NaN 
# [1] "chr2"
#           Reference
# Prediction    CT    TT
#         CT   554     0
#         TT     0 10058
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.996524e-01   1.000000e+00   9.477949e-01  7.818193e-248            NaN 
# [1] "chr20"
#           Reference
# Prediction    CT    TT
#         CT   564     0
#         TT     0 11114
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.996842e-01   1.000000e+00   9.517041e-01  8.819413e-252            NaN 
# [1] "chr3"
#           Reference
# Prediction    CT    TT
#         CT   517     0
#         TT     0 10070
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.996516e-01   1.000000e+00   9.511665e-01  6.341237e-231            NaN 
# [1] "chr4"
#           Reference
# Prediction    CT    TT
#         CT   656     0
#         TT     0 12301
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.997153e-01   1.000000e+00   9.493710e-01  4.342404e-293            NaN 
# [1] "chr5"
#           Reference
# Prediction   CT   TT
#         CT  446    0
#         TT    0 8844
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.996030e-01   1.000000e+00   9.519914e-01  3.167780e-199            NaN 
# [1] "chr6"
#           Reference
# Prediction    CT    TT
#         CT   551     0
#         TT     0 10560
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.996681e-01   1.000000e+00   9.504095e-01  3.686143e-246            NaN 
# [1] "chr7"
#           Reference
# Prediction   CT   TT
#         CT  552    0
#         TT    0 9959
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.996491e-01   1.000000e+00   9.474836e-01  5.560318e-247            NaN 
# [1] "chr8"
#           Reference
# Prediction   CT   TT
#         CT  464    0
#         TT    0 9040
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.996119e-01   1.000000e+00   9.511785e-01  2.523916e-207            NaN 
# [1] "chr9"
#           Reference
# Prediction    CT    TT
#         CT   592     0
#         TT     0 11389
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.996922e-01   1.000000e+00   9.505884e-01  2.130294e-264            NaN 

for (item in names(G14A_CHH_grl)) {
    assign(
        paste0("DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr", item),
        selectDIMP(
            eval(str2expression(
                paste0("PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
            )),
            div.col = 9,
            cutpoint = eval(str2expression(
                paste0(
                    "cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr",
                    item,
                    '$cutpoint'
                )
            )),
        )
    )
    
}

ls(pattern = "DIMPsYI")
#  [1] "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr1"  "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr10" "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr11"
#  [4] "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr12" "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr13" "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr14"
#  [7] "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr15" "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr16" "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr17"
# [10] "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr18" "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr19" "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr2" 
# [13] "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr20" "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr3"  "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr4" 
# [16] "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr5"  "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr6"  "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr7" 
# [19] "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr8"  "DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr9" 

unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr1, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2511      627      883    16000    19356    32053 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr2, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1724      462      723    11898    12948    20684 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr3, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1486      429      575    10079    12698    20423 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr4, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1920      544      717    12113    15452    24380 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr5, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1305      456      603    10284    11712    18372 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr6, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1820      472      771    12360    13616    21550 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr7, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1695      437      644    10984    12946    20250 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr8, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1403      386      536    10043    11018    17295 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr9, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1787      498      769    12779    14380    23895 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr10, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1894      504      800    13106    15035    24274 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr11, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1206      318      451     7551     8962    14459 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr12, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1306      387      619    10017    11116    17412 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr13, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1369      356      480     9004    10094    15182 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr14, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1956      463      615    12595    15828    25911 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr15, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1714      589      948    12476    13627    20363 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr16, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1341      387      614     9645    10901    17223 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr17, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1430      408      659    10525    11384    17190 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr18, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2120      649     1060    15465    16524    26572 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr19, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1843      584      870    13823    15354    25085 
unlist(lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr20, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    1878      419      838    12634    14662    23583 


DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    data.frame(
        "G1A_CHH" = 0,
        "G2A_CHH" = 0,
        "G3A_CHH" = 0,
        "G13A_CHH" = 0,
        "G14A_CHH" = 0,
        "G15A_CHH" = 0
    )
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr1, length
        )))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr2, length
        )))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr3, length
        )))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr4, length
        )))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr5, length
        )))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr6, length
        )))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr7, length
        )))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr8, length
        )))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr9, length
        )))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr10, length)
        ))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr11, length)
        ))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr12, length)
        ))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr13, length)
        ))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr14, length)
        ))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr15, length)
        ))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr16, length)
        ))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr17, length)
        ))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr18, length)
        ))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr19, length)
        ))
DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr20, length)
        ))

colSums(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS)
 # G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
 #   33708     9375    14175   233381   267613   426156

for (item in names(G14A_CHH_grl)) {
    assign(
        paste0("DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr", item),
        selectDIMP(
            eval(str2expression(
                paste0("PS_G123refctrl_vs_G131415_tvcut02_CHH_sum_chr", item)
            )),
            div.col = 9,
            cutpoint = eval(str2expression(
                paste0(
                    "cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr",
                    item,
                    '$cutpoint'
                )
            )),
        )
    )
    
}

ls(pattern = "DIMPsML")
#  [1] "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr1"  "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr10" "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr11"
#  [4] "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr12" "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr13" "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr14"
#  [7] "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr15" "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr16" "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr17"
# [10] "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr18" "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr19" "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr2" 
# [13] "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr20" "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr3"  "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr4" 
# [16] "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr5"  "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr6"  "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr7" 
# [19] "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr8"  "DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr9" 

unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr1, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#   10760     7113     9113    29131    19356    32053 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr2, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6013     3512     5754    20275    12948    20684 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr3, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6028     3988     6028    18734    12698    20423 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr4, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    7977     5362     7288    22122    15452    24380 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr5, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    4978     3179     5018    18012    11712    18372 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr6, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6166     3633     5856    21100    13616    21550 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr7, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6059     3726     5750    19654    12946    20250 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr8, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    5871     3772     5089    16439    11018    17295 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr9, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6634     4073     6631    22990    14380    23895 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr10, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    3970     1702     2944    22857    15035    24274 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr11, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    2509     1157     1938    13872     8962    14459 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr12, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    5398     3593     4812    17047    11116    17412 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr13, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    4812     2714     4163    14988    10094    15182 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr14, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    7614     5047     7770    24032    15828    25911 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr15, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    5470     3261     4992    18947    13627    20363 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr16, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    5585     3667     4703    16598    10901    17223 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr17, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    4677     2617     4287    16690    11384    17190 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr18, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6103     3391     5555    24632    16524    26604 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr19, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6607     4161     6976    24012    15354    25085 
unlist(lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr20, length))
# G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
#    6515     3727     6503    22607    14662    23583 


DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    data.frame(
        "G1A_CHH" = 0,
        "G2A_CHH" = 0,
        "G3A_CHH" = 0,
        "G13A_CHH" = 0,
        "G14A_CHH" = 0,
        "G15A_CHH" = 0
    )
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr1, length
        )))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr2, length
        )))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr3, length
        )))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr4, length
        )))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr5, length
        )))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr6, length
        )))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr7, length
        )))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr8, length
        )))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(lapply(
            DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr9, length
        )))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr10, length)
        ))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr11, length)
        ))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr12, length)
        ))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr13, length)
        ))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr14, length)
        ))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr15, length)
        ))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr16, length)
        ))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr17, length)
        ))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr18, length)
        ))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr19, length)
        ))
DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS <-
    rbind(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS,
        unlist(
            lapply(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr20, length)
        ))

colSums(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS)
 # G1A_CHH  G2A_CHH  G3A_CHH G13A_CHH G14A_CHH G15A_CHH 
 #  119746    73395   111170   404739   267613   426188

pryr::mem_used()
# 46.4 GB

length(ls(pattern = "best"))
# 20
length(ls(pattern = "critical"))
# 20
length(ls(pattern = "cutpointsML"))
# 20
length(ls(pattern = "cutpointsYI"))
# 20
length(ls(pattern = "DIMPsML"))
# 21
length(ls(pattern = "DIMPsYI"))
# 21
length(ls(pattern = "HD"))
# 20
length(ls(pattern = "PS"))
# 60
length(ls(pattern = "^PS_"))
# 20

## 20200110_1530pm saved:
#
# save(control_names_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/control_names_CHH.RData")
# save(treatment_names_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/treatment_names_CHH.RData")
# save(covr, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/covr.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr1.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr10.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr11.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr12.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr13.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr14.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr15.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr16.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr17.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr18.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr19.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr2.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr20.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr3.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr4.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr5.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr6.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr7.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr8.RData")
# save(critical.val_G123refctrl_vs_G131415_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/critical.val_G123refctrl_vs_G131415_sum_chr9.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr1.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr10.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr11.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr12.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr13.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr14.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr15.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr16.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr17.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr18.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr19.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr2.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr20.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr3.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr4.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr5.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr6.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr7.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr8.RData")
# save(cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsML_PS_G123refctrl_vs_G131415_CHH_sum_chr9.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr1.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr10.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr11.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr12.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr13.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr14.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr15.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr16.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr17.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr18.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr19.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr2.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr20.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr3.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr4.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr5.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr6.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr7.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr8.RData")
# save(cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/cutpointsYI_PS_G123refctrl_vs_G131415_CHH_sum_chr9.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr1.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr10.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr11.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr12.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr13.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr14.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr15.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr16.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr17.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr18.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr19.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr2.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr20.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr3.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr4.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr5.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr6.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr7.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr8.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_chr9.RData")
# save(DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsML_G123refctrl_vs_G131415_CHH_sum_COUNTS.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr1.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr10.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr11.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr12.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr13.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr14.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr15.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr16.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr17.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr18.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr19.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr2.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr20.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr3.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr4.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr5.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr6.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr7.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr8.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_chr9.RData")
# save(DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/DIMPsYI_G123refctrl_vs_G131415_CHH_sum_COUNTS.RData")
# save(F4_P37_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/F4_P37_CHH.RData")
# save(F4_WT_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G123refctrl_vs_G131415_sepChroms_CHH_20200110_R/F4_WT_CHH.RData")





date()
# [1] "Fri Jan 10 16:20:07 2020"

