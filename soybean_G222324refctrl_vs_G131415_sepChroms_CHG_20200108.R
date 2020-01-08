
.libPaths()

library(BiocManager)
library(MethylIT)
library(MethylIT.utils)

rm(list=ls())
length(ls()) 
# [1] 0

# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F4_P37_CG_4_29_2019.RData")
# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F6_P37_CG_4_29_2019.RData")

load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F4_P37_CHG_4_29_2019.RData")
load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F6_P37_CHG_4_29_2019.RData")

names(F4_P37_CHG)
# [1] "G13A_CHG" "G14A_CHG" "G15A_CHG"

names(F6_P37_CHG)
# [1] "G22A_CHG" "G23A_CHG" "G24A_CHG"

ls()
# [1] "F4_P37_CHG" "F6_P37_CHG"

seqnames(F4_P37_CHG$G14A_CHG)
# factor-Rle of length 38351000 with 20 runs
#   Lengths: 2280891 2117733 1429711 1618695 1865889 1965855 ... 2158782 1742307 2066290 1805914 1952798 2050573
#   Values :       1      10      11      12      13      14 ...       4       5       6       7       8       9
# Levels(20): 1 10 11 12 13 14 15 16 17 18 19 2 20 3 4 5 6 7 8 9

G22A_CHG_grl <- split(F6_P37_CHG$G22A_CHG,seqnames(F4_P37_CHG$G14A_CHG))
G23A_CHG_grl <- split(F6_P37_CHG$G23A_CHG,seqnames(F4_P37_CHG$G14A_CHG))
G24A_CHG_grl <- split(F6_P37_CHG$G24A_CHG,seqnames(F4_P37_CHG$G14A_CHG))

G13A_CHG_grl <- split(F4_P37_CHG$G13A_CHG,seqnames(F4_P37_CHG$G13A_CHG))
G14A_CHG_grl <- split(F4_P37_CHG$G14A_CHG,seqnames(F4_P37_CHG$G13A_CHG))
G15A_CHG_grl <- split(F4_P37_CHG$G15A_CHG,seqnames(F4_P37_CHG$G13A_CHG))

names(G22A_CHG_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

names(G14A_CHG_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 


date()
# [1] "Wed Jan  8 06:53:32 2020"

for(item in names(G14A_CHG_grl)) {
    assign(
        paste0("soy_WT_G222324_Ref_CHG_sum_chr", item),
        poolFromGRlist(
            list(G22A_CHG_grl[[item]], G23A_CHG_grl[[item]], G24A_CHG_grl[[item]]),
            stat = "sum",
            num.cores = 12L
        ))}

date()
# [1] "Wed Jan  8 06:57:55 2020"

ls()
#  [1] "F4_P37_CHG"                       "F6_P37_CHG"                       "G13A_CHG_grl"                    
#  [4] "G14A_CHG_grl"                     "G15A_CHG_grl"                     "G22A_CHG_grl"                    
#  [7] "G23A_CHG_grl"                     "G24A_CHG_grl"                     "item"                            
# [10] "soy_WT_G222324_Ref_CHG_sum_chr1"  "soy_WT_G222324_Ref_CHG_sum_chr10" "soy_WT_G222324_Ref_CHG_sum_chr11"
# [13] "soy_WT_G222324_Ref_CHG_sum_chr12" "soy_WT_G222324_Ref_CHG_sum_chr13" "soy_WT_G222324_Ref_CHG_sum_chr14"
# [16] "soy_WT_G222324_Ref_CHG_sum_chr15" "soy_WT_G222324_Ref_CHG_sum_chr16" "soy_WT_G222324_Ref_CHG_sum_chr17"
# [19] "soy_WT_G222324_Ref_CHG_sum_chr18" "soy_WT_G222324_Ref_CHG_sum_chr19" "soy_WT_G222324_Ref_CHG_sum_chr2" 
# [22] "soy_WT_G222324_Ref_CHG_sum_chr20" "soy_WT_G222324_Ref_CHG_sum_chr3"  "soy_WT_G222324_Ref_CHG_sum_chr4" 
# [25] "soy_WT_G222324_Ref_CHG_sum_chr5"  "soy_WT_G222324_Ref_CHG_sum_chr6"  "soy_WT_G222324_Ref_CHG_sum_chr7" 
# [28] "soy_WT_G222324_Ref_CHG_sum_chr8"  "soy_WT_G222324_Ref_CHG_sum_chr9" 

length(ls(pattern="soy_WT_G222324_Ref_CHG_sum_chr"))
# [1] 20

## 20190108 saved:
#
# save(F4_P37_CHG,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/F4_P37_CHG.RData")
# save(F6_P37_CHG,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/F6_P37_CHG.RData")
# save(G13A_CHG_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G13A_CHG_grl.RData")
# save(G14A_CHG_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G14A_CHG_grl.RData")
# save(G15A_CHG_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G15A_CHG_grl.RData")
# save(G22A_CHG_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G22A_CHG_grl.RData")
# save(G23A_CHG_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G23A_CHG_grl.RData")
# save(G24A_CHG_grl,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G24A_CHG_grl.RData")
# save(item,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/item.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr1.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr10.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr11.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr12.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr13.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr14.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr15.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr16.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr17.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr18.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr19.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr2.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr20.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr3.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr4.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr5.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr6.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr7.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr8.RData")
# save(soy_WT_G222324_Ref_CHG_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/soy_WT_G222324_Ref_CHG_sum_chr9.RData")

date()
# [1] "Wed Jan  8 07:10:24 2020"

## 20190107: 61 minute runtime:
#
for(item in names(G14A_CHG_grl)) {
    assign(
        paste0(
            "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr",
            item
        ),
        estimateDivergence(
            ref = eval(str2expression(paste0("soy_WT_G222324_Ref_CHG_sum_chr",item))),
            indiv = list(G22A_CHG_grl[[item]], G23A_CHG_grl[[item]], G24A_CHG_grl[[item]],G13A_CHG_grl[[item]], G14A_CHG_grl[[item]], G15A_CHG_grl[[item]]),
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
# [1] "Wed Jan  8 07:53:34 2020"


length(ls(pattern="HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr"))
# [1] 20

names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr1)
# NULL

names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr1) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr2) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr3) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr4) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr5) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr6) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr7) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr8) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr9) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr10) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr11) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr12) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr13) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr14) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr15) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr16) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr17) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr18) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr19) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr20) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")

names(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr1)
# [1] "G22A_CHG" "G23A_CHG" "G24A_CHG" "G13A_CHG" "G14A_CHG" "G15A_CHG"


for(item in names(G14A_CHG_grl)) {
    assign(paste0("critical.val_G222324refctrl_vs_G131415_sum_chr", item),
    do.call(rbind, lapply(eval(str2expression(paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr",item))), function(x) {
    hd.95 = quantile(x$hdiv, 0.95)
    baytv.95 = quantile(abs(x$bay.TV), 0.95)
    tv.95 = quantile(abs(x$TV), 0.95)
    return(c(tv = tv.95, baytv = baytv.95, hd = hd.95, num.sites.hd95 = sum(x$hdiv > hd.95), num.sites.tv95 = sum(x$TV > tv.95), num.sites.baytv95 = sum(x$bay.TV > baytv.95)))}))
)
}
    
for(item in names(G1A_CHG_grl)) {
 print(eval(str2expression(paste0("critical.val_G222324refctrl_vs_G131415_sum_chr",item))))
}
## 20190107: 
#

date()

## 20190107: almost 3 hours runtime:
#


names(G14A_CHG_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

names(G14A_CHG_grl[2:20])
# [1] "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 

date()
# [1] "Wed Jan  8 07:54:00 2020"

# for(item in names(G14A_CHG_grl[2:20])) {
for(item in names(G14A_CHG_grl)) {
    assign(
        paste0("bestFits_G222324refctrl_vs_G131415_CHG_sum_chr", item),
        gofReport(
            eval(str2expression(
                paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr", item)
            )),
            model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
            column = 9,
            absolute = FALSE,
            output = c("best.model"),
            num.cores = 64L,
            verbose = FALSE
        )
    )
    print(eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHG_sum_chr", item,"$bestModel"))))
}

#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val   g2p_AIC g2p_R.Cross.val   g3p_AIC g3p_R.Cross.val
# G22A_CHG -7584962       0.9991041      NA              NA -10565043       0.9999027 -10621414       0.9999052
# G23A_CHG -7658858       0.9993426      NA              NA -12359669       0.9999836 -12374905       0.9999837
# G24A_CHG -5948426       0.9993763      NA              NA  -9892953       0.9999891  -9909137       0.9999892
# G13A_CHG -7060058       0.9995098      NA              NA -13625797       0.9999986 -13625480       0.9999986
# G14A_CHG -9135621       0.9996111      NA              NA -13533507       0.9999831 -13533728       0.9999830
# G15A_CHG -7703284       0.9996106      NA              NA -11337459       0.9999819 -11347299       0.9999818
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -14777629        0.9999970      gg3p
# G23A_CHG -13507179        0.9999948      gg3p
# G24A_CHG -11686867        0.9999983      gg3p
# G13A_CHG -14061637        0.9999992      gg3p
# G14A_CHG -17005264        0.9999989      gg3p
# G15A_CHG -13147132        0.9999967      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val   g2p_AIC g2p_R.Cross.val   g3p_AIC g3p_R.Cross.val
# G22A_CHG -6207278       0.9990687      NA              NA  -8572331       0.9998904  -8615457       0.9998927
# G23A_CHG -6549004       0.9993885      NA              NA -11098214       0.9999911 -11108088       0.9999912
# G24A_CHG -5807955       0.9994318      NA              NA -10378117       0.9999956 -10395774       0.9999957
# G13A_CHG -6299351       0.9995235      NA              NA -11664468       0.9999978 -11664799       0.9999978
# G14A_CHG -7384598       0.9996041      NA              NA -10941902       0.9999825 -10942073       0.9999824
# G15A_CHG -6635636       0.9995966      NA              NA  -9987025       0.9999848  -9990042       0.9999847
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -11553935        0.9999946      gg3p
# G23A_CHG -11592814        0.9999959      gg3p
# G24A_CHG -11086541        0.9999980      gg3p
# G13A_CHG -12694628        0.9999993      gg3p
# G14A_CHG -12550528        0.9999964      gg3p
# G15A_CHG -11948079        0.9999978      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -3253615       0.9990527      NA              NA -4482697       0.9998855 -4505044       0.9998878
# G23A_CHG -3272248       0.9993538      NA              NA -5285612       0.9999840 -5295986       0.9999841
# G24A_CHG -2499663       0.9993662      NA              NA -4079062       0.9999865 -4086666       0.9999866
# G13A_CHG -2899690       0.9995044      NA              NA -5561219       0.9999985 -5560277       0.9999985
# G14A_CHG -3974583       0.9995855      NA              NA -6144410       0.9999876 -6145293       0.9999876
# G15A_CHG -3183358       0.9996069      NA              NA -4650002       0.9999805 -4653996       0.9999805
#          gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -6636108        0.9999974      gg3p
# G23A_CHG -6403137        0.9999981      gg3p
# G24A_CHG -4829581        0.9999978      gg3p
# G13A_CHG -5717369        0.9999990      gg3p
# G14A_CHG -8271638        0.9999996      gg3p
# G15A_CHG -5627857        0.9999975      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -4510032       0.9990798      NA              NA -6246649       0.9998955 -6275134       0.9998974
# G23A_CHG -4758658       0.9993588      NA              NA -7843940       0.9999873 -7854389       0.9999874
# G24A_CHG -4799678       0.9994792      NA              NA -8645162       0.9999965 -8646614       0.9999965
# G13A_CHG -4861159       0.9995153      NA              NA -9179803       0.9999983 -9179114       0.9999983
# G14A_CHG -5302931       0.9995947      NA              NA -8031562       0.9999852 -8030556       0.9999851
# G15A_CHG -5106301       0.9995988      NA              NA -7701280       0.9999852 -7703810       0.9999851
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG  -9042008        0.9999971      gg3p
# G23A_CHG  -9473596        0.9999985      gg3p
# G24A_CHG  -8579668        0.9999963       g3p
# G13A_CHG  -9697130        0.9999992      gg3p
# G14A_CHG -10188590        0.9999990      gg3p
# G15A_CHG  -9365278        0.9999983      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P"  "Gamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -3366827       0.9990684      NA              NA -4625980       0.9998875 -4653944       0.9998905
# G23A_CHG -3474090       0.9993463      NA              NA -5601215       0.9999839 -5611444       0.9999840
# G24A_CHG -2899186       0.9994257      NA              NA -4979381       0.9999931 -4987640       0.9999932
# G13A_CHG -3357144       0.9995661      NA              NA -5577605       0.9999937 -5579175       0.9999937
# G14A_CHG -4094377       0.9996261      NA              NA -5865945       0.9999775 -5864843       0.9999774
# G15A_CHG -3569335       0.9996207      NA              NA -5088294       0.9999764 -5091644       0.9999763
#          gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -6532193        0.9999956      gg3p
# G23A_CHG -6530187        0.9999968      gg3p
# G24A_CHG -5275230        0.9999963      gg3p
# G13A_CHG -6791857        0.9999994      gg3p
# G14A_CHG -7618375        0.9999987      gg3p
# G15A_CHG -6138223        0.9999966      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val   g2p_AIC g2p_R.Cross.val   g3p_AIC g3p_R.Cross.val
# G22A_CHG -6320572       0.9990690      NA              NA  -8735045       0.9998911  -8775279       0.9998932
# G23A_CHG -6351149       0.9993438      NA              NA -10201874       0.9999834 -10215505       0.9999834
# G24A_CHG -4393337       0.9993210      NA              NA  -6934862       0.9999788  -6953037       0.9999790
# G13A_CHG -5530905       0.9995130      NA              NA -10667713       0.9999986 -10668258       0.9999986
# G14A_CHG -7712341       0.9995999      NA              NA -11675898       0.9999856 -11675191       0.9999855
# G15A_CHG -6075825       0.9996024      NA              NA  -9148909       0.9999849  -9159673       0.9999850
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -11982930        0.9999953      gg3p
# G23A_CHG -11392519        0.9999959      gg3p
# G24A_CHG  -8760876        0.9999983      gg3p
# G13A_CHG -11204066        0.9999993      gg3p
# G14A_CHG -13935470        0.9999982      gg3p
# G15A_CHG -11898226        0.9999993      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -5071649       0.9991199      NA              NA -7171930       0.9999158 -7204331       0.9999176
# G23A_CHG -5334611       0.9993987      NA              NA -9026684       0.9999911 -9041908       0.9999912
# G24A_CHG -5310699       0.9994798      NA              NA -9374926       0.9999956 -9377533       0.9999956
# G13A_CHG -5445786       0.9995504      NA              NA -9288470       0.9999950 -9286821       0.9999950
# G14A_CHG -6069662       0.9996199      NA              NA -8747953       0.9999779 -8746280       0.9999778
# G15A_CHG -5702014       0.9996316      NA              NA -8045058       0.9999752 -8052954       0.9999751
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -10365992        0.9999979      gg3p
# G23A_CHG -10274857        0.9999979      gg3p
# G24A_CHG  -9325074        0.9999954       g3p
# G13A_CHG -10432936        0.9999988      gg3p
# G14A_CHG -11069088        0.9999983      gg3p
# G15A_CHG  -9858501        0.9999970      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P"  "Gamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -4065653       0.9990776      NA              NA -5631954       0.9998942 -5659429       0.9998966
# G23A_CHG -4159060       0.9993614      NA              NA -6823126       0.9999872 -6828817       0.9999871
# G24A_CHG -3612630       0.9994035      NA              NA -6211672       0.9999928 -6216556       0.9999928
# G13A_CHG -4080413       0.9995276      NA              NA -7340206       0.9999973 -7339221       0.9999973
# G14A_CHG -5072745       0.9996780      NA              NA -6636697       0.9999571 -6637467       0.9999569
# G15A_CHG -4221884       0.9996434      NA              NA -5816562       0.9999695 -5823858       0.9999695
#          gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -8214393        0.9999974      gg3p
# G23A_CHG -7997876        0.9999977      gg3p
# G24A_CHG -6900998        0.9999977      gg3p
# G13A_CHG -7868509        0.9999987      gg3p
# G14A_CHG -8519494        0.9999967      gg3p
# G15A_CHG -7235163        0.9999969      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -4594919       0.9990796      NA              NA -6376621       0.9998973 -6404360       0.9998990
# G23A_CHG -4909600       0.9993832      NA              NA -8176378       0.9999895 -8184770       0.9999895
# G24A_CHG -5098074       0.9994981      NA              NA -8882763       0.9999953 -8882144       0.9999953
# G13A_CHG -5116288       0.9995491      NA              NA -8580192       0.9999943 -8579835       0.9999942
# G14A_CHG -5435042       0.9996226      NA              NA -7885315       0.9999797 -7887394       0.9999796
# G15A_CHG -5285375       0.9996164      NA              NA -7695882       0.9999803 -7699481       0.9999802
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG  -8889645        0.9999956      gg3p
# G23A_CHG  -9212617        0.9999972      gg3p
# G24A_CHG  -8919610        0.9999955      gg3p
# G13A_CHG  -9400780        0.9999980      gg3p
# G14A_CHG -10870812        0.9999995      gg3p
# G15A_CHG  -9290131        0.9999973      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val   g2p_AIC g2p_R.Cross.val   g3p_AIC g3p_R.Cross.val
# G22A_CHG -6604636       0.9991521      NA              NA  -9385766       0.9999237  -9423526       0.9999251
# G23A_CHG -6831606       0.9993840      NA              NA -11457914       0.9999901 -11470707       0.9999902
# G24A_CHG -6384099       0.9994475      NA              NA -11316665       0.9999954 -11331020       0.9999954
# G13A_CHG -6930848       0.9995688      NA              NA -11371787       0.9999931 -11372033       0.9999930
# G14A_CHG -7910458       0.9996358      NA              NA -11162411       0.9999747 -11161221       0.9999746
# G15A_CHG -7146322       0.9996045      NA              NA -10628082       0.9999834 -10638754       0.9999834
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -11961071        0.9999942      gg3p
# G23A_CHG -12026371        0.9999955      gg3p
# G24A_CHG -10472177        0.9999917       g3p
# G13A_CHG -12330868        0.9999975      gg3p
# G14A_CHG -13748941        0.9999979      gg3p
# G15A_CHG -12205130        0.9999964      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P"  "Gamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val   g2p_AIC g2p_R.Cross.val   g3p_AIC g3p_R.Cross.val
# G22A_CHG -6379235       0.9990795      NA              NA  -8855527       0.9998960  -8897220       0.9998982
# G23A_CHG -6603577       0.9993604      NA              NA -10864417       0.9999874 -10876315       0.9999874
# G24A_CHG -5644896       0.9993955      NA              NA  -9571852       0.9999913  -9587606       0.9999913
# G13A_CHG -6350628       0.9995173      NA              NA -12100896       0.9999984 -12100430       0.9999984
# G14A_CHG -7606258       0.9995959      NA              NA -11668792       0.9999871 -11667927       0.9999871
# G15A_CHG -6821112       0.9996147      NA              NA -10044760       0.9999819 -10051574       0.9999819
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -12300825        0.9999963      gg3p
# G23A_CHG -11631331        0.9999953      gg3p
# G24A_CHG -10718191        0.9999975      gg3p
# G13A_CHG -12723805        0.9999992      gg3p
# G14A_CHG -14420193        0.9999991      gg3p
# G15A_CHG -12639134        0.9999989      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -5047953       0.9990376      NA              NA -6923402       0.9998803 -6958202       0.9998827
# G23A_CHG -5343060       0.9993866      NA              NA -8936940       0.9999899 -8952518       0.9999899
# G24A_CHG -4632113       0.9994430      NA              NA -8307278       0.9999960 -8313948       0.9999960
# G13A_CHG -5210122       0.9995401      NA              NA -9314369       0.9999971 -9314915       0.9999970
# G14A_CHG -6008883       0.9995840      NA              NA -9300267       0.9999877 -9300206       0.9999876
# G15A_CHG -5558764       0.9996101      NA              NA -8229604       0.9999828 -8232571       0.9999827
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG  -9984295        0.9999963      gg3p
# G23A_CHG -10203068        0.9999977      gg3p
# G24A_CHG  -8620804        0.9999973      gg3p
# G13A_CHG -10365705        0.9999992      gg3p
# G14A_CHG -12616250        0.9999997      gg3p
# G15A_CHG -10061287        0.9999980      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val   g2p_AIC g2p_R.Cross.val   g3p_AIC g3p_R.Cross.val
# G22A_CHG -6302137       0.9990833      NA              NA  -8745568       0.9998967  -8791754       0.9998992
# G23A_CHG -6743984       0.9994175      NA              NA -11687383       0.9999935 -11703364       0.9999936
# G24A_CHG -6040808       0.9994360      NA              NA -10590252       0.9999946 -10605176       0.9999946
# G13A_CHG -6717068       0.9995463      NA              NA -11759625       0.9999963 -11760562       0.9999963
# G14A_CHG -7572012       0.9996464      NA              NA -10555022       0.9999730 -10558363       0.9999728
# G15A_CHG -7119087       0.9996101      NA              NA -10511747       0.9999824 -10517373       0.9999823
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -12064918        0.9999962      gg3p
# G23A_CHG -11750518        0.9999956      gg3p
# G24A_CHG -11050505        0.9999966      gg3p
# G13A_CHG -12388091        0.9999983      gg3p
# G14A_CHG -12220131        0.9999953      gg3p
# G15A_CHG -12042535        0.9999966      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -4778130       0.9990473      NA              NA -6573759       0.9998836 -6607866       0.9998861
# G23A_CHG -4856946       0.9993761      NA              NA -8007488       0.9999876 -8028202       0.9999879
# G24A_CHG -3534133       0.9993571      NA              NA -5732327       0.9999847 -5745318       0.9999849
# G13A_CHG -4320385       0.9995039      NA              NA -8506818       0.9999989 -8503550       0.9999989
# G14A_CHG -5861048       0.9996004      NA              NA -8843176       0.9999848 -8842244       0.9999848
# G15A_CHG -4690291       0.9996063      NA              NA -6971742       0.9999837 -6980759       0.9999836
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG  -9825359        0.9999976      gg3p
# G23A_CHG  -9601751        0.9999984      gg3p
# G24A_CHG  -6901183        0.9999980      gg3p
# G13A_CHG  -8710279        0.9999992      gg3p
# G14A_CHG -11522745        0.9999993      gg3p
# G15A_CHG  -8564124        0.9999982      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val   g2p_AIC g2p_R.Cross.val   g3p_AIC g3p_R.Cross.val
# G22A_CHG -5660203       0.9991001      NA              NA  -7885770       0.9999022  -7924783       0.9999044
# G23A_CHG -5591552       0.9993558      NA              NA  -9091370       0.9999850  -9104793       0.9999851
# G24A_CHG -3789241       0.9993050      NA              NA  -5901792       0.9999748  -5922160       0.9999753
# G13A_CHG -4850763       0.9995219      NA              NA  -8880606       0.9999976  -8877631       0.9999975
# G14A_CHG -7010975       0.9996264      NA              NA -10023330       0.9999768 -10022643       0.9999767
# G15A_CHG -5190313       0.9996089      NA              NA  -7657160       0.9999823  -7662949       0.9999822
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -11591387        0.9999978      gg3p
# G23A_CHG -10699219        0.9999975      gg3p
# G24A_CHG  -7382341        0.9999976      gg3p
# G13A_CHG  -9532666        0.9999990      gg3p
# G14A_CHG -11771458        0.9999967      gg3p
# G15A_CHG  -9308077        0.9999978      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -4892227       0.9990476      NA              NA -6718856       0.9998827 -6752674       0.9998851
# G23A_CHG -5199249       0.9993920      NA              NA -8896409       0.9999922 -8908066       0.9999922
# G24A_CHG -4769788       0.9994727      NA              NA -9123809       0.9999983 -9134257       0.9999983
# G13A_CHG -5107991       0.9994996      NA              NA -9866711       0.9999986 -9866673       0.9999986
# G14A_CHG -5815419       0.9995985      NA              NA -8747623       0.9999845 -8746854       0.9999845
# G15A_CHG -5453888       0.9995973      NA              NA -8232813       0.9999853 -8235463       0.9999853
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG  -9877556        0.9999971      gg3p
# G23A_CHG -10610140        0.9999990      gg3p
# G24A_CHG  -9259227        0.9999985      gg3p
# G13A_CHG -10026890        0.9999989      gg3p
# G14A_CHG -11314926        0.9999992      gg3p
# G15A_CHG -10164176        0.9999986      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -5580865       0.9990976      NA              NA -7794076       0.9999039 -7833186       0.9999061
# G23A_CHG -5854334       0.9993757      NA              NA -9749860       0.9999892 -9761694       0.9999892
# G24A_CHG -5324302       0.9994569      NA              NA -9651441       0.9999966 -9660867       0.9999966
# G13A_CHG -5870100       0.9995577      NA              NA -9877566       0.9999944 -9877052       0.9999944
# G14A_CHG -6607645       0.9996149      NA              NA -9748710       0.9999821 -9747543       0.9999820
# G15A_CHG -6208664       0.9996190      NA              NA -9091442       0.9999814 -9101456       0.9999814
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -11459334        0.9999978      gg3p
# G23A_CHG -11328254        0.9999979      gg3p
# G24A_CHG  -9914203        0.9999974      gg3p
# G13A_CHG -11352405        0.9999989      gg3p
# G14A_CHG -12594717        0.9999990      gg3p
# G15A_CHG -11554250        0.9999987      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -4856758       0.9990634      NA              NA -6690193       0.9998874 -6724285       0.9998899
# G23A_CHG -5012348       0.9993573      NA              NA -8210261       0.9999864 -8221775       0.9999864
# G24A_CHG -3784358       0.9993824      NA              NA -6269899       0.9999884 -6280759       0.9999885
# G13A_CHG -4528626       0.9994931      NA              NA -8966465       0.9999990 -8966440       0.9999990
# G14A_CHG -5838589       0.9996012      NA              NA -8811090       0.9999850 -8811082       0.9999850
# G15A_CHG -5023135       0.9995880      NA              NA -7703313       0.9999867 -7711393       0.9999867
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG  -9849378        0.9999973      gg3p
# G23A_CHG  -9877908        0.9999982      gg3p
# G24A_CHG  -7315315        0.9999979      gg3p
# G13A_CHG  -9064043        0.9999991      gg3p
# G14A_CHG -11863796        0.9999995      gg3p
# G15A_CHG  -9534108        0.9999988      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val  g2p_AIC g2p_R.Cross.val  g3p_AIC g3p_R.Cross.val
# G22A_CHG -3776252       0.9990819      NA              NA -5229628       0.9998955 -5254887       0.9998975
# G23A_CHG -3793377       0.9993634      NA              NA -6133882       0.9999847 -6145515       0.9999849
# G24A_CHG -2917365       0.9993669      NA              NA -4735627       0.9999855 -4742363       0.9999856
# G13A_CHG -3569620       0.9995536      NA              NA -5908362       0.9999934 -5910048       0.9999934
# G14A_CHG -4788458       0.9997015      NA              NA -5971034       0.9999403 -5971959       0.9999399
# G15A_CHG -3748434       0.9996164      NA              NA -5462400       0.9999801 -5468780       0.9999801
#          gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -7659342        0.9999975      gg3p
# G23A_CHG -7153329        0.9999970      gg3p
# G24A_CHG -5212400        0.9999947      gg3p
# G13A_CHG -6698610        0.9999984      gg3p
# G14A_CHG -7535924        0.9999939      gg3p
# G15A_CHG -7056799        0.9999988      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 
#   |=====================================================================================================| 100%
# 
#  *** Creating report ... 
#           w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val   g2p_AIC g2p_R.Cross.val   g3p_AIC g3p_R.Cross.val
# G22A_CHG -5929155       0.9990692      NA              NA  -8192131       0.9998906  -8236452       0.9998933
# G23A_CHG -6127103       0.9993710      NA              NA -10106772       0.9999882 -10114831       0.9999882
# G24A_CHG -5049916       0.9994114      NA              NA  -8608076       0.9999921  -8623015       0.9999922
# G13A_CHG -5613624       0.9995091      NA              NA -11057270       0.9999989 -11055849       0.9999989
# G14A_CHG -7157810       0.9995875      NA              NA -11064719       0.9999876 -11063022       0.9999876
# G15A_CHG -6164888       0.9996019      NA              NA  -9227821       0.9999838  -9241123       0.9999839
#           gg3p_AIC gg3p_R.Cross.val bestModel
# G22A_CHG -11442308        0.9999964      gg3p
# G23A_CHG -11937765        0.9999980      gg3p
# G24A_CHG  -9368499        0.9999969      gg3p
# G13A_CHG -11788008        0.9999995      gg3p
# G14A_CHG -13692575        0.9999992      gg3p
# G15A_CHG -11682033        0.9999989      gg3p
# 
#   G22A_CHG   G23A_CHG   G24A_CHG   G13A_CHG   G14A_CHG   G15A_CHG 
# "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" "GGamma3P" 

date()
# [1] "Wed Jan  8 12:08:31 2020"


# TODO 20190108

ls()
#  [1] "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr1"                        "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr10"                      
#  [3] "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr11"                       "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr12"                      
#  [5] "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr13"                       "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr14"                      
#  [7] "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr15"                       "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr16"                      
#  [9] "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr17"                       "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr18"                      
# [11] "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr19"                       "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr2"                       
# [13] "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr20"                       "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr3"                       
# [15] "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr4"                        "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr5"                       
# [17] "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr6"                        "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr7"                       
# [19] "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr8"                        "bestFits_G222324refctrl_vs_G131415_CHG_sum_chr9"                       
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
# [41] "F4_P37_CHG"                                                             "F6_P37_CHG"                                                            
# [43] "G13A_CHG_grl"                                                           "G14A_CHG_grl"                                                          
# [45] "G15A_CHG_grl"                                                           "G22A_CHG_grl"                                                          
# [47] "G23A_CHG_grl"                                                           "G24A_CHG_grl"                                                          
# [49] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr1"  "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr10"
# [51] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr11" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr12"
# [53] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr13" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr14"
# [55] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr15" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr16"
# [57] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr17" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr18"
# [59] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr19" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr2" 
# [61] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr20" "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr3" 
# [63] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr4"  "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr5" 
# [65] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr6"  "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr7" 
# [67] "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr8"  "HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr9" 
# [69] "item"                                                                   "soy_WT_G222324_Ref_CHG_sum_chr1"                                       
# [71] "soy_WT_G222324_Ref_CHG_sum_chr10"                                       "soy_WT_G222324_Ref_CHG_sum_chr11"                                      
# [73] "soy_WT_G222324_Ref_CHG_sum_chr12"                                       "soy_WT_G222324_Ref_CHG_sum_chr13"                                      
# [75] "soy_WT_G222324_Ref_CHG_sum_chr14"                                       "soy_WT_G222324_Ref_CHG_sum_chr15"                                      
# [77] "soy_WT_G222324_Ref_CHG_sum_chr16"                                       "soy_WT_G222324_Ref_CHG_sum_chr17"                                      
# [79] "soy_WT_G222324_Ref_CHG_sum_chr18"                                       "soy_WT_G222324_Ref_CHG_sum_chr19"                                      
# [81] "soy_WT_G222324_Ref_CHG_sum_chr2"                                        "soy_WT_G222324_Ref_CHG_sum_chr20"                                      
# [83] "soy_WT_G222324_Ref_CHG_sum_chr3"                                        "soy_WT_G222324_Ref_CHG_sum_chr4"                                       
# [85] "soy_WT_G222324_Ref_CHG_sum_chr5"                                        "soy_WT_G222324_Ref_CHG_sum_chr6"                                       
# [87] "soy_WT_G222324_Ref_CHG_sum_chr7"                                        "soy_WT_G222324_Ref_CHG_sum_chr8"                                       
# [89] "soy_WT_G222324_Ref_CHG_sum_chr9" 

## 20200107 and 20200108 saved:
#
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr1.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr10.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr11.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr12.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr13.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr14.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr15.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr16.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr17.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr18.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr19.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr2.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr20.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr3.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr4.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr5.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr6.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr7.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr8.RData")
save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr9.RData")

save(critical.val_G222324refctrl_vs_G131415_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr1.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr10.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr11.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr12.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr13.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr14.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr15.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr16.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr17.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr18.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr19.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr2.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr20.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr3.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr4.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr5.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr6.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr7.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr8.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr9.RData")

save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr1.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr10.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr11.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr12.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr13.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr14.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr15.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr16.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr17.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr18.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr19.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr2.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr20.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr3.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr4.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr5.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr6.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr7.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr8.RData")
save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr9.RData")

















































date()

covr <-
    lapply(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr1, function(x) {
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

## 20190107: we see:
#   - num.siteGreater_500 are all 0, confirming that our "high.coverage=300" filter in estimateDivergence worked
#   - 1st quartile values of about 10, showing that our data has good coverage, >10, for 75% of sites
#   - our data is extremely similar for all 6 samples in terms of coverage
# See also soybean_G222324refctrl_vs_G1415_sepChroms_20200102.R for more discussion.
#
#          Min. 1st Qu. Median Mean 3rd Qu. Max. 60% 95% 99% 99.9% 99.99% num.siteGreater_8 q60_to_500 num.siteGreater_500
# G22A_CHG    0      11     19   21      29   74  22  44  53    61     66           1121143     551216                   0
# G23A_CHG    0      10     16   17      24   61  19  35  41    46     51           1054091     518053                   0
# G24A_CHG    0       9     12   12      16   38  13  21  24    29     32            803591     451824                   0
# G13A_CHG    0       9     14   14      19   54  15  26  32    39     44            889076     492700                   0
# G14A_CHG    0      18     29   30      40  144  33  60  75    92    107           1303672     560101                   0
# G15A_CHG    0      10     15   15      20   54  17  28  34    40     46            974495     474239                   0

length(mcols(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr1$G13A_CHG)$hdiv)
# [1] 1095092

critical.val_G222324refctrl_vs_G131415_sum_chr1["G13A_CHG",][4]
# num.sites.hd95 
#          54755 

54755/1095092
# [1] 0.05000037

critical.val_G222324refctrl_vs_G131415_sum_chr1
#             tv.95% baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHG 0.2002639 0.1562063 0.6946715          64826          27829             38239
# G23A_CHG 0.2084204 0.1586989 0.6554687          62184          25241             35120
# G24A_CHG 0.2195122 0.1692335 0.6498147          47906          22157             29061
# G13A_CHG 0.3174603 0.2401692 1.4702534          54755          16204             20504
# G14A_CHG 0.2364672 0.1942815 1.5440043          68421          34795             39090
# G15A_CHG 0.2945736 0.2281038 1.3465868          57705          26326             34640

critical.val_G222324refctrl_vs_G131415_sum_chr10
#             tv.95% baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHG 0.2000000 0.1551442 0.7015498          53389          22639             31978
# G23A_CHG 0.2051282 0.1570055 0.6706583          52566          21919             30841
# G24A_CHG 0.2125000 0.1668865 0.7041412          46098          21765             28547
# G13A_CHG 0.3074074 0.2354406 1.5052272          48653          14544             18576
# G14A_CHG 0.2355863 0.1941635 1.5686334          55461          27660             31227
# G15A_CHG 0.2857143 0.2234985 1.3603745          49979          22999             30155
critical.val_G222324refctrl_vs_G131415_sum_chr13
#             tv.95% baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G22A_CHG 0.2030075 0.1605769 0.7333922          28964          11878             15726
# G23A_CHG 0.2096774 0.1636141 0.7038220          28185          11752             16119
# G24A_CHG 0.2177778 0.1728692 0.7184290          23055          11353             14598
# G13A_CHG 0.3162393 0.2457333 1.5886873          25554           7352              9064
# G14A_CHG 0.2402516 0.2010938 1.6533418          30490          15866             17747
# G15A_CHG 0.2919604 0.2324837 1.4366278          26643          12089             15585

date()
# [1] "Wed Jan  8 13:43:25 2020"

## 20200107:
#
# We call getPotentialDIMP with tv.cut values shown here:
#   - tv.col= 7L,
#   - tv.cut=0.2,
# ..normally, per a best practice, we would use a calculated value like this for chr1:
#
max(critical.val_G222324refctrl_vs_G131415_sum_chr1[,"tv.95%"])
# [1] 0.3174603
#
# ..but that will cause too many DIMPs to be filtered out.  Therefore, the best 
# ..practice may be to use a hardcoded 0.2 in cases where tv.95% values are 
# ..lower than, say 0.25. 
#
# UPDATE 20200108, we will use 0.2 anyway, even though maybe 0.317 is better
#

for(item in names(G14A_CHG_grl)) {
    assign(
        paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr", item),
        getPotentialDIMP(
            eval(str2expression(
                paste0("HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr", item)
            )),
            nlms = eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHG_sum_chr", item,"$nlms"))),
            div.col = 9L,
            tv.col= 7L,
            tv.cut=0.2,
            dist.name = eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHG_sum_chr", item,"$bestModel")))
        )
    )
    #print(eval(str2expression(paste0("bestFits_G222324refctrl_vs_G131415_CHG_sum_chr", item,"$bestModel"))))
}

## 20190108 saved:
#
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr1,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr1.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr10,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr10.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr11,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr11.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr12,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr12.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr13,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr13.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr14,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr14.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr15,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr15.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr16,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr16.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr17,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr17.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr18,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr18.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr19,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr19.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr2,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr2.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr20,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr20.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr3,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr3.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr4,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr4.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr5,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr5.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr6,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr6.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr7,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr7.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr8,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr8.RData")
# save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr9,file="/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr9.RData")

head(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr1$G14A_CHG,n=2)
# GRanges object with 2 ranges and 10 metadata columns:
#       seqnames    ranges strand |        c1        t1        c2        t2                p1                p2                TV            bay.TV
#          <Rle> <IRanges>  <Rle> | <numeric> <numeric> <numeric> <numeric>         <numeric>         <numeric>         <numeric>         <numeric>
#   [1]        1      1729      + |        12        31        14        10 0.314345246489285 0.597570390891131 0.304263565891473 0.283225144401846
#   [2]        1      1911      + |        16        43        19        21 0.298025133951155 0.493932202170905 0.203813559322034  0.19590706821975
#                   hdiv              wprob
#              <numeric>          <numeric>
#   [1] 2.63379308163417 0.0117679679825894
#   [2] 1.97794617126656 0.0277744305447254
#   -------
#   seqinfo: 20 sequences from an unspecified genome; no seqlengths

date()
# [1] "Wed Jan  8 13:45:42 2020"

names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr1)
# [1] "G22A_CHG" "G23A_CHG" "G24A_CHG" "G13A_CHG" "G14A_CHG" "G15A_CHG"

ls(pattern="^PS")
#  [1] "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr1"  "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr10"
#  [3] "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr11" "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr12"
#  [5] "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr13" "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr14"
#  [7] "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr15" "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr16"
#  [9] "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr17" "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr18"
# [11] "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr19" "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr2" 
# [13] "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr20" "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr3" 
# [15] "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr4"  "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr5" 
# [17] "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr6"  "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr7" 
# [19] "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr8"  "PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr9" 

mcols(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr1$G14A_CHG)
# DataFrame with 36691 rows and 10 columns
#              c1        t1        c2        t2                p1                p2                 TV             bay.TV             hdiv
#       <numeric> <numeric> <numeric> <numeric>         <numeric>         <numeric>          <numeric>          <numeric>        <numeric>
# 1            12        31        14        10 0.314345246489285 0.597570390891131  0.304263565891473  0.283225144401846 2.63379308163417
# 2            16        43        19        21 0.298025133951155 0.493932202170905  0.203813559322034   0.19590706821975 1.97794617126656
# 3            17        28         4        24  0.40340453463338 0.210496169526808 -0.234920634920635 -0.192908365106572 1.58608229393587
# 4            46         7        17        14 0.854899734955187 0.563780341862891 -0.319537431527693 -0.291119393092296 4.35369065639937
# 5            46        18        14        20 0.716776007311982 0.440356778848721 -0.306985294117647 -0.276419228463261 3.64527292184683
# ...         ...       ...       ...       ...               ...               ...                ...                ...              ...
# 36687        30        19         8        20 0.617905445176546  0.33545675880109 -0.326530612244898 -0.282448686375456 2.99662757484538
# 36688        76        60        21        52 0.562525663347411 0.308245336470489 -0.271152296535052 -0.254280326876922 6.43414305612511
# 36689        34        18        24         3 0.656175122391656 0.862235479294681  0.235042735042735  0.206060356903025 2.21162771733591
# 36690        24         0         9         3 0.954322104216674 0.733162656460435              -0.25 -0.221159447756239 1.81772055260489
# 36691         9        32        10        10  0.26174665846363 0.530526982041973  0.280487804878049  0.268780323578343  2.1636651764288
#                      wprob
#                  <numeric>
# 1       0.0117679679825894
# 2       0.0277744305447254
# 3       0.0475097294690465
# 4      0.00142767099435564
# 5      0.00333961802964503
# ...                    ...
# 36687  0.00743396778355233
# 36688 0.000130194484070888
# 36689    0.020353497570476
# 36690   0.0345019485195349
# 36691   0.0216838823583459

treatment_names_CHG <- names(F4_P37_CHG)

control_names_CHG <- names(F6_P37_CHG)

control_names_CHG
# [1] "G22A_CHG" "G23A_CHG" "G24A_CHG"

treatment_names_CHG
# [1] "G13A_CHG" "G14A_CHG" "G15A_CHG"

names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr1) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr2) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr3) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr4) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr5) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr6) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr7) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr8) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr9) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr10) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr11) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr12) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr13) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr14) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr15) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr16) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr17) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr18) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr19) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")
names(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr20) <- c("G22A_CHG","G23A_CHG","G24A_CHG","G13A_CHG","G14A_CHG","G15A_CHG")


















date()
# [1] "Wed Jan  8 13:46:44 2020"

#
# cutpoint with machine learning:
#
for(item in names(G14A_CHG_grl)) {
    assign(
        paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr", item),
        estimateCutPoint(
            eval(str2expression(
                paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr", item)
            )),
            simple = FALSE,
            control.names = control_names_CHG,
            treatment.names = treatment_names_CHG,
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
# [1] "Wed Jan  8 05:29:29 2020"

cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr1$cutpoint
# 1.345082


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr2$cutpoint
# 1.352122


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr3$cutpoint
# 1.30541


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr4$cutpoint
# 1.2782


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr5$cutpoint
# 1.364309


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr6$cutpoint
# 1.395088


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr7$cutpoint
# 1.329742


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr8$cutpoint
# 1.342863


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr9$cutpoint
# 1.359704


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr10$cutpoint
# 1.350709


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr11$cutpoint
# 1.304544


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr12$cutpoint
# 1.38042


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr13$cutpoint
# 1.422605


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr14$cutpoint
# 1.302386


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr15$cutpoint
# 1.361717


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr16$cutpoint
# 1.426917


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr17$cutpoint
# 1.395653


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr18$cutpoint
# 1.365431


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr19$cutpoint
# 1.364383


cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr20$cutpoint
# 1.385449



for(item in names(G14A_CHG_grl)) {
    print(paste0("chr",item))
    print(eval(str2expression(
        paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr", item,"$testSetPerformance$table")
    )))
    print(eval(str2expression(
        paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr", item,"$testSetPerformance$overall")
    )))
}
# [1] "chr1"
#           Reference
# Prediction    CT    TT
#         CT  3646     0
#         TT     0 33426
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9999005      1.0000000      0.9016508      0.0000000            NaN 
# [1] "chr10"
#           Reference
# Prediction    CT    TT
#         CT  3413     0
#         TT     0 27470
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998806      1.0000000      0.8894861      0.0000000            NaN 
# [1] "chr11"
#           Reference
# Prediction    CT    TT
#         CT  1784     0
#         TT     0 14456
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9997729      1.0000000      0.8901478      0.0000000            NaN 
# [1] "chr12"
#           Reference
# Prediction    CT    TT
#         CT  2385     0
#         TT     0 19672
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998328      1.0000000      0.8918711      0.0000000            NaN 
# [1] "chr13"
#           Reference
# Prediction    CT    TT
#         CT  1963     0
#         TT     0 15542
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9997893      1.0000000      0.8878606      0.0000000            NaN 
# [1] "chr14"
#           Reference
# Prediction    CT    TT
#         CT  3133     0
#         TT     0 28563
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998836      1.0000000      0.9011547      0.0000000            NaN 
# [1] "chr15"
#           Reference
# Prediction    CT    TT
#         CT  2883     0
#         TT     0 23281
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998590      1.0000000      0.8898104      0.0000000            NaN 
# [1] "chr16"
#           Reference
# Prediction    CT    TT
#         CT  2025     0
#         TT     0 19724
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998304      1.0000000      0.9068923      0.0000000            NaN 
# [1] "chr17"
#           Reference
# Prediction    CT    TT
#         CT  2686     0
#         TT     0 20084
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998380      1.0000000      0.8820378      0.0000000            NaN 
# [1] "chr18"
#           Reference
# Prediction    CT    TT
#         CT  3816     0
#         TT     0 30116
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998913      1.0000000      0.8875398      0.0000000            NaN 
# [1] "chr19"
#           Reference
# Prediction    CT    TT
#         CT  3278     0
#         TT     0 27962
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998819      1.0000000      0.8950704      0.0000000            NaN 
# [1] "chr2"
#           Reference
# Prediction    CT    TT
#         CT  2740     0
#         TT     0 21882
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998502      1.0000000      0.8887174      0.0000000            NaN 
# [1] "chr20"
#           Reference
# Prediction    CT    TT
#         CT  3191     0
#         TT     0 27401
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998794      1.0000000      0.8956917      0.0000000            NaN 
# [1] "chr3"
#           Reference
# Prediction    CT    TT
#         CT  2461     0
#         TT     0 21581
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998466      1.0000000      0.8976375      0.0000000            NaN 
# [1] "chr4"
#           Reference
# Prediction    CT    TT
#         CT  3120     0
#         TT     0 26786
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998767      1.0000000      0.8956731      0.0000000            NaN 
# [1] "chr5"
#           Reference
# Prediction    CT    TT
#         CT  2556     0
#         TT     0 20820
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998422      1.0000000      0.8906571      0.0000000            NaN 
# [1] "chr6"
#           Reference
# Prediction    CT    TT
#         CT  3132     0
#         TT     0 24761
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998678      1.0000000      0.8877138      0.0000000            NaN 
# [1] "chr7"
#           Reference
# Prediction    CT    TT
#         CT  2415     0
#         TT     0 21506
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998458      1.0000000      0.8990427      0.0000000            NaN 
# [1] "chr8"
#           Reference
# Prediction    CT    TT
#         CT  2312    10
#         TT     6 19279
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.9992595      0.9961369      0.9987978      0.9995767      0.8927200      0.0000000      0.4532547 
# [1] "chr9"
#           Reference
# Prediction    CT    TT
#         CT  3225     0
#         TT     0 26581
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998762      1.0000000      0.8918003      0.0000000            NaN 


# 
# cutpoint simple
#
for(item in names(G14A_CHG_grl)) {
    assign(
        paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr", item),
        estimateCutPoint(
            eval(str2expression(
                paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr", item)
            )),
            simple = TRUE,
            control.names = control_names_CHG,
            treatment.names = treatment_names_CHG,
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
cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr1$cutpoint
# 1.394336


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr2$cutpoint
# 1.488507


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr3$cutpoint
# 1.382339


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr4$cutpoint
# 1.393688


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr5$cutpoint
# 1.482074


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr6$cutpoint
# 1.544343


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr7$cutpoint
# 1.390631


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr8$cutpoint
# 1.484061


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr9$cutpoint
# 1.490388


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr10$cutpoint
# 1.486571


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr11$cutpoint
# 1.434369


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr12$cutpoint
# 1.516209


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr13$cutpoint
# 1.55992


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr14$cutpoint
# 1.428777


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr15$cutpoint
# 1.449031


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr16$cutpoint
# 1.481279


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr17$cutpoint
# 1.565575


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr18$cutpoint
# 1.52719


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr19$cutpoint
# 1.483298


cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr20$cutpoint
# 1.5243


for(item in names(G14A_CHG_grl)) {
    print(paste0("chr",item))
    print(eval(str2expression(
        paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr", item,"$testSetPerformance$table")
    )))
    print(eval(str2expression(
        paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr", item,"$testSetPerformance$overall")
    )))
}
# [1] "chr1"
#           Reference
# Prediction    CT    TT
#         CT  3270     0
#         TT     0 32800
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998977      1.0000000      0.9093429      0.0000000            NaN 
# [1] "chr10"
#           Reference
# Prediction    CT    TT
#         CT  2614     0
#         TT     0 26172
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998719      1.0000000      0.9091920      0.0000000            NaN 
# [1] "chr11"
#           Reference
# Prediction    CT    TT
#         CT  1350     0
#         TT     0 13764
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9997560      1.0000000      0.9106788      0.0000000            NaN 
# [1] "chr12"
#           Reference
# Prediction    CT    TT
#         CT  1853     0
#         TT     0 18716
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998207      1.0000000      0.9099130      0.0000000            NaN 
# [1] "chr13"
#           Reference
# Prediction    CT    TT
#         CT  1517     0
#         TT     0 14842
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9997745      1.0000000      0.9072682      0.0000000            NaN 
# [1] "chr14"
#           Reference
# Prediction    CT    TT
#         CT  2336     0
#         TT     0 27211
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998752      1.0000000      0.9209395      0.0000000            NaN 
# [1] "chr15"
#           Reference
# Prediction    CT    TT
#         CT  2421     0
#         TT     0 22517
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998521      1.0000000      0.9029192      0.0000000            NaN 
# [1] "chr16"
#           Reference
# Prediction    CT    TT
#         CT  1820     0
#         TT     0 19366
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998259      1.0000000      0.9140942      0.0000000            NaN 
# [1] "chr17"
#           Reference
# Prediction    CT    TT
#         CT  1980     0
#         TT     0 18901
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998234      1.0000000      0.9051770      0.0000000            NaN 
# [1] "chr18"
#           Reference
# Prediction    CT    TT
#         CT  2794     0
#         TT     0 28364
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998816      1.0000000      0.9103280      0.0000000            NaN 
# [1] "chr19"
#           Reference
# Prediction    CT    TT
#         CT  2574     0
#         TT     0 26760
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998743      1.0000000      0.9122520      0.0000000            NaN 
# [1] "chr2"
#           Reference
# Prediction    CT    TT
#         CT  2090     0
#         TT     0 20849
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998392      1.0000000      0.9088888      0.0000000            NaN 
# [1] "chr20"
#           Reference
# Prediction    CT    TT
#         CT  2445     0
#         TT     0 26090
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998707      1.0000000      0.9143158      0.0000000            NaN 
# [1] "chr3"
#           Reference
# Prediction    CT    TT
#         CT  2083     0
#         TT     0 20949
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998398      1.0000000      0.9095606      0.0000000            NaN 
# [1] "chr4"
#           Reference
# Prediction    CT    TT
#         CT  2434     0
#         TT     0 25586
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998684      1.0000000      0.9131335      0.0000000            NaN 
# [1] "chr5"
#           Reference
# Prediction    CT    TT
#         CT  2025     0
#         TT     0 19971
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998323      1.0000000      0.9079378      0.0000000            NaN 
# [1] "chr6"
#           Reference
# Prediction    CT    TT
#         CT  2368     0
#         TT     0 23503
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998574      1.0000000      0.9084689      0.0000000            NaN 
# [1] "chr7"
#           Reference
# Prediction    CT    TT
#         CT  2103     0
#         TT     0 21004
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998404      1.0000000      0.9089886      0.0000000            NaN 
# [1] "chr8"
#           Reference
# Prediction    CT    TT
#         CT  1777     0
#         TT     0 18393
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998171      1.0000000      0.9118989      0.0000000            NaN 
# [1] "chr9"
#           Reference
# Prediction    CT    TT
#         CT  2498     0
#         TT     0 25320
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998674      1.0000000      0.9102020      0.0000000            NaN 


for(item in names(G14A_CHG_grl)) {
    assign(
        paste0("DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr", item),
        selectDIMP(
            eval(str2expression(
                paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr", item)
            )),
            div.col = 9,
            cutpoint = eval(str2expression(
                paste0("cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr", item, '$cutpoint')
            )),
            )
        )

}

ls(pattern="DIMPsYI")
#  [1] "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr1"  "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr10" "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr11"
#  [4] "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr12" "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr13" "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr14"
#  [7] "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr15" "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr16" "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr17"
# [10] "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr18" "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr19" "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr2" 
# [13] "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr20" "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr3"  "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr4" 
# [16] "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr5"  "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr6"  "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr7" 
# [19] "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr8"  "DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr9" 

unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr1, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    3706     4242     4203    41705    36691    42726 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr2, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2017     2510     3278    29184    23917    26261 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr3, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2636     2790     2113    26020    23957    25814 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr4, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2906     3650     2081    30692    30242    28157 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr5, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    1835     2244     3591    27966    22628    25751 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr6, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2543     2611     3725    32628    26844    28576 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr7, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2298     2740     2711    26965    23728    27226 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr8, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    1831     2175     2173    21997    22207    19105 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr9, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2416     3057     3633    33308    28966    31046 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr10, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2694     2913     4218    35516    29938    32313 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr11, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    1664     1710     1477    17659    16007    16316 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr12, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    1526     1865     3794    25931    21586    23592 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr13, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    1613     1811     2083    19402    16911    17865 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr14, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2933     3223     2327    33960    32139    31520 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr15, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2201     2657     4597    29610    25359    28636 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr16, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    1882     2264     2625    24047    22156    24324 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr17, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    1830     2026     3826    26833    21933    22944 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr18, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2912     2988     4372    38485    32306    33062 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr19, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2547     3040     3931    36547    30356    33263 



unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr20, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2302     2697     4137    36329    30201    32247 

date()
# [1] "Wed Jan  8 14:10:30 2020"

DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- data.frame("G4A_CHG"=0, "G5A_CHG"=0, "G6A_CHG"=0, "G16A_CHG"=0, "G17A_CHG"=0, "G18A_CHG"=0)
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr1, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr2, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr3, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr4, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr5, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr6, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr7, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr8, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr9, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr10, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr11, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr12, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr13, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr14, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr15, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr16, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr17, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr18, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr19, length)))
DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr20, length)))

colSums(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS)
# G4A_CHG  G5A_CHG  G6A_CHG G16A_CHG G17A_CHG G18A_CHG 
#   46292    53213    64895   594784   518072   550744 

for(item in names(G14A_CHG_grl)) {
    assign(
        paste0("DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr", item),
        selectDIMP(
            eval(str2expression(
                paste0("PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr", item)
            )),
            div.col = 9,
            cutpoint = eval(str2expression(
                paste0("cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr", item, '$cutpoint')
            )),
            )
        )

}

ls(pattern="DIMPsML")
#  [1] "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr1"  "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr10" "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr11"
#  [4] "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr12" "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr13" "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr14"
#  [7] "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr15" "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr16" "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr17"
# [10] "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr18" "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr19" "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr2" 
# [13] "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr20" "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr3"  "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr4" 
# [16] "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr5"  "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr6"  "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr7" 
# [19] "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr8"  "DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr9" 

unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr1, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#     4176     4782     4786    41705    36691    45464 


unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr2, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2730     3397     4535    29184    23917    30852 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr3, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    3123     3328     2672    26020    23957    28534 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr4, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    3840     4699     2898    30692    30242    32859 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr5, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2428     2907     4693    27966    22628    29558 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr6, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    3498     3550     5226    32628    26844    34076 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr7, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2653     3177     3214    26965    23728    29455 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr8, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2383     2951     3028    21997    22207    22774 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr9, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    3288     4082     4916    33308    28966    36417 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr10, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    3597     3944     5758    35516    29938    37929 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr11, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2204     2323     2108    17659    16007    19312 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr12, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2106     2613     4967    25931    21586    27678 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr13, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2102     2401     2820    19402    16911    20768 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr14, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    3994     4371     3354    33963    32139    37236 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr15, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2679     3243     5562    29610    25359    31717 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr16, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2119     2571     2990    24047    22156    25824 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr17, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    2629     2878     5343    26833    21933    27953 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr18, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    4086     4302     6293    38485    32306    40192 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr19, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    3393     3989     5235    36547    30356    38557 
unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr20, length))
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    3143     3663     5723    36329    30201    38039 



date()
# [1] "Wed Jan  8 14:13:19 2020"

DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- data.frame("G22A_CHG"=0, "G23A_CHG"=0, "G24A_CHG"=0, "G13A_CHG"=0, "G14A_CHG"=0, "G15A_CHG"=0)
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr1, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr2, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr3, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr4, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr5, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr6, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr7, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr8, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr9, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr10, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr11, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr12, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr13, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr14, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr15, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr16, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr17, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr18, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr19, length)))
DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS <- rbind(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS,unlist(lapply(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr20, length)))

colSums(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS)
# G22A_CHG G23A_CHG G24A_CHG G13A_CHG G14A_CHG G15A_CHG 
#    60171    69171    86121   594787   518072   635194 

date()
# [1] "Wed Jan  8 14:13:44 2020"

pryr::mem_used()
# 17 GB

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

## 20190108: successfully saved:
# #
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr1.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr10.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr11.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr12.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr13.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr14.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr15.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr16.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr17.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr18.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr19.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr2.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr20.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr3.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr4.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr5.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr6.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr7.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr8.RData")
# # save(bestFits_G222324refctrl_vs_G131415_CHG_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/bestFits_G222324refctrl_vs_G131415_CHG_sum_chr9.RData")
save(control_names_CHG, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/control_names_CHG.RData")
save(treatment_names_CHG, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/treatment_names_CHG.RData")
save(covr, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/covr.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr1.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr10.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr11.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr12.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr13.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr14.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr15.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr16.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr17.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr18.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr19.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr2.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr20.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr3.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr4.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr5.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr6.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr7.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr8.RData")
save(critical.val_G222324refctrl_vs_G131415_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/critical.val_G222324refctrl_vs_G131415_sum_chr9.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr1.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr10.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr11.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr12.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr13.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr14.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr15.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr16.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr17.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr18.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr19.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr2.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr20.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr3.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr4.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr5.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr6.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr7.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr8.RData")
save(cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsML_PS_G222324refctrl_vs_G131415_CHG_sum_chr9.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr1.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr10.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr11.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr12.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr13.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr14.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr15.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr16.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr17.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr18.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr19.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr2.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr20.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr3.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr4.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr5.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr6.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr7.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr8.RData")
save(cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/cutpointsYI_PS_G222324refctrl_vs_G131415_CHG_sum_chr9.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr1.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr10.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr11.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr12.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr13.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr14.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr15.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr16.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr17.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr18.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr19.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr2.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr20.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr3.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr4.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr5.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr6.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr7.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr8.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_chr9.RData")
save(DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsML_G222324refctrl_vs_G131415_CHG_sum_COUNTS.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr1.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr10.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr11.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr12.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr13.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr14.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr15.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr16.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr17.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr18.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr19.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr2.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr20.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr3.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr4.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr5.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr6.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr7.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr8.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_chr9.RData")
save(DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/DIMPsYI_G222324refctrl_vs_G131415_CHG_sum_COUNTS.RData")
save(F4_P37_CHG, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/F4_P37_CHG.RData")
save(F6_P37_CHG, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/F6_P37_CHG.RData")

# TODO 20190108_1414pm

save(G22A_CHG_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G22A_CHG_grl.RData")
save(G23A_CHG_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G23A_CHG_grl.RData")
save(G24A_CHG_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G24A_CHG_grl.RData")
save(G13A_CHG_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G13A_CHG_grl.RData")
save(G14A_CHG_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G14A_CHG_grl.RData")
save(G15A_CHG_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/G15A_CHG_grl.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr1.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr10.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr11.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr12.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr13.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr14.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr15.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr16.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr17.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr18.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr19.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr2.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr20.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr3.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr4.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr5.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr6.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr7.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr8.RData")
# # save(HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/HD_mc12_4_mm3_hc300_p999_WT_G222324refctrl_vs_G131415_CHG_sum_ID_chr9.RData")
# save(item, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/item.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr1.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr10.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr11.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr12.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr13.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr14.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr15.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr16.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr17.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr18.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr19.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr2.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr20.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr3.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr4.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr5.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr6.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr7.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr8.RData")
# # save(PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G222324refctrl_vs_G131415_sepChroms_CHG_20200108_R/PS_G222324refctrl_vs_G131415_tvcut02_CHG_sum_chr9.RData")



date()
# [1] "Wed Jan  8 06:11:24 2020"

