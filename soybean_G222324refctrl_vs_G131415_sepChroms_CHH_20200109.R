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






