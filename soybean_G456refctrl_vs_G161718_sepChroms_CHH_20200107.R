# soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107.R


.libPaths()

library(BiocManager)
library(MethylIT)
library(MethylIT.utils)

rm(list=ls())
length(ls())
# [1] 0

# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F6_R10_CG_4_29_2019.RData")
# load(file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/methyl_deduplicate/F4_R10_CG_4_29_2019.RData")
load(file = "/data/users/twm118/soybean/F6_R10_CHH_4_29_2019.RData")
load(file = "/data/users/twm118/soybean/F4_R10_CHH_4_29_2019.RData")

names(F6_R10_CHH)
# [1] "G4A_CHH" "G5A_CHH" "G6A_CHH"

names(F4_R10_CHH)
# [1] "G16A_CHH" "G17A_CHH" "G18A_CHH"

ls()
# [1] "F4_R10_CHH" "F6_R10_CHH"

seqnames(F4_R10_CHH$G16A_CHH)
# factor-Rle of length 254667767 with 20 runs
# Lengths: 15619475 13855537  9164311 10652937 11910229 13405530 ... 11297345 13530428 11878478 12636861 13507325
# Values :        1       10       11       12       13       14 ...        5        6        7        8        9
# Levels(20): 1 10 11 12 13 14 15 16 17 18 19 2 20 3 4 5 6 7 8 9

G4A_CHH_grl <- split(F6_R10_CHH$G4A_CHH,seqnames(F6_R10_CHH$G4A_CHH))
G5A_CHH_grl <- split(F6_R10_CHH$G5A_CHH,seqnames(F6_R10_CHH$G4A_CHH))
G6A_CHH_grl <- split(F6_R10_CHH$G6A_CHH,seqnames(F6_R10_CHH$G4A_CHH))

G16A_CHH_grl <- split(F4_R10_CHH$G16A_CHH,seqnames(F4_R10_CHH$G16A_CHH))
G17A_CHH_grl <- split(F4_R10_CHH$G17A_CHH,seqnames(F4_R10_CHH$G16A_CHH))
G18A_CHH_grl <- split(F4_R10_CHH$G18A_CHH,seqnames(F4_R10_CHH$G16A_CHH))

names(G4A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9"

names(G16A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9"


date()
# [1] "Wed Jan  8 04:30:06 2020"

for(item in names(G4A_CHH_grl)) {
  assign(
    paste0("soy_WT_G456_Ref_CHH_sum_chr", item),
    poolFromGRlist(
      list(G4A_CHH_grl[[item]], G5A_CHH_grl[[item]], G6A_CHH_grl[[item]]),
      stat = "sum",
      num.cores = 12L
    ))}

date()
# [1] "Wed Jan  8 05:02:41 2020"

length(ls(pattern="soy_WT_G456_Ref_CHH_sum_chr"))
# [1] 20

date()
# [1] "Wed Jan  8 05:07:23 2020"

## 20190104: 61 minute runtime:
#
for(item in names(G4A_CHH_grl)) {
  assign(
    paste0(
      "HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr",
      item
    ),
    estimateDivergence(
      ref = eval(str2expression(paste0("soy_WT_G456_Ref_CHH_sum_chr",item))),
      indiv = list(G4A_CHH_grl[[item]], G5A_CHH_grl[[item]], G6A_CHH_grl[[item]],G16A_CHH_grl[[item]], G17A_CHH_grl[[item]], G18A_CHH_grl[[item]]),
      Bayesian = TRUE,
      min.coverage = c(12,4),
      min.meth = 3,
      high.coverage = 300,
      percentile = 0.999,
      num.cores = 24L,
      tasks = 20L,
      verbose = FALSE
    )
  )
}

## 20200108: sadly we get this error from R:
#
# Error in result[[njob]] <- value :
#   attempt to select less than one element in OneIndex
#
# ..and we see, in linux, that the OOM killer killed it:
#
# dmesg -T|tail -2
# [Wed Jan  8 06:06:26 2020] Out of memory: Kill process 14073 (rsession) score 409 or sacrifice child
# [Wed Jan  8 06:06:26 2020] Killed process 14073 (rsession), UID 1002, total-vm:55926780kB, anon-rss:53824996kB, file-rss:0kB, shmem-rss:0kB
#
# ..and so we check which ones completed successfully,
# ..and we notice that it probably died on chr18:
#
ls(pattern="^HD")
# [1] "HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr1"
# [2] "HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr10"
# [3] "HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr11"
# [4] "HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr12"
# [5] "HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr13"
# [6] "HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr14"
# [7] "HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr15"
# [8] "HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr16"
# [9] "HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr17"
#
# ..we do not have 'pryr::mem_used()' here on cornerPC46
# ..but linux 'htop' reveals that this R process is using 40GB
# ..and cornerPC46 has 128GB of RAM
# ..rstudio environment pane shows that these 6 big objects:
#
ls(pattern="_grl")
# [1] "G16A_CHH_grl" "G17A_CHH_grl" "G18A_CHH_grl" "G4A_CHH_grl"  "G5A_CHH_grl"  "G6A_CHH_grl"
#
# ..are consuming about 4.5GB each, but estimateDivergence needs them
# ..therefore there is not much RAM to save here with rm()
# ..and so we will start the loop again but with fewer "num.cores"
# ..and we will start only for the remaining chromosomes.
#
# Notice the 2 changes from the failed run above:
#   - we only process chromosomes: 18 19 2 20 3 4 5 6 7 8 9
#   - we change num.cores from 24L to 12L
#
names(G4A_CHH_grl[10:20])
# [1] "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9"

date()
# [1] "Wed Jan  8 06:23:40 2020"
# [1] "Wed Jan  8 06:30:48 2020"
# [1] "Wed Jan  8 06:49:52 2020"

# for(item in names(G4A_CHH_grl)) {
for(item in names(G4A_CHH_grl[10:20])) {
  assign(
    paste0(
      "HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr",
      item
    ),
    estimateDivergence(
      ref = eval(str2expression(paste0("soy_WT_G456_Ref_CHH_sum_chr",item))),
      indiv = list(G4A_CHH_grl[[item]], G5A_CHH_grl[[item]], G6A_CHH_grl[[item]],G16A_CHH_grl[[item]], G17A_CHH_grl[[item]], G18A_CHH_grl[[item]]),
      Bayesian = TRUE,
      min.coverage = c(12,4),
      min.meth = 3,
      high.coverage = 300,
      percentile = 0.999,
      num.cores = 2L,
      tasks = 20L,
      verbose = FALSE
    )
  )
}
date()
# [1] "Wed Jan  8 06:28:15 2020"
# [1] "Wed Jan  8 09:18:18 2020" GOOD With num.cores=2
#
# 20200108, we failed again with 12L, so we change to 6L, above, and run again.
# ..again we failed, change to num.cores=2L
#  ..this error was obtained quickly with num.cores=12L
#
# Error in result[[njob]] <- value :
#   attempt to select less than one element in OneIndex
#
# ..and fresh OOM killing as shown in linux:
# dmesg -T|tail -2
# [Wed Jan  8 06:27:47 2020] Out of memory: Kill process 15343 (rsession) score 409 or sacrifice child
# [Wed Jan  8 06:27:47 2020] Killed process 15343 (rsession), UID 1002, total-vm:55898532kB, anon-rss:53801464kB, file-rss:248kB, shmem-rss:0kB






date()

## TODO 20200108

length(ls(pattern="HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr"))
# [1] 20

names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr1)
# NULL

names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr1) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr2) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr3) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr4) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr5) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr6) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr7) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr8) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr9) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr10) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr11) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr12) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr13) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr14) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr15) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr16) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr17) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr18) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr19) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr20) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")

names(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr1)
# [1] "G4A_CHH"  "G5A_CHH"  "G6A_CHH"  "G16A_CHH" "G17A_CHH" "G18A_CHH"

for(item in names(G4A_CHH_grl)) {
  assign(paste0("critical.val_G456ref_v_G161718_sum_chr", item),
         do.call(rbind, lapply(eval(str2expression(paste0("HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr",item))), function(x) {
           hd.95 = quantile(x$hdiv, 0.95)
           baytv.95 = quantile(abs(x$bay.TV), 0.95)
           tv.95 = quantile(abs(x$TV), 0.95)
           return(c(tv = tv.95, baytv = baytv.95, hd = hd.95, num.sites.hd95 = sum(x$hdiv > hd.95), num.sites.tv95 = sum(x$TV > tv.95), num.sites.baytv95 = sum(x$bay.TV > baytv.95)))}))
  )
}

for(item in names(G4A_CHH_grl)) {
  print(eval(str2expression(paste0("critical.val_G456ref_v_G161718_sum_chr",item))))
}
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1904762 0.09163842 0.3245104          58895          21134             11388
# G5A_CHH  0.1904762 0.09212633 0.3205800          54854          23287             21481
# G6A_CHH  0.2105263 0.11685444 0.5165605          58329          46709             55481
# G16A_CHH 0.2333333 0.13013719 0.7991224          66307          21326             13375
# G17A_CHH 0.2333333 0.12981531 0.8790658          72864          38137             42220
# G18A_CHH 0.2697368 0.15963263 1.1611405          72033          54488             66141
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1825397 0.09084005 0.3591027          52468          18108              8098
# G5A_CHH  0.1785714 0.08994970 0.3529417          50180          21048             16879
# G6A_CHH  0.2000000 0.11270666 0.5489035          52018          41336             49246
# G16A_CHH 0.2272727 0.12894136 0.8508447          56860          18274             11281
# G17A_CHH 0.2285714 0.12824847 0.9151370          60165          30971             33876
# G18A_CHH 0.2666667 0.15896677 1.2117620          58836          44331             53738
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.2000000 0.09441736 0.3264068          25806           8804              4963
# G5A_CHH  0.1904762 0.09401921 0.3274650          24629          10714              9489
# G6A_CHH  0.2142857 0.11799028 0.5159587          25879          20312             24485
# G16A_CHH 0.2400000 0.13445766 0.8241023          29504           9654              6093
# G17A_CHH 0.2377014 0.13252119 0.8897738          32577          16622             17976
# G18A_CHH 0.2727273 0.16114272 1.1677970          32372          24032             29291
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1875000 0.09306963 0.3596789          38345          13444              6034
# G5A_CHH  0.1820652 0.09177561 0.3527098          36429          15890             12729
# G6A_CHH  0.2033493 0.11322826 0.5410908          38033          30006             35227
# G16A_CHH 0.2307692 0.13172224 0.8614453          42176          13234              7665
# G17A_CHH 0.2325581 0.13058698 0.9177720          44176          22393             23633
# G18A_CHH 0.2662338 0.15611351 1.1943054          44614          32306             39393
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1948052 0.09728317 0.3721606          31597          11028              4496
# G5A_CHH  0.1867322 0.09480905 0.3597007          29758          12670             10132
# G6A_CHH  0.2142857 0.11800592 0.5596819          31489          24432             29390
# G16A_CHH 0.2380952 0.13713419 0.8870122          34842          11346              6491
# G17A_CHH 0.2407407 0.13536473 0.9360303          36937          18572             19387
# G18A_CHH 0.2767361 0.16209719 1.2213054          36966          27037             32527
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1958525 0.09003267 0.3003177          45056          16159             10311
# G5A_CHH  0.1916667 0.09091757 0.3010949          42449          18577             18416
# G6A_CHH  0.2142857 0.11843451 0.4899590          43810          35379             42230
# G16A_CHH 0.2380952 0.12949975 0.7603280          50886          16520             10868
# G17A_CHH 0.2347779 0.12870773 0.8435333          57887          30017             33983
# G18A_CHH 0.2733333 0.16071091 1.1110523          55891          42470             51795
# tv.95% baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1869369 0.1024305 0.4214902          46248          15843              4734
# G5A_CHH  0.1805843 0.1001418 0.4200619          45828          19269              9993
# G6A_CHH  0.2022027 0.1173702 0.5798857          45959          35903             39174
# G16A_CHH 0.2337398 0.1405235 0.9231150          49079          15809              8016
# G17A_CHH 0.2380952 0.1377909 0.9598436          50401          25368             25295
# G18A_CHH 0.2727273 0.1644449 1.2253115          50122          36412             42797
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1888889 0.09174362 0.3357559          34109          12122              6091
# G5A_CHH  0.1826923 0.09090123 0.3438261          33879          14512             11988
# G6A_CHH  0.2083333 0.11604070 0.5177373          32476          25874             30885
# G16A_CHH 0.2307692 0.13076781 0.8285219          38299          12126              7312
# G17A_CHH 0.2321429 0.12977780 0.8809419          40273          20802             22566
# G18A_CHH 0.2747253 0.16129931 1.1223304          37093          27630             33769
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1777778 0.09188569 0.3989029          43635          14702              5667
# G5A_CHH  0.1711495 0.08966572 0.4028436          43730          18379             12865
# G6A_CHH  0.1947368 0.11014255 0.5836680          43691          34715             40603
# G16A_CHH 0.2250000 0.13054240 0.9085465          45546          14313              8043
# G17A_CHH 0.2285714 0.12919515 0.9494686          46497          23466             25160
# G18A_CHH 0.2625995 0.15622276 1.2413028          46640          34602             41761
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1875000 0.09686220 0.3813952          57488          19901              7435
# G5A_CHH  0.1818182 0.09514367 0.3805992          56342          24096             16377
# G6A_CHH  0.2051282 0.11567727 0.5521157          56852          44630             51317
# G16A_CHH 0.2323232 0.13504576 0.8797712          62414          20132             11211
# G17A_CHH 0.2361111 0.13423549 0.9300259          64858          33448             34696
# G18A_CHH 0.2727273 0.16238964 1.1827803          63015          46779             55942
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1857101 0.08965861 0.3405134          52927          18230              9212
# G5A_CHH  0.1818182 0.08927676 0.3378186          50138          21308             18452
# G6A_CHH  0.2000000 0.11256979 0.5330619          52652          42574             50152
# G16A_CHH 0.2280702 0.12776736 0.8279601          58201          18488             11219
# G17A_CHH 0.2294118 0.12723336 0.8980458          62230          32270             35728
# G18A_CHH 0.2631579 0.15579295 1.1898503          62175          46598             56754
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1875000 0.09176020 0.3496657          43984          14974              6935
# G5A_CHH  0.1818182 0.09075362 0.3471527          41899          17992             15370
# G6A_CHH  0.2035714 0.11361016 0.5396564          43791          34683             41429
# G16A_CHH 0.2307692 0.12998353 0.8505238          48936          15194              9466
# G17A_CHH 0.2307692 0.12939886 0.9079249          51412          26286             28516
# G18A_CHH 0.2666667 0.15753718 1.1965969          51293          38089             46373
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1773836 0.08865199 0.3658618          55856          19173              8476
# G5A_CHH  0.1753247 0.08773561 0.3566975          52922          22010             18193
# G6A_CHH  0.1969697 0.11017474 0.5542195          55443          44080             52717
# G16A_CHH 0.2207792 0.12589130 0.8626460          60401          19037             11468
# G17A_CHH 0.2230769 0.12538094 0.9185580          62350          31956             34780
# G18A_CHH 0.2571429 0.15335873 1.2298194          63117          47526             57459
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.2000000 0.09485560 0.3097057          33681          12366              7323
# G5A_CHH  0.1969697 0.09469539 0.3137925          32834          14261             13237
# G6A_CHH  0.2198276 0.12030538 0.4960994          33769          26750             31837
# G16A_CHH 0.2500000 0.13524590 0.7798451          38501          11934              8149
# G17A_CHH 0.2442396 0.13431571 0.8624315          43579          22390             24923
# G18A_CHH 0.2840909 0.16525358 1.1143956          41347          30465             37779
# tv.95% baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.2012987 0.1006009 0.3321660          40499          14437              6862
# G5A_CHH  0.2000000 0.1000599 0.3331217          38877          15710             12819
# G6A_CHH  0.2230769 0.1251319 0.5076638          39455          31394             36274
# G16A_CHH 0.2500000 0.1405969 0.8011721          45656          14242              8712
# G17A_CHH 0.2500000 0.1380512 0.8777175          51815          24912             28809
# G18A_CHH 0.2857143 0.1689910 1.1373727          50047          36735             44932
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1800910 0.08990919 0.3609436          43125          14660              6435
# G5A_CHH  0.1765841 0.08887716 0.3546354          41013          17582             14245
# G6A_CHH  0.2000000 0.11088734 0.5486005          42847          33643             40464
# G16A_CHH 0.2228070 0.12700695 0.8591029          46960          14976              8777
# G17A_CHH 0.2272727 0.12674724 0.9163694          48722          24875             27139
# G18A_CHH 0.2592593 0.15360216 1.2129037          49309          36746             44520
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1846154 0.09127318 0.3634446          48897          16718              7387
# G5A_CHH  0.1772358 0.08970745 0.3639027          47768          20524             16055
# G6A_CHH  0.2000000 0.11220593 0.5509863          48389          38207             45539
# G16A_CHH 0.2272727 0.12920007 0.8674127          53403          17103             10144
# G17A_CHH 0.2307692 0.12868570 0.9292160          55459          27800             30701
# G18A_CHH 0.2637363 0.15620262 1.2210198          55809          41355             50225
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1923077 0.09393117 0.3401528          39690          14191              7046
# G5A_CHH  0.1894737 0.09308234 0.3341323          37123          16290             14219
# G6A_CHH  0.2083333 0.11576401 0.5316350          40447          31728             38129
# G16A_CHH 0.2352941 0.13338832 0.8344913          44529          14383              8676
# G17A_CHH 0.2352941 0.13171611 0.8993802          48130          24489             26291
# G18A_CHH 0.2706767 0.15902302 1.1812207          48277          35478             43241
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.2017544 0.10115517 0.3456403          29501          10571              4958
# G5A_CHH  0.1971667 0.09989524 0.3436020          28242          11920              9060
# G6A_CHH  0.2222222 0.12387943 0.5212351          28899          22441             26642
# G16A_CHH 0.2500000 0.14170298 0.8428042          33379          10391              6414
# G17A_CHH 0.2500000 0.13938735 0.8936298          36376          17531             19525
# G18A_CHH 0.2857143 0.16803722 1.1451498          34886          25144             30742
# tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1904762 0.09098733 0.3307379          46469          16169              8695
# G5A_CHH  0.1861538 0.09090463 0.3337897          44753          18949             17065
# G6A_CHH  0.2058824 0.11416338 0.5273649          47434          37668             45364
# G16A_CHH 0.2311213 0.12977463 0.8190595          52449          17092             10625
# G17A_CHH 0.2307692 0.12901199 0.8905877          56942          29403             32084
# G18A_CHH 0.2666667 0.15661681 1.1818432          57382          42688             52259

## 20190208_0900am saved on cornerPc46
#
# save(critical.val_G456ref_v_G161718_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr1.RData")
# save(critical.val_G456ref_v_G161718_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr10.RData")
# save(critical.val_G456ref_v_G161718_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr11.RData")
# save(critical.val_G456ref_v_G161718_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr12.RData")
# save(critical.val_G456ref_v_G161718_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr13.RData")
# save(critical.val_G456ref_v_G161718_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr14.RData")
# save(critical.val_G456ref_v_G161718_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr15.RData")
# save(critical.val_G456ref_v_G161718_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr16.RData")
# save(critical.val_G456ref_v_G161718_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr17.RData")
# save(critical.val_G456ref_v_G161718_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr18.RData")
# save(critical.val_G456ref_v_G161718_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr19.RData")
# save(critical.val_G456ref_v_G161718_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr2.RData")
# save(critical.val_G456ref_v_G161718_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr20.RData")
# save(critical.val_G456ref_v_G161718_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr3.RData")
# save(critical.val_G456ref_v_G161718_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr4.RData")
# save(critical.val_G456ref_v_G161718_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr5.RData")
# save(critical.val_G456ref_v_G161718_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr6.RData")
# save(critical.val_G456ref_v_G161718_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr7.RData")
# save(critical.val_G456ref_v_G161718_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr8.RData")
# save(critical.val_G456ref_v_G161718_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr9.RData")
#
# save(F4_R10_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/F4_R10_CHH.RData")
# save(F6_R10_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/F6_R10_CHH.RData")
# save(G16A_CHH_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G16A_CHH_grl.RData")
# save(G17A_CHH_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G17A_CHH_grl.RData")
# save(G18A_CHH_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G18A_CHH_grl.RData")
# save(G4A_CHH_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G4A_CHH_grl.RData")
# save(G5A_CHH_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G5A_CHH_grl.RData")
# save(G6A_CHH_grl, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G6A_CHH_grl.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr1.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr10.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr11.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr12.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr13.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr14.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr15.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr16.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr17.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr18.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr19.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr2.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr20.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr3.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr4.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr5.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr6.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr7.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr8.RData")
# save(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr9.RData")

date()
# [1] "Wed Jan  8 09:21:31 2020"

## 20190104: almost 3 hours runtime:
#

names(G4A_CHH_grl)
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9"

for(item in names(G4A_CHH_grl)) {
  assign(
    paste0("bestFits_G456ref_v_G161718_CHH_sum_chr", item),
    gofReport(
      eval(str2expression(
        paste0("HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr", item)
      )),
      model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P", "GGamma3P"),
      column = 9,
      absolute = FALSE,
      output = c("best.model"),
      num.cores = 64L,
      verbose = FALSE
    )
  )
  print(eval(str2expression(paste0("bestFits_G456ref_v_G161718_CHH_sum_chr", item,"$bestModel"))))
}

## 20200108: these were run on another machine and copied to here, so we load():

load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr10.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr11.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr12.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr13.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr14.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr15.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr16.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr17.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr18.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr19.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr1.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr20.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr2.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr3.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr4.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr5.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr6.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr7.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr8.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr9.RData")

load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr10.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr11.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr12.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr13.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr14.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr15.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr16.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr17.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr18.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr19.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr1.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr20.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr2.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr3.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr4.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr5.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr6.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr7.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr8.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr9.RData")

load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr10.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr11.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr12.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr13.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr14.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr15.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr16.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr17.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr18.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr19.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr1.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr20.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr2.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr3.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr4.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr5.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr6.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr7.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr8.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/critical.val_G456ref_v_G161718_sum_chr9.RData")

load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G16A_CHH_grl.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G17A_CHH_grl.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G18A_CHH_grl.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G4A_CHH_grl.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G5A_CHH_grl.RData")
load("/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G6A_CHH_grl.RData")
















date()
# [1] "Wed Jan  8 14:26:48 2020"


covr <-
  lapply(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr1, function(x) {
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

## 20190108: we see:
#   - num.siteGreater_500 are all 0, confirming that our "high.coverage=300" filter in estimateDivergence worked
#   - 1st quartile values of about 10, showing that our data has good coverage, >10, for 75% of sites
#   - our data is extremely similar for all 6 samples in terms of coverage
#

## 20200104: these 6 samples look ok, not great:
#    - the tv.95% values are low, we prefer >0.25 but we see lower values, ~0.13
#    - the 3 control/reference samples G4/G5/G6 are different than treatment G16/G17/G18
#    - the hd.95% values for G4/G5 look more like treatment samples
#    - the hd.95% value for G16 is exceptionally high
#    -
#
# See also soybean_G222324refctrl_vs_G1415_sepChroms_20200102.R for more discussion.
#
#          Min. 1st Qu. Median Mean 3rd Qu. Max. 60% 95% 99% 99.9% 99.99% num.siteGreater_8 q60_to_500 num.siteGreater_500
# G4A_CHH     0      10     15   15      20   43  17  27  31    35     38           1010357     507869                   0
# G5A_CHH     0      10     14   14      18   41  16  24  28    32     36            945145     452744                   0
# G6A_CHH     0      10     15   15      20   44  16  26  30    35     38           1000789     525153                   0
# G16A_CHH    0      12     18   19      24   61  20  34  40    47     53           1203190     584059                   0
# G17A_CHH    0      12     19   21      27   81  22  40  48    58     66           1327122     613538                   0
# G18A_CHH    0      10     16   17      23   76  19  35  43    53     62           1230376     584295                   0

length(mcols(HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr1$G4A_CHH)$hdiv)
# [1] 1178246

critical.val_G456ref_v_G161718_sum_chr1["G4A_CHH",][4]
# num.sites.hd95
#          58895

58895/1178246
# [1] 0.04998532

critical.val_G456ref_v_G161718_sum_chr1
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1904762 0.09163842 0.3245104          58895          21134             11388
# G5A_CHH  0.1904762 0.09212633 0.3205800          54854          23287             21481
# G6A_CHH  0.2105263 0.11685444 0.5165605          58329          46709             55481
# G16A_CHH 0.2333333 0.13013719 0.7991224          66307          21326             13375
# G17A_CHH 0.2333333 0.12981531 0.8790658          72864          38137             42220
# G18A_CHH 0.2697368 0.15963263 1.1611405          72033          54488             66141

critical.val_G456ref_v_G161718_sum_chr10
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1825397 0.09084005 0.3591027          52468          18108              8098
# G5A_CHH  0.1785714 0.08994970 0.3529417          50180          21048             16879
# G6A_CHH  0.2000000 0.11270666 0.5489035          52018          41336             49246
# G16A_CHH 0.2272727 0.12894136 0.8508447          56860          18274             11281
# G17A_CHH 0.2285714 0.12824847 0.9151370          60165          30971             33876
# G18A_CHH 0.2666667 0.15896677 1.2117620          58836          44331             53738

critical.val_G456ref_v_G161718_sum_chr13
#             tv.95%  baytv.95%    hd.95% num.sites.hd95 num.sites.tv95 num.sites.baytv95
# G4A_CHH  0.1948052 0.09728317 0.3721606          31597          11028              4496
# G5A_CHH  0.1867322 0.09480905 0.3597007          29758          12670             10132
# G6A_CHH  0.2142857 0.11800592 0.5596819          31489          24432             29390
# G16A_CHH 0.2380952 0.13713419 0.8870122          34842          11346              6491
# G17A_CHH 0.2407407 0.13536473 0.9360303          36937          18572             19387
# G18A_CHH 0.2767361 0.16209719 1.2213054          36966          27037             32527

date()
# [1] "Wed Jan  8 14:33:56 2020"

## 20200104:
#
# We call getPotentialDIMP with tv.cut values shown here:
#   - tv.col= 7L,
#   - tv.cut=0.2,
# ..normally, per a best practice, we would use a calculated value like this for chr1:
#
max(critical.val_G456ref_v_G161718_sum_chr1[,"tv.95%"])
# [1] 0.2697368
#
# ..but that will cause too many DIMPs to be filtered out.  Therefore, the best
# ..practice may be to use a hardcoded 0.2 in cases where tv.95% values are
# ..lower than, say 0.25.
# UPDATE: we use 0.2 anyway!
#
date()
# [1] "Wed Jan  8 14:34:17 2020"

for(item in names(G4A_CHH_grl)) {
  assign(
    paste0("PS_G456ref_v_G161718_tvcut02_CHH_sum_chr", item),
    getPotentialDIMP(
      eval(str2expression(
        paste0("HD_mc12_4_mm3_hc300_p999_WT_G456ref_v_G161718_CHH_sum_ID_chr", item)
      )),
      nlms = eval(str2expression(paste0("bestFits_G456ref_v_G161718_CHH_sum_chr", item,"$nlms"))),
      div.col = 9L,
      tv.col= 7L,
      tv.cut=0.2,
      dist.name = eval(str2expression(paste0("bestFits_G456ref_v_G161718_CHH_sum_chr", item,"$bestModel")))
    )
  )
  #print(eval(str2expression(paste0("bestFits_G456ref_v_G161718_CHH_sum_chr", item,"$bestModel"))))
}

head(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr1$G4A_CHH,n=2)
# GRanges object with 2 ranges and 10 metadata columns:
#       seqnames    ranges strand |        c1        t1        c2        t2                p1                 p2
#          <Rle> <IRanges>  <Rle> | <numeric> <numeric> <numeric> <numeric>         <numeric>          <numeric>
#   [1]        1      7167      + |        29        24         5        11 0.455038056185482  0.250187898758817
#   [2]        1     10045      + |         8        31         0        10  0.19163233380284 0.0965654850840401
#                       TV              bay.TV              hdiv                wprob
#                <numeric>           <numeric>         <numeric>            <numeric>
#   [1] -0.234669811320755  -0.204850157426665  1.20801296343337 0.000510417614253018
#   [2] -0.205128205128205 -0.0950668487188002 0.323926304633803    0.049737730827657
#   -------
#   seqinfo: 20 sequences from an unspecified genome; no seqlengths

date()
# [1] "Wed Jan  8 14:40:04 2020"

names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr1)
# [1] "G4A_CHH"  "G5A_CHH"  "G6A_CHH"  "G16A_CHH" "G17A_CHH" "G18A_CHH"

ls(pattern="^PS")

mcols(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr1$G4A_CHH)
# DataFrame with 10889 rows and 10 columns
#              c1        t1        c2        t2                p1                 p2                 TV              bay.TV              hdiv                wprob
#       <numeric> <numeric> <numeric> <numeric>         <numeric>          <numeric>          <numeric>           <numeric>         <numeric>            <numeric>
# 1            29        24         5        11 0.455038056185482  0.250187898758817 -0.234669811320755  -0.204850157426665  1.20801296343337 0.000510417614253018
# 2             8        31         0        10  0.19163233380284 0.0965654850840401 -0.205128205128205 -0.0950668487188002 0.323926304633803    0.049737730827657
# 3             9        22         0         8 0.244915383652962  0.105861445589706 -0.290322580645161  -0.139053938063256 0.486250398404176   0.0188183362442925
# 4            25        36         4        16 0.356503920932964  0.189144327044461  -0.20983606557377  -0.167359593888502  1.12823415568325 0.000734989702940755
# 5            14        28         9         6 0.284306291263774  0.403205585294809  0.266666666666667   0.118899294031034    0.367462852428    0.037917838449712
# ...         ...       ...       ...       ...               ...                ...                ...                 ...               ...                  ...
# 10885        24        17         4         9 0.463040477204245  0.240510746091598 -0.277673545966229  -0.222529731112647  1.16257262244745 0.000627768742315357
# 10886        22        21         4        11  0.41389541315836  0.223192736238561 -0.244961240310078  -0.190702676919799  1.00013702256138  0.00133784697427037
# 10887         8        16         0         7 0.262537817787698  0.111214527242007 -0.333333333333333  -0.151323290545691  0.47344791134106   0.0202463582493898
# 10888         8        32         6         7 0.188242990891311  0.318102921109904  0.261538461538462   0.129859930218592 0.471541212656636   0.0204690120899908
# 10889         5        29         4         7   0.1512316291189   0.26074230896336  0.216577540106952    0.10951067984446 0.332462993855164   0.0471258276686812

# treatment_names_CHH <- names(F4_R10_CHH)

# control_names_CHH <- names(F6_R10_CHH)

treatment_names_CHH <- c("G16A_CHH","G17A_CHH","G18A_CHH")
control_names_CHH <-c( "G4A_CHH","G5A_CHH","G6A_CHH")

control_names_CHH
# [1] "G4A_CHH" "G5A_CHH" "G6A_CHH"

treatment_names_CHH
# [1] "G16A_CHH" "G17A_CHH" "G18A_CHH"

names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr1) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr2) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr3) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr4) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr5) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr6) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr7) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr8) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr9) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr10) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr11) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr12) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr13) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr14) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr15) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr16) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr17) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr18) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr19) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")
names(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr20) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")

date()
[1] "Wed Jan  8 14:42:35 2020"

#
# cutpoint with machine learning:
#
for(item in names(G4A_CHH_grl)) {
  assign(
    paste0("cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr", item),
    estimateCutPoint(
      eval(str2expression(
        paste0("PS_G456ref_v_G161718_tvcut02_CHH_sum_chr", item)
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

cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr1$cutpoint
# 0.5177557
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr2$cutpoint
# 0.3802809
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr3$cutpoint
# 0.3929554
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr4$cutpoint
# 0.3345746
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr5$cutpoint
# 0.5533773
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr6$cutpoint
# 0.3645171
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr7$cutpoint
# 0.5346509
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr8$cutpoint
# 0.5224205
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr9$cutpoint
# 0.4227848
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr10$cutpoint
# 0.5522671
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr11$cutpoint
# 0.5172369
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr12$cutpoint
# 0.3534049
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr13$cutpoint
# 0.3726247
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr14$cutpoint
# 0.4869896
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr15$cutpoint
# 0.4663416
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr16$cutpoint
# 0.5189578
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr17$cutpoint
# 0.4093049
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr18$cutpoint
# 0.3967037
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr19$cutpoint
# 0.5335728
cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr20$cutpoint
# 0.55874

for(item in names(G4A_CHH_grl)) {
  print(paste0("chr",item))
  print(eval(str2expression(
    paste0("cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr", item,"$testSetPerformance$table")
  )))
  print(eval(str2expression(
    paste0("cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr", item,"$testSetPerformance$overall")
  )))
}
# [1] "chr1"
#           Reference
# Prediction    CT    TT
#         CT  3386     0
#         TT   821 18950
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.645464e-01   8.709669e-01   9.620840e-01   9.668913e-01   8.183271e-01   0.000000e+00  3.989462e-180 
# [1] "chr10"
#           Reference
# Prediction    CT    TT
#         CT  2579     0
#         TT   580 14481
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.671202e-01   8.795255e-01   9.643821e-01   9.697028e-01   8.209184e-01   0.000000e+00  1.020347e-127 
# [1] "chr11"
#           Reference
# Prediction   CT   TT
#         CT 1659    0
#         TT  351 8941
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.679481e-01   8.852940e-01   9.644766e-01   9.711678e-01   8.164551e-01   0.000000e+00   6.986895e-78 
# [1] "chr12"
#           Reference
# Prediction    CT    TT
#         CT  2666     0
#         TT   532 10802
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.620000e-01   8.854934e-01   9.586997e-01   9.651065e-01   7.715714e-01   0.000000e+00  2.821889e-117 
# [1] "chr13"
#           Reference
# Prediction   CT   TT
#         CT 2611    0
#         TT  434 9900
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.664735e-01   9.019796e-01   9.632273e-01   9.695078e-01   7.647740e-01   0.000000e+00   5.957843e-96 
# [1] "chr14"
#           Reference
# Prediction    CT    TT
#         CT  2938     0
#         TT   654 15999
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.666173e-01   8.800581e-01   9.640061e-01   9.690888e-01   8.166505e-01   0.000000e+00  8.200098e-144 
# [1] "chr15"
#           Reference
# Prediction    CT    TT
#         CT  2804     0
#         TT   562 12642
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.648926e-01   8.873930e-01   9.619244e-01   9.676902e-01   7.897301e-01   0.000000e+00  8.399085e-124 
# [1] "chr16"
#           Reference
# Prediction    CT    TT
#         CT  2033     0
#         TT   395 10522
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.694981e-01   8.932047e-01   9.663902e-01   9.723925e-01   8.125097e-01   0.000000e+00   1.837506e-87 
# [1] "chr17"
#           Reference
# Prediction    CT    TT
#         CT  2620     0
#         TT   481 10763
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.653058e-01   8.942614e-01   9.621242e-01   9.682901e-01   7.763272e-01   0.000000e+00  3.522812e-106 
# [1] "chr18"
#           Reference
# Prediction    CT    TT
#         CT  3985     0
#         TT   597 16206
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.712815e-01   9.123385e-01   9.689201e-01   9.735098e-01   7.795844e-01   0.000000e+00  2.046360e-131 
# [1] "chr19"
#           Reference
# Prediction    CT    TT
#         CT  2621     0
#         TT   631 14560
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.645744e-01   8.716423e-01   9.617547e-01   9.672411e-01   8.174265e-01   0.000000e+00  8.240770e-139 
# [1] "chr2"
#           Reference
# Prediction    CT    TT
#         CT  2995     0
#         TT   560 12677
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.655002e-01   8.930914e-01   9.625771e-01   9.682550e-01   7.809882e-01   0.000000e+00  2.287174e-123 
# [1] "chr20"
#           Reference
# Prediction    CT    TT
#         CT  2453     0
#         TT   571 13717
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.658921e-01   8.756210e-01   9.630306e-01   9.685903e-01   8.193656e-01   0.000000e+00  9.256844e-126 
# [1] "chr3"
#           Reference
# Prediction    CT    TT
#         CT  2975     0
#         TT   436 12752
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.730248e-01   9.150152e-01   9.704095e-01   9.754680e-01   7.889624e-01   0.000000e+00   2.186747e-96 
# [1] "chr4"
#           Reference
# Prediction    CT    TT
#         CT  4571     0
#         TT   418 15642
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.797392e-01   9.431237e-01   9.777235e-01   9.816182e-01   7.581794e-01   0.000000e+00   1.809600e-92 
# [1] "chr5"
#           Reference
# Prediction    CT    TT
#         CT  1868     0
#         TT   498 11195
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.632770e-01   8.609784e-01   9.599736e-01   9.663797e-01   8.255291e-01   0.000000e+00  7.044652e-110 
# [1] "chr6"
#           Reference
# Prediction    CT    TT
#         CT  3400     0
#         TT   478 12946
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.715882e-01   9.162956e-01   9.689649e-01   9.740467e-01   7.694960e-01   0.000000e+00  1.583751e-105 
# [1] "chr7"
#           Reference
# Prediction    CT    TT
#         CT  2160     0
#         TT   516 12588
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.661950e-01   8.734873e-01   9.632049e-01   9.690056e-01   8.246855e-01   0.000000e+00  8.541107e-114 
# [1] "chr8"
#           Reference
# Prediction    CT    TT
#         CT  2119     0
#         TT   455 11213
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.669979e-01   8.833864e-01   9.638794e-01   9.699172e-01   8.133024e-01   0.000000e+00  1.602356e-100 
# [1] "chr9"
#           Reference
# Prediction    CT    TT
#         CT  3117     0
#         TT   611 14009
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   9.655522e-01   8.896061e-01   9.627626e-01   9.681878e-01   7.898179e-01   0.000000e+00  1.844571e-134 

#
# cutpoint simple
#
for(item in names(G4A_CHH_grl)) {
  assign(
    paste0("cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr", item),
    estimateCutPoint(
      eval(str2expression(
        paste0("PS_G456ref_v_G161718_tvcut02_CHH_sum_chr", item)
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
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr1$cutpoint
# 0.8632062
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr2$cutpoint
# 0.8913822
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr3$cutpoint
# 0.8318985
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr4$cutpoint
# 0.8634961
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr5$cutpoint
# 0.8835118
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr6$cutpoint
# 0.9109304
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr7$cutpoint
# 0.8677201
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr8$cutpoint
# 0.860343
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr9$cutpoint
# 0.873138
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr10$cutpoint
# 0.9001492
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr11$cutpoint
# 0.8614706
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr12$cutpoint
# 0.8823884
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr13$cutpoint
# 0.9003125
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr14$cutpoint
# 0.8289525
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr15$cutpoint
# 0.9400275
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr16$cutpoint
# 0.8529181
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr17$cutpoint
# 0.9125234
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr18$cutpoint
# 0.9128672
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr19$cutpoint
# 0.8855917
cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr20$cutpoint
# 0.902225

for(item in names(G4A_CHH_grl)) {
  print(paste0("chr",item))
  print(eval(str2expression(
    paste0("cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr", item,"$testSetPerformance$table")
  )))
  print(eval(str2expression(
    paste0("cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr", item,"$testSetPerformance$overall")
  )))
}
# [1] "chr1"
#           Reference
# Prediction    CT    TT
#         CT  1392     0
#         TT     0 18254
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9998123      1.0000000      0.9291459      0.0000000            NaN 
# [1] "chr10"
#           Reference
# Prediction    CT    TT
#         CT  1138     0
#         TT     0 14034
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9997569      1.0000000      0.9249934      0.0000000            NaN 
# [1] "chr11"
#           Reference
# Prediction   CT   TT
#         CT  696    0
#         TT    0 8609
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#   1.000000e+00   1.000000e+00   9.996036e-01   1.000000e+00   9.252015e-01  6.746118e-315            NaN 
# [1] "chr12"
#           Reference
# Prediction    CT    TT
#         CT   803     0
#         TT     6 10456
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#     0.99946738     0.99599109     0.99884107     0.99980451     0.92818464     0.00000000     0.04122683 
# [1] "chr13"
#           Reference
# Prediction   CT   TT
#         CT  818    0
#         TT    8 9597
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#     0.99923247     0.99471721     0.99848822     0.99966858     0.92075218     0.00000000     0.01332833 
# [1] "chr14"
#           Reference
# Prediction    CT    TT
#         CT  1110     0
#         TT     0 15375
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9997763      1.0000000      0.9326661      0.0000000            NaN 
# [1] "chr15"
#           Reference
# Prediction    CT    TT
#         CT   981     0
#         TT     0 12276
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9997218      1.0000000      0.9260014      0.0000000            NaN 
# [1] "chr16"
#           Reference
# Prediction    CT    TT
#         CT   851     0
#         TT     2 10177
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.9998187      0.9987280      0.9993452      0.9999780      0.9226655      0.0000000      0.4795001 
# [1] "chr17"
#           Reference
# Prediction    CT    TT
#         CT   912     0
#         TT    10 10487
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#    0.999123499    0.994070899    0.998388676    0.999579606    0.919186607    0.000000000    0.004426526 
# [1] "chr18"
#           Reference
# Prediction    CT    TT
#         CT  1213     0
#         TT     0 15811
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9997833      1.0000000      0.9287477      0.0000000            NaN 
# [1] "chr19"
#           Reference
# Prediction    CT    TT
#         CT  1110     0
#         TT     2 14047
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.9998681      0.9990287      0.9995235      0.9999840      0.9266442      0.0000000      0.4795001 
# [1] "chr2"
#           Reference
# Prediction    CT    TT
#         CT   941     0
#         TT     3 12175
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.9997713      0.9982853      0.9993319      0.9999528      0.9280433      0.0000000      0.2482131 
# [1] "chr20"
#           Reference
# Prediction    CT    TT
#         CT  1115     0
#         TT     0 13337
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9997448      1.0000000      0.9228480      0.0000000            NaN 
# [1] "chr3"
#           Reference
# Prediction    CT    TT
#         CT   902     0
#         TT     3 12150
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.9997702      0.9982163      0.9993286      0.9999526      0.9306779      0.0000000      0.2482131 
# [1] "chr4"
#           Reference
# Prediction    CT    TT
#         CT  1146     0
#         TT     0 14888
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      1.0000000      1.0000000      0.9997700      1.0000000      0.9285269      0.0000000            NaN 
# [1] "chr5"
#           Reference
# Prediction    CT    TT
#         CT   882     0
#         TT    10 10831
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#    0.999146976    0.993901612    0.998431820    0.999590868    0.923910262    0.000000000    0.004426526 
# [1] "chr6"
#           Reference
# Prediction    CT    TT
#         CT   964     0
#         TT     4 12596
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#      0.9997051      0.9977709      0.9992451      0.9999196      0.9286346      0.0000000      0.1336144 
# [1] "chr7"
#           Reference
# Prediction    CT    TT
#         CT   932     0
#         TT     6 12111
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#     0.99954019     0.99654380     0.99899947     0.99983124     0.92811710     0.00000000     0.04122683 
# [1] "chr8"
#           Reference
# Prediction    CT    TT
#         CT   879     0
#         TT     6 10892
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#     0.99949053     0.99632329     0.99889144     0.99981301     0.92485353     0.00000000     0.04122683 
# [1] "chr9"
#           Reference
# Prediction    CT    TT
#         CT  1064     0
#         TT     6 13531
#       Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
#     0.99958907     0.99696672     0.99910579     0.99984918     0.92671735     0.00000000     0.04122683 

for(item in names(G4A_CHH_grl)) {
  assign(
    paste0("DIMPsYI_G456ref_v_G161718_CHH_sum_chr", item),
    selectDIMP(
      eval(str2expression(
        paste0("PS_G456ref_v_G161718_tvcut02_CHH_sum_chr", item)
      )),
      div.col = 9,
      cutpoint = eval(str2expression(
        paste0("cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr", item, '$cutpoint')
      )),
    )
  )

}

ls(pattern="DIMPsYI")

unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr1, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 2282     1545     3316    23088    24700    31794 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr2, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1732     1119     2247    16352    16163    20730 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr3, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1308      945     2215    15311    16498    20215 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr4, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1657     1176     2571    18033    19695    24777 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr5, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1702     1168     2027    14927    14835    18279 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr6, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1807     1326     2199    16722    16702    21686 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr7, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1609     1055     2352    16372    16408    19912 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr8, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1332      927     2027    13913    14135    17627 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr9, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1894     1290     2539    17889    18194    23286 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr10, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 2019     1480     2561    18704    18467    24226 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr11, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1124      788     1593    11384    11297    14613 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr12, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1489      978     1905    14128    14186    17609 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr13, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1515      919     1978    12610    12442    15462 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr14, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1696     1093     2742    19287    20462    26635 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr15, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1733     1335     2211    15812    15715    20769 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr16, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1403     1079     1927    13775    13623    17236 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr17, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1786     1340     2025    14741    13967    17491 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr18, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 2211     1527     2807    20187    20565    27172 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr19, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 1978     1367     2557    18985    19012    24455 
unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr20, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 2127     1418     2501    18258    17821    23406 

date()

DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- data.frame("G4A_CHH"=0, "G5A_CHH"=0, "G6A_CHH"=0, "G16A_CHH"=0, "G17A_CHH"=0, "G18A_CHH"=0)
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr1, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr2, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr3, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr4, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr5, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr6, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr7, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr8, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr9, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr10, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr11, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr12, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr13, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr14, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr15, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr16, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr17, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr18, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr19, length)))
DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsYI_G456ref_v_G161718_CHH_sum_chr20, length)))

colSums(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS)
 # G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
 #   34404    23875    46300   330478   334887   427380 

for(item in names(G4A_CHH_grl)) {
  assign(
    paste0("DIMPsML_G456ref_v_G161718_CHH_sum_chr", item),
    selectDIMP(
      eval(str2expression(
        paste0("PS_G456ref_v_G161718_tvcut02_CHH_sum_chr", item)
      )),
      div.col = 9,
      cutpoint = eval(str2expression(
        paste0("cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr", item, '$cutpoint')
      )),
    )
  )

}

ls(pattern="DIMPsML_G456ref_v_G161718_CHH")
#  [1] "DIMPsML_G456ref_v_G161718_CHH_sum_chr1"  "DIMPsML_G456ref_v_G161718_CHH_sum_chr10" "DIMPsML_G456ref_v_G161718_CHH_sum_chr11"
#  [4] "DIMPsML_G456ref_v_G161718_CHH_sum_chr12" "DIMPsML_G456ref_v_G161718_CHH_sum_chr13" "DIMPsML_G456ref_v_G161718_CHH_sum_chr14"
#  [7] "DIMPsML_G456ref_v_G161718_CHH_sum_chr15" "DIMPsML_G456ref_v_G161718_CHH_sum_chr16" "DIMPsML_G456ref_v_G161718_CHH_sum_chr17"
# [10] "DIMPsML_G456ref_v_G161718_CHH_sum_chr18" "DIMPsML_G456ref_v_G161718_CHH_sum_chr19" "DIMPsML_G456ref_v_G161718_CHH_sum_chr2" 
# [13] "DIMPsML_G456ref_v_G161718_CHH_sum_chr20" "DIMPsML_G456ref_v_G161718_CHH_sum_chr3"  "DIMPsML_G456ref_v_G161718_CHH_sum_chr4" 
# [16] "DIMPsML_G456ref_v_G161718_CHH_sum_chr5"  "DIMPsML_G456ref_v_G161718_CHH_sum_chr6"  "DIMPsML_G456ref_v_G161718_CHH_sum_chr7" 
# [19] "DIMPsML_G456ref_v_G161718_CHH_sum_chr8"  "DIMPsML_G456ref_v_G161718_CHH_sum_chr9" 

unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr1, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 5843     5006    12319    26660    24701    31794 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr2, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 6756     5738     7681    18782    16165    20730 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr3, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 5280     4815     8012    18109    16498    20215 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr4, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 8463     7613    10281    21614    19755    24777 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr5, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 3783     3075     6594    16665    14835    18279 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr6, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 7745     6689     7780    18517    16702    21686 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr7, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 4007     3195     7609    18601    16424    19912 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr8, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 3527     2825     7037    15459    14135    17627 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr9, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 6275     5533     8936    20280    18194    23286 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr10, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 4837     4020     8823    20972    18470    24226 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr11, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 2931     2392     5670    12957    11318    14613 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr12, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 6291     5461     6399    15804    14198    17609 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr13, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 5847     4692     6047    13991    12462    15462 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr14, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 4899     4035    10616    22502    20467    26635 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr15, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 5933     5168     7197    17541    15715    20769 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr16, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 3480     3034     6668    15349    13684    17236 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr17, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 6126     5311     6007    16061    13967    17491 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr18, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 8592     7364     9547    22085    20569    27172 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr19, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 4926     3962     9272    21496    19130    24455 
unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr20, length))
# G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
# 4911     3853     8356    20145    17850    23406

date()
[1] "Wed Jan  8 14:54:17 2020"

DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- data.frame("G4A_CHH"=0, "G5A_CHH"=0, "G6A_CHH"=0, "G16A_CHH"=0, "G17A_CHH"=0, "G18A_CHH"=0)
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr1, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr2, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr3, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr4, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr5, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr6, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr7, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr8, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr9, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr10, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr11, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr12, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr13, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr14, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr15, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr16, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr17, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr18, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr19, length)))
DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS <- rbind(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS,unlist(lapply(DIMPsML_G456ref_v_G161718_CHH_sum_chr20, length)))

colSums(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS)
 # G4A_CHH  G5A_CHH  G6A_CHH G16A_CHH G17A_CHH G18A_CHH 
 #  110452    93781   160851   373590   335239   427380

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

## 20190108: successfully saved:
#
# save(bestFits_G456ref_v_G161718_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr1.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr10.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr11.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr12.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr13.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr14.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr15.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr16.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr17.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr18.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr19.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr2.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr20.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr3.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr4.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr5.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr6.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr7.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr8.RData")
# save(bestFits_G456ref_v_G161718_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/bestFits_G456ref_v_G161718_CHH_sum_chr9.RData")
save(control_names_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/control_names_CHH.RData")
save(covr, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/covr.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr1.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr10.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr11.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr12.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr13.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr14.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr15.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr16.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr17.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr18.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr19.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr2.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr20.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr3.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr4.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr5.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr6.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr7.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr8.RData")
save(cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsML_PS_G456ref_v_G161718_CHH_sum_chr9.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr1.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr10.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr11.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr12.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr13.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr14.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr15.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr16.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr17.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr18.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr19.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr2.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr20.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr3.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr4.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr5.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr6.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr7.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr8.RData")
save(cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/cutpointsYI_PS_G456ref_v_G161718_CHH_sum_chr9.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr1.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr10.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr11.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr12.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr13.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr14.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr15.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr16.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr17.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr18.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr19.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr2.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr20.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr3.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr4.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr5.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr6.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr7.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr8.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_chr9.RData")
save(DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsML_G456ref_v_G161718_CHH_sum_COUNTS.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr1.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr10.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr11.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr12.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr13.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr14.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr15.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr16.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr17.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr18.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr19.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr2.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr20.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr3.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr4.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr5.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr6.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr7.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr8.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_chr9.RData")
save(DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DIMPsYI_G456ref_v_G161718_CHH_sum_COUNTS.RData")
save(item, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/item.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr1.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr10.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr11.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr12.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr13.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr14.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr15.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr16.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr17.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr18.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr19.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr2.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr20.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr3.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr4.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr5.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr6.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepCMon Jan  6 08:25:07 EST 2020Mon Jan  6 08:25:07 EST 2020hroms_20200103_R/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr7.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr8.RData")
save(PS_G456ref_v_G161718_tvcut02_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/PS_G456ref_v_G161718_tvcut02_CHH_sum_chr9.RData")
save(soy_WT_G456_Ref_CHH_sum_chr1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr1.RData")
save(soy_WT_G456_Ref_CHH_sum_chr10, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr10.RData")
save(soy_WT_G456_Ref_CHH_sum_chr11, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr11.RData")
save(soy_WT_G456_Ref_CHH_sum_chr12, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr12.RData")
save(soy_WT_G456_Ref_CHH_sum_chr13, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr13.RData")
save(soy_WT_G456_Ref_CHH_sum_chr14, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr14.RData")
save(soy_WT_G456_Ref_CHH_sum_chr15, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr15.RData")
save(soy_WT_G456_Ref_CHH_sum_chr16, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr16.RData")
save(soy_WT_G456_Ref_CHH_sum_chr17, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr17.RData")
save(soy_WT_G456_Ref_CHH_sum_chr18, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr18.RData")
save(soy_WT_G456_Ref_CHH_sum_chr19, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr19.RData")
save(soy_WT_G456_Ref_CHH_sum_chr2, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr2.RData")
save(soy_WT_G456_Ref_CHH_sum_chr20, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr20.RData")
save(soy_WT_G456_Ref_CHH_sum_chr3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr3.RData")
save(soy_WT_G456_Ref_CHH_sum_chr4, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr4.RData")
save(soy_WT_G456_Ref_CHH_sum_chr5, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr5.RData")
save(soy_WT_G456_Ref_CHH_sum_chr6, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr6.RData")
save(soy_WT_G456_Ref_CHH_sum_chr7, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr7.RData")
save(soy_WT_G456_Ref_CHH_sum_chr8, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr8.RData")
save(soy_WT_G456_Ref_CHH_sum_chr9, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soy_WT_G456_Ref_CHH_sum_chr9.RData")
save(treatment_names_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/treatment_names_CHH.RData")


date()
# [1] "Sun Jan  5 07:28:34 2020"

soybean_gff3 = import("/data5/soybean/Glycine_max.Glycine_max_v2.1/Glycine_max.Glycine_max_v2.1.44_chrom1thru10only.gff3")

genes <-soybean_gff3[soybean_gff3$type=="gene"]

Genes2kb = GeneUpDownStream(genes, upstream = 2000, downstream = 2000)

G4A_YI_CHH_allChrom <- unlist(as(list(DIMPsYI_G456ref_v_G161718_CHH_sum_chr1$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr10$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr11$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr12$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr13$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr14$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr15$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr16$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr17$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr18$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr19$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr2$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr20$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr3$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr4$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr5$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr6$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr7$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr8$G4A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr9$G4A_CHH),"GRangesList"))
G5A_YI_CHH_allChrom <- unlist(as(list(DIMPsYI_G456ref_v_G161718_CHH_sum_chr1$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr10$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr11$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr12$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr13$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr14$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr15$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr16$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr17$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr18$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr19$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr2$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr20$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr3$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr4$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr5$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr6$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr7$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr8$G5A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr9$G5A_CHH),"GRangesList"))
G6A_YI_CHH_allChrom <- unlist(as(list(DIMPsYI_G456ref_v_G161718_CHH_sum_chr1$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr10$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr11$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr12$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr13$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr14$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr15$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr16$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr17$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr18$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr19$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr2$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr20$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr3$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr4$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr5$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr6$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr7$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr8$G6A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr9$G6A_CHH),"GRangesList"))
G16A_YI_CHH_allChrom <- unlist(as(list(DIMPsYI_G456ref_v_G161718_CHH_sum_chr1$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr10$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr11$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr12$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr13$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr14$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr15$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr16$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr17$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr18$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr19$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr2$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr20$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr3$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr4$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr5$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr6$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr7$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr8$G16A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr9$G16A_CHH),"GRangesList"))
G17A_YI_CHH_allChrom <- unlist(as(list(DIMPsYI_G456ref_v_G161718_CHH_sum_chr1$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr10$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr11$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr12$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr13$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr14$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr15$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr16$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr17$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr18$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr19$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr2$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr20$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr3$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr4$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr5$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr6$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr7$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr8$G17A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr9$G17A_CHH),"GRangesList"))
G18A_YI_CHH_allChrom <- unlist(as(list(DIMPsYI_G456ref_v_G161718_CHH_sum_chr1$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr10$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr11$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr12$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr13$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr14$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr15$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr16$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr17$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr18$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr19$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr2$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr20$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr3$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr4$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr5$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr6$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr7$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr8$G18A_CHH,DIMPsYI_G456ref_v_G161718_CHH_sum_chr9$G18A_CHH),"GRangesList"))

Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH <- getDIMPatGenes(GR = list(G4A_YI_CHH_allChrom, G5A_YI_CHH_allChrom, G6A_YI_CHH_allChrom, G16A_YI_CHH_allChrom, G17A_YI_CHH_allChrom, G18A_YI_CHH_allChrom),GENES=Genes2kb)
names(Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH) <- c("G4A_CHH","G5A_CHH","G6A_CHH","G16A_CHH","G17A_CHH","G18A_CHH")

## Local MethylIT experts prefer to see the "treatment" DIMPs when listed, bad thing but we do it here:
Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH <- Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH[c("G16A_CHH","G17A_CHH","G18A_CHH","G4A_CHH","G5A_CHH","G6A_CHH")]

condition = data.frame(condition = factor(c("TT","TT","TT","CT","CT","CT"), levels = c("CT", "TT")))
rownames(condition)
# [1] "1" "2" "3" "4" "5" "6"

names(Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH)
# [1] "G16A_CHH" "G17A_CHH" "G18A_CHH" "G4A_CHH"  "G5A_CHH"  "G6A_CHH"

rownames(condition) <- names(Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH)
rownames(condition)
# [1] "G16A_CHH" "G17A_CHH" "G18A_CHH" "G4A_CHH"  "G5A_CHH"  "G6A_CHH"

Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr = uniqueGRanges(Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH, type = "equal", verbose = TRUE, ignore.strand = FALSE, columns=2 )
head(Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr)
# GRanges object with 6 ranges and 6 metadata columns:
#       seqnames        ranges strand |     DIMPs   DIMPs.1   DIMPs.2   DIMPs.3   DIMPs.4   DIMPs.5
#          <Rle>     <IRanges>  <Rle> | <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
#   [1]        1   25355-30320      - |         1         0         0         0         0         0
#   [2]        1   56975-69527      - |         1         1         2         0         1         0
#   [3]        1   65770-71968      + |         1         0         1         0         0         0
#   [4]        1   88152-97947      - |         1         1         1         0         0         1
#   [5]        1   88289-93197      + |         1         1         1         0         0         1
#   [6]        1 114094-129845      + |         2         0         0         0         0         0

dim(mcols(Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr))
# [1] 29104     6

colnames(mcols(Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr)) <- rownames(condition)
colnames(mcols(Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr))
# [1] "G16A_CHH" "G17A_CHH" "G18A_CHH" "G4A_CHH"  "G5A_CHH"  "G6A_CHH"

DIMR_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr <- glmDataSet(GR = Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr, colData = condition)
class(DIMR_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr)
# [1] "RangedGlmDataSet"

head(DIMR_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr)
# $GR
# GRanges object with 6 ranges and 0 metadata columns:
#       seqnames        ranges strand
#          <Rle>     <IRanges>  <Rle>
#   [1]        1   25355-30320      -
#   [2]        1   56975-69527      -
#   [3]        1   65770-71968      +
#   [4]        1   88152-97947      -
#   [5]        1   88289-93197      +
#   [6]        1 114094-129845      +
#   -------
#   seqinfo: 20 sequences from an unspecified genome; no seqlengths
#
# $counts
#      G16A_CHH G17A_CHH G18A_CHH G4A_CHH G5A_CHH G6A_CHH
# [1,]       1       0       0      0      0      0
# [2,]       1       1       2      0      1      0
# [3,]       1       0       1      0      0      0
# [4,]       1       1       1      0      0      1
# [5,]       1       1       1      0      0      1
# [6,]       2       0       0      0      0      0
#
# $colData
#         x.colData.j...
# G16A_CHH             TT
# G17A_CHH             TT
# G18A_CHH             TT
# G4A_CHH              CT
# G5A_CHH              CT
# G6A_CHH              CT
#
# $sampleNames
# [1] "G16A_CHH" "G17A_CHH" "G18A_CHH" "G4A_CHH"  "G5A_CHH"  "G6A_CHH"
#
# $levels
# [1] "CT" "TT"
#
# $optionData
# NULL
#
# attr(,"class")
# [1] "RangedGlmDataSet"

dim(DIMR_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr$counts)
# [1] 29104     6

DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT = countTest2(
  DIMR_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr,
  num.cores = 24L,
  countFilter = T
)
# *** Number of GR after filtering counts 2346
# *** GLM...

class(DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT)
# [1] "GRanges"

dim(mcols(DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT))
# [1] 1208   14

min(abs(mcols(DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT)$log2FC))
# [1] 0.6061358

head(DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT)

DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1 = countTest2(
  DIMR_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr,
  num.cores = 24L,
  countFilter = T,
  Minlog2FC=1.0
)
# *** Number of GR after filtering counts 2346
# *** GLM...

dim(mcols(DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1))
# [1] 1009   14

min(abs(mcols(DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1)$log2FC))
# [1] 1.003302

hitsDMG_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1_ALL <- findOverlaps(Genes2kb,DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1, type = "equal",ignore.strand = FALSE)
hitsDMG_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1_ALL
# Hits object with 1009 hits and 0 metadata columns:
#          queryHits subjectHits
#          <integer>   <integer>
#      [1]        35           1
#      [2]       216           2
#      [3]       236           3
#      [4]       451           4
#      [5]       504           5
#      ...       ...         ...
#   [1005]     55386        1005
#   [1006]     55409        1006
#   [1007]     55519        1007
#   [1008]     55520        1008
#   [1009]     55561        1009
#   -------
#   queryLength: 55589 / subjectLength: 1009

DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1_ID = Genes2kb[queryHits(hitsDMG_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1_ALL),]
length(DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1_ID$gene_id)
# [1] 1009

head(DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1_ID$gene_id)
# [1] "GLYMA_01G003500" "GLYMA_01G021600" "GLYMA_01G023600" "GLYMA_01G045100" "GLYMA_01G050400" "GLYMA_01G068300"

save(soybean_gff3, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/soybean_gff3.RData")
save(genes, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/genes.RData")
save(Genes2kb, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/Genes2kb.RData")
save(G4A_YI_CHH_allChrom, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G4A_YI_CHH_allChrom.RData")
save(G5A_YI_CHH_allChrom, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G5A_YI_CHH_allChrom.RData")
save(G6A_YI_CHH_allChrom, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G6A_YI_CHH_allChrom.RData")
save(G16A_YI_CHH_allChrom, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G16A_YI_CHH_allChrom.RData")
save(G17A_YI_CHH_allChrom, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G17A_YI_CHH_allChrom.RData")
save(G18A_YI_CHH_allChrom, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/G18A_YI_CHH_allChrom.RData")
save(Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH.RData")
save(condition, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/condition.RData")
save(Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_isF_c2_ugr.RData")
save(DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT.RData")
save(DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1.RData")
save(hitsDMG_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1_ALL, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/hitsDMG_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1_ALL.RData")
save(DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1_ID, file = "/data5/soybean/F17FTSUSAT0456_SOYvvkM/soybean_G456refctrl_vs_G161718_sepChroms_CHH_20200107/DMGs_Genes_DIMPs_YI_soybean_G456ref_v_G161718_CHH_countTest2_cfT_ml2fc1_ID.RData")


