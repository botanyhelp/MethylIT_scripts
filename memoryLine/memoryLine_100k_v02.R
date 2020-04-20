# source("https://bioconductor.org/biocLite.R")
# biocLite()
.libPaths(c("/usr/lib64/R/library","/usr/share/R/library","~/R/x86_64-redhat-linux-gnu-library/3.6"))
library(BiocManager)
library(MethylIT)
library(MethylIT.utils)
library(devtools)
library(ggplot2) # graphic
library(reshape2) # To reshape the data frame
library(grid) # For multiple plots
library(gridExtra) # For multiple plots
pryr::mem_used()
# 498 MB

date()
# [1] "Sat Dec  8 09:18:16 2018"
# [1] "Tue Jul 16 13:58:25 2019"
# [1] "Tue Jul 30 08:25:30 2019"
# [1] "Thu Jan 30 10:59:09 2020"
# [1] "Wed Feb  5 12:33:19 2020"
# [1] "Fri Feb  7 16:47:05 2020"
# [1] "Wed Feb 19 11:32:53 2020"
# [1] "Thu Feb 20 11:57:15 2020"
# [1] "Fri Feb 28 13:18:58 2020"
# [1] "Mon Mar  2 05:28:52 2020"
# [1] "Tue Mar  3 08:51:43 2020"
# [1] "Wed Apr 15 08:53:05 2020"
# [1] "Fri Apr 17 07:22:58 2020"
# [1] "Mon Apr 20 08:40:29 2020"


# samples_files <- c(
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/M-2-1.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/M-2-2.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/M-2-3.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/M-2-4.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/M-2-5.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/W-2-1.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/W-2-2.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/W-2-3.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/W-2-4.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/W-2-5.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/M-3-1.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/M-3-2.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/M-3-3.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/M-3-4.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/M-3-5.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/W-3-1.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/W-3-2.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/W-3-3.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/W-3-4.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/W-3-5.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/M-4-1.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/M-4-2.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/M-4-3.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/M-4-4.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/M-4-5.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/W-4-1.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/W-4-2.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/W-4-3.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/W-4-4.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/W-4-5.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/M-5-1.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/M-5-2.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/M-5-3.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/M-5-4.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/M-5-5.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/W-5-1.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/W-5-2.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/W-5-3.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/W-5-4.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/W-5-5.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/M-6-1.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/M-6-2.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/M-6-3.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/M-6-4.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/M-6-5.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/W-6-1.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/W-6-2.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/W-6-3.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/W-6-4.cov.CX_report_CG.txt",
# "/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/W-6-5.cov.CX_report_CG.txt"
# )
#
# memory_M_2_1 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/M-2-1.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_2_2 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/M-2-2.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_2_3 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/M-2-3.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_2_4 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/M-2-4.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_2_5 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/M-2-5.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_2_1 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/W-2-1.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_2_2 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/W-2-2.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_2_3 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/W-2-3.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_2_4 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/W-2-4.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_2_5 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_2/W-2-5.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_3_1 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/M-3-1.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_3_2 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/M-3-2.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_3_3 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/M-3-3.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_3_4 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/M-3-4.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_3_5 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/M-3-5.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_3_1 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/W-3-1.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_3_2 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/W-3-2.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_3_3 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/W-3-3.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_3_4 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/W-3-4.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_3_5 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_3/W-3-5.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_4_1 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/M-4-1.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_4_2 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/M-4-2.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_4_3 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/M-4-3.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_4_4 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/M-4-4.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_4_5 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/M-4-5.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_4_1 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/W-4-1.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_4_2 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/W-4-2.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_4_3 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/W-4-3.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_4_4 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/W-4-4.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_4_5 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_4/W-4-5.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_5_1 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/M-5-1.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_5_2 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/M-5-2.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_5_3 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/M-5-3.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_5_4 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/M-5-4.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_5_5 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/M-5-5.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_5_1 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/W-5-1.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_5_2 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/W-5-2.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_5_3 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/W-5-3.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_5_4 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/W-5-4.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_5_5 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_5/W-5-5.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_6_1 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/M-6-1.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_6_2 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/M-6-2.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_6_3 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/M-6-3.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_6_4 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/M-6-4.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_M_6_5 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/M-6-5.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_6_1 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/W-6-1.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_6_2 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/W-6-2.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_6_3 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/W-6-3.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_6_4 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/W-6-4.cov.CX_report_CG.txt",sep='\t',header=F)
# memory_W_6_5 <- read.delim("/data2/FCHLFGMCCXY_L4_ARAkxqHAABCAAA/cdts-wh.genomics.cn/Upload_06132018/cov/Gen_6/W-6-5.cov.CX_report_CG.txt",sep='\t',header=F)
#
# # 20190719 saved:
# # save(list=c("memory_M_2_1","memory_M_2_2","memory_M_2_3","memory_M_2_4","memory_M_2_5","memory_W_2_1","memory_W_2_2","memory_W_2_3","memory_W_2_4","memory_W_2_5","memory_M_3_1","memory_M_3_2","memory_M_3_3","memory_M_3_4","memory_M_3_5","memory_W_3_1","memory_W_3_2","memory_W_3_3","memory_W_3_4","memory_W_3_5","memory_M_4_1","memory_M_4_2","memory_M_4_3","memory_M_4_4","memory_M_4_5","memory_W_4_1","memory_W_4_2","memory_W_4_3","memory_W_4_4","memory_W_4_5","memory_M_5_1","memory_M_5_2","memory_M_5_3","memory_M_5_4","memory_M_5_5","memory_W_5_1","memory_W_5_2","memory_W_5_3","memory_W_5_4","memory_W_5_5","memory_M_6_1","memory_M_6_2","memory_M_6_3","memory_M_6_4","memory_M_6_5","memory_W_6_1","memory_W_6_2","memory_W_6_3","memory_W_6_4","memory_W_6_5"), file="/data/users/twm118/memoryLine/memoryLine_all50samples_CG.RData")
#
# #NO
# #save(list=c(coral_3_CG_2,coral_4_CG_2,coral_5_CG_2,coral_6_CG_2,coral_7_CG_2,coral_8_CG_2,coral_9_CG_2,coral_10_CG_2,coral_11_CG_2,coral_12_CG_2,coral_13_CG_2,coral_14_CG_2,coral_15_CG_2,coral_16_CG_2,coral_17_CG_2,coral_18_CG_2,coral_19_CG_2,coral_20_CG_2,coral_21_CG_2,coral_22_CG_2,coral_23_CG_2,coral_24_CG_2,coral_25_CG_2,coral_26_CG_2,coral_27_CG_2,coral_28_CG_2,coral_29_CG_2,coral_30_CG_2,coral_31_CG_2,coral_32_CG_2),file="/data/users/twm118/coral/rdata20181208/coral_3_32_CG_2_sampleDFs.RData")
# memory_M_2_1_100k <- memory_M_2_1 %>% filter(V2 < 100000)
# memory_M_2_2_100k <- memory_M_2_2 %>% filter(V2 < 100000)
# memory_M_2_3_100k <- memory_M_2_3 %>% filter(V2 < 100000)
# memory_M_2_4_100k <- memory_M_2_4 %>% filter(V2 < 100000)
# memory_M_2_5_100k <- memory_M_2_5 %>% filter(V2 < 100000)
# memory_W_2_1_100k <- memory_W_2_1 %>% filter(V2 < 100000)
# memory_W_2_2_100k <- memory_W_2_2 %>% filter(V2 < 100000)
# memory_W_2_3_100k <- memory_W_2_3 %>% filter(V2 < 100000)
# memory_W_2_4_100k <- memory_W_2_4 %>% filter(V2 < 100000)
# memory_W_2_5_100k <- memory_W_2_5 %>% filter(V2 < 100000)
# memory_M_3_1_100k <- memory_M_3_1 %>% filter(V2 < 100000)
# memory_M_3_2_100k <- memory_M_3_2 %>% filter(V2 < 100000)
# memory_M_3_3_100k <- memory_M_3_3 %>% filter(V2 < 100000)
# memory_M_3_4_100k <- memory_M_3_4 %>% filter(V2 < 100000)
# memory_M_3_5_100k <- memory_M_3_5 %>% filter(V2 < 100000)
# memory_W_3_1_100k <- memory_W_3_1 %>% filter(V2 < 100000)
# memory_W_3_2_100k <- memory_W_3_2 %>% filter(V2 < 100000)
# memory_W_3_3_100k <- memory_W_3_3 %>% filter(V2 < 100000)
# memory_W_3_4_100k <- memory_W_3_4 %>% filter(V2 < 100000)
# memory_W_3_5_100k <- memory_W_3_5 %>% filter(V2 < 100000)
# memory_M_4_1_100k <- memory_M_4_1 %>% filter(V2 < 100000)
# memory_M_4_2_100k <- memory_M_4_2 %>% filter(V2 < 100000)
# memory_M_4_3_100k <- memory_M_4_3 %>% filter(V2 < 100000)
# memory_M_4_4_100k <- memory_M_4_4 %>% filter(V2 < 100000)
# memory_M_4_5_100k <- memory_M_4_5 %>% filter(V2 < 100000)
# memory_W_4_1_100k <- memory_W_4_1 %>% filter(V2 < 100000)
# memory_W_4_2_100k <- memory_W_4_2 %>% filter(V2 < 100000)
# memory_W_4_3_100k <- memory_W_4_3 %>% filter(V2 < 100000)
# memory_W_4_4_100k <- memory_W_4_4 %>% filter(V2 < 100000)
# memory_W_4_5_100k <- memory_W_4_5 %>% filter(V2 < 100000)
# memory_M_5_1_100k <- memory_M_5_1 %>% filter(V2 < 100000)
# memory_M_5_2_100k <- memory_M_5_2 %>% filter(V2 < 100000)
# memory_M_5_3_100k <- memory_M_5_3 %>% filter(V2 < 100000)
# memory_M_5_4_100k <- memory_M_5_4 %>% filter(V2 < 100000)
# memory_M_5_5_100k <- memory_M_5_5 %>% filter(V2 < 100000)
# memory_W_5_1_100k <- memory_W_5_1 %>% filter(V2 < 100000)
# memory_W_5_2_100k <- memory_W_5_2 %>% filter(V2 < 100000)
# memory_W_5_3_100k <- memory_W_5_3 %>% filter(V2 < 100000)
# memory_W_5_4_100k <- memory_W_5_4 %>% filter(V2 < 100000)
# memory_W_5_5_100k <- memory_W_5_5 %>% filter(V2 < 100000)
# memory_M_6_1_100k <- memory_M_6_1 %>% filter(V2 < 100000)
# memory_M_6_2_100k <- memory_M_6_2 %>% filter(V2 < 100000)
# memory_M_6_3_100k <- memory_M_6_3 %>% filter(V2 < 100000)
# memory_M_6_4_100k <- memory_M_6_4 %>% filter(V2 < 100000)
# memory_M_6_5_100k <- memory_M_6_5 %>% filter(V2 < 100000)
# memory_W_6_1_100k <- memory_W_6_1 %>% filter(V2 < 100000)
# memory_W_6_2_100k <- memory_W_6_2 %>% filter(V2 < 100000)
# memory_W_6_3_100k <- memory_W_6_3 %>% filter(V2 < 100000)
# memory_W_6_4_100k <- memory_W_6_4 %>% filter(V2 < 100000)
# memory_W_6_5_100k <- memory_W_6_5 %>% filter(V2 < 100000)
#
# memory_DFs <- c(memory_M_2_1_100k,memory_M_2_2_100k,memory_M_2_3_100k,memory_M_2_4_100k,memory_M_2_5_100k,memory_W_2_1_100k,memory_W_2_2_100k,memory_W_2_3_100k,memory_W_2_4_100k,memory_W_2_5_100k,memory_M_3_1_100k,memory_M_3_2_100k,memory_M_3_3_100k,memory_M_3_4_100k,memory_M_3_5_100k,memory_W_3_1_100k,memory_W_3_2_100k,memory_W_3_3_100k,memory_W_3_4_100k,memory_W_3_5_100k,memory_M_4_1_100k,memory_M_4_2_100k,memory_M_4_3_100k,memory_M_4_4_100k,memory_M_4_5_100k,memory_W_4_1_100k,memory_W_4_2_100k,memory_W_4_3_100k,memory_W_4_4_100k,memory_W_4_5_100k,memory_M_5_1_100k,memory_M_5_2_100k,memory_M_5_3_100k,memory_M_5_4_100k,memory_M_5_5_100k,memory_W_5_1_100k,memory_W_5_2_100k,memory_W_5_3_100k,memory_W_5_4_100k,memory_W_5_5_100k,memory_M_6_1_100k,memory_M_6_2_100k,memory_M_6_3_100k,memory_M_6_4_100k,memory_M_6_5_100k,memory_W_6_1_100k,memory_W_6_2_100k,memory_W_6_3_100k,memory_W_6_4_100k,memory_W_6_5_100k)
#
# # 20190719 saved, but with 17 warnings like this:
# # > warnings()
# # Warning messages:
# # 1: In get(results[[i]], pos = which(search() == packages[[i]])) :
# #   internal error -3 in R_decompress1
# # save(list=c("memory_M_2_1_100k","memory_M_2_2_100k","memory_M_2_3_100k","memory_M_2_4_100k","memory_M_2_5_100k","memory_W_2_1_100k","memory_W_2_2_100k","memory_W_2_3_100k","memory_W_2_4_100k","memory_W_2_5_100k","memory_M_3_1_100k","memory_M_3_2_100k","memory_M_3_3_100k","memory_M_3_4_100k","memory_M_3_5_100k","memory_W_3_1_100k","memory_W_3_2_100k","memory_W_3_3_100k","memory_W_3_4_100k","memory_W_3_5_100k","memory_M_4_1_100k","memory_M_4_2_100k","memory_M_4_3_100k","memory_M_4_4_100k","memory_M_4_5_100k","memory_W_4_1_100k","memory_W_4_2_100k","memory_W_4_3_100k","memory_W_4_4_100k","memory_W_4_5_100k","memory_M_5_1_100k","memory_M_5_2_100k","memory_M_5_3_100k","memory_M_5_4_100k","memory_M_5_5_100k","memory_W_5_1_100k","memory_W_5_2_100k","memory_W_5_3_100k","memory_W_5_4_100k","memory_W_5_5_100k","memory_M_6_1_100k","memory_M_6_2_100k","memory_M_6_3_100k","memory_M_6_4_100k","memory_M_6_5_100k","memory_W_6_1_100k","memory_W_6_2_100k","memory_W_6_3_100k","memory_W_6_4_100k","memory_W_6_5_100k"), file="/data/users/twm118/memoryLine/memoryLine100k_all50samples_CG.RData")
#
# #coral_3_32_CG_2_gene706_DFs <- c(coral_3_CG_2_gene706,coral_4_CG_2_gene706,coral_5_CG_2_gene706,coral_6_CG_2_gene706,coral_7_CG_2_gene706,coral_8_CG_2_gene706,coral_9_CG_2_gene706,coral_10_CG_2_gene706,coral_11_CG_2_gene706,coral_12_CG_2_gene706,coral_13_CG_2_gene706,coral_14_CG_2_gene706,coral_15_CG_2_gene706,coral_16_CG_2_gene706,coral_17_CG_2_gene706,coral_18_CG_2_gene706,coral_19_CG_2_gene706,coral_20_CG_2_gene706,coral_21_CG_2_gene706,coral_22_CG_2_gene706,coral_23_CG_2_gene706,coral_24_CG_2_gene706,coral_25_CG_2_gene706,coral_26_CG_2_gene706,coral_27_CG_2_gene706,coral_28_CG_2_gene706,coral_29_CG_2_gene706,coral_30_CG_2_gene706,coral_31_CG_2_gene706,coral_32_CG_2_gene706)
# #NO
# #save(list=coral_3_32_CG_2_gene706_DFs,file="/data/users/twm118/coral/rdata20181208/coral_3_32_CG_2_gene706_sampleDFs.RData")
#
# #write.table(coral_32_CG_2_gene706,file="/var/tmp/coral_32_CG_2_gene706.csv",sep='\t',quote=F, row.names=F, col.names=F)
#
# write.table(memory_M_2_1_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_2_1_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_2_2_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_2_2_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_2_3_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_2_3_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_2_4_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_2_4_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_2_5_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_2_5_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_2_1_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_2_1_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_2_2_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_2_2_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_2_3_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_2_3_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_2_4_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_2_4_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_2_5_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_2_5_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_3_1_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_3_1_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_3_2_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_3_2_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_3_3_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_3_3_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_3_4_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_3_4_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_3_5_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_3_5_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_3_1_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_3_1_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_3_2_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_3_2_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_3_3_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_3_3_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_3_4_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_3_4_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_3_5_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_3_5_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_4_1_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_4_1_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_4_2_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_4_2_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_4_3_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_4_3_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_4_4_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_4_4_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_4_5_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_4_5_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_4_1_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_4_1_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_4_2_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_4_2_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_4_3_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_4_3_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_4_4_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_4_4_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_4_5_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_4_5_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_5_1_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_5_1_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_5_2_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_5_2_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_5_3_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_5_3_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_5_4_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_5_4_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_5_5_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_5_5_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_5_1_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_5_1_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_5_2_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_5_2_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_5_3_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_5_3_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_5_4_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_5_4_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_5_5_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_5_5_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_6_1_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_6_1_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_6_2_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_6_2_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_6_3_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_6_3_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_6_4_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_6_4_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_M_6_5_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_M_6_5_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_6_1_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_6_1_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_6_2_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_6_2_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_6_3_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_6_3_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_6_4_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_6_4_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
# write.table(memory_W_6_5_100k,file="/data/users/twm118/memoryLine/memoryLine_100k/memory_W_6_5_100k.txt",sep='\t',quote=F, row.names=F, col.names=F)
#
#
# samples_files <- c("/data/users/twm118/memoryLine/memoryLine_100k/memory_M_2_1_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_2_2_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_2_3_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_2_4_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_2_5_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_2_1_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_2_2_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_2_3_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_2_4_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_2_5_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_3_1_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_3_2_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_3_3_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_3_4_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_3_5_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_3_1_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_3_2_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_3_3_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_3_4_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_3_5_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_4_1_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_4_2_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_4_3_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_4_4_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_4_5_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_4_1_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_4_2_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_4_3_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_4_4_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_4_5_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_5_1_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_5_2_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_5_3_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_5_4_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_5_5_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_5_1_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_5_2_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_5_3_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_5_4_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_5_5_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_6_1_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_6_2_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_6_3_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_6_4_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_M_6_5_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_6_1_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_6_2_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_6_3_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_6_4_100k.txt",
# "/data/users/twm118/memoryLine/memoryLine_100k/memory_W_6_5_100k.txt"
# )
#
#
# #experiment.path <- '/data2/coral_data_20181016/'
# #pattern = "*_CG.txt"
# #samples_id <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32")
# #samples_id <- c("TG1up03","TG1up04","TG2up05","TG2up06","TG1up07","TG1up08","CG2un09","CG2un10","CG1un11","CG1un12","TG2up13","TG2up14","CG1un15","CG1un16","CG1up17","CG1up18","CG1up19","CG1up20","TG1un21","TG1un22","TG1un23","TG1un24","CG2up25","CG2up26","CG2up27","CG2up28","TG2un29","TG2un30","TG2un31","TG2un32")
#
# samples_id <- c("M_2_1","M_2_2","M_2_3","M_2_4","M_2_5","W_2_1","W_2_2","W_2_3","W_2_4","W_2_5","M_3_1","M_3_2","M_3_3","M_3_4","M_3_5","W_3_1","W_3_2","W_3_3","W_3_4","W_3_5","M_4_1","M_4_2","M_4_3","M_4_4","M_4_5","W_4_1","W_4_2","W_4_3","W_4_4","W_4_5","M_5_1","M_5_2","M_5_3","M_5_4","M_5_5","W_5_1","W_5_2","W_5_3","W_5_4","W_5_5","M_6_1","M_6_2","M_6_3","M_6_4","M_6_5","W_6_1","W_6_2","W_6_3","W_6_4","W_6_5")
# length(samples_files)
# length(samples_id)
#
# memoryLine_100k_samples_GR_CG <-readCounts2GRangesList(filenames = samples_files,
#                                        sample.id = samples_id,
#                                        columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
#                                        verbose = TRUE)

## 20190729 1729pm saved:
# save(memoryLine_100k_samples_GR_CG, file="/data/users/twm118/memoryLine/memoryLine_100k/memoryLine_100k_samples_GR_CG.RData")
#load("/data/users/twm118/memoryLine/memoryLine_100k/memoryLine_100k_samples_GR_CG.RData")
#load("http://li1077-179.members.linode.com/memoryLine_100k/memoryLine_100k_samples_GR_CG.RData")


## These next two lines will:
#  1. attempt to create a new directory in your home directory named "memoryLine_100k"
#  2. set the current working directory to be that new directory at "~/memoryLine_100k"

dir.create(file.path("~", "memoryLine_100k"))
setwd(file.path("~", "memoryLine_100k"))

## These next 50 lines will download and save 50 text files, each file is almost 1 megabyte, almost 50MB total
#  After downloading, you will have 50 files, one each for 50 arabidopsis samples, in "~/memoryLine_100k"
#  Each file is a subset of the bismark counts data--only the first 100,000 base pairs of each chromosome are represented

library(RCurl)
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_2_1_100k.txt", destfile="memory_M_2_1_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_2_2_100k.txt", destfile="memory_M_2_2_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_2_3_100k.txt", destfile="memory_M_2_3_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_2_4_100k.txt", destfile="memory_M_2_4_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_2_5_100k.txt", destfile="memory_M_2_5_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_2_1_100k.txt", destfile="memory_W_2_1_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_2_2_100k.txt", destfile="memory_W_2_2_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_2_3_100k.txt", destfile="memory_W_2_3_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_2_4_100k.txt", destfile="memory_W_2_4_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_2_5_100k.txt", destfile="memory_W_2_5_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_3_1_100k.txt", destfile="memory_M_3_1_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_3_2_100k.txt", destfile="memory_M_3_2_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_3_3_100k.txt", destfile="memory_M_3_3_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_3_4_100k.txt", destfile="memory_M_3_4_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_3_5_100k.txt", destfile="memory_M_3_5_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_3_1_100k.txt", destfile="memory_W_3_1_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_3_2_100k.txt", destfile="memory_W_3_2_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_3_3_100k.txt", destfile="memory_W_3_3_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_3_4_100k.txt", destfile="memory_W_3_4_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_3_5_100k.txt", destfile="memory_W_3_5_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_4_1_100k.txt", destfile="memory_M_4_1_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_4_2_100k.txt", destfile="memory_M_4_2_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_4_3_100k.txt", destfile="memory_M_4_3_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_4_4_100k.txt", destfile="memory_M_4_4_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_4_5_100k.txt", destfile="memory_M_4_5_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_4_1_100k.txt", destfile="memory_W_4_1_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_4_2_100k.txt", destfile="memory_W_4_2_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_4_3_100k.txt", destfile="memory_W_4_3_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_4_4_100k.txt", destfile="memory_W_4_4_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_4_5_100k.txt", destfile="memory_W_4_5_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_5_1_100k.txt", destfile="memory_M_5_1_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_5_2_100k.txt", destfile="memory_M_5_2_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_5_3_100k.txt", destfile="memory_M_5_3_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_5_4_100k.txt", destfile="memory_M_5_4_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_5_5_100k.txt", destfile="memory_M_5_5_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_5_1_100k.txt", destfile="memory_W_5_1_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_5_2_100k.txt", destfile="memory_W_5_2_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_5_3_100k.txt", destfile="memory_W_5_3_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_5_4_100k.txt", destfile="memory_W_5_4_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_5_5_100k.txt", destfile="memory_W_5_5_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_6_1_100k.txt", destfile="memory_M_6_1_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_6_2_100k.txt", destfile="memory_M_6_2_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_6_3_100k.txt", destfile="memory_M_6_3_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_6_4_100k.txt", destfile="memory_M_6_4_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_M_6_5_100k.txt", destfile="memory_M_6_5_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_6_1_100k.txt", destfile="memory_W_6_1_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_6_2_100k.txt", destfile="memory_W_6_2_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_6_3_100k.txt", destfile="memory_W_6_3_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_6_4_100k.txt", destfile="memory_W_6_4_100k.txt", method="libcurl")
download.file("http://li1077-179.members.linode.com/memoryLine_100k/memory_W_6_5_100k.txt", destfile="memory_W_6_5_100k.txt", method="libcurl")

## We create a character vector containing the file path to all 50 files:

samples_files <- c("~/memoryLine_100k/memory_M_2_1_100k.txt",
"~/memoryLine_100k/memory_M_2_2_100k.txt",
"~/memoryLine_100k/memory_M_2_3_100k.txt",
"~/memoryLine_100k/memory_M_2_4_100k.txt",
"~/memoryLine_100k/memory_M_2_5_100k.txt",
"~/memoryLine_100k/memory_W_2_1_100k.txt",
"~/memoryLine_100k/memory_W_2_2_100k.txt",
"~/memoryLine_100k/memory_W_2_3_100k.txt",
"~/memoryLine_100k/memory_W_2_4_100k.txt",
"~/memoryLine_100k/memory_W_2_5_100k.txt",
"~/memoryLine_100k/memory_M_3_1_100k.txt",
"~/memoryLine_100k/memory_M_3_2_100k.txt",
"~/memoryLine_100k/memory_M_3_3_100k.txt",
"~/memoryLine_100k/memory_M_3_4_100k.txt",
"~/memoryLine_100k/memory_M_3_5_100k.txt",
"~/memoryLine_100k/memory_W_3_1_100k.txt",
"~/memoryLine_100k/memory_W_3_2_100k.txt",
"~/memoryLine_100k/memory_W_3_3_100k.txt",
"~/memoryLine_100k/memory_W_3_4_100k.txt",
"~/memoryLine_100k/memory_W_3_5_100k.txt",
"~/memoryLine_100k/memory_M_4_1_100k.txt",
"~/memoryLine_100k/memory_M_4_2_100k.txt",
"~/memoryLine_100k/memory_M_4_3_100k.txt",
"~/memoryLine_100k/memory_M_4_4_100k.txt",
"~/memoryLine_100k/memory_M_4_5_100k.txt",
"~/memoryLine_100k/memory_W_4_1_100k.txt",
"~/memoryLine_100k/memory_W_4_2_100k.txt",
"~/memoryLine_100k/memory_W_4_3_100k.txt",
"~/memoryLine_100k/memory_W_4_4_100k.txt",
"~/memoryLine_100k/memory_W_4_5_100k.txt",
"~/memoryLine_100k/memory_M_5_1_100k.txt",
"~/memoryLine_100k/memory_M_5_2_100k.txt",
"~/memoryLine_100k/memory_M_5_3_100k.txt",
"~/memoryLine_100k/memory_M_5_4_100k.txt",
"~/memoryLine_100k/memory_M_5_5_100k.txt",
"~/memoryLine_100k/memory_W_5_1_100k.txt",
"~/memoryLine_100k/memory_W_5_2_100k.txt",
"~/memoryLine_100k/memory_W_5_3_100k.txt",
"~/memoryLine_100k/memory_W_5_4_100k.txt",
"~/memoryLine_100k/memory_W_5_5_100k.txt",
"~/memoryLine_100k/memory_M_6_1_100k.txt",
"~/memoryLine_100k/memory_M_6_2_100k.txt",
"~/memoryLine_100k/memory_M_6_3_100k.txt",
"~/memoryLine_100k/memory_M_6_4_100k.txt",
"~/memoryLine_100k/memory_M_6_5_100k.txt",
"~/memoryLine_100k/memory_W_6_1_100k.txt",
"~/memoryLine_100k/memory_W_6_2_100k.txt",
"~/memoryLine_100k/memory_W_6_3_100k.txt",
"~/memoryLine_100k/memory_W_6_4_100k.txt",
"~/memoryLine_100k/memory_W_6_5_100k.txt"
)


## We create a character vector to give a name/ID to all 50 samples:

samples_id <- c("M_2_1","M_2_2","M_2_3","M_2_4","M_2_5","W_2_1","W_2_2","W_2_3","W_2_4","W_2_5","M_3_1","M_3_2","M_3_3","M_3_4","M_3_5","W_3_1","W_3_2","W_3_3","W_3_4","W_3_5","M_4_1","M_4_2","M_4_3","M_4_4","M_4_5","W_4_1","W_4_2","W_4_3","W_4_4","W_4_5","M_5_1","M_5_2","M_5_3","M_5_4","M_5_5","W_5_1","W_5_2","W_5_3","W_5_4","W_5_5","M_6_1","M_6_2","M_6_3","M_6_4","M_6_5","W_6_1","W_6_2","W_6_3","W_6_4","W_6_5")

## We have 50 files and 50 sample IDs:
length(samples_files)
length(samples_id)



memoryLine_100k_samples_GR_CG <-readCounts2GRangesList(filenames = samples_files,
                                       sample.id = samples_id,
                                       columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                       verbose = TRUE)

names(memoryLine_100k_samples_GR_CG)
#  [1] "M_2_1" "M_2_2" "M_2_3" "M_2_4" "M_2_5" "W_2_1" "W_2_2" "W_2_3" "W_2_4" "W_2_5" "M_3_1" "M_3_2" "M_3_3" "M_3_4" "M_3_5" "W_3_1" "W_3_2" "W_3_3" "W_3_4" "W_3_5" "M_4_1" "M_4_2" "M_4_3" "M_4_4" "M_4_5" "W_4_1" "W_4_2" "W_4_3" "W_4_4"
# [30] "W_4_5" "M_5_1" "M_5_2" "M_5_3" "M_5_4" "M_5_5" "W_5_1" "W_5_2" "W_5_3" "W_5_4" "W_5_5" "M_6_1" "M_6_2" "M_6_3" "M_6_4" "M_6_5" "W_6_1" "W_6_2" "W_6_3" "W_6_4" "W_6_5"

# ref_gen2_MM_CG <- poolFromGRlist(list(ml_CG_gen2$`M-2-1`, ml_CG_gen2$`M-2-2`, ml_CG_gen2$`M-2-3`, ml_CG_gen2$`M-2-4`, ml_CG_gen2$`M-2-5` ), stat = "sum", num.cores = 12L)
# ml_samples_W31W32W33_Ref_CG_sum <- poolFromGRlist(list(memoryLine_100k_samples_GR_CG$W_3_1, memoryLine_100k_samples_GR_CG$W_3_2, memoryLine_100k_samples_GR_CG$W_3_3), stat = "sum", num.cores = 12L)
ml_samples_W31W32_Ref_CG_mean <- poolFromGRlist(list(memoryLine_100k_samples_GR_CG$W_3_1, memoryLine_100k_samples_GR_CG$W_3_2), stat = "mean", num.cores = 2L)
## 20190730 1023am saved:
# save(ml_samples_W31W32W33_Ref_CG_sum, file="/data/users/twm118/memoryLine/memoryLine_100k/ml_samples_W31W32W33_Ref_CG_sum.RData")
# ml_samples_W34W35M31M32M33M34M35_LR <- memoryLine_100k_samples_GR_CG[c("W_3_4","W_3_5","M_3_1","M_3_2","M_3_3","M_3_4","M_3_5")]
ml_samples_W33W34W35M31M32M33M34M35_LR <- memoryLine_100k_samples_GR_CG[c("W_3_3","W_3_4","W_3_5","M_3_1","M_3_2","M_3_3","M_3_4","M_3_5")]
## 20190730 1025am saved:
# save(ml_samples_W34W35M31M32M33M34M35_LR, file="/data/users/twm118/memoryLine/memoryLine_100k/ml_samples_W34W35M31M32M33M34M35_LR.RData")

# names(ml_samples_W34W35M31M32M33M34M35_LR)
# [1] "W_3_4" "W_3_5" "M_3_1" "M_3_2" "M_3_3" "M_3_4" "M_3_5"
names(ml_samples_W33W34W35M31M32M33M34M35_LR)
# [1] "W_3_3" "W_3_4" "W_3_5" "M_3_1" "M_3_2" "M_3_3" "M_3_4" "M_3_5"

# ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR = estimateDivergence(ref = ml_samples_W31W32W33_Ref_CG_sum,
#                                                      indiv = ml_samples_W34W35M31M32M33M34M35_LR,
#                                                      Bayesian = TRUE,
#                                                      min.coverage = 4,
#                                                      high.coverage = 300,
#                                                      percentile = 0.999,
#                                                      num.cores = 64L, tasks = 20L, verbose = TRUE )

ml_HD_CG_sum_mincov4_W31W32_Ref_W33W34W35M31M32M33M34M35_LR = estimateDivergence(ref = ml_samples_W31W32_Ref_CG_mean,
                                                     indiv = ml_samples_W33W34W35M31M32M33M34M35_LR,
                                                     Bayesian = TRUE,
                                                     min.coverage = 4,
                                                     high.coverage = 300,
                                                     percentile = 0.999,
                                                     num.cores = 64L, tasks = 20L, verbose = TRUE )
## 20190730 1029am saved:
# save(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, file="/data/users/twm118/memoryLine/memoryLine_100k/ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR.RData")

critical.val <- do.call(rbind, lapply(ml_HD_CG_sum_mincov4_W31W32_Ref_W33W34W35M31M32M33M34M35_LR, function(x) {
    hd.95 = quantile(x$hdiv, 0.95)
    tv.95 = quantile(abs(x$TV), 0.95)
    btv.95 = quantile(abs(x$bay.TV), 0.95)
    return(c(tv = tv.95, hd = hd.95, BTV = btv.95,
            num.sites.hd95 = sum(x$hdiv > hd.95),
            num.sites.tv95 = sum(x$bay.TV > tv.95),
            num.sites.btv.95 = sum(x$TV > btv.95)
    ))})
)
critical.val
#          tv.95%   hd.95%   BTV.95% num.sites.hd95 num.sites.tv95 num.sites.btv.95
# W_3_3 0.3333333 1.196889 0.2709579            532            134              414
# W_3_4 0.3333333 1.218544 0.2751964            515            116              380
# W_3_5 0.3333333 1.263522 0.2620380            537            115              394
# M_3_1 0.4000000 1.353473 0.3201460            522             61              265
# M_3_2 0.3750000 1.638139 0.2981684            522             95              328
# M_3_3 0.4166667 1.606706 0.3308608            514             72              306
# M_3_4 0.3431788 1.974152 0.2846830            533            141              413
# M_3_5 0.3750000 1.515590 0.2947317            521            104              368

covr <- lapply(ml_HD_CG_sum_mincov4_W31W32_Ref_W33W34W35M31M32M33M34M35_LR, function(x) {
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
        q60, quantile(x, c(0.95, 0.99, 0.999, 0.9999)),
        num.siteGreater_8 = sum(x >= 8),
        q60_to_500 = sum((x >= q60) & (x <= 500)),
        num.siteGreater_500 = sum(x > 500)
))})
)
#       Min. 1st Qu. Median Mean 3rd Qu. Max. 60%    95%    99%   99.9%   99.99% num.siteGreater_8 q60_to_500 num.siteGreater_500
# W_3_3    0       8     13   30      43  288  17 100.00 154.70 259.370 284.8740              8147       4334                   0
# W_3_4    0       7     12   29      47  300  16  95.00 145.00 261.000 296.7459              7451       4218                   0
# W_3_5    0       8     13   31      47  300  18 103.55 159.00 290.000 298.8542              8445       4420                   0
# M_3_1    0       7     11   27      42  294  14  86.00 124.68 229.000 274.7408              7361       4333                   0
# M_3_2    0       7     12   23      33  295  16  71.00 107.61 232.561 278.7805              7688       4223                   0
# M_3_3    0       6     10   21      33  291  14  57.70  93.34 225.936 264.9468              6610       4123                   0
# M_3_4    0       8     12   24      34  298  16  73.00 122.54 264.770 297.9354              8087       4463                   0
# M_3_5    0       7     11   23      34  300  14  66.00 111.80 258.580 295.9580              7133       4345                   0

gof <- gofReport(HD = ml_HD_CG_sum_mincov4_W31W32_Ref_W33W34W35M31M32M33M34M35_LR,
                model = c("Weibull2P", "Weibull3P",
                        "Gamma2P", "Gamma3P"),
                column = 9,
                output = "all",
                confl_model = TRUE,
                num.cores = 4L,
                task = 2L,
                verbose = FALSE
)

gof$stats
#         w2p_AIC w2p_R.Cross.val w3p_AIC w3p_R.Cross.val   g2p_AIC g2p_R.Cross.val   g3p_AIC g3p_R.Cross.val bestModel
# W_3_3 -60164.48       0.9989067      NA              NA -81213.20       0.9998383       Inf       0.0000000       g2p
# W_3_4 -58118.26       0.9988738      NA              NA -76425.95       0.9997881       Inf       0.0000000       g2p
# W_3_5 -62350.26       0.9990648      NA              NA -86914.39       0.9998960       Inf       0.0000000       g2p
# M_3_1 -59883.82       0.9989681      NA              NA -84100.23       0.9998911 -84040.91       0.9998898       g2p
# M_3_2 -59804.17       0.9989874      NA              NA -82165.02       0.9998699       Inf       0.0000000       g2p
# M_3_3 -61468.69       0.9992200      NA              NA -84696.15       0.9999083       Inf       0.0000000       g2p
# M_3_4 -63788.20       0.9992086      NA              NA -85602.87       0.9998886       Inf       0.0000000       g2p
# M_3_5 -56967.74       0.9986595      NA              NA -75315.16       0.9997463       Inf       0.0000000       g2p

gof$bestModel
#     W_3_3     W_3_4     W_3_5     M_3_1     M_3_2     M_3_3     M_3_4     M_3_5
# "Gamma2P" "Gamma2P" "Gamma2P" "Gamma2P" "Gamma2P" "Gamma2P" "Gamma2P" "Gamma2P"

ps <- getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32_Ref_W33W34W35M31M32M33M34M35_LR,
                        nlms = gof$nlms,
                        div.col = 9L,
                        tv.col = 0.25,
                        dist.name = gof$bestModel)
ps$M_3_1
# GRanges object with 490 ranges and 10 metadata columns:
#         seqnames    ranges strand |        c1        t1        c2        t2                p1                 p2                 TV             bay.TV             hdiv                wprob
#            <Rle> <IRanges>  <Rle> | <numeric> <numeric> <numeric> <numeric>         <numeric>          <numeric>          <numeric>          <numeric>        <numeric>            <numeric>
#     [1]        1       511      + |        12         6         4         8 0.657498626527702  0.358218900525954 -0.333333333333333 -0.299279726001747 1.41553604183559   0.0480885184127172
#     [2]        1     24196      + |         6        14        16         0 0.318307935551991  0.955511076808764                0.7  0.637203141256773 10.3100227360395 1.66670919262768e-07
#     [3]        1     28159      + |        14         0        12         9 0.954940809833813   0.56869348945189 -0.428571428571429 -0.386247320381923  4.4116920565169 0.000576094301569358
#     [4]        1     28925      - |         0         7         6         5 0.100938111911038  0.544026606696251  0.545454545454545  0.443088494785213 2.40753943816567   0.0104970003179963
#     [5]        1     29383      + |        12         2         0         6 0.826455912966856  0.117611179697775 -0.857142857142857 -0.708844733269081 5.66823493712532 9.84316359254483e-05
#     ...      ...       ...    ... .       ...       ...       ...       ...               ...                ...                ...                ...              ...                  ...
#   [486]       Mt     42218      - |         6        39        15        34 0.147417308446873  0.313752083571642  0.172789115646259  0.166334775124769 1.91943032948207   0.0219325755354997
#   [487]       Mt     43386      + |         6        24         2        36 0.217469180519584 0.0731798004298523 -0.147368421052632 -0.144289380089731 1.53527590730327   0.0397879005939534
#   [488]       Mt     43413      + |         8        23         2        36 0.272205073009019 0.0731798004298523 -0.205432937181664 -0.199025272579167 2.64092667745193  0.00742253709498072
#   [489]       Mt     43416      + |         8        22         2        35 0.280828427192174 0.0750711561337632 -0.212612612612613 -0.205757271058411 2.67805190457703  0.00702641287216851
#   [490]       Mt     43428      + |         8        22         3        30 0.280828427192174  0.112552281928813 -0.175757575757576 -0.168276145263361 1.51274842879368   0.0412252524937753
#   -------
#   seqinfo: 7 sequences from an unspecified genome; no seqlengths

### Cutpoint estimated using a model classifier
cutpoint <- estimateCutPoint(LR = ps,
                            simple = FALSE,
                            div.col = 9,
                            control.names = c("W_3_3", "W_3_4", "W_3_5"),
                            treatment.names = c("M_3_1", "M_3_2", "M_3_3",
                                                "M_3_4", "M_3_5"),
                            column = c(hdiv = TRUE, bay.TV = TRUE,
                                        wprob = TRUE, pos = TRUE),
                            classifier1 = "pca.qda", n.pc = 4,
                            classifier2 = "pca.logistic",
                            center = TRUE,
                            scale = TRUE,
                            clas.perf = TRUE,
                            verbose = FALSE
)
cutpoint$cutpoint
# [1] 1.290194

cutpoint$testSetPerformance
# Confusion Matrix and Statistics
#
#           Reference
# Prediction  CT  TT
#         CT 244   0
#         TT   1 493
#
#                Accuracy : 0.9986
#                  95% CI : (0.9925, 1)
#     No Information Rate : 0.668
#     P-Value [Acc > NIR] : <2e-16
#
#                   Kappa : 0.9969
#
#  Mcnemar's Test P-Value : 1
#
#             Sensitivity : 1.0000
#             Specificity : 0.9959
#          Pos Pred Value : 0.9980
#          Neg Pred Value : 1.0000
#              Prevalence : 0.6680
#          Detection Rate : 0.6680
#    Detection Prevalence : 0.6694
#       Balanced Accuracy : 0.9980
#
#        'Positive' Class : TT

cutpoint$testSetModel.FDR
# [1] 0.002024291

dmps <- selectDIMP(ps, div.col = 9, cutpoint = cutpoint$cutpoint)
data.frame(dmps =unlist(lapply(dmps, length)))
#       dmps
# W_3_3  451
# W_3_4  447
# W_3_5  502
# M_3_1  490
# M_3_2  525
# M_3_3  558
# M_3_4  627
# M_3_5  551


# ftp://ftp.ensemblgenomes.org/pub/plants/current/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3.gz
# AG_gff3_ftp = import("ftp://ftp.ensemblgenomes.org/pub/plants/current/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.44.gff3.gz")

AG_gff3_ftp = import("ftp://ftp.ensemblgenomes.org/pub/plants/current/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3.gz")
# AG_gff3 = import("/data/TAIR10_gff3/Arabidopsis_thaliana.TAIR10.38.gff3.gz")
# dim(as.data.frame(AG_gff3[AG_gff3$type=="gene"]))
# [1] 27655    25

dim(as.data.frame(AG_gff3_ftp[AG_gff3_ftp$type=="gene"]))
# [1] 27655    25

genes <-AG_gff3_ftp[AG_gff3_ftp$type=="gene"]
# genes <-AG_gff3_ftp[AG_gff3_ftp$type=="gene", AG_gff3_ftp$biotype=="protein_coding"]

#AG_gff3_ftp_omitNA <- AG_gff3_ftp[!is.na(AG_gff3_ftp$biotype),]
#NO genes <- AG_gff3_ftp_omitNA[AG_gff3_ftp_omitNA$type=="gene", AG_gff3_ftp_omitNA$biotype=="protein_coding"]


nams <- names(dmps)
dmps_at_genes <- getDIMPatGenes(GR = dmps, GENES = genes, ignore.strand = TRUE)
dmps_at_genes <- uniqueGRanges(dmps_at_genes, columns = 2L,
                                ignore.strand = TRUE, type = "equal")


colnames(mcols(dmps_at_genes)) <- nams
dmps_at_genes

colData <- data.frame(condition = factor(c("WT", "WT", "WT",
                                            "ML", "ML", "ML", "ML", "ML"),
                                        levels = c("WT", "ML")),
                    nams,
                    row.names = 2)
## A RangedGlmDataSet is created
ds <- glmDataSet(GR = dmps_at_genes, colData = colData)

dmgs <- countTest2(DS = ds, num.cores = 4L,
                    tasks = 2L,
                    minCountPerIndv = 7,
                    maxGrpCV = c(1.2, 1.2),
                    Minlog2FC = 0.5,
                    CountPerBp = 0.001,
                    test = "LRT",
                    verbose = TRUE)

dmgs
# GRanges object with 2 ranges and 16 metadata columns:
#       seqnames      ranges strand |     W_3_3     W_3_4     W_3_5     M_3_1     M_3_2     M_3_3     M_3_4     M_3_5            log2FC  scaled.deviance              pvalue          model            adj.pval     CT.SignalDensity    TT.SignalDensity SignalDensityVariation
#          <Rle>   <IRanges>  <Rle> | <numeric> <numeric> <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>         <numeric>        <numeric>           <numeric>       <factor>           <numeric>            <numeric>           <numeric>              <numeric>
#   [1]        1 23121-31227      * |         4         2         7         9         9        11         9         8 0.752866415257533 7.62903400703822 0.00574360569975566 Neg.Binomial.W 0.00574360569975566 0.000534517495168784 0.00113482175897373   0.000600304263804942
#   [2]       Mt 23663-24235      * |         1         0         5         2        18        10         9         9  1.26224171244991 7.78372913797753 0.00527188625687915   Neg.Binomial 0.00574360569975566  0.00349040139616056  0.0167539267015707     0.0132635253054101
#   -------
#   seqinfo: 7 sequences from an unspecified genome; no seqlengths


# hits = findOverlaps(Genes2kb,Genes_DIMPs_autism_CG_sum_GGamma3P, type = "equal",ignore.strand = FALSE)
hits = findOverlaps(genes, dmgs, type = "equal",ignore.strand = FALSE)

genes[queryHits(hits)]$ID
# [1] "gene:AT1G01040" "gene:ATMG00070"











Genes2kb = GeneUpDownStream(genes, upstream = 2000, downstream = 2000)
# for(i in names(ml_samples_W34W35M31M32M33M34M35_LR)){assign(paste0("DIMPs_ml_",i,"_gene"),getDIMPatGenes(GR = DIMPs_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR[[i]],GENES=Genes2kb))}
# head(DIMPs_ml_M_3_1_gene)

genesOrig <- genes

genes <- Genes2kb

dmps_at_genes <- getDIMPatGenes(GR = dmps, GENES = genes, ignore.strand = TRUE)
dmps_at_genes <- uniqueGRanges(dmps_at_genes, columns = 2L,
                                ignore.strand = TRUE, type = "equal")


colnames(mcols(dmps_at_genes)) <- nams
dmps_at_genes

colData <- data.frame(condition = factor(c("WT", "WT", "WT",
                                            "ML", "ML", "ML", "ML", "ML"),
                                        levels = c("WT", "ML")),
                    nams,
                    row.names = 2)
## A RangedGlmDataSet is created
ds <- glmDataSet(GR = dmps_at_genes, colData = colData)

dmgs <- countTest2(DS = ds, num.cores = 4L,
                    tasks = 2L,
                    minCountPerIndv = 7,
                    maxGrpCV = c(1.2, 1.2),
                    Minlog2FC = 0.5,
                    CountPerBp = 0.001,
                    test = "LRT",
                    verbose = TRUE)

dmgs

# GRanges object with 2 ranges and 16 metadata columns:
#       seqnames      ranges strand |     W_3_3     W_3_4     W_3_5     M_3_1     M_3_2     M_3_3     M_3_4     M_3_5            log2FC  scaled.deviance             pvalue        model           adj.pval     CT.SignalDensity    TT.SignalDensity SignalDensityVariation
#          <Rle>   <IRanges>  <Rle> | <numeric> <numeric> <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>         <numeric>        <numeric>          <numeric>  <character>          <numeric>            <numeric>           <numeric>              <numeric>
#   [1]        3 29145-36670      * |         4         5         4         8         7        10         7         8 0.613104472886409 4.04994770455367 0.0441727132899051      Poisson 0.0441727132899051 0.000575781734431748  0.0010629816635663   0.000487199929134556
#   [2]       Mt  9918-14241      * |        21        21        27        19        41        38        68        41  0.58778666490212 5.58300673736071 0.0181355609011105 Neg.Binomial 0.0362711218022211  0.00531914893617021 0.00957446808510638    0.00425531914893617

hits = findOverlaps(Genes2kb, dmgs, type = "equal",ignore.strand = FALSE)

Genes2kb[queryHits(hits)]$ID

# [1] "gene:AT3G01090" "gene:ATMG00030"








library(org.At.tair.db)
# [1] "gene:AT1G01050" "gene:AT3G01090" "gene:AT4G00080" "gene:AT4G00190" "gene:ATMG00030"

Genes2kb = GeneUpDownStream(genes, upstream = 2000, downstream = 2000)
# for(i in names(ml_samples_W34W35M31M32M33M34M35_LR)){assign(paste0("DIMPs_ml_",i,"_gene"),getDIMPatGenes(GR = DIMPs_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR[[i]],GENES=Genes2kb))}
# head(DIMPs_ml_M_3_1_gene)

genesOrig <- genes

genes <- Genes2kb
hits = findOverlaps(Genes2kb, dmgs, type = "equal",ignore.strand = FALSE)
Genes2kb[queryHits(hits)]$ID



# library(TxDb)
library(BiocManager)
#NO BiocManager::install("BSgenome.Athaliana.TAIR.TAIR10")
#NO BiocManager("BSgenome.Athaliana.TAIR.TAIR10")
#NO library(BSgenome.Athaliana.TAIR.TAIR10)

BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9")
BiocManager("BSgenome.Athaliana.TAIR.TAIR9")
library(BSgenome.Athaliana.TAIR.TAIR9)

devtools::install_github("lileiting/BSgenome.Athaliana10.TAIR.TAIR10")
library(BSgenome.Athaliana10.TAIR.TAIR10)


BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")
library("TxDb.Athaliana.BioMart.plantsmart28")

library(OrganismDbi)
gd<-list(join1=c(GO.db="GOID",org.At.tair.db="GO"),join2=c(org.At.tair.db="ENTREZID",TxDb.Athaliana.BioMart.plantsmart28="GENEID"))
makeOrganismPackage("AraTha.10",gd,"A.thaliana","1.0.0","me <me@abc.com>","me <me@abc.com>",".",license="Artistic-2.0")
install.packages("AraTha.10",repos=NULL,type="source")
library(AraTha.10)


keytypes(AraTha.10)
#  [1] "ARACYC"       "ARACYCENZYME" "CDSID"        "CDSNAME"      "DEFINITION"   "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "EXONID"       "EXONNAME"     "GENEID"       "GENENAME"
# [14] "GO"           "GOALL"        "GOID"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PMID"         "REFSEQ"       "SYMBOL"       "TAIR"         "TERM"         "TXID"         "TXNAME"

exonRanges <- exonsBy(TxDb.Athaliana.BioMart.plantsmart28, by = "gene")
head(exonRanges)

## Here we see the first few:
##
select(org.At.tair.db, head(names(exonRanges),n=22), "SYMBOL","TAIR")
#         TAIR  SYMBOL
# 1  AT1G01010 ANAC001
# 2  AT1G01010  NAC001
# 3  AT1G01020    ARV1
# 4  AT1G01030    NGA3
# 5  AT1G01040    ASU1
# 6  AT1G01040  ATDCL1
# 7  AT1G01040     CAF
# 8  AT1G01040    DCL1
# 9  AT1G01040   EMB60
# 10 AT1G01040   EMB76
# 11 AT1G01040    SIN1
# 12 AT1G01040    SUS1
# 13 AT1G01046    <NA>
# 14 AT1G01050  AtPPa1
# 15 AT1G01050    PPa1


dim(select(org.At.tair.db, names(exonRanges), "SYMBOL","TAIR"))
# 'select()' returned 1:many mapping between keys and columns
# [1] 42610     2


map_ENSEMBL_TAIR <- select(org.At.tair.db, names(exonRanges), "SYMBOL","TAIR")

map_ENSEMBL_TAIR[which(map_ENSEMBL_TAIR$TAIR=="AT3G01090")]

map_ENSEMBL_TAIR[which(map_ENSEMBL_TAIR$TAIR %in% c("AT3G01090", "ATMG00030")),]

Genes2kb[queryHits(hits)]$ID
# [1] "gene:AT3G01090" "gene:ATMG00030"

gsub("gene:","",Genes2kb[queryHits(hits)]$ID)
# [1] "AT3G01090" "ATMG00030"

map_ENSEMBL_TAIR[which(map_ENSEMBL_TAIR$TAIR %in% gsub("gene:","",Genes2kb[queryHits(hits)]$ID)),]
#            TAIR  SYMBOL
# 17757 AT3G01090  AKIN10
# 17758 AT3G01090   KIN10
# 17759 AT3G01090 SNRK1.1
# 42449 ATMG00030 ORF107A

dmgs_map_ENSEMBLE_TAIR <- map_ENSEMBL_TAIR[which(map_ENSEMBL_TAIR$TAIR %in% gsub("gene:","",Genes2kb[queryHits(hits)]$ID)),]

dmgs_map_ENSEMBLE_TAIR$SYMBOL
# [1] "AKIN10"  "KIN10"   "SNRK1.1" "ORF107A"

# irfunc <- select(AraTha.10,key="AT1G01050",keytype = "SYMBOL",columns = c("GO","TERM"))
# irfunc <- select(AraTha.10,key="AT3G01090",keytype = "SYMBOL",columns = c("GO","TERM"))
# irfunc <- select(AraTha.10,key="AT4G00080",keytype = "SYMBOL",columns = c("GO","TERM"))

irfunc <- select(AraTha.10,key="PPa1",keytype = "SYMBOL",columns = c("GO","TERM"))

irfunc[,c("ONTOLOGY","TERM","EVIDENCE","SYMBOL")]
#    ONTOLOGY                                                    TERM EVIDENCE SYMBOL
# 1        BP RNA splicing, via endonucleolytic cleavage and ligation      RCA   PPa1
# 2        MF                        inorganic diphosphatase activity      IDA   PPa1
# 3        CC                                                 nucleus      IDA   PPa1
# 4        CC                                               cytoplasm      IDA   PPa1
# 5        CC                                               cytoplasm      ISM   PPa1
# 6        CC                                                 cytosol      IDA   PPa1
# 7        BP                                       metabolic process      ISS   PPa1
# 8        BP                         methionine biosynthetic process      RCA   PPa1
# 9        CC                                                membrane      ISS   PPa1
# 10       BP                                           lipid storage      IMP   PPa1

irfunc <- select(AraTha.10, key=dmgs_map_ENSEMBLE_TAIR$SYMBOL,keytype = "SYMBOL",columns = c("GO","TERM"))

irfunc
#     SYMBOL         GO EVIDENCE ONTOLOGY                                           TERM
# 1   AKIN10 GO:0000152      IPI       CC               nuclear ubiquitin ligase complex
# 2   AKIN10 GO:0003006      IMP       BP developmental process involved in reproduction
# 3   AKIN10 GO:0004672      ISS       MF                        protein kinase activity
# 4   AKIN10 GO:0004672      IDA       MF                        protein kinase activity
# 5   AKIN10 GO:0004674      IDA       MF       protein serine/threonine kinase activity
# 6   AKIN10 GO:0005515      IPI       MF                                protein binding
# 7   AKIN10 GO:0005634      ISM       CC                                        nucleus
# 8   AKIN10 GO:0006007      RCA       BP                      glucose catabolic process
# 9   AKIN10 GO:0006486      RCA       BP                          protein glycosylation
# 10  AKIN10 GO:0009594      IDA       BP                          detection of nutrient
# 11  AKIN10 GO:0009738      IMP       BP      abscisic acid-activated signaling pathway
# 12  AKIN10 GO:0010050      IMP       BP                        vegetative phase change
# 13  AKIN10 GO:0010182      IMP       BP               sugar mediated signaling pathway
# 14  AKIN10 GO:0010260      IMP       BP                        animal organ senescence
# 15  AKIN10 GO:0080022      IMP       BP                       primary root development
# 16   KIN10 GO:0000152      IPI       CC               nuclear ubiquitin ligase complex
# 17   KIN10 GO:0003006      IMP       BP developmental process involved in reproduction
# 18   KIN10 GO:0004672      IDA       MF                        protein kinase activity
# 19   KIN10 GO:0004672      ISS       MF                        protein kinase activity
# 20   KIN10 GO:0004674      IDA       MF       protein serine/threonine kinase activity
# 21   KIN10 GO:0005515      IPI       MF                                protein binding
# 22   KIN10 GO:0005634      ISM       CC                                        nucleus
# 23   KIN10 GO:0006007      RCA       BP                      glucose catabolic process
# 24   KIN10 GO:0006486      RCA       BP                          protein glycosylation
# 25   KIN10 GO:0009594      IDA       BP                          detection of nutrient
# 26   KIN10 GO:0009738      IMP       BP      abscisic acid-activated signaling pathway
# 27   KIN10 GO:0010050      IMP       BP                        vegetative phase change
# 28   KIN10 GO:0010182      IMP       BP               sugar mediated signaling pathway
# 29   KIN10 GO:0010260      IMP       BP                        animal organ senescence
# 30   KIN10 GO:0080022      IMP       BP                       primary root development
# 31 SNRK1.1 GO:0000152      IPI       CC               nuclear ubiquitin ligase complex
# 32 SNRK1.1 GO:0003006      IMP       BP developmental process involved in reproduction
# 33 SNRK1.1 GO:0004672      IDA       MF                        protein kinase activity
# 34 SNRK1.1 GO:0004672      ISS       MF                        protein kinase activity
# 35 SNRK1.1 GO:0004674      IDA       MF       protein serine/threonine kinase activity
# 36 SNRK1.1 GO:0005515      IPI       MF                                protein binding
# 37 SNRK1.1 GO:0005634      ISM       CC                                        nucleus
# 38 SNRK1.1 GO:0006007      RCA       BP                      glucose catabolic process
# 39 SNRK1.1 GO:0006486      RCA       BP                          protein glycosylation
# 40 SNRK1.1 GO:0009594      IDA       BP                          detection of nutrient
# 41 SNRK1.1 GO:0009738      IMP       BP      abscisic acid-activated signaling pathway
# 42 SNRK1.1 GO:0010050      IMP       BP                        vegetative phase change
# 43 SNRK1.1 GO:0010182      IMP       BP               sugar mediated signaling pathway
# 44 SNRK1.1 GO:0010260      IMP       BP                        animal organ senescence
# 45 SNRK1.1 GO:0080022      IMP       BP                       primary root development
# 46 ORF107A GO:0003674       ND       MF                             molecular_function
# 47 ORF107A GO:0005739      ISM       CC                                  mitochondrion
# 48 ORF107A GO:0008150       ND       BP                             biological_process

exonRanges[exonRanges$TAIR=="AT1G01050",]

exonRanges[exonRanges$TAIR=="AT3G01090",]



# library(org.Hs.eg.db)

# biocLite("TxDb.Celegans.UCSC.ce6.ensGene")
# library(TxDb.Celegans.UCSC.ce6.ensGene)
# biocLite("celegans.db")
# library(org.Ce.eg.db)
# library(OrganismDbi)
# gd<-list(join1=c(GO.db="GOID",org.Ce.eg.db="GO"),join2=c(org.Ce.eg.db="ENTREZID",TxDb.Celegans.UCSC.ce6.ensGene="GENEID"))
# makeOrganismPackage("Cen.ele6",gd,"C.elegans","1.0.0","me <me@abc.com>","me <me@abc.com>",".",license="Artistic-2.0")
# install.packages("Cen.ele6",repos=NULL,type="source")






nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, num.cores = 64L, verbose = TRUE)
## 201907230 1030am saved:
# save(nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, file="/data/users/twm118/memoryLine/memoryLine_100k/nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR.RData")
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull2P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="Weibull2P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull3P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="Weibull3P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma2P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="Gamma2P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma3P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="Gamma3P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma3P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="GGamma3P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma4P = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="GGamma4P", num.cores = 64L, verbose = TRUE)
nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_ECDF = nonlinearFitDist(ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, column = 9, dist.name="ECDF", num.cores = 64L, verbose = TRUE)


PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)

# coral3_32_HD_CHG_sum_mincov4_Genotype2_Underside_PS_CHG_sum_GGamma3P <- getPotentialDIMP(LR = coral3_32_HD_CHG_sum_mincov4_Genotype2_Underside, nlms = coral3_32_HD_CHG_sum_mincov4_Genotype2_Underside_GGamma3P_nlms, dist.name="Weibull2P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)

PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull2P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull2P, dist.name="Weibull2P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)

PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull3P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull3P, dist.name="Weibull3P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
## 20200205_1246pm works today, not sure why??
#
#NO: PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull3P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Weibull3P, dist.name="Weibull3P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
# Error in pweibull(q - m[3], shape = m[1], scale = m[2], lower.tail = FALSE) :
#   Non-numeric argument to mathematical function

PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma2P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma2P, dist.name="Gamma2P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma3P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_Gamma3P, dist.name="Gamma3P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma3P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma3P, dist.name="GGamma3P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma4P = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_GGamma4P, dist.name="GGamma4P", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)
PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_ECDF = getPotentialDIMP(LR = ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, nlms = nlms_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_ECDF, dist.name="ECDF", div.col = 9, alpha = 0.05, tv.col = 7, tv.cut = 0.2)

## 20190730 1031am saved:
# save(PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, file="/data/users/twm118/memoryLine/memoryLine_100k/PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR.RData")

cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR = estimateCutPoint(LR = PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, control.names = c( "W_3_4", "W_3_5"),treatment.names = c("M_3_1","M_3_2","M_3_3","M_3_4","M_3_5"),div.col = 9, verbose = FALSE)

class(cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR)
# [1] "CutPoint" "list"

library(sloop)
# MethylIT:::print.CutPoint
sloop::s3_dispatch(print(cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR))
## WRONG: only shows one defined method even though it looked for 3:
#    print.CutPoint
#    print.list
# => print.default

## GOOD: shows two methods defined and it chose the right one:
# => print.CutPoint
#    print.list
#  * print.default

cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR$cutpoint
# [1] 2.514052
## 20200130: something changed!
# [1] 2.669684


## 20190730 1031am saved:
# save(cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, file="/data/users/twm118/memoryLine/memoryLine_100k/cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR.RData")
## 20200130: something changed!, saved:
# save(cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, file="/data/users/twm118/memoryLine/memoryLine_100k/cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR_20200130somethingChanged.RData")


DIMPs_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR <- selectDIMP(PS_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, div.col = 9, cutpoint=cutpoints_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR$cutpoint)

## 20190730 1031am saved:
# save(DIMPs_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR, file="/data/users/twm118/memoryLine/memoryLine_100k/DIMPs_ml_HD_CG_sum_mincov4_W31W32W33_Ref_W34W35M31M32M33M34M35_LR.RData")
# ftp://ftp.ensemblgenomes.org/pub/plants/current/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3.gz
AG_gff3_ftp = import("ftp://ftp.ensemblgenomes.org/pub/plants/current/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.44.gff3.gz")
AG_gff3 = import("/data/TAIR10_gff3/Arabidopsis_thaliana.TAIR10.38.gff3.gz")
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

## 20200205_1250pm
# GRanges object with 6 ranges and 2 metadata columns:
#       seqnames      ranges strand |    GeneID     DIMPs
#          <Rle>   <IRanges>  <Rle> |  <factor> <integer>
#   [1]        1   1631-7899      + | AT1G01010         1
#   [2]        1  4788-11130      - | AT1G01020         1
#   [3]        1 21121-33227      + | AT1G01040         9
#   [4]        1 29170-35171      - | AT1G01050         5
#   [5]        1 31365-39871      - | AT1G01060         1
#   [6]        1 49953-56737      + | AT1G01110         2
#   -------
#   seqinfo: 7 sequences from an unspecified genome; no seqlengths



# Genes_DIMPs_coral3_32_mincov4_CG_sum_GGamma3P <- uniqueGRanges(list(DIMPs_coral3_32_mincov4_CG_sum_GGamma3P_Genotype1_Upperside_TG1up03_gene[,2],DIMPs_coral3_32_mincov4_CG_sum_GGamma3P_Genotype1_Upperside_TG1up04_gene[,2],DIMPs_coral3_32_mincov4_CG_sum_GGamma3P_Genotype1_Upperside_TG1up07_gene[,2],DIMPs_coral3_32_mincov4_CG_sum_GGamma3P_Genotype1_Upperside_TG1up08_gene[,2],DIMPs_coral3_32_mincov4_CG_sum_GGamma3P_Genotype1_Upperside_CG1up17_gene[,2],DIMPs_coral3_32_mincov4_CG_sum_GGamma3P_Genotype1_Upperside_CG1up18_gene[,2],DIMPs_coral3_32_mincov4_CG_sum_GGamma3P_Genotype1_Upperside_CG1up19_gene[,2],DIMPs_coral3_32_mincov4_CG_sum_GGamma3P_Genotype1_Upperside_CG1up20_gene[,2]), type = "equal", verbose = TRUE, ignore.strand = TRUE)

Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35 = uniqueGRanges(list(DIMPs_ml_M_3_1_gene[,2],DIMPs_ml_M_3_2_gene[,2],DIMPs_ml_M_3_3_gene[,2],DIMPs_ml_M_3_4_gene[,2],DIMPs_ml_M_3_5_gene[,2],DIMPs_ml_W_3_4_gene[,2],DIMPs_ml_W_3_5_gene[,2]), type = "equal", verbose = TRUE, chromosomes = c("1", "2", "3", "4", "5"), ignore.strand = TRUE )

colnames( mcols(Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35)) <- c("M_3_1","M_3_2","M_3_3","M_3_4","M_3_5","W_3_4","W_3_5")

## 20190730 1031am saved:
# save(Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35, file="/data/users/twm118/memoryLine/memoryLine_100k/Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35.RData"

Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35_ID = subsetByOverlaps(Genes2kb, Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35, type = "equal", ignore.strand = FALSE)

dmps = data.frame( mcols( Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35 ) )
dmps = apply( dmps, 2, as.numeric )
rownames(dmps) <- Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35_ID$gene_id

condition = data.frame(condition = factor(c("TT","TT","TT","TT","TT","CT","CT"), levels = c("CT", "TT")))
rownames(condition)
# [1] "1" "2" "3" "4" "5" "6" "7"
names(mcols(Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35))
# [1] "M_3_1" "M_3_2" "M_3_3" "M_3_4" "M_3_5" "W_3_4" "W_3_5"

rownames(condition) <- names(mcols(Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35))

DIMR_ml_samples_W34W35M31M32M33M34M35 <- glmDataSet(GR = Genes_DIMPs_ml_samples_W34W35M31M32M33M34M35, colData = condition)

class(DIMR_ml_samples_W34W35M31M32M33M34M35)
# [1] "RangedGlmDataSet"

sloop::s3_dispatch(print(DIMR_ml_samples_W34W35M31M32M33M34M35))
## WRONG:
#    print.RangedGlmDataSet
# => print.default
## GOOD:
# => print.RangedGlmDataSet
#  * print.default

## 20190730 1031am saved:
# save(DIMR_ml_samples_W34W35M31M32M33M34M35, file="/data/users/twm118/memoryLine/memoryLine_100k/DIMR_ml_samples_W34W35M31M32M33M34M35.RData")

DMGs_ml_samples_W34W35M31M32M33M34M35_countTest2 = countTest2(DIMR_ml_samples_W34W35M31M32M33M34M35, num.cores = 24L, countFilter = T)
dim(mcols(DMGs_ml_samples_W34W35M31M32M33M34M35_countTest2))
# [1] 29 14
# [1] 31 14    ## 20190730 1657pm found 2 more DMGs because it used the latest version of MethylIT:
## 20200130: something changed:
# [1] 21 15












# ## 20190730 1031am saved:
# save(, file="/data/users/twm118/memoryLine/memoryLine_100k/.RData")
#
# ## 20190730 1031am saved:
# save(, file="/data/users/twm118/memoryLine/memoryLine_100k/.RData")
#
# ## 20190730 1031am saved:
# save(, file="/data/users/twm118/memoryLine/memoryLine_100k/.RData")









