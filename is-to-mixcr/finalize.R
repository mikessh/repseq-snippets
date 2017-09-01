#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
df = fread(args[1])[nSeqCDR3 != ""]  # filter cases with no CDR3

colnames(df) = c("count", "cdr3nt", "cdr3aa", "v", "d", "j", "cdr3Start", "vEnd", "dStart", "dEnd", "jStart")
df$freq = df$count / sum(df$count)
df$vEnd = df$vEnd - df$cdr3Start
df$dStart = df$dStart - df$cdr3Start
df$dEnd = df$dEnd - df$cdr3Start
df$jStart = df$jStart - df$cdr3Start
df$cdr3Start = NULL

setcolorder(df, c("count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "vEnd", "dStart", "dEnd", "jStart"))

fwrite(df, file = args[1], quote = F, sep = "\t", row.names = F)