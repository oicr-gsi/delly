library(DNAcopy)
cmd_args=commandArgs(trailingOnly = TRUE)

varscanResults<-cmd_args[1]
varscanSmoothed = paste(varscanResults, "segmented", sep=".")

cn <- read.table(varscanResults,header=T)
CNA.object <-CNA( genomdat = cn[,6], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
segs2 = segs$output
write.table(segs2[,2:6], file=varscanSmoothed, row.names=F, col.names=F, quote=F, sep="\t")