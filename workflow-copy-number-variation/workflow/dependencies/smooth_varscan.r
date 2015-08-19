# This script relies on wrapper which is supposed to set environment variable
# R_LIBS that points to a directory with needed modules

library(DNAcopy)
cmd_args=commandArgs(trailingOnly = TRUE)

varscanResults<-cmd_args[1]
sampleID<-cmd_args[2]
varscanSmoothed = paste(varscanResults, "segmented", sep=".")

cn <- read.table(varscanResults,header=T)
CNA.object <-CNA( genomdat = cn[,6], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio', sampleid = sampleID)
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
segs2 = segs$output
write.table(segs2[,2:6], file=varscanSmoothed, row.names=F, col.names=F, quote=F, sep="\t")

options(bitmapType="cairo")
# Visualizaton 1 of 2
png_plot_w = paste(varscanResults, "w_plot", "png", sep=".")

png(filename=png_plot_w, bg="white", width=600, height=600, units="px")
plot(segs, plot.type="w")
blah<-dev.off()

# Visualizaton 2 of 2
png_plot_s = paste(varscanResults, "s_plot", "png", sep=".")

png(filename=png_plot_s, bg="white", width=600, height=600, units="px")
plot(segs, plot.type="s")
blah<-dev.off()
