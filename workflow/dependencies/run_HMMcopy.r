# Below commands are for creating segmentation files and producing some QC images
# outputBasename used to construct names for files produced by this script

library(HMMcopy)

cmd_args=commandArgs(trailingOnly = TRUE)
# Arguments should be passed as: normal.wig, tumor.wig, refGC.wig, ref_mappable.wig, outputBasename
normalReads<-cmd_args[1]
tumorReads<-cmd_args[2]
gcContent<-cmd_args[3]
refMappable<-cmd_args[4]
outputBasename<-cmd_args[5]

tum_uncorrected_reads <- wigsToRangedData(tumorReads, gcContent, refMappable)
norm_uncorrected_reads <- wigsToRangedData(normalReads, gcContent, refMappable)
tum_corrected_copy <- correctReadcount(tum_uncorrected_reads)
norm_corrected_copy <- correctReadcount(norm_uncorrected_reads)

# Should take no longer than a few minutes on a human genome.
# The correctReadcount requires at least about 1000 bins to work properly.
# 2. Segmentation

# Below commands in R
# Normalizing Tumour by Normal
tum_corrected_copy$copy <- tum_corrected_copy$copy - norm_corrected_copy$copy

param <- HMMsegment(tum_corrected_copy, getparam = TRUE) # retrieve converged parameters via EM
param$mu <- log(c(1, 1.4, 2, 2.7, 3, 4.5) / 2, 2)
param$m <- param$mu
segmented_copy <- HMMsegment(tum_corrected_copy, param) # perform segmentation via Viterbi

# 3. Export
# Export to SEG format for CNAseq segmentation
segFile<-paste(outputBasename, "seg", sep = ".")
tsvFile<-paste(outputBasename, "tsv", sep = ".")

rangedDataToSeg(tum_corrected_copy, file = segFile)
write.table(segmented_copy$segs, file = tsvFile, quote = FALSE, sep = "\t")

# 4. Visualization - produce some images with hard-coded dimensions

# Bias plots:
print("Producing CG Bias plot...")
png(filename = paste(outputBasename,"bias_plot","png", sep="."),width = 1200, height = 580, units="px", pointsize=15, bg="white")
plotBias(tum_corrected_copy)  # May be one plot per comparison  1200x580
dev.off()

# Segmentation plots:
# need to do it one plot per chromosome 1200x450
print("Producing Segmentation plots...")
par(mfrow = c(1, 1))
for (c in 1:length(chroms)) {
 png(filename = paste(outputBasename,"s_plot", chroms[c], "png", sep="."),width = 1200, height = 450, units="px", pointsize=15, bg="white")
 plotSegments(tum_corrected_copy, segmented_copy, pch = ".", ylab = "Tumour Copy Number", xlab = "Chromosome Position",chr = chroms[c], main = paste("Segmentation for Chromosome",chroms[c], sep=" "))
 cols <- stateCols()
 legend("topleft", c("HOMD", "HETD", "NEUT", "GAIN", "AMPL", "HLAMP"), fill = cols, horiz = TRUE, bty = "n", cex = 0.9)
 dev.off()
}


# Correction plots:
# need to do it one plot per chromosome 1200x680

print("Producing Correction plots...")
chroms<-unique(segmented_copy$segs$chr)

slice_valid<-function (correctOutput) {
    copy <- correctOutput$reads/median(correctOutput$reads, na.rm = TRUE)
    if (length(copy[is.na(copy)])> 0) {
      FALSE 
    } else {
      TRUE
    }
}

for (c in 1:length(chroms)) {
 
 if (!grepl("_",chroms[c]) && slice_valid(tum_corrected_copy[paste(chroms[c])]) > 0 ) {
  print(paste("Plotting for chromosome",c,sep=" "))
  png(filename = paste(outputBasename,"c_plot",chroms[c],"png", sep="."),width = 1200, height = 680, units="px", pointsize=15, bg="white")
  plotCorrection(tum_corrected_copy, pch = ".", chr= chroms[c], na.rm = TRUE)
  dev.off()
 } else {
  print(paste("Chromosome",c,"cannot be plotted with  plotCorrection",sep = " "))
 }
}

