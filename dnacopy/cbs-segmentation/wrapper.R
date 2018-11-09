# __author__ = "Dohoon Lee"
# __copyright__ = "Copyright 2018, Dohoon Lee"
# __email__ = "dohlee.bioinfo@gmail.com"
# __license__ = "MIT"

# The code below is adopted from
# http://varscan.sourceforge.net/copy-number-calling.html#copy-number-segmentation

library(DNAcopy)

cn <- read.table(snakemake@input[["copynumber_calls"]], header=TRUE)
CNA.object <-CNA(genomdat=cn[,7], chrom=cn[,1], maploc=cn[,2], data.type='logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, alpha=0.01, nperm=10000, min.width=2, undo.splits="sdundo", undo.SD=3)
psegs <- segments.p(segs)
write.table(psegs, file=snakemake@output[["segments"]], row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
