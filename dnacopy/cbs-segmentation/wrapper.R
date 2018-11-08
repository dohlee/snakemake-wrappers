# __author__ = "Dohoon Lee"
# __copyright__ = "Copyright 2018, Dohoon Lee"
# __email__ = "dohlee.bioinfo@gmail.com"
# __license__ = "MIT"

# The code below is adopted from
# http://varscan.sourceforge.net/copy-number-calling.html#copy-number-segmentation

library(DNAcopy)

cn <- read.table(snakemake@input[["copynumber_calls"]])
CNA.object <-CNA(genomdat = cn[,6], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
segs2 = segs$output
write.table(segs2[,2:6], file=snakemake@output[["segments"]], row.names=F, col.names=F, quote=F, sep="\t")
