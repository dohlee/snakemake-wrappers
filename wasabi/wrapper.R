library(wasabi)

sfdirs = dirname(unlist(snakemake@input))
prepare_fish_for_sleuth(sfdirs)
