library(sleuth)
library(dplyr)

directories = dirname(unlist(snakemake@input))

manifest = read.csv(snakemake@params[['manifest']])
manifest = dplyr::select(manifest, sample=name, condition)
manifest = dplyr::mutate(manifest, path=directories)

so = sleuth_prep(manifest, extra_bootstrap_summary=TRUE)
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)
sleuth_table = dplyr::filter(sleuth_table, qval < snakemake@params[['qval']])

write.table(sleuth_table, snakemake@output[[1]])
