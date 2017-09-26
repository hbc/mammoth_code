library(readr)
library(dplyr)
library(tidyr)

lynch = read.csv("../table/Lynch_table.csv")
variants = read_delim("../../final/2017-09-06_elephants/split/merged-parsed-flank-wheader.tsv", delim="\t")
variants_seq = read_delim("../../final/2017-09-06_elephants/split/merged.fa", col_names = FALSE, delim="\t")
variants$african_region = variants_seq$X3
variants = variants %>% distinct()

lynch$id = paste0(lynch$loxodonta.scaffold, lynch$position.in.scaffold)
variants$pos = variants$pos - 1
variants$id = paste0(variants$chrom, variants$pos)

table = variants

table$is_lynch_table = table$id %in% lynch$id

write.table(table %>% select(-id), "../table/all_23genomes.xls", sep="\t", row.names=F)

lynch$is_in_bcbio_analysis = lynch$id %in% variants$id
write.table(lynch %>% select(-id), "../table/lynch_w_23bcbio.xls", sep="\t", row.names=F)

