library(readr)
library(dplyr)
library(tidyr)


# Parse functional table: remove duplicates, score again genomes, split column
full = read_tsv("../table/all_23genomes_fnc.xls")
variants_genotype = read_delim("../../final/2017-09-06_elephants/b1-gatk-haplotype-joint-annotated-decomposed-parsed-genotype.tsv", col_names = TRUE, delim="\t")
variants_genotype$pos = variants_genotype$pos - 1

full_clean = full %>%
    separate(change, into=c("vc_aa_ref","vc_aa_alt"),
             sep = "[0-9]+", remove = F,extra = "merge")

# full_clean$mammoth_genome = rowSums(data.frame(wr=!grepl("None", full_clean$Wrangel) & !is.na(full_clean$Wrangel),
#                                                oi=!grepl("None", full_clean$oimyakon) & !is.na(full_clean$oimyakon),
#                                                m4=!grepl("None", full_clean$M4) & !is.na(full_clean$M4),
#                                                m25=!grepl("None", full_clean$M25) & !is.na(full_clean$M25)))
# 
# full_clean$asian_genome = rowSums(data.frame(  uno=!grepl("None", full_clean$Uno) & !is.na(full_clean$Uno),
#                                                asha=!grepl("None", full_clean$Asha) & !is.na(full_clean$Asha),
#                                                pa=!grepl("None", full_clean$Parvathy) & !is.na(full_clean$Parvathy)))

full_clean_gen = left_join(full_clean %>% mutate(id=paste0(chrom, pos)),
                            variants_genotype %>% mutate(id=paste0(chrom, pos)) %>%
                                select(-chrom, -pos),
                            by="id", suffix = c("", "_GT"))

write.table(full_clean_gen %>% select(-id) %>% distinct(),
            "../table/all_23genomes_fnc_clean.xls", sep="\t", row.names=F)
