library(readr)
library(dplyr)
library(tidyr)


# Parse functional table: remove duplicates, score again genomes, split column
full = read_tsv("../table/all_genomes_fnc.xls")
variants_genotype = read_delim("../../mammoth_vc/work/joint/gatk-haplotype-joint/batch1/split/batch1-joint-effects-filterSNP-filterINDEL-gatkclean-parsed-genotype.tsv", col_names = TRUE, delim="\t")
variants_genotype$pos = variants_genotype$pos - 1

full_clean = full %>%
    distinct(genome_pos, aa_ref, aa_pos, tx, exon, gene, wrangel_new_aa, oimyako_new_aa, .keep_all = TRUE) %>%
    filter(wrangel_new_aa == oimyako_new_aa | is.na(oimyako_new_aa) | is.na(wrangel_new_aa)) %>%
    separate(vc_change, into=c("vc_aa_ref","vc_aa_alt"),sep = "[0-9]+", remove = F,extra = "merge")

full_clean$mammoth_genome = rowSums(data.frame(wr=!grepl("None", full_clean$Wrangel) & !is.na(full_clean$Wrangel),
                                               oi=!grepl("None", full_clean$oimyakon) & !is.na(full_clean$oimyakon),
                                               m4=!grepl("None", full_clean$M4) & !is.na(full_clean$M4),
                                               m25=!grepl("None", full_clean$M25) & !is.na(full_clean$M25)))

full_clean$asian_genome = rowSums(data.frame(  uno=!grepl("None", full_clean$Uno) & !is.na(full_clean$Uno),
                                               asha=!grepl("None", full_clean$Asha) & !is.na(full_clean$Asha),
                                               pa=!grepl("None", full_clean$Parvathy) & !is.na(full_clean$Parvathy)))

full_clean_gen = inner_join(full_clean %>% mutate(id=paste0(chrom, genome_pos)),
                            variants_genotype %>% mutate(id=paste0(chrom, pos)) %>%
                                select(-chrom, -pos) %>%
                                rename(Asha_G=Asha, M25_G=M25, M4_G=M4,
                                       Parvathy_G=Parvathy, Uno_G=Uno,
                                       Wrangel_G=Wrangel, oimyako_G=oimyakon),
                            by="id")

write.table(full_clean_gen %>% select(-id) %>% distinct(),
            "../table/all_genomes_fnc_clean.xls", sep="\t", row.names=F)
