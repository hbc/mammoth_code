library(readr)
library(dplyr)
library(tidyr)


# Parse functional table: remove duplicates, score again genomes, split column
full = read_tsv("../table/all_23genomes_fnc.xls")
variants_genotype = read_delim("../../elephants/2017-09-06_elephants/b1-gatk-haplotype-joint-annotated-decomposed-parsed-genotype.tsv", col_names = TRUE, delim="\t")
variants_genotype$pos = variants_genotype$pos - 1

full_clean = full %>%
    separate(change, into=c("vc_aa_ref","vc_aa_alt"),
             sep = "[0-9]+", remove = F,extra = "merge")


full_clean_gen = left_join(full_clean %>% mutate(id=paste0(chrom, pos)),
                            variants_genotype %>% mutate(id=paste0(chrom, pos)) %>%
                                select(-chrom, -pos),
                            by="id", suffix = c("", "_GT"))

write.table(full_clean_gen %>% select(-id) %>% distinct(),
            "../table/all_23genomes_fnc_clean.xls", sep="\t", row.names=F)


# full_clean_gen = read_tsv("../table/all_23genomes_fnc_clean.xls.gz")

# elephant_counter should count to max 12 (Asha, E_maximus_D, E_maximus_E, EmaximusZ, L_africana_B, L_africana_C, L_cyclotis_A, L_cyclotis_F, M_columbi_U, P_antiquus_N, P_antiquus_O, Parvathy)
# 
# mammoth_counter should only count to max of 9. (M25, M4, M_columbi_U, M_primigenius_G	M_primigenius_H	M_primigenius_S	M_primigenius_V, Uno,  Wrangel)


species = names(full_clean_gen)[9:31]
e = species[grepl("L_", species) | grepl("E", species) | grepl("P", species) | grepl("Asha", species) | grepl("Uno", species)]
m = species[grepl("M_am", species)]
mm = setdiff(species, c(e, m))

full_clean_gen[["elephant_counter"]] = full_clean_gen[,e] %>% 
    apply(., 1, function(r){ sum(r!="None") })

full_clean_gen[["mastodont_counter"]] = full_clean_gen[,m] %>% 
    apply(., 1, function(r){ sum(r!="None") })

full_clean_gen[["mammoth_counter"]] = full_clean_gen[,mm] %>% 
    apply(., 1, function(r){ sum(r!="None") })

write_tsv(full_clean_gen, "../table/all_23genomes_fnc_clean.xls.gz")
write_csv(full_clean_gen %>% select(-contains("_1"), -african_region),
          "../table/all_23genomes_fnc_clean_no_flank.csv.gz")

#this is to kept only very well annotated variants, with
# no warnings and in the HIGH or missense category
snpeff = read_tsv("../../elephants/2017-09-06_elephants/b1-gatk-haplotype-joint-annotated-decomposed-snpeff.tsv")
snpeff$pos = snpeff$pos - 1 

snpeff %>% filter(impact  %in%  c("HIGH") | type == "missense_variant") %>% 
    filter(!grepl("WAR", warning)) %>% 
    inner_join(full_clean_gen %>% select(-ref, -alt, -gene, -change), by = c("pos", "chrom")) %>% write_csv("../table/all_23genomes_fnc_clean_high.csv.gz")

