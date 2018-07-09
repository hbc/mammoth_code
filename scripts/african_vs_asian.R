library(tidyverse)
df = read_csv("../table/all_23genomes_fnc_clean_high.csv.gz")

maximus = df %>% select(chrom:vc_aa_alt, contains("maximus"), contains("L_africana"), -contains("_1")) %>% filter(E_maximus_D=="Hom", E_maximus_E=="Hom", EmaximusZ=="Hom",
                L_africana_B=="None", L_africana_C=="None",
                L_africana_B_GT!="./.", L_africana_C_GT!="./.") %>% distinct()
write_csv(maximus, "../table/african_vs_asian.csv.gz")
