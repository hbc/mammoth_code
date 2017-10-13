library(rdrop2)
token <- readRDS("~/.droptoken.rds")
# Then pass the token to each drop_ function
d = drop_acc(dtoken = token)
dropdir = "HBC Team Folder (1)/Consults/george_church/churhc_variants_mammoth"
drop_upload("../table/all_23genomes_fnc_clean.xls.zip", dest=dropdir)
