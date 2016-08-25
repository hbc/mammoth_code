lynch = read.csv("../test/Lynch_table.csv")
wragel = read.table("~/orch/scratch/church_mammoth/res_wrangel/changes.tsv", sep=" ", header=T)
oimyako = read.table("~/orch/scratch/church_mammoth/res_oimyako/changes.tsv", sep=" ", header=T)

wragel$id = paste0(wragel$tx, wragel$aa_pos, wragel$ref_aa)
oimyako$id = paste0(oimyako$tx, oimyako$aa_pos, oimyako$ref_aa)
lynch$id = paste0(lynch$ensembl.transcript, lynch$pepetide.position, lynch$elephant.AA)

library(dplyr)
rows = rbind(wragel[,c(1:6,8,9)],oimyako[,c(1:6,8,9)]) %>% distinct()
rows$id = paste0(rows$tx, rows$aa_pos, rows$ref_aa)
t = wragel[,c(7,10,12:16)]
names(t)[1:6] = paste0(names(t)[1:6],"_wrangel")
all = left_join(rows, t, by="id")

t = oimyako[,c(7,10,12:16)]
names(t)[1:6] = paste0(names(t)[1:6],"_oimyako")
all = left_join(all, t, by="id")
all$is_lynch = all$id %in% lynch$id

all_pos=wragel$genome_pos[match(all$id,wragel$id)]
all_pos[is.na(all_pos)]=(oimyako$genome_pos[match(all$id,oimyako$id)])[is.na(all_pos)]
all_gene=as.character(wragel$gene[match(all$id,wragel$id)])
all_gene[is.na(all_gene)]=as.character((oimyako$gene[match(all$id,oimyako$id)])[is.na(all_gene)])
all_id = paste0(all_pos, all_gene)


all = all %>% select(-id) %>% distinct()
write.table(all, "../test/wrangel_oimyako.xls", sep="\t", row.names=F)


variants = read.table("~/orch/scratch/church_mammoth/mammoth_vc/work/joint/gatk-haplotype-joint/batch1/batch1-joint-effects-filterSNP-filterINDEL-gatkclean-parsed.tsv", header=T)
variants$id = paste0(variants$pos, variants$gene)

all$id = all_id
all2 = dplyr::full_join(all, variants, by="id")

all2$mammoth_genome = rowSums(data.frame(wr=!is.na(all2$new_aa_wrangel), 
                                         oi=!is.na(all2$new_aa_oimyako),
                                         m4=!grepl("None", all2$M4) & !is.na(all2$M4),
                                         m25=!grepl("None", all2$M25) & !is.na(all2$M25)))

all2$asian_genome = rowSums(data.frame(  uno=!grepl("None", all2$Uno) & !is.na(all2$Uno),
                                         asha=!grepl("None", all2$Asha) & !is.na(all2$Asha),
                                         pa=!grepl("None", all2$Parvathy) & !is.na(all2$Parvathy)))

all_genomes = all2 %>% mutate(vc_ref=ref, vc_alt=alt, gene=gene.x) %>% 
         select(tx,exon,number,tx_pos, gene, ref_nc, aa_pos, ref_aa, ref_new_wrangel, new_aa_wrangel,
                mapped.annotated_wrangel, found_gap_wrangel, missing_exons_wrangel, flank_wrangel,
                ref_new_oimyako, new_aa_oimyako,
                mapped.annotated_oimyako, found_gap_oimyako, missing_exons_oimyako, flank_oimyako,
                is_lynch, vc_ref, vc_alt, change, M25, M4, Asha, Parvathy, Uno,
                mammoth_genome, asian_genome)

write.table(all_genomes, "../test/all_genomes.xls", sep="\t", row.names=F)




