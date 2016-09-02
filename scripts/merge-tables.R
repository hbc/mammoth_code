lynch = read.csv("../test/Lynch_table.csv")
wragel = read.table("~/orch/scratch/church_mammoth/res_wrangel/changes.tsv", sep=" ", header=T)
oimyako = read.table("~/orch/scratch/church_mammoth/res_oimyako/changes.tsv", sep=" ", header=T)
variants = read.table("~/orch/scratch/church_mammoth/mammoth_vc/work/joint/gatk-haplotype-joint/batch1/split/merged-parsed-flank-wheader.tsv", header=T,sep="\t")

wragel$id = paste0(wragel$chrom, wragel$genome_pos)
oimyako$id = paste0(oimyako$chrom, oimyako$genome_pos)
lynch$id = paste0(lynch$loxodonta.scaffold, lynch$position.in.scaffold+1)
variants$id = paste0(variants$chrom, variants$pos)

library(dplyr)
all = full_join(wragel, oimyako, by="id")
all_w_vc = full_join(all, variants, by="id")

combined_gene = as.character(all_w_vc$gene.x)
combined_gene[is.na(combined_gene)] = as.character(all_w_vc$gene.y[is.na(combined_gene)])
combined_gene[is.na(combined_gene)] = as.character(all_w_vc$gene[is.na(combined_gene)])

combined_tx = as.character(all_w_vc$tx.x)
combined_tx[is.na(combined_tx)] = as.character(all_w_vc$tx.y[is.na(combined_tx)])

combined_exon = as.character(all_w_vc$exon.x)
combined_exon[is.na(combined_exon)] = as.character(all_w_vc$exon.y[is.na(combined_exon)])

combined_txpos = as.character(all_w_vc$tx_pos.x)
combined_txpos[is.na(combined_txpos)] = as.character(all_w_vc$tx_pos.y[is.na(combined_txpos)])

combined_refnc = as.character(all_w_vc$ref_nc.x)
combined_refnc[is.na(combined_refnc)] = as.character(all_w_vc$ref_nc.y[is.na(combined_refnc)])

combined_aaref = as.character(all_w_vc$ref_aa.x)
combined_aaref[is.na(combined_aaref)] = as.character(all_w_vc$ref_aa.y[is.na(combined_aaref)])

combined_aapos = as.character(all_w_vc$aa_pos.x)
combined_aapos[is.na(combined_aapos)] = as.character(all_w_vc$aa_pos.y[is.na(combined_aapos)])

combined_enum = as.character(all_w_vc$number.x)
combined_enum[is.na(combined_enum)] = as.character(all_w_vc$number.y[is.na(combined_enum)])

combined_chrom = as.character(all_w_vc$chrom.x)
combined_chrom[is.na(combined_chrom)] = as.character(all_w_vc$chrom.y[is.na(combined_chrom)])
combined_chrom[is.na(combined_chrom)] = as.character(all_w_vc$chrom[is.na(combined_chrom)])

combined_pos = (all_w_vc$genome_pos.x)
combined_pos[is.na(combined_pos)] = (all_w_vc$genome_pos.y[is.na(combined_pos)])
combined_pos[is.na(combined_pos)] = (all_w_vc$pos[is.na(combined_pos)])

table = all_w_vc %>%  
    mutate(gene.1=combined_gene, tx=combined_tx, exon=combined_exon, 
           exon_number=number.x,
           tx_pos=combined_txpos, ref_nc=combined_refnc, 
           aa_pos = combined_aapos, aa_ref = combined_aaref,
           chrom = combined_chrom, genome_pos = combined_pos) %>%
    select(id, gene=gene.1, tx, exon, exon_number, tx_pos, ref_nc, aa_pos, aa_ref,
           chrom, genome_pos,
           wrangel_new_nc=ref_new.x,wrangel_new_aa = new_aa.x, 
           wrangel_annotated = mapped.annotated.x, wrangel_found_gap = found_gap.x,
           wrangel_missing_exons = missing_exons.x, wrangle_flank = flank.x,
           oimyako_new_nc = ref_new.y, oimyako_new_aa = new_aa.y,
           oimyako_annotated = mapped.annotated.y, oimyako_found_gap = found_gap.y,
           oimyako_missing_exons = missing_exons.y, oimyako_flank = flank.y,
           vc_ref=ref, vc_alt=alt, 
           vc_gene = gene,vc_change=change,
           M25, M4, Asha, Parvathy, Uno,
           M25_flank = M25.1, M4_flank = M4.1, Asha_flank = Asha.1, 
           Parvathy_flank = Parvathy.1, Uno_flank = Uno.1
           )
table$wrangel_annotated = gsub("/"," out of ", table$wrangel_annotated)
table$oimyako_annotated = gsub("/"," out of ", table$oimyako_annotated)

table$is_lynch_table = table$id %in% lynch$id 


table$mammoth_genome = rowSums(data.frame(wr=!is.na(table$wrangel_new_aa), 
                                          oi=!is.na(table$oimyako_new_aa),
                                          m4=!grepl("None", table$M4) & !is.na(table$M4),
                                          m25=!grepl("None", table$M25) & !is.na(table$M25)))

table$asian_genome = rowSums(data.frame(  uno=!grepl("None", table$Uno) & !is.na(table$Uno),
                                          asha=!grepl("None", table$Asha) & !is.na(table$Asha),
                                          pa=!grepl("None", table$Parvathy) & !is.na(table$Parvathy)))

write.table(table %>% select(-id), "../test/all_genomes.xls", sep="\t", row.names=F)

lynch$is_in_bcbio_analysis = lynch$id %in% variants$id
write.table(lynch %>% select(-id), "../test/lynch_w_bcbio.xls", sep="\t", row.names=F)



