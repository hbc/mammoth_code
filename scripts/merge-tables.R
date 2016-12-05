
lynch = read.csv("../test/Lynch_table.csv")
wragel = read_delim("~/orch/scratch/church_mammoth/res_wrangel/changes.tsv", delim =" ")
oimyako = read_delim("~/orch/scratch/church_mammoth/res_oimyako/changes.tsv", delim=" ")
variants = read_delim("~/orch/scratch/church_mammoth/mammoth_vc/work/joint/gatk-haplotype-joint/batch1/split/merged-parsed-flank-wheader.tsv", delim="\t")
variants_seq = read_delim("~/orch/scratch/church_mammoth/mammoth_vc/work/joint/gatk-haplotype-joint/batch1/split/merged.fa", col_names = FALSE, delim="\t")
variants$african_region = variants_seq$X3
variants = variants %>% distinct()

wragel$id = paste0(wragel$chrom, wragel$genome_pos)
oimyako$id = paste0(oimyako$chrom, oimyako$genome_pos)
lynch$id = paste0(lynch$loxodonta.scaffold, lynch$position.in.scaffold)
variants$pos = variants$pos - 1
variants$id = paste0(variants$chrom, variants$pos)

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
           wrangel_annotated = `mapped/annotated.x`, wrangel_found_gap = found_gap.x,
           wrangel_missing_exons = missing_exons.x,
           wrangel_flank = flank_mammoth.x, wrangle_ref_flank = flank_african.x,
           oimyako_new_nc = ref_new.y, oimyako_new_aa = new_aa.y,
           oimyako_annotated = `mapped/annotated.y`, oimyako_found_gap = found_gap.y,
           oimyako_missing_exons = missing_exons.y,
           oimyako_flank = flank_mammoth.y, oimyako_ref_flank = flank_african.y,
           vc_ref=ref, vc_alt=alt,
           vc_gene = gene,vc_change=change,
           M25, M4, Wrangel, oimyakon, Asha, Parvathy, Uno,
           M25_flank = M25_1, M4_flank = M4_1,
           Wrangel_flank = Wrangel_1, Oimyako_flank = oimyakon_1,
           Asha_flank = Asha_1,
           Parvathy_flank = Parvathy_1, Uno_flank = Uno_1,
           vc_african_flank = african_region
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


# remove not exact variants
table_filtered = table %>%
    distinct(genome_pos, aa_ref, aa_pos, tx, exon, gene, wrangel_new_aa, oimyako_new_aa, .keep_all = TRUE) %>%
    filter(wrangel_new_aa == oimyako_new_aa | is.na(oimyako_new_aa) | is.na(wrangel_new_aa)) #%>%
#filter(gene=="A2ML1") %>%
#select(genome_pos, aa_ref, aa_pos, tx, exon, gene, wrangel_new_aa, oimyako_new_aa, vc_change)

write.table(table_filtered %>% select(-id), "../test/all_genomes.xls", sep="\t", row.names=F)

lynch$is_in_bcbio_analysis = lynch$id %in% variants$id
write.table(lynch %>% select(-id), "../test/lynch_w_bcbio.xls", sep="\t", row.names=F)

