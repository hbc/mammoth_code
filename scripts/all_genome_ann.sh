cd ~/scratch/church_mammoth/mammoth_code/table/ann
cut -f 9-10 ../all_genomes.xls | sed 's/"//g' | sed 1d | awk '{print $1"\t"$2"\t"$2+1"\t"$1":""}' | sort -u > african_pos.bed
CrossMap.py bed ~/scratch/church_mammoth/tools/liftover/loxAfr3ToHg19.over.chain.gz african_pos.bed  hg19_pos.bed
awk '{print $1"\t"$2+1}' hg19_pos.bed | sort -k1,1 -k2n > hg19_dbnsfp.bed
java -classpath ~/scratch/church_mammoth/tools/snpEff/db/hg19/dbNSFP search_dbNSFP32a -i hg19_dbnsfp.bed -v hg19 -o hg19_ann.txt

cd ..
~/scratch/church_mammoth/conda/bin/python ~/scratch/church_mammoth/mammoth_code/scripts/merge_func_table.py  > all_genomes_fnc.xls

