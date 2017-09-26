cd $PATHROOT/mammoth_code/table/ann

cut -f 1-2 ../all_23genomes.xls | sed 's/"//g' | sed 1d | awk '{print $1"\t"$2"\t"$2+1"\t"$1":"$2}' | sort -u > african_pos.bed

$PATHROOT/conda/bin/CrossMap.py bed $PATHROOT/tools/liftover/loxAfr3ToHg19.over.chain.gz african_pos.bed  hg19_pos.bed

awk '{print $1"\t"$2+1}' hg19_pos.bed | sort -k1,1 -k2n > hg19_dbnsfp.bed

java -classpath $PATHROOT/tools/dbNSFP search_dbNSFP35a -i hg19_dbnsfp.bed -v hg19 -o hg19_ann.txt

cd ..

$PATHROOT/conda/bin/python $PATHROOT/mammoth_code/scripts/merge_func_table.py  > all_23genomes_fnc.xls

