===============
mammoth
===============

Package created to help Church's Lab to find changes in Mammoth genomes.

## blast based analysis

We developed a python package to blast all African genes to the assembled genomes.
The command line used for both mammoth genomes are inside `analysis/blast`.

## Variant calling analysis

All genomes were analyzed with bcbio-nextgen framework using the variant calling pipeline to detect the difference against
the African genome. The config files to run bcbio are inside: `analysis/bcbio` folder. All commands used are at `bcbio-nextgen-commands.log`.

After this was done, we used the bash script `vcf_parsing.sh` to parse the final VCF in order to get all the information
showed in the final table.

 * VCF files were splited in multiple small files to run in parallel
  * get mutation affecting protein coding genes
  * get sequence from the specific genome to show the NT that changed with 200 flank regions 
  * all small files are merged together to have the full list of variants in one file again
 * script to get the African sequences for all the variants with flank regions (`get_african_sequence.py`)
 * script to get the genotype right, when there are multiple alleles, or missing information (`parse_vcf_genotpye.py`)
 * create table with R script `merge-tables2.R` merging the bcbio analysis with flank regions and genotype output from previous scripts
 * `all_genome_ann.sh` will annotate mutation impact with dbNSFP mapping these variants to human variants
 * clean tables with `clean-table2.R` script, producing the final output



Note:
Some scripts have hardcode full path, so this set of scripts are not designed to be automatically run from scratch.

