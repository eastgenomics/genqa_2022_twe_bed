# genqa_2022_twe_bed
Program to create a bed file for a TWE germline NGS submission in the 2022 GenQA. Running genqa.py generates the output file 'genqa_bed.bed', which describes 44930 genomic intervals.

## PART 1
Create a bed file of clinical transcript exons for all our panel genes

a. Get the following files from DNAnexus:
- The latest file in 001_reference > dynamic_files > gene_panels
- The latest file in 001_reference > dynamic_files > nirvana_genes2transcripts
- The latest genomic.symbols.exon file in 001_reference > annotation > b37

b. Extract all HGNC IDs from the gene_panels file
- Viewing the original .tsv in excel shows that it has 16297 rows;
    removing duplicates leaves 3122 unique HGNC IDs

c. For each HGNC ID from (b), identify the clinical transcript
- Each gene should have exactly 1 clinical transcript

d. For each gene/clinical transcript from (c), identify exon locations

e. Do some sanity tests on the output
- Do all our HGNC IDs have exactly 1 clinical transcript
- Do all their clinical transcripts have at least 1 exon

## PART 2
Create (using bedtools intersect) a bed file which is the intersect of:
- The bed file from part 1
- The TWE (whole-exome) capture bed

## Reference files

dx_200900_TWE.bed
- Original TWE whole-exome capture bed
- Identical to 001_reference > bed_files > b37 > kits > twist_exome > Twist_ComprehensiveExome_targets_hg19_noChr_2020-09.bed

dx_221027_panels2genes.tsv
- Describes all HGNC IDs for genes in all EGLH panels
- Identical to 001_reference > dynamic_files > gene_panels > 221027_genepanels.tsv

dx_221025_genes2transcripts.tsv
- Describes clinical and non-clinical EGLH transcripts for all HGNC IDs
- Identical to 001_reference > dynamic_files > nirvana_genes2transcripts > 221025_g2t.tsv

dx_220711_transcripts2exons.tsv
- Describes positions of exons for gene transcripts
- Identical to 001_reference > annotation > b37 > GCF_000001405.25_GRCh37.p13_genomic.symbols.exon_5bp_v2.0.0.tsv

hgnc_20221103_dump
- Downloaded from https://www.genenames.org/download/custom/ on 3 Nov 2022
- Contains HGNC ID, current gene symbol, previous symbols, and alias symbols
