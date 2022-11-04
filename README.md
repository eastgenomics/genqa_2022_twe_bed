# genqa_2022_twe_bed

PART 1
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

PART 2
Create (using bedtools intersect) a bed file which is the intersect of:
- The bed file from part 1
- The TWE (whole-exome) capture bed
