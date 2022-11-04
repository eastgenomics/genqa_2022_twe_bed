#!usr/bin/env python

"""
Author: Jay Miles
Date: 4 Nov 2022

Create a bed file of clinical transcript exons for all our panel genes,
then intersect this with the TWE (whole-exome) capture bed to create a
bed file for submission to GenQA 2022.
"""


import hgnc
import pandas as pd
import subprocess


def get_panel_genes(fp):
    """ Identify all unique HGNC IDs from the gene_panels file.

    args:
        fp [str]: map of panels to associated genes

    returns:
        genes [list]: list of HGNC IDs [str]
    """

    # read in the col of hgnc ids as a series

    with open(fp, 'r') as reader:

        series = pd.read_csv(
            reader, sep='\t', header=None, usecols=[2]).squeeze('columns')

    # get unique values, convert to list, strip whitespace

    hgnc_series = series.unique()
    hgnc_list = hgnc_series.tolist()
    genes = [hgnc_id.strip() for hgnc_id in hgnc_list if hgnc_id]

    # confirm list is the right length

    assert len(genes) == 3122, \
        f'{len(genes)} unique non-empty genes identified, should be 3122'

    return genes


def get_gene_transcripts(fp, genes):
    """ Identify the clinical transcript for each gene in a list of
    unique HGNC IDs.

    args:
        fp [str]: map of HGNC IDs to associated transcripts
        genes [list]: list of HGNC IDs [str]

    returns:
        gene_transcripts [dict]: in the form {hgnc_id: transcript, ...}
    """

    # read in genes2transcripts as df

    gene_transcripts = {}

    with open(fp, 'r') as reader:
        df = pd.read_csv(reader, sep='\t', header=None)

    # get clinical transcript for all our hgnc ids

    for index, row in df.iterrows():
        if (row[0] in genes) and (row[2] == 'clinical_transcript'):

            gene_transcripts[row[0].strip()] = row[1].strip()

    # confirm right numbers of unique hgnc ids and transcripts

    keys = []
    values = []

    for key, value in gene_transcripts.items():
        if key and (key not in keys):
            keys.append(key)
        if value and (value not in values):
            values.append(value)

    assert len(keys) == 3122, \
        f'{len(keys)} unique non-empty HGNC IDs instead of 3122'

    assert len(values) == 3122, \
        f'{len(values)} unique non-empty clinical transcripts instead of 3122'

    return gene_transcripts


def create_eglh_bed(fp, transcripts, hgnc_df):
    """ Get the genomic positions for each exon in the specified
    clinical transcripts, and use these to create a bed file.

    args:
        fp [str]: containing positions of each transcript's exons
        transcripts [dict]: in the form {hgnc_id: transcript, ...}
        hgnc_df [pandas df]: maps HGNC IDs to gene symbols
    """

    # read in genes2exons file as df, initialise new df fields

    with open(fp, 'r') as reader:
        df = pd.read_csv(reader, sep='\t', header=None)

    chrom_field = []
    start_field = []
    end_field = []
    exon_field = []
    transcript_field = []
    symbol_field = []
    hgnc_field = []

    # find rows with our transcripts in

    for index, row in df.iterrows():
        if row[4] in transcripts.values():

            # get hgnc for row's transcript using hgnc-transcript dict

            hgnc_id = list(transcripts.keys())\
                [list(transcripts.values()).index(row[4])]

            # get hgnc for row's gene symbol using hgnc dump file

            symbol_hgnc = hgnc.get_hgnc_from_symbol(hgnc_df, row[3])

            # confirm the two hgnc ids match

            assert symbol_hgnc == hgnc_id, \
                f'g2e line {index+1}: ' \
                f'row transcript {row[4]} corresponds to {hgnc_id}, but ' \
                f'row gene {row[3]} corresponds to {symbol_hgnc}.'

            # append row data to new df fields

            chrom_field.append(row[0])
            start_field.append(row[1])
            end_field.append(row[2])
            symbol_field.append(row[3])
            transcript_field.append(row[4])
            exon_field.append(row[5])
            hgnc_field.append(hgnc_id)

    # confirm all 3122 hgnc ids have at least one exon

    for hgnc_id in transcripts.keys():
        assert hgnc_id in hgnc_field, f'{hgnc_id} has no exons'

    # create df for data verification
    # (also asserts that all fields are the same length)

    test_df = pd.DataFrame(
        list(zip(chrom_field, start_field, end_field,
        symbol_field, transcript_field, exon_field, hgnc_field)))

    test_df.columns = [
        'chrom', 'start', 'end', 'symbol', 'transcript', 'exon', 'hgnc']

    # confirm hgnc, symbol and transcript fields match in every row
    # (can't use get_symbol_from_hgnc, bc g2e uses some previous/alias symbols)

    exon_check = {}

    for index, row in test_df.iterrows():

        row_hgnc = row[6]
        row_symbol = row[3]
        rows_trans = row[4]

        # get the hgnc id of that row's gene symbol
        test_hgnc = hgnc.get_hgnc_from_symbol(hgnc_df, row_symbol)

        # confirm the hgnc id of that row is the same
        assert test_hgnc == row_hgnc, \
            f'{row_symbol} fetches {test_hgnc}, not {row_hgnc}'

        # confirm the transcript of that row is correct for the hgnc id
        assert transcripts[row_hgnc] == rows_trans, \
            f'{row_hgnc} should have {transcripts[row_hgnc]}, not {rows_trans}'

        # confirm no hgnc ids have the same exon number in multiple rows
        if row_hgnc not in exon_check.keys():
            exon_check[row_hgnc] = [row[5]]

        else:
            assert row[5] not in exon_check[row_hgnc], \
                f'{row_hgnc} has multiple exon {row[5]}'

            exon_check[row_hgnc].append(row[5])

    # create bed file df and write out

    sorted_df = test_df.sort_values(['chrom', 'start'])
    bed_df = sorted_df[['chrom', 'start', 'end']]

    bed_df.to_csv('eglh_bed.bed', sep='\t', header=False, index=False)


def intersect_beds(twe_bed, eglh_bed, fp):
    """ Find the intersect between the TWE capture bed file, and the bed
    file of all exons for all clinical transcripts of all our panel genes.

    args:
        twe_bed [str]
        eglh_bed [str]
        fp [str]: path to output bed file

    returns:
        intersect_bed [str]
    """

    output = subprocess.run([
        'bedtools', 'intersect',
        '-a', twe_bed,
        '-b', eglh_bed],
        capture_output=True,
        text=True)

    with open(fp, 'w') as writer:
        writer.write(output.stdout)


def main():
    """ reference files """

    p2g = 'dx_221027_panels2genes.tsv'
    g2t = 'dx_221025_genes2transcripts.tsv'
    t2e = 'dx_220711_transcripts2exons.tsv'

    twe_bed = 'dx_200900_TWE.bed'
    eglh_bed = 'eglh_bed.bed'
    final_bed = 'genqa_bed.bed'

    hgnc_dump = 'hgnc_20221103_dump.tsv'
    hgnc_df = hgnc.import_hgnc_dump(hgnc_dump)

    """ function calls """

    gene_list = get_panel_genes(p2g)
    transcript_dict = get_gene_transcripts(g2t, gene_list)
    create_eglh_bed(t2e, transcript_dict, hgnc_df)
    intersect_beds(twe_bed, eglh_bed, final_bed)


if __name__ == '__main__':
    main()
