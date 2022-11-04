#!usr/bin/env python

"""
Author: Jay Miles
Date: 4 Nov 2022

Functions for dealing with data from a tab-delimited text file
('filename') containing a dump of the HGNC database.

Note that depending on when the dump file was created, the columns
containing current, alias and previous gene symbols may have slightly
different titles, and the script may need to be modified accordingly.
"""


import pandas as pd


def import_hgnc_dump(filename):
    """ Read in a tab-separated dump of the HGNC database and create a
    pandas dataframe from it.

    args:
        filename [path]: containing dump of HGNC database

    returns:
        hgnc_df [pandas dataframe]
    """

    hgnc_df = pd.read_csv(filename, sep='\t')

    hgnc_df = rename_columns(hgnc_df)

    return hgnc_df


def rename_columns(hgnc_df):
    """ Renames dataframe columns to standardised values.

    Args:
        df [pandas dataframe]

    Returns:
        df [pandas dataframe]: columns have been renamed
    """

    renamed_columns = [
        'hgnc_id',
        'symbol',
        'prev_symbol',
        'alias_symbol']

    hgnc_df.columns = renamed_columns

    return hgnc_df


def get_hgnc_from_symbol(hgnc_df, gene_symbol):
    """ Get the HGNC ID for a supplied gene symbol, if one exists. A
    gene symbol may appear in the 'symbol', 'alias_symbol' or
    'prev_symbol' columns depending on whether it is current or not, so
    the function looks through all three in turn.

    args:
        hgnc_df: pandas df containing dump of HGNC site
        gene_symbol [str]: query gene symbol

    returns:
        hgnc_id [str], or None: HGNC ID of query gene
    """

    hgnc_id = None

    # if a row exists where this gene is the official gene symbol,
    # get that row's index and hgnc id

    try:

        target_index = hgnc_df.index[hgnc_df['symbol'] == gene_symbol]

        hgnc_id = hgnc_df.loc[target_index[0], 'hgnc_id']

    except IndexError:

        # or see if it's in the 'previous symbols' field

        try:
            i = 0

            for row in hgnc_df['prev_symbol']:

                if type(row) == str:
                    symbols = [ele.strip() for ele in row.split(',')]

                    for symbol in symbols:
                        if symbol == gene_symbol:

                            hgnc_id = hgnc_df.iloc[i].loc['hgnc_id']
                            break

                i += 1

        # or see if it's in the 'alias symbols' field

        except IndexError:

            j = 0

            for row in hgnc_df['alias_symbol']:

                if type(row) == str:
                    symbols = [ele.strip() for ele in row.split(',')]

                    for symbol in symbols:
                        if symbol == gene_symbol:

                            hgnc_id = hgnc_df.iloc[j].loc['hgnc_id']
                            break

                j += 1

    return hgnc_id


def get_symbol_from_hgnc(hgnc_df, hgnc_id):
    """ Get the current gene symbol associated with a supplied HGNC id,
    if one exists.

    args:
        hgnc_df: pandas df containing dump of HGNC site
        hgnc_id [str]: query HGNC id

    returns:
        gene_symbol [str], or None: current symbol of query gene
    """

    try:
        # if a row exists with this HGNC id, get the symbol for that row

        target_index = hgnc_df.index[hgnc_df['hgnc_id'] == hgnc_id]

        gene_symbol = hgnc_df.loc[target_index[0], 'symbol']

    except IndexError:

        gene_symbol = None

    return gene_symbol
