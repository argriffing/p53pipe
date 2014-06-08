"""
Reduce the p53 data.

This reduction requires the genetic code to compute the column
corresponding to a minimum number of possible changes from a given codon
to a given amino acid.

Make a file of tab separated values (tsv) with columns:
 - codon position starting at 1 and ending at 393
 - wild type human codon
 - corresponding wild type amino acid
 - query amino acid
 - min possible changes from the wild type human codon to the query amino acid
 - number of disease observations with codon distance 0 from wild type
 - number of disease observations with codon distance 1 from wild type
 - number of disease observations with codon distance 2 from wild type
 - number of disease observations with codon distance 3 from wild type

"""
from __future__ import division, print_function, absolute_import

import argparse
import csv

"""
 
 File # Mutation position   Exon    Codon   WT codon    Mutant codon    WT AA   Mutant AA   Protein variant cDNA variant    Name    ATCC    Variation type  Event   Type    Complexity

"""

header_out = (
        'codonpos', 'wtcodon', 'wtaa', 'mtaa', 'minrequired',
        'n0', 'n1', 'n2', 'n3')

important_cols_of_input_header = (
        'Codon', 'WT codon', 'Mutant codon', 'WT AA', 'Mutant AA')


def submain(fin_reader, fout_writer):
    # args are tsv reader and writer
    header_in = None
    for i, row in enumerate(fin_reader):

        # extract the header from the first row
        if not i:
            header_in = row
            continue

        # extract relevant cols
        print(row[:10])
        if i > 10:
            return


def main(args):
    # read in universal newline mode
    with open(args.filename_in, 'rU') as fin:
        fin_reader = csv.reader(fin, delimiter='\t')
        with open(args.filename_out, 'w') as fout:
            fout_writer = csv.writer(fout, delimiter='\t')
            return submain(fin_reader, fout_writer)


if __name__ == '__main__':
    raise NotImplementedError
    parser = argparse.ArgumentParser()
    parser.add_argument('--filename-in', required=True,
            help='name of tab separated data input file')
    parser.add_argument('--filename-out', required=True,
            help='name of tab separated data output file')
    args = parser.parse_args()
    main(args)

