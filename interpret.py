"""
Interpret molecular cancer disease data in light of the genetic code.

Not all interpretations are equally reasonable.
They are described below.

"""
from __future__ import division, print_function, absolute_import

import contextlib
import sys
import argparse
import csv

UNKNOWN = 'UNKNOWN'
BENIGN = 'BENIGN'
LETHAL = 'LETHAL'


def interpret_1(wildaa, mutantaa, dmin, dcounts, threshold):
    """
    First interpretation.

    At each codon position, each amino acid is considered either lethal
    or benign.  Every non-reference disease-associated amino acid at the
    position is considered lethal and every amino acid that is not
    disease-associated is considered benign.
    """
    if sum(dcounts) >= threshold:
        return LETHAL
    else:
        return BENIGN


def interpret_2(wildaa, mutantaa, dmin, dcounts, threshold):
    """
    Second interpretation.

    As in the previous interpretation, at each codon position each
    amino acid is considered either lethal or benign.  Synonymous
    mutations and mutations with more than one nucleotide change per codon
    are removed from the data set.  From the remaining data, every
    disease-associated amino acid at the position is considered lethal and
    every amino acid that is not disease-associated is considered benign.
    In particular, if for a given codon position the original data
    contains evidence of disease association for a given amino acid but
    the only such evidence is through >1 nucleotide change from the
    reference, then this amino acid will be interpreted as benign at that
    position.
    """
    if sum(dcounts[:2]) >= threshold:
        return LETHAL
    else:
        return BENIGN


def interpret_3(wildaa, mutantaa, dmin, dcounts, threshold):
    """
    Third interpretation.

    At each codon position, the disease state of each possible amino
    acid is considered either lethal or benign or unknown.  Every
    disease-associated amino acid is considered to be lethal.  Among the
    amino acids without evidence of disease association, the amino acids
    that are reachable by a single nucleotide point mutation from the
    reference codon are considered to be benign.  The remaining amino
    acids (those without evidence of disease-association and which are not
    reachable by a single nucleotide point mutation from the reference
    codon) are assumed to have an unknown (missing) disease state.
    """
    if sum(dcounts) >= threshold:
        return LETHAL
    elif dmin > 1:
        return UNKNOWN
    else:
        return BENIGN


def interpret_4(wildaa, mutantaa, dmin, dcounts, threshold):
    """
    Fourth interpretation.

    All amino acids except the reference amino acid
    are considered to be lethal in the reference process.
    """
    if wildaa == mutantaa:
        return BENIGN
    else:
        return LETHAL


@contextlib.contextmanager
def open_in(filename):
    if filename == '-':
        yield sys.stdin
    else:
        with open(filename, 'rU') as fin:
            yield fin


@contextlib.contextmanager
def open_out(filename):
    if filename == '-':
        yield sys.stdout
    else:
        with open(filename, 'w') as fout:
            yield fout


def submain(f_interp, threshold, fin, fout):
    # filter the tsv file
    header_out = 'position', 'residue', 'status'
    for i, row_in in enumerate(fin):
        if not i:
            fout_writer.write_row(header_out)
            continue
        pos, wildcodon, waa, maa, dmin, d0, d1, d2, d3 = row_in
        pos = int(pos)
        dmin = int(dmin)
        dcounts = np.array([int(x) for x in (d0, d1, d2, d3)])
        status = f_interp(waa, maa, dmin, dcounts, threshold)
        row_out = pos, maa, status
        fout.writerow(row_out)


def main(args):
    # define the interpretation function
    f_interp = {
            1 : interpret_1,
            2 : interpret_2,
            3 : interpret_3,
            4 : interpret_4,
            }[args.interpretation]
    
    # open the tsv reader and writer and call a submain function
    with open_in(args.infile) as fin:
        fin_reader = csv.reader(fin, delimiter='\t')
        with open_out(args.outfile) as fout:
            fout_writer = csv.writer(fout, delimiter='\t')
            return submain(f_interp, threshold, fin_reader, fout_writer)
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', default='-',
            help='reduced but not yet interpreted disease tsv file')
    parser.add_argument('-o', '--outfile', default='-',
            help='interpreted disease tsv file')
    parser.add_argument('--interpretation',
            type=int, default=1, choices=(1, 2, 3, 4),
            help='interpretation number (default=1)')
    parser.add_argument('--threshold', type=int, default=1,
            help='threshold number of disease observations (default=1)')
    main(parser.parse_args())

