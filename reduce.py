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

from collections import defaultdict
import argparse
import csv


header_out = (
        'codonpos', 'wildcodon', 'wildaa', 'mutantaa', 'minrequired',
        'tumor0', 'tumor1', 'tumor2', 'tumor3')


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def get_codon_aa_distance_table(codon_aa_pairs):
    """
    Get a table mapping (codon, aa) to a distance between 0 and 3

    The distance is the minimum number of nucleotide sites within the
    input codon that must be changed to reach a codon that codes for the
    output amino acid.

    """
    codon_to_aa = dict(codon_aa_pairs)
    codons = set(codon_to_aa)
    codon_aa_distance = {}
    for ca in codons:
        for cb in codons:
            edge = ca, codon_to_aa[cb]
            d = hamming_distance(ca, cb)
            best = codon_aa_distance.get(edge, d)
            codon_aa_distance[edge] = min(best, d)
    return codon_aa_distance


def submain(human_codons, codon_aa_pairs, fin_reader, fout_writer):
    nsites = len(human_codons)
    codon_to_aa = dict(codon_aa_pairs)
    aa_set = set(aa for c, aa in codon_aa_pairs if aa != 'STOP')
    codon_set = set(c for c, aa in codon_aa_pairs if aa != 'STOP')
    codon_aa_dist_table = get_codon_aa_distance_table(codon_aa_pairs)
    # for each (codonpos, mutantaa, distance) count observations in tumors
    counts = defaultdict(int)
    codonpos_to_wildcodon = {}
    # args are tsv reader and writer
    header_in = None
    to_idx = None
    for i, row in enumerate(fin_reader):

        # extract the header from the first row
        if not i:
            header_in = row
            to_idx = dict((name, i) for i, name in enumerate(header_in))
            fout_writer.writerow(header_out)
            continue

        # extract useful entries from the data row
        codonpos_string = row[to_idx['Codon']]
        try:
            codonpos = int(codonpos_string)
        except ValueError as e:
            continue
        wildcodon = row[to_idx['WT codon']].upper()
        mutantcodon = row[to_idx['Mutant codon']].upper()
        wildaa = row[to_idx['WT AA']].upper()
        mutantaa = row[to_idx['Mutant AA']].upper()

        # skip data rows that do not represent nonsynonymous mutations
        if mutantaa in ('STOP', 'FS.'):
            continue
        if wildaa == mutantaa:
            continue

        # report weird codons
        for codon in wildcodon, mutantcodon:
            if codon not in codon_set:
                raise Exception('unrecognized codon:', codon)

        # report weird amino acids
        for aa in wildaa, mutantaa:
            if aa not in aa_set:
                raise Exception('unrecognized amino acid:', aa)

        # report weird codon translations
        for codon, aa in (mutantcodon, mutantaa), (wildcodon, wildaa):
            if codon_to_aa[codon] != aa:
                raise Exception('bad translation %s -> %s', (codon, aa))

        # check the wild type codon at this position vs the human reference
        if wildcodon != human_codons[codonpos-1]:
            raise Exception('wild codon does not match human sequence codon')

        # report out of bound sites
        if not (1 <= codonpos <= nsites):
            raise Exception('unrecognized codon site:', codonpos)


        # increment the count
        d = hamming_distance(wildcodon, mutantcodon)
        if d < codon_aa_dist_table[wildcodon, mutantaa]:
            raise Exception('impossibly small number of nt changes')
        counts[(codonpos, mutantaa, d)] += 1

    # write the output data rows
    for codonpos in range(1, nsites+1):
        wildcodon = human_codons[codonpos-1]
        wildaa = codon_to_aa[wildcodon]
        for mutantaa in sorted(aa_set):
            tumor = {}
            for i in range(4):
                tumor[i] = counts.get((codonpos, mutantaa, i), 0)
            row = (codonpos, wildcodon, wildaa, mutantaa,
                    codon_aa_dist_table[wildcodon, mutantaa],
                    tumor[0], tumor[1], tumor[2], tumor[3])
            fout_writer.writerow(row)


def main(args):
    # read the genetic code
    codon_aa_pairs = []
    with open(args.code) as fin:
        for line in fin:
            state, aa, codon = line.split()
            codon_aa_pairs.append((codon.upper(), aa.upper()))

    # read the human reference sequence
    with open(args.ref) as fin:
        human_codons = fin.read().split()

    # check the sequence length
    if len(human_codons) != args.nsites:
        raise Exception('sequence length discrepancy')

    # read the input tsv data file in universal newline mode
    with open(args.filename_in, 'rU') as fin:
        fin_reader = csv.reader(fin, delimiter='\t')
        with open(args.filename_out, 'w') as fout:
            fout_writer = csv.writer(fout, delimiter='\t')
            return submain(human_codons, codon_aa_pairs,
                    fin_reader, fout_writer)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--filename-in', required=True,
            help='name of tab separated data input file')
    parser.add_argument('--filename-out', required=True,
            help='name of tab separated data output file')
    parser.add_argument('--ref', required=True,
            help='human p53 reference codon sequence file name')
    parser.add_argument('--code', required=True,
            help='universal genetic code filename')
    parser.add_argument('--nsites', type=int, default=393,
            help='number of aligned codon sites')
    args = parser.parse_args()
    main(args)

