"""
Reduce the p53 data.

This reduction requires the genetic code to compute the column
corresponding to a minimum number of possible changes from a given codon
to a given amino acid.

Make a tsv file with columns:
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
