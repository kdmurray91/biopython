"""Demonstration use-case of the Bio.Webtools.TAIR module.
See comments in source code for pointers on the use of this module.
"""

from Bio.Webtools import TAIR
from Bio import SeqIO
from optparse import OptionParser


# Get the commandline options, this section is unrelated to the function of the
# TAIR module.
parser = OptionParser()
parser.add_option('-f', '--file', dest='filename', default=None)
parser.add_option('-m', '--mode', dest='mode', default='direct',
        help="'direct', 'ncbi_rna' or 'ncbi_protein'")
parser.add_option('-d', '--dataset', dest='dataset', default="transcript",
        help="one of: upstream_500, downstream_3000, intergenic, 5prime_utr, \
        upstream_1000, intron, downstream_500, cds, 3prime_utr, genomic, \
        protein, gene, transcript, upstream_3000, downstream_1000")
parser.add_option('-a', '--agis', dest='agis',
        default="AT5G63980.1,AT5G63280.1,AT5G23980.1",
        help='Comma seperated list of AGIs to fetch')
(options, args) = parser.parse_args()


# The TAIR module get functions require a AGIs to be given as a python list
# of strings, containing AGIs. This gets a list of strings from the CSV
# string given.
agis = options.agis.split(",")

# Get the sequences
if options.mode == "direct":
    seqs = TAIR.get(agis, options.dataset, "representative")
elif options.mode == "ncbi_protein":
    # Note this is how you specfiy the NCBI node
    seqs = TAIR.get_protein_from_ncbi(agis)
elif options.mode == "ncbi_rna":
    # Note this is how you specfiy the NCBI node
    seqs = TAIR.get_rna_from_ncbi(agis)
else:
    raise ValueError("Invalid mode: %s, see --help\n" % options.mode)


if options.filename is None:
    # If a filename is not specified, print a fasta to stdout
    for seq in seqs:
        print seq.format("fasta")
else:
    # If we have a filename, write a fasta to the file specified
    file_handle = open(options.filename, "wb")
    SeqIO.write(seqs, file_handle, "fasta")
    file_handle.close()
