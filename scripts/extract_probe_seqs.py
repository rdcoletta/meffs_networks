#!/usr/bin/python3

'''
created by Rafael Della Coletta
2019-09-19
'''

import argparse as ap
import glob
import natsort
from Bio import SeqIO

# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script gets the positions of each SNP of a SNPchip and
             extracts the probes sequences around that SNP from the reference
             genome used to build the chip. Each chromosome should be in
             separate fasta file.''')
# add positional arguments
parser.add_argument("SNPchip_data", type=str,
                    help="file with SNP chip data to extract SNP positions")
parser.add_argument("ref_gen", type=str,
                    help="relative path to folder containing ref gen fasta"
                    " files. This script requires one chromosome per fasta.")
parser.add_argument("probe_size", type=int,
                    help="number of bases to be extracted around the SNP")
parser.add_argument("output_name", type=str,
                    help="name of the output")
# pass arguments into variables
args = parser.parse_args()
SNPchip_data = args.SNPchip_data
ref_gen = args.ref_gen
probe_size = args.probe_size
output_name = args.output_name


# open SNPchip file
SNPchip = open(SNPchip_data, "r")

# create dictionary to store genomic coordinates of SNPs used in chip
# (chromosome number as keys and positions as values)
SNPs_coord = {}

# skip header
SNPchip.readline()

# retrieve SNP coordinates
for line in SNPchip:
    line = line.strip()
    line = line.split("\t")
    chr = line[2]
    pos = line[3]
    if chr not in SNPs_coord:
        # create a key-value pair, where the value is a list of positions
        SNPs_coord[chr] = [pos]
    else:
        # if chr is already in the dict, just append the position to the list
        SNPs_coord[chr].append(pos)

# close SNPchip file
SNPchip.close()

# create new file to write output
outfile = open(output_name, "w")

# for each chromosome (for each key in dict)...
for chr in SNPs_coord.keys():
    # get full name of all fasta files containing ref seq for a chromosome
    ref_gen_chr = glob.glob(ref_gen + "*.fa")
    # natural sort list (because name has letters and numbers)
    ref_gen_chr = natsort.natsorted(ref_gen_chr)
    # transform variable into string (glob output is a list), and subtract 1
    # because python is 0-index and first chromosome is number 1
    ref_gen_chr = str(ref_gen_chr[int(chr) - 1])
    print("Analyzing file '", ref_gen_chr, "'", sep="")
    # read reference genome chromosome sequence
    ref_seq = SeqIO.read(ref_gen_chr, "fasta")
    # store only sequence
    my_seq = ref_seq.seq
    # for each SNP in chromosome, get sequence around it according to probe size
    for SNP in SNPs_coord[chr]:
        # correct for python being 0-index and SNP positions 1-index
        SNP = int(SNP) - 1
        # calculate number of bases that will flank each SNP
        up_flank = round(probe_size / 2)
        down_flank = probe_size - up_flank
        # get probe sequence (and make sure it's a string)
        start_pos = SNP - up_flank
        end_pos = SNP + down_flank
        probe = str(my_seq[start_pos:end_pos])
        # print results to table:
        #   add 1 to start and SNP position to return to 1-index again, but
        #   end_pos is already in 1-index
        print(">probe-B73v2_chr-", chr, "_start-", start_pos + 1, "_SNP-",
              SNP + 1, "_end-", end_pos, "\n", probe, sep="", file=outfile)
