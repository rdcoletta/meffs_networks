#!/usr/bin/python3

'''
created by Rafael Della Coletta
2019-09-24
'''


import argparse as ap
import pandas as pd


# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script first creates a file with SNP positions of refgen B73v2
             and their respective positions in v4, and then converts a hapmap
             file(s) with v2 coordinates into a hapmap with v4 coordinates. It
             also removes reads that to a different chromosome in v4.

note: the name of the reads used in for alignment shoud have the chromosome
      number and SNP position''')
# add positional arguments
parser.add_argument("alignment", type=str,
                    help="alignment file in SAM format")
parser.add_argument("output_SNP_pos", type=str,
                    help="path to where you want to save the coordinates of a "
                    "SNP in v2 and v4 ref gen")
parser.add_argument("hmp_files", type=str,
                    help="comma-separated list of hapmap files with v2 "
                    "coordinates to be converted to v4")
# pass arguments into variables
args = parser.parse_args()
alignment = args.alignment
output_SNP_pos = args.output_SNP_pos
hmp_files = args.hmp_files
hmp_files = hmp_files.split(",")


# create a dictionary with the chromosome IDs used in v4, so I can convert
# them to numerical form (easier to work with)
chr_IDs = {"CM007647.1": 1, "CM007648.1": 2, "CM007649.1": 3, "CM000780.4": 4,
           "CM000781.4": 5, "CM000782.4": 6, "CM007650.1": 7, "CM000784.4": 8,
           "CM000785.4": 9, "CM000786.4": 10}


# open sam file
sam_file = open(alignment, "r")
# open output file for SNP positions
outfile = open(output_SNP_pos, "w")
outfile_discard = output_SNP_pos.rstrip("txt")
outfile_discard = open(outfile_discard + "discarded.txt", "w")
# write header of table with v2 and v4 SNP positions
print("chr_v2", "pos_v2", "chr_v4", "pos_v4", "SNP", sep="\t", file=outfile)

# create a dictionary relating v2 to v4 coordinates
# (this will be useful to modify hapmap files)
v2_v4_coord = {}

# keep track of number of discarded reads due to mapping into diff chromosome
count_discard = 0

# parse sam file to get v4 positions for a particular SNP
for line in sam_file:
    line = line.strip()
    # if line starts with @ (header) -- skip
    if line[0] == "@":
        continue
    # otherwise, extract info from line
    else:
        line = line.split("\t")
        # only proceed if read is mapped (i.e., FLAG is not 4)
        flag = int(line[1])
        if flag != 4:
            # get v2 info from SAM file
            probe_ID = line[0].split("_")
            chr_v2 = probe_ID[1].split("-")[1]
            start_v2 = probe_ID[2].split("-")[1]
            SNP_v2 = probe_ID[3].split("-")[1]
            end_v2 = probe_ID[4].split("-")[1]
            # get v4 info from SAM file
            chr_v4 = chr_IDs[line[2]]
            start_v4 = line[3]
            probe_seq = line[9]
            # don't use probes that mapped in different chromosomes
            if str(chr_v2) == str(chr_v4):
                # get SNP index in seq -- it's rev comp if SAM flag is 16
                if flag == 16:
                    # seq in sam file is already the reverse complement
                    # so need to use reverse coordinates from v2
                    start_to_SNP = int(end_v2) - int(SNP_v2)
                else:
                    start_to_SNP = int(SNP_v2) - int(start_v2)
                # get SNP v4 position based on start and SNP v2 positions
                SNP_v4 = int(start_v4) + start_to_SNP
                SNP_nt = probe_seq[start_to_SNP]
                # write results
                print(chr_v2, SNP_v2, chr_v4, SNP_v4, SNP_nt,
                      sep="\t", file=outfile)
                # write v2 coord as dictionary keys and v4 as values
                dict_key = str(chr_v2) + "," + str(SNP_v2)
                dict_val = str(chr_v4) + "," + str(SNP_v4)
                v2_v4_coord[dict_key] = dict_val
            else:
                # keep track of probes that mapped to different chromosome
                print("\t".join(line), file=outfile_discard)
                count_discard += 1

print(count_discard, "reads discarded due to mapping into different chromosome"
      " in v4")

# close files
sam_file.close()
outfile.close()
outfile_discard.close()


# parse hapmap files
for hmp in hmp_files:
    # open (and close automatically) one hapmap file
    with open(hmp, "r") as hmp_file:
        # open new output file
        hmp_out = hmp.rstrip(".hmp.txt")
        hmp_out = open(hmp_out + ".v4.hmp.txt", "w")
        # print header
        header = hmp_file.readline()
        header = header.strip()
        print(header, file=hmp_out)
        # parse file line by line
        for line in hmp_file:
            line = line.strip()
            line = line.split()
            # get coordinates in hapmap file
            hmp_chr = line[2]
            hmp_pos = line[3]
            # only print lines that has v2 coordinates on dictionary
            if hmp_chr + "," + hmp_pos in v2_v4_coord.keys():
                v4_coord = v2_v4_coord[hmp_chr + "," + hmp_pos]
                chr_v4 = v4_coord.split(",")[0]
                pos_v4 = v4_coord.split(",")[1]
                line[2] = chr_v4
                line[3] = pos_v4
                print("\t".join(line), file=hmp_out)
        # close output
        hmp_out.close()


# after everything is done, sort the data frame by chr and pos
for hmp in hmp_files:
    hmp = hmp.rstrip(".hmp.txt")
    hmp = hmp + ".v4.hmp.txt"
    # open file as data frame
    hmp_df = pd.read_table(hmp, sep="\t", keep_default_na=False)
    # sort by chromosome and then by position
    hmp_sorted = hmp_df.sort_values(["chrom", "pos"])
    # write sorted hapmap
    hmp_sorted.to_csv(hmp, sep="\t", index=False)
