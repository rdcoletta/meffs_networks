#!/usr/bin/python3

'''
created by Rafael Della Coletta
2019-07-02
'''

import argparse as ap

# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script reads in a table with two columns containing the
             genotype names and source IDs of RILs genotyped for the USDA
             project, and outputs a new table containing the names of lines
             that make up a biparental cross.''')
# add positional arguments
parser.add_argument("id_table", type=str,
                    help="table with genotype names and source IDs")
parser.add_argument("output_name", type=str,
                    help="name of the output table")
# pass arguments into variables
args = parser.parse_args()
id_table = args.id_table
output_name = args.output_name


# open input file
infile = open(id_table, "r")

# skip header
infile.readline()

# create dictionary to store biparental crosses as keys and RIL IDs as values
biparental_crosses = {}

# parse file
for line in infile:
    line = line.strip()
    line = line.split("\t")
    cross = line[0].split("-")[0]
    ril_name = line[0]
    if cross not in biparental_crosses:
        # introduce new key-value pair (where value is an element from a list)
        biparental_crosses[cross] = [ril_name]
    else:
        # update the correspondent list by adding more values
        biparental_crosses[cross].append(ril_name)

# remove duplicated ril names
for cross in biparental_crosses.keys():
    biparental_crosses[cross] = set(biparental_crosses[cross])

# open output file
outfile = open(output_name, "w")

# print header
print("cross", "RILs", sep="\t", file=outfile)

# write table
for key, val in biparental_crosses.items():
    print(key, ",".join(val), sep="\t", file=outfile)


# close files
infile.close()
outfile.close()
