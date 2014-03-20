# -*- coding: utf-8 -*-
import os
import sys
import re


__author__ = 'Sergey Aganezov Jr.'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'

def sb_to_adj(sb):
    if sb.startswith("-"):
        return "{sb}h,{sb}t".format(sb=sb[1:])  # slice from one to omit "-" sign
    else:
        return "{sb}t,{sb}h".format(sb=sb)

def main(file, destination):
    with open(file, "r") as source:
        with open(destination, "w") as dest:
            data = source.read()
            _, *genomes_raw_data = re.split('>', data)  # first element doesn't matter
            for genomes_raw_data in genomes_raw_data:
                genome_name, *chrs_raw_data = re.split("\\n", genomes_raw_data)
                chrs_raw_data = filter(lambda x: x.strip(), chrs_raw_data)  # get rid of empty lines
                print(">{gn}".format(gn=genome_name), file=dest)
                for chromosome_raw_data in chrs_raw_data:
                    chromosome_raw_data = chromosome_raw_data.strip()
                    *sb_data, _ = re.split(r" ", chromosome_raw_data)  # as last shall be "$" char
                    result = "inf "  # as we have linear chromosomes
                    result += " ".join(map(sb_to_adj, sb_data))
                    result += " inf"  # as we have linear chromosomes
                    print(result, file=dest)
                print(file=dest)

if __name__ == "__main__":
    cmd_args = sys.argv[1:]
    source_files = cmd_args
    for file in source_files:
        dest = os.path.basename(file)
        dest = os.path.sep.join(["adj", dest])
        main(file, dest)