# -*- coding: utf-8 -*-
from collections import defaultdict
import os
import re
import sys
import itertools

__author__ = 'Sergey Aganezov Jr.'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'


def count_bp(adj_list_1, adj_list_2):
    result = 0
    set_1, set_2 = set(adj_list_1), set(adj_list_2)
    for adj_l, adj_r in set_1:
        if (adj_l, adj_r) not in set_2 or (adj_r, adj_l) not in set_2:
            result += 1
    for adj_l, adj_r in set_2:
        if (adj_l, adj_r) not in set_1 or (adj_r, adj_l) not in set_1:
            result += 1
    return result


def main(source_file, destination=None):
    adj_dict = defaultdict(list)
    with open(source_file, "r") as source:
        with open(destination, "w") as dest:
            data = source.read()
            _, *genomes_raw_data = re.split('>', data)  # first element doesn't matter
            for genomes_raw_data in genomes_raw_data:
                genome_name, *chrs_raw_data = re.split("\\n", genomes_raw_data)
                chrs_raw_data = filter(lambda x: x.strip(), chrs_raw_data)  # get rid of empty lines
                for chromosome_raw_data in chrs_raw_data:
                    chromosome_raw_data = chromosome_raw_data.strip()
                    adj_data = re.split(r",", chromosome_raw_data)
                    adj_data = list(map(lambda x: tuple(x.split()), adj_data))
                    adj_dict[genome_name].extend(adj_data)
            for genome_a, genome_b in itertools.combinations(adj_dict, 2):
                print("{g1},{g2}:{bpc}".format(g1=genome_a,
                                               g2=genome_b,
                                               bpc=count_bp(adj_dict[genome_a], adj_dict[genome_b])), file=dest)


if __name__ == "__main__":
    cmd_args = sys.argv[1:]
    source_files = cmd_args
    for file in source_files:
        dest = os.path.basename(file)
        dest = os.path.sep.join(["bp", dest])
        main(file, dest)