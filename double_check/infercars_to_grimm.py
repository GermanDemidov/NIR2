# -*- coding: utf-8 -*-
from collections import defaultdict
import os
import re
import sys
import itertools

__author__ = 'Sergey Aganezov Jr.'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'


def main(file=None, dest=None):
    genomes = defaultdict(lambda: defaultdict(list))
    with open(file, "r") as source:
        data = source.read()
        _, *synteny_blocks = re.split('>', data)
        result = []
        for block in synteny_blocks:
            sb_split_data = block.split('\n')
            sb_split_data = list(filter(lambda x: x.strip(), sb_split_data))
            synteny_id, *sb_split_data = sb_split_data
            result.append((synteny_id, sb_split_data))
        for sb_id, sb_raw_data in result:
            for sb_raw_line in sb_raw_data:
                genome, *rest = sb_raw_line.split(".")
                rest = "".join(rest)
                chromosome, *rest = rest.split(":")
                rest = "".join(rest)
                start_end, strand = rest.split(" ")
                start, end = map(int, start_end.split("-"))
                genomes[genome][chromosome].append((sb_id, int(start), int(end), strand))
        for genome in genomes:
            for chromosome in genomes[genome]:
                genomes[genome][chromosome] = sorted(genomes[genome][chromosome], key=lambda x: (x[1], x[2]))
        with open(dest, "w") as dest:
            for genome in genomes:
                print(">" + genome, file=dest)
                for chromosome in genomes[genome]:
                    print(" ".join(map(lambda x: x[0] if x[3] == "+" else "-{}".format(x[0]), genomes[genome][chromosome])),
                          "$", file=dest)
                print(file=dest)

if __name__ == "__main__":
    cmd_args = sys.argv[1:]
    source_files = cmd_args
    for file in source_files:
        dest = os.path.basename(file)
        dest = os.path.sep.join(["grimm", dest])
        main(file, dest)