# -*- coding: utf-8 -*-
from collections import defaultdict
import os
import sys
import itertools

__author__ = 'Sergey Aganezov Jr.'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'


def read_distances(source_file):
    pw_distances = defaultdict(dict)
    with open(source_file, "r") as source:
        header, *rows = list(map(lambda x: x.strip(), source.readlines()))
        column_genomes = header.split("\t")
        for row in rows:
            genome, *distances = row.split("\t")
            for value, pair_genome in zip(distances, column_genomes):
                pw_distances[genome][pair_genome] = int(value)
    return pw_distances


def read_breakpoints(pw_bp_file):
    pw_breakpoints = defaultdict(dict)
    with open(pw_bp_file, "r") as source:
        for line in source:
            line = line.strip()
            genomes, bp_count = line.split(":")
            bp_count = int(bp_count)
            genome_a, genome_b = genomes.split(",")
            pw_breakpoints[genome_a][genome_b] = bp_count
            pw_breakpoints[genome_b][genome_a] = bp_count
    return pw_breakpoints


def compute_bp_reuse_rate(pw_distance_file, pw_bp_file, destination_file):
    pw_distances = read_distances(pw_distance_file)
    pw_breakpoints = read_breakpoints(pw_bp_file)
    genomes = sorted(pw_distances.keys())
    with open(destination_file, "w") as dest:
        print(os.path.basename(pw_distance_file), file=dest)
        header = "{:^10}".format("") + "|".join("{:^10}".format(x) for x in genomes) + "|"
        print(header, file=dest)
        for genome in genomes:
            print("-"*len(header), file=dest)
            row = "{:^9}".format(genome)
            row += "|"
            row += "|".join(map(lambda x: "{:^10.4}".format(4 * pw_distances[genome][x] / pw_breakpoints[genome][x] if genome != x else ""),
                                                           genomes))
            row += "|"
            print(row, file=dest)

if __name__ == "__main__":
    cmd_args = sys.argv[1:]
    source_files = cmd_args
    assert len(source_files) % 2 == 0  # files have to be in even numbers, as we look for pairs (distance:break_point) of files
    for distance_file, break_point_file in zip(source_files[::2], source_files[1::2]):
        dest = os.path.basename(distance_file)
        dest = os.path.sep.join(["reuse", dest])
        compute_bp_reuse_rate(distance_file, break_point_file, dest)