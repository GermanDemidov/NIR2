# -*- coding: utf-8 -*-
import sys
import re
import os


__author__ = 'Sergey Aganezov Jr.'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'


def different_genome_check(sb_raw_data_list):
    return len(set(x.split(".")[0] for x in sb_raw_data_list)) == len(sb_raw_data_list)


def filtration(source_file, output_file):
    with open(source_file, "r") as source:
        data = source.read()
        _, *synteny_blocks = re.split('>', data)
        result = []
        for block in synteny_blocks:
            sb_split_data = block.split('\n')
            sb_split_data = list(filter(lambda x: x.strip(), sb_split_data))
            synteny_id, *sb_split_data = sb_split_data
            result.append((synteny_id, sb_split_data))
        result = list(filter(lambda x: len(x[1]) == 6, result))
        result = list(filter(lambda x: different_genome_check(x[1]), result))
        with open(output_file, "w") as destination:
            for sb_id, sb_list in result:
                print(">" + sb_id, file=destination)
                for sb in sb_list:
                    print(sb, file=destination)
                print(file=destination)

if __name__ == "__main__":
    cmd_args = sys.argv[1:]
    source_files = cmd_args
    for file in source_files:
        output_file = os.path.basename(file)
        output_file = os.path.sep.join(["filtered", output_file])
        filtration(file, output_file)