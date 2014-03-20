# -*- coding: utf-8 -*-

__author__ = 'Sergey Aganezov Jr.'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'


import sys


def rename_sb_ids(source_file):
    cnt = 1
    with open(source_file, "r") as source:
        data = source.readlines()
        with open(source_file, "w") as dest:
            for line in data:
                if line.startswith(">"):
                    print(">" + str(cnt), file=dest)
                    cnt += 1
                else:
                    print(line, end="", file=dest)

if __name__ == "__main__":
    cmd_args = sys.argv[1:]
    source_files = cmd_args
    for file in source_files:
        rename_sb_ids(file)