# -*- coding: utf-8 -*-

__author__ = 'Sergey Aganezov Jr.'
__email__ = 'aganezov@gwu.edu'
__status__ = 'develop'

import sys


def main(source_file):
    result = set()
    with open(source_file, "r") as source:
        for line in source:
            if line.startswith(">"):
                result.add(line.strip()[1:])
    return len(result)

if __name__ == "__main__":
    cmd_args = sys.argv[1:]
    source_files = cmd_args
    for file in source_files:
        print(main(file))