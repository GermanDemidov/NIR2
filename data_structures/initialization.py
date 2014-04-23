#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import sys
import re
from collections import defaultdict
from optparse import OptionParser
import optparse
from sets import Set
import copy
from data_structures.manager import Manager

def readBlocks( f , genomes ):
    # generator function
    # in: file
    # out: syntheny blocks in form of dictionary
    #      { ">name" : [[genome, chromosome, start, end, sign],...] }
    name = ""
    while True:
        line = name or f.readline()
        if not line:
            break
        blocks = []
        while True:
            name = f.readline()
            name.strip()
            if name == "":
                break
            if not name or name.startswith(">"):
                break
            elif name != "" and name.split() != []:
                name = name.replace(":"," ")
                name = name.replace("."," ")
                sign = name[-2:-1]
                name = name.replace("-"," ")
                name = name.replace("+"," ")
                block = name.split()
                block.append(sign)

                if len(block) > 3:
                    block[2] = int(block[2])
                    block[3] = int(block[3])
                blocks.append(block)

            for element in blocks:
                m = re.search("([^\.]*)", element[0])
                genome = m.group(0)
                if genome not in genomes:
                    genomes.append(genome)
        yield (line.split(), blocks)



def clusterization(synthenyBlocks, genomes, withoutDupDel=True):
    # in: dict of syntheny blocks and list of genomes
    # out: cluster of syntheny blocks without deletions and duplications
    #      [ names ]
    noDuplications = []
    noDeletions = []
    noDupDel = []
    
    for name, synBlock in synthenyBlocks.items():
        counterOfOccur = {}
        for genome in genomes:
            counterOfOccur[genome] = 0
        for location in synBlock:
            # organism means genome
            organism = location[0]
            for genome in genomes:
                if organism[:len(genome)] == genome:
                    counterOfOccur[genome] += 1
                    
        noDupl = True
        noDel = True
        noDD = True
        for value in counterOfOccur.itervalues():
            if value > 1:
                noDupl = False
            if value == 0:
                noDel = False
            if value == 0 or value > 1:
                noDD = False
        if noDupl:
            noDuplications.append(name)
        if noDel:
            noDeletions.append(name)
        if noDD:
            noDupDel.append(name)
    if withoutDupDel == True:
        return sorted(noDupDel)
    else:
        return sorted(noDel)


def clusterByChr(genomes, synthenyBlocks):
    # in: genomes and dict of syntheny blocks
    # out: dict of chromosomes with syntheny blocks in each genome
    #      chromosomes[ genome ] = {chromosome : [[ name, start, end, sign ],...]}
    chromosomes = {}
    for genome in genomes:
        chromosomes[genome] = defaultdict(list)
        
    noDuplAndDelBeforeMap = clusterization(synthenyBlocks, genomes)
    mapping = {}
    counter = 1
    for name in noDuplAndDelBeforeMap:
        mapping[name] = counter
        counter += 1
    for name in noDuplAndDelBeforeMap:
        infoList = synthenyBlocks[name]
        for element in infoList:
            chromosomes[element[0]][element[1]].append((
                mapping[name], element[2], element[3], element[4]))
            
    return chromosomes


def main(filename, fileToMerge):
    synthenyBlocks = {}
    synthenyBlocksToMerge = {}
    genomes = []
    
    lowResolutionName = ''.join(c for c in filename if c.isdigit()) + "K"
    highResolutionName = ''.join(c for c in fileToMerge if c.isdigit()) + "K"
    print lowResolutionName
    print highResolutionName


    with open(filename) as f:
        for line, blocks in readBlocks(f, genomes):
            synthenyBlocks[line[0]] = blocks
    with open(fileToMerge) as f:
        for line, blocks in readBlocks(f, genomes):
            synthenyBlocksToMerge[line[0]] = blocks
    chrLow = clusterByChr(genomes, synthenyBlocks)
    chrHigh = clusterByChr(genomes, synthenyBlocksToMerge)


if __name__ == "__main__":

    parser = OptionParser(version = "%prog 1.0", conflict_handler="error", description =
                            """This script merges one file with high-resolution syntheny
                            blocks to another file with low-resolution syntheny blocks.""")

    parser.add_option("-f", "--file", dest="filename", help="file with syntheny blocks of low resolution", metavar="FILE")
    parser.add_option("-m", "--merge", dest="fileToMerge", help="file with syntheny blocks of high resolution", metavar="FILE")

    options, args = parser.parse_args()
    
    main(options.filename, options.fileToMerge)
