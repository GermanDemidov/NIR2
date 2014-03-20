#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import sys
import re
from collections import defaultdict



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
                blocks.append(block)
            for element in blocks:
                m = re.search("([^\.]*)", element[0])
                genome = m.group(0)
                if genome not in genomes:
                    genomes.append(genome)
        yield (line.split(), blocks)

def sortWithKey(lst):
    lst.sort(key=lambda x: x[0])

def merge(synthenyBlocks, synthenyBlocksToMerge, genomes):
    '''
    in: synthenyBlocks and synthenyBlocksToMerge in format:
        { name of synBlock : [[genome1, chr, start, finish, sign], ...] }
    out: statistics in table format: num of synblocks + num of merged synblocks
        and merged synBlocks in file
    '''
    tuplesFromSynBlocks = {}
    for genome in genomes:
        tuplesFromSynBlocks[genome] = defaultdict(list)
    for names in synthenyBlocks:
        for lists in synthenyBlocks[names]:
            tuplesFromSynBlocks[lists[0]][lists[1]].append((int(lists[2]), int(lists[3]), lists[4]))

    tuplesFromSynBlocksToMerge = {}
    for genome in genomes:
        tuplesFromSynBlocksToMerge[genome] = defaultdict(list)
    for names in synthenyBlocksToMerge:
        for lists in synthenyBlocksToMerge[names]:
            tuplesFromSynBlocksToMerge[lists[0]][lists[1]].append((int(lists[2]), int(lists[3]), lists[4]))

    for genome in genomes:
        for chromosome in tuplesFromSynBlocks[genome]:
            sortWithKey(tuplesFromSynBlocks[genome][chromosome])
            
    for genome in genomes:
        for chromosome in tuplesFromSynBlocksToMerge[genome]:
            sortWithKey(tuplesFromSynBlocksToMerge[genome][chromosome])
    # in merged synBlocks we want to write chromosome cart (to convert it in
    # GRIMM format)
    mergedSynBlocks = {}
    counterOfAddedSynBlocks = {}
    for genome in genomes:
        mergedSynBlocks[genome] = defaultdict(list)
        counterOfAddedSynBlocks[genome] = 0

    for genome in genomes:
        print genome
        for chromosome in tuplesFromSynBlocks[genome]:
            print chromosome
            print tuplesFromSynBlocks[genome][chromosome]
            print tuplesFromSynBlocksToMerge[genome][chromosome]
            it1 = iter(tuplesFromSynBlocks[genome][chromosome])

            it2 = iter(tuplesFromSynBlocksToMerge[genome][chromosome])

            counter = 0
            mergedSynBlocks[genome][chromosome] = list(mergeIters(it1, it2))
            
            counterOfAddedSynBlocks[genome] += len(mergedSynBlocks[genome][chromosome]) - len(tuplesFromSynBlocks[genome][chromosome])

            prevEnd = 0
            for elem in mergedSynBlocks[genome][chromosome]:
                if elem[0] < prevEnd:
                    print
                    print elem, "WOW"
                    print
                prevEnd = elem[1]
            print mergedSynBlocks[genome][chromosome]
            #print counterOfAddedSynBlocks[genome]

    
def mergeIters(b, a):
    """Merges two iterators a and b, returning a single iterator that yields
    the elements of a and b in non-decreasing order.  a and b are assumed to each
    yield their elements in non-decreasing order."""

    done = object()
    aNext = next(a, done)
    startOfFree = 0
    bNext = next(b, done)
    endOfFree = bNext[0]

    while (aNext is not done) or (bNext is not done):
        #print aNext, bNext
        if ((aNext is not done) and (aNext[0] > startOfFree and aNext[1] < endOfFree)):
            if aNext[0] < endOfFree:
                aNext = next(a, done)
            else:
                yield aNext
                aNext = next(a, done)
        elif bNext is not done:
            while aNext is not done and aNext[0] <= endOfFree:
                aNext = next(a, done)
            #print bNext, "second"
            yield bNext
            startOfFree = bNext[1]
            bNext = next(b, done)
            if bNext is not done:
                endOfFree = bNext[0]
            else:
                endOfFree = startOfFree
                startOfFree = 1000000000
        else:
            if aNext is not done:
                if aNext[0] > endOfFree:
                    print "third", aNext
                    yield aNext
                    aNext = next(a, done)
                else:
                    aNext = next(a, done)

def main(filename, fileToMerge):
    synthenyBlocks = {}
    synthenyBlocksToMerge = {}
    genomes = []
    
    with open(filename) as f:
        for line, blocks in readBlocks(f, genomes):
            synthenyBlocks[line[0]] = blocks
    with open(fileToMerge) as f:
        for line, blocks in readBlocks(f, genomes):
            synthenyBlocksToMerge[line[0]] = blocks

    merge(synthenyBlocks, synthenyBlocksToMerge, genomes)



print '''Please, enter 2 files which you want to merge.
Second file should be in higher resolution than first.
Please, remember: your data should be filtered
(i.e. contains no duplications or deletions).
'''
cmdParams = sys.argv[1:]
filename = cmdParams[0]
fileToMerge = cmdParams[1]
main(filename, fileToMerge)
