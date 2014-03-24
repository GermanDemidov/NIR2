#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import sys
import re
from collections import defaultdict
from optparse import OptionParser
import optparse


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
            tuplesFromSynBlocks[lists[0]][lists[1]].append(((lists[2]), (lists[3]), lists[4]))

    tuplesFromSynBlocksToMerge = {}
    for genome in genomes:
        tuplesFromSynBlocksToMerge[genome] = defaultdict(list)
    for names in synthenyBlocksToMerge:
        for lists in synthenyBlocksToMerge[names]:
            tuplesFromSynBlocksToMerge[lists[0]][lists[1]].append(((lists[2]), (lists[3]), lists[4]))

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
        for chromosome in tuplesFromSynBlocks[genome]:
            mergedSynBlocks[genome][chromosome] = formMap(tuplesFromSynBlocks[genome][chromosome],
                                                                  tuplesFromSynBlocksToMerge[genome][chromosome])
    return mergedSynBlocks


def countStatistics(synthenyBlocks, synthenyBlocksToMerge, genomes):
    mergedSynBlocks = merge(synthenyBlocks, synthenyBlocksToMerge, genomes)
    for genome in mergedSynBlocks:
        M1, M2, M3, N1, N2, N3 = (0, 0, 0, 0, 0, 0)
        for chromosome in mergedSynBlocks[genome]:
            if len(mergedSynBlocks[genome][chromosome]) > 0:
                aDict = mergedSynBlocks[genome][chromosome][0]
                for key in aDict:
                    countAdj = 0
                    countBlocks = 0
    
                    for values in aDict[key]:
                        if len(values) == 3:
                            countBlocks += 1
                        elif len(values) == 2:
                            countAdj += 1
                        else:
                            print "error in len of value"
                    if countBlocks == 1:
                        M2 += 1
                    elif countBlocks > 1:
                        M1 += 1
                    N1 += countAdj
                    
            if len(mergedSynBlocks[genome][chromosome]) > 0:

                aAdjDict = mergedSynBlocks[genome][chromosome][1]
                for key in aDict:
                    countAdj = 0
                    countBlocks = 0
    
                    for values in aDict[key]:
                        if len(values) == 3:
                            countBlocks += 1
                        elif len(values) == 2:
                            countAdj += 1
                        else:
                            print "error in len of value"
                    if countBlocks > 1:
                        M3 += 1
                    if countAdj == 1:
                        N2 += 1
                    elif countAdj > 1:
                        N3 += countAdj
        print genome
        print 'M1 = {0}, M2 = {1}, M3 = {2}, N1 = {3}, N2 = {4}, N3 = {5}'.format(M1, M2, M3, N1, N2, N3)


            
def formAdj(chromosome, maximum):
    """make an adjacency list from a list of syntheny blocks"""
    adjList = []
    start = 0
    for elem in chromosome:
        adjList.append((start, elem[0]))
        start = elem[1]
    adjList.append((start, maximum))
    return adjList



def formMap(lowResGenome, highResGenome):
    """Returns elements of class that lays between the borders of elements of class b. Class b:
    elements of syntheny blocks and adj list"""
    maximum = 3 * (10 ** 9) # 3 billions of bp length - no one chromosome can reach this limit! =)
    # of course we can compute this value, but...why we need this if we can get it for free? =)
    a = iter(lowResGenome)
    b = iter(highResGenome)
    aAdj = iter(formAdj(lowResGenome, maximum))

    bAdj = iter(formAdj(highResGenome, maximum))
    aDict = defaultdict(list)
    aAdjDict = defaultdict(list)

    currLow = next(aAdj)
    currHigh = next(bAdj)
    done = object()
    while (currLow is not done) or (currHigh is not done):
        if currLow[1] == currHigh[1] == maximum and currHigh[0] >= currLow[0]:
            aAdjDict[currLow].append(currHigh)
            break
        elif currLow[1] == currHigh[1] == maximum:
            break
        elif currLow[1] == maximum:
            if len(currHigh) == 3:
                if currHigh[0] >= currLow[0]:
                    aAdjDict[currLow].append(currHigh)
                currHigh = next(bAdj, done)

            elif len(currHigh) == 2:
                if currHigh[0] >= currLow[0]:
                    aAdjDict[currLow].append(currHigh)
                currHigh = next(b, done)
        elif currHigh[1] == maximum:
            if len(currLow) == 3:
                currLow = next(aAdj, done)
            elif len(currLow) == 2:
                currLow = next(a, done)
        elif (currHigh[0] == currLow[0] and currHigh[1] == currLow[1]):
            if len(currLow) == 2 and len(currHigh) == 2:
                aAdjDict[currLow].append(currHigh)
                currHigh = next(b, done)
                currLow = next(a, done)
                continue
            elif len(currHigh) == 3 and len(currLow) == 3:
                aDict[currLow].append(currHigh)
                currHigh = next(bAdj, done)
                currLow = next(aAdj, done)
                continue
            elif len(currHigh) == 3 and len(currLow) == 2:
                aAdjDict[currLow].append(currHigh)
                print currHigh, currLow
                print "ERROR"
                currHigh = next(bAdj, done)
                currLow = next(a, done)
        elif (currHigh[0] >= currLow[0] and currHigh[1] <= currLow[1]):
            if len(currHigh) == 3 and len(currLow) == 3:
                aDict[currLow].append(currHigh)
                currHigh = next(bAdj, done)
            elif len(currHigh) == 2 and len(currLow) == 2:
                aAdjDict[currLow].append(currHigh)
                currHigh = next(b, done)
            elif len(currHigh) == 3 and len(currLow) == 2:
                aAdjDict[currLow].append(currHigh)
                currHigh = next(bAdj, done)
            elif (len(currHigh) == 2 or len(currHigh) == 1) and len(currLow) == 3:
                aDict[currLow].append(currHigh)
                currHigh = next(b, done)
        elif (currHigh[0] < currLow[1] and currHigh[1] > currLow[1]):
            if len(currHigh) == 2:
                """if len(currLow) == 2:
                    aAdjDict[currLow].append(currHigh)
                elif len(currLow) == 3:
                    aDict[currLow].append(currHigh)"""
                currHigh = next(b, done)

            elif len(currHigh) == 3:
                """if len(currLow) == 2:
                    aAdjDict[currLow].append(currHigh)
                elif len(currLow) == 3:
                    aDict[currLow].append(currHigh)"""
                currHigh = next(bAdj, done)

        elif (currHigh[0] < currLow[0] and currHigh[1] > currLow[0]):
            if len(currHigh) == 2:
                """if len(currLow) == 2:
                    aAdjDict[currLow].append(currHigh)
                elif len(currLow) == 3:
                    aDict[currLow].append(currHigh)"""
                currHigh = next(b, done)

            elif len(currHigh) == 3:
                """if len(currLow) == 2:
                    aAdjDict[currLow].append(currHigh)
                elif len(currLow) == 3:
                    aDict[currLow].append(currHigh)"""
                currHigh = next(bAdj, done)

        elif currHigh[0] >= currLow[1]:
            if len(currLow) == 2:
                currLow = next(a, done)
            elif len(currLow) == 3:
                currLow = next(aAdj, done)

    return [aDict, aAdjDict]

            


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

    countStatistics(synthenyBlocks, synthenyBlocksToMerge, genomes)


if __name__ == "__main__":

    parser = OptionParser(version = "%prog 1.0", conflict_handler="error", description =
                            """This script merges one file with high-resolution syntheny
                            blocks to another file with low-resolution syntheny blocks.""")

    parser.add_option("-f", "--file", dest="filename", help="file with syntheny blocks of low resolution", metavar="FILE")
    parser.add_option("-m", "--merge", dest="fileToMerge", help="file with syntheny blocks of high resolution", metavar="FILE")

    options, args = parser.parse_args()
    
    main(options.filename, options.fileToMerge)
