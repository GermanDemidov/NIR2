#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import sys
import re
from collections import defaultdict
from optparse import OptionParser
import optparse
from mapofblocks import *

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


    sixDimVect = {}
    tuplesFromSynBlocksToMerge = {}
    for genome in genomes:
        tuplesFromSynBlocksToMerge[genome] = defaultdict(list)
    for names in synthenyBlocksToMerge:
        for lists in synthenyBlocksToMerge[names]:
            tuplesFromSynBlocksToMerge[lists[0]][lists[1]].append(((lists[2]), (lists[3]), lists[4]))
            sixDimVect[lists[0]] = ["" for i in range(6)]
    for genome in genomes:
        for chromosome in tuplesFromSynBlocks[genome]:
            sortWithKey(tuplesFromSynBlocks[genome][chromosome])
            
    for genome in genomes:
        for chromosome in tuplesFromSynBlocksToMerge[genome]:
            sortWithKey(tuplesFromSynBlocksToMerge[genome][chromosome])
    # in merged synBlocks we want to write chromosome cart (to convert it in
    # GRIMM format)
    mergedSynBlocks = {}
    softMergedSynBlocks = {}
    counterOfAddedSynBlocks = {}
    for genome in genomes:
        mergedSynBlocks[genome] = defaultdict(list)
        softMergedSynBlocks[genome] = defaultdict(list)
        counterOfAddedSynBlocks[genome] = 0

    for genome in genomes:
        # print genome
        for chromosome in tuplesFromSynBlocks[genome]:
            """
            
            print "Synteny blocks: "
            print tuplesFromSynBlocks[genome][chromosome]
            print "Adjacencies: "
            print formAdj(tuplesFromSynBlocks[genome][chromosome], "infty")
            print "Synteny blocks of high resolution "
            print tuplesFromSynBlocksToMerge[genome][chromosome]
            print "Adjacencies of high resolution: "
            print formAdj(tuplesFromSynBlocksToMerge[genome][chromosome], "infty")"""
            softMergedSynBlocks[genome][chromosome] = formMapSmooth(tuplesFromSynBlocks[genome][chromosome],
                                                                  tuplesFromSynBlocksToMerge[genome][chromosome],
                                                                    sixDimVect)
            # mergedSynBlocks[genome][chromosome] = formMapStrict(tuplesFromSynBlocks[genome][chromosome],
            #                                                      tuplesFromSynBlocksToMerge[genome][chromosome])


            """print "Each synteny block of low res contains (of high res): "
            for k, v in softMergedSynBlocks[genome][chromosome][0].iteritems():
                print k, v
            print "Each adjacencie of low res contains (of high res): "
            for k, v in softMergedSynBlocks[genome][chromosome][1].iteritems():
                print k, v"""
            
    return softMergedSynBlocks


def countStatistics(synthenyBlocks, synthenyBlocksToMerge, genomes):
    mergedSynBlocks = merge(synthenyBlocks, synthenyBlocksToMerge, genomes)
    for genome in mergedSynBlocks:
        overlapStats = 0
        length = 0
        M1, M2, M3, N1, N2, N3 = (0, 0, 0, 0, 0, 0)
        fake = 0
        totlen = 0
        for chromosome in mergedSynBlocks[genome]:
                #M1, M2, M3, N1, N2, N3 = (0, 0, 0, 0, 0, 0)
                #fake = 0
                aDict, aAdjDict, tmpOverlap, fakesNum = mergedSynBlocks[genome][chromosome]
                """print
                print aDict
                print aAdjDict
                print"""
                
                totlen += len(aDict)
                fake += fakesNum
                for key in aDict:
                    countAdj = 0
                    countInnerAdj = 0
                    countBlocks = 0
                    tmpAdj = (-100,0)
                    
                    for values in aDict[key]:
                        if len(values) == 3:
                            countBlocks += 1
                        elif len(values) == 2:
                            if values[1] - values[0] != 0:
                                tmpAdj = values
                            if values[1] - values[0] != 0:
                                countAdj += 1
                            if values[0] > key[0] and values[1] < key[1]:# and values[1] - values[0] != 0:
                                countInnerAdj += 1
                            elif values[0] >= key[0] and values[1] <= key[1] and values[1] - values[0] != 0:
                                countInnerAdj += 1
                        else:
                            print "error in len of value"


                    if countBlocks > 1 or countAdj > 2 or tmpAdj[0] > key[0]:
                        M1 += 1
                    elif countBlocks == 1:
                        M2 += 1
                    elif countBlocks == 0:
                        print "Hmmm..."
                        print key
                        print aDict[key]
                    else:
                        print key
                    N1 += countInnerAdj
                    
                for key in aAdjDict:
                    countAdj = 0
                    countBlocks = 0
                    flagOfDivisionByBlock = False
                    tmpAdj = ()

                    for values in aAdjDict[key]:
                        if len(values) == 3:
                            countBlocks += 1
                            if values[0] >= key[0] and values[1] <= key[1] and values[1] - values[0] != 0:
                                flagOfDivisionByBlock = True
                        elif len(values) == 2:
                            if values[1] < 3000000000:
                                tmpAdj = values
                            if values[1] - values[0] != 0:
                                countAdj += 1
                        else:
                            continue
                            print "error in len of value"
                    if countBlocks >= 1 and countAdj == 0:
                        countAdj += 2
                    if flagOfDivisionByBlock == True:
                        M3 += 1
                    if countAdj == 0 or countAdj == 1:# and countBlocks == 0:
                        N2 += 1
                    elif countBlocks >= 1:
                        """(countAdj > 1 or countBlocks > 1 or
                          (countAdj >= 1 and countBlocks >= 1)):"""
                        N3 += countAdj
                """print "Fake: ", fake
                # 18, 15, 11, 1, 5,
                print "Number of blocks in high: ", M1 + N1 + M2 + N3 - M3 - fake
                print 'M1 = {0}, M2 = {1}, M3 = {2}, N1 = {3}, N2 = {4}, N3 = {5}'.format(M1, M2, M3, N1, N2, N3)"""
                
                        
                overlapStats += len(tmpOverlap)
        print genome
        print totlen
        print "Fake: ", fake
        print "Number of blocks in high: ", M1 + N1 + M2 + N3 - M3 - fake
        print "Num of synteny blocks of low res that has overlaps =", overlapStats
        print 'M1 = {0}, M2 = {1}, M3 = {2}, N1 = {3}, N2 = {4}, N3 = {5}'.format(M1, M2, M3, N1, N2, N3)






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
