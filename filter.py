#!/usr/bin/env python 
# -*- coding: utf-8 -*-


import sys
import re
from collections import defaultdict
import operator
from sets import Set

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
                sign = name[-5:-4]
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


def clusterization(synthenyBlocks, genomes):
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
    return sorted(noDupDel)


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


def outGRIMM(genomes, synthenyBlocks, outFile):
    # in: genomes and synthenyBlocks (for clusterByChr function) and name of output file
    # out: GRIMM-like file with sorted order of syntheny blocks
    #      > name of first genome
    #      1 -4 2 $ // chr1
    #      -3 5 6 $ // chr2
    #      > name of second genome
    #      ...
    # plus this function returns list of adjacency
    chromosomes = clusterByChr(genomes, synthenyBlocks)

    formAdjList = {}
    
    with open(outFile, "wb") as f:
        for genome in chromosomes:
            formAdjList[genome] = {}
            stringToWrite = ""
            stringToWrite += (">" + genome + "\n")
            for chrName in sorted(chromosomes[genome].iterkeys()):
                representation = []
                rightOrder = chromosomes[genome][chrName]
                rightOrder.sort(key = operator.itemgetter(1))
                for synBlock in rightOrder:
                    strSign = synBlock[3] + "1"
                    signedName = (synBlock[0] * int(strSign))
                    head = str(synBlock[0]) + "h"
                    tail = str(synBlock[0]) + "t"
                    if (signedName > 0):
                        representation.append((tail, head))
                    else:
                        representation.append((head, tail))
                    stringToWrite += (str(signedName) + " ")
                formAdjList[genome][chrName] = representation
                stringToWrite += "$\n"
            f.write(stringToWrite)
    return formAdjList


def adjacList(formAdjList):
    # in: dict of synblocks on chromosomes of genomes
    #     dict[genome][chr] = [()]
    # out: dict of adjacencies on chromosomes of genomes
    #      [(infty, head),(tail, head),...,(tail,infty)]
    adjList = {}
    for genome in formAdjList:
        genAdjList = []
        for chromosome in formAdjList[genome]:
            previousTail = "0"
            for tup in formAdjList[genome][chromosome]:
                genAdjList.append((previousTail, tup[0]))
                previousTail = tup[1]
            genAdjList.append((previousTail, "0"))
        adjList[genome] = Set(genAdjList)
    return adjList


def pairwiseNumberOfBreakpoints(adjList, adjOutFile):
    # in: dict of adjacencies in genome
    #     adjList[genome] = Set(adj)
    # out: matrix of pairwise number of breakpoints

    dictOfSpecies = {}
    topString = "\t"
    rows = []
    for genome1 in adjList:
        topString += (str(genome1) + "\t")
        row = (str(genome1) + "\t")
        for genome2 in adjList:
            sumOfUniques = 0
            for elem in adjList[genome1]:
                if elem not in adjList[genome2] and (elem[1], elem[0]) not in adjList[genome2]:
                    sumOfUniques += 1
            for elem in adjList[genome2]:
                if elem not in adjList[genome1] and (elem[1], elem[0]) not in adjList[genome1]:
                    sumOfUniques += 1
            row += (str(sumOfUniques) + "\t")
            dictOfSpecies[(genome1, genome2)] = sumOfUniques
        row += "\n"
        rows.append(row)
    topString += "\n"
    adjOutFile.write(topString)
    for elem in rows:
        adjOutFile.write(elem)
    return dictOfSpecies


def outputFilteredData(synthenyBlocks, genomes, outFiltFile):
    # in: dict of syntheny blocks, list of genomes and output file
    # out: file of sorted list of syntheny blocks without duplications and deletions
    namesOfFilteredSynBlocks = clusterization(synthenyBlocks, genomes)
    listOfNames = (namesOfFilteredSynBlocks)
    dictForSorting = {} # just tmp dict! useless
    listForSorting = [] # tmp list! useless
    for name in listOfNames:
        tmpName = re.search(r'\d+', name).group(0)
        dictForSorting[int(tmpName)] = name
        listForSorting.append(int(tmpName))
    
    for name in sorted(listForSorting):
        outString = dictForSorting[name] + "\n\n\n"
        for lists in synthenyBlocks[dictForSorting[name]]:
            newFormattedString = (lists[0] + "." +
                                  lists[1] + ":" + lists[2] + "-" +
                                  lists[3] + " " + lists[4] + "\n\n\n")
            outString += newFormattedString
        outFiltFile.write(outString)


def parseDistMatrix(synthenyBlocks, genomes, f):
    # in: file with matrices of distances in format
    #     # xK
    #     species
    #     specie dist1 dist2 ...
    # out: dict with pairwise distances between two species

    mapOfDists = {}
    
    trueNumOfK = re.search(r'\d+', filename).group(0)
    listOfLines = f.readlines()
    
    for ind in range(len(listOfLines)):
        line = listOfLines[ind]
        if line[0] == "#":
            numberOfK = re.search(r'\d+', line).group(0)
            if int(numberOfK) == int(trueNumOfK):
                ind += 1
                line = listOfLines[ind]
                line = line.strip()
                namesOfSpecies = line.split("\t")
                while line[0] != "#" and ind < len(listOfLines) - 1 and line != "\n":
                    counter = 1
                    ind += 1
                    line = listOfLines[ind]
                    line = line.strip()
                    listOfDists = line.split("\t")
                    if listOfDists != ['']:
                        for genome in namesOfSpecies:
                            mapOfDists[(genome,listOfDists[0])] = int(listOfDists[counter])
                            counter += 1
                    else:
                        break
                break
            else:
                continue
        else:
            continue
    return mapOfDists

def calculateBreakpointReuse(dictOfPairwiseDist, dictOfPairwiseBreakpoints, genomes, filename):
    # in: dict with pairwise distances between two species
    # out: file with breakpoint reuses for different sizes of syn blocks
    #     # xK
    #     species
    #     specie breakReUse1 ...
    
    numOfK = re.search(r'\d+', filename).group(0)
    line = "# Size of Block, *1.000 = " + str(numOfK) + "\n\n"
    with open(filename + ".reuse.txt", "wb") as f:
        for genome1 in genomes:
            line += ("breakpoint reuse rate between " + genome1 + " and:\n")
            for genome2 in genomes:
                if genome2 != genome1:
                    line += (genome2 + "\t")
                    line += str(4.0 * dictOfPairwiseDist[(genome1, genome2)]
                                / dictOfPairwiseBreakpoints[(genome1, genome2)])
                    line += "\n"
            line += "\n"
        f.write(line)
            
            

    

###############
# entry point #
###############
def main(filename, outFilteredFiles = False, distanceMatrixFile = ""):
    synthenyBlocks = {}
    genomes = []

    with open(filename) as f:
        for line, blocks in readBlocks(f, genomes):
            synthenyBlocks[line[0]] = blocks
    outFile = filename + ".output.txt"
    adjOutFile = filename + ".adjOutput.txt"
    formAdjList = outGRIMM(genomes, synthenyBlocks, outFile)               
    adjList = adjacList(formAdjList)
    outInFile = True
    
    with open(adjOutFile, "wb") as f:
        dictOfPairwiseBreakpoints = pairwiseNumberOfBreakpoints(adjList, f)
        
    if (outFilteredFiles == True):
        filteredDataOutFile = filename + ".filtered.txt"
        with open(filteredDataOutFile, "wb") as f:
            outputFilteredData(synthenyBlocks, genomes, f)
            
    if distanceMatrixFile:
        outInFile = False 
        with open(distanceMatrixFile) as f:
            dictOfPairwiseDist = parseDistMatrix(synthenyBlocks, genomes, f)
        calculateBreakpointReuse(dictOfPairwiseDist, dictOfPairwiseBreakpoints, genomes, filename)

        

    
cmdParams = sys.argv[1:]
filename = cmdParams[0]
outFilteredFiles = False
if len(cmdParams) >= 2 and cmdParams[1] == "True":
    outFilteredFiles = True
if len(cmdParams) == 3 and cmdParams[2] != "":
    distanceMatrixFile = cmdParams[2]
main(filename, outFilteredFiles, distanceMatrixFile)
