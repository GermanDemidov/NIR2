#!/usr/bin/env python 
# -*- coding: utf-8 -*-
from collections import defaultdict
import sys

            
def formAdj(chromosome, maximum):
    """make an adjacency list from a list of syntheny blocks"""
    adjList = []
    start = 0
    for elem in chromosome:
        adjList.append((start, elem[0]))
        start = elem[1]
    adjList.append((start, maximum))
    return adjList

def formFullMapAdj(chromosome, maximum):
    """make an adjacency list from a list of syntheny blocks"""
    adjList = []
    start = 0
    for i in range(len(chromosome)):
        adjList.append((start, chromosome[i][0]))
        start = chromosome[i][1]
        adjList.append(chromosome[i])
    adjList.append((start, maximum))
    return adjList

def appending(aDict, aAdjDict, currentLow, currentHigh):
    if len(currentLow) == 3 and currentHigh not in aDict[currentLow]:
        aDict[currentLow].append(currentHigh)
    elif len(currentLow) == 2 and currentHigh not in aAdjDict[currentLow]:
        aAdjDict[currentLow].append(currentHigh)


def isInside(block, med):
    if len(block) == 3 and block[0] <= med and block[1] >= med:
        return True
    elif len(block) == 2 and block[0] < med and block[1] > med:
        return True
    else:
        return False


def formMapSmooth(lowResGenome, highResGenome, sixDimVect):
    """returns statistics for M1, M2, M3, N1, N2, N3 for different genomes/chromosomes"""
    maximum = 30000000000
    lowRes = formFullMapAdj(lowResGenome, maximum)
    highRes = formFullMapAdj(highResGenome, maximum)

    length1 = len(lowRes)
    length2 = len(highRes)

    dictOfFakeSynBlocks = defaultdict(bool)
    listOfFakes = []
    k = 0
    l = 0
    while k < len(lowResGenome):
        if k < len(lowResGenome):
            currentLow = lowResGenome[k]

        dictOfFakeSynBlocks[currentLow]
            
        for currentHigh in highResGenome:
            med = (currentHigh[0] + currentHigh[1]) / 2
            if isInside(currentLow, med):

                dictOfFakeSynBlocks[currentLow] = True

        l += 1
        k += 1
    for key in sorted(dictOfFakeSynBlocks.iterkeys()):
        if dictOfFakeSynBlocks[key] == False:
            listOfFakes.append(key)

    m = 0
    newHighRes = []
    filteredHighRes = []
    if len(listOfFakes) > 0:
        newHighRes.extend(listOfFakes)
        newHighRes.extend(highResGenome)
        newHighRes.sort(key=lambda x: x[0])
        prevElem = (0,0)
        nextElem = (maximum, maximum)
        for i in range(len(newHighRes)):
            elem = newHighRes[i]
            if i > 0:
                prevElem = newHighRes[i-1]
            else:
                prevElem = (0, 0, "+")
            if i < len(newHighRes) - 1:
                nextElem = newHighRes[i+1]
            else:
                nextElem = (maximum, maximum, "+")
            if elem[0] < prevElem[1] and elem[1] > nextElem[0]:

                continue
            elif elem[1] < nextElem[0]:
                filteredHighRes.append(elem)
            elif elem in listOfFakes:
                filteredHighRes.append(elem)
            elif elem[1] >= nextElem[0]: #and nextElem in listOfFakes:
                filteredHighRes.append((elem[0], nextElem[0], elem[2]))
            elif elem[0] <= prevElem[1]: #and prevElem in listOfFakes:
                filteredHighRes.append((prevElem[1], elem[1], elem[2]))
            elif elem[0] >= prevElem[1] and elem[1] <= nextElem[0]:
                filteredHighRes.append(elem)
            else:
                print elem
                print nextElem, "NEXT"
                elem = list(elem)
                elem[0] = nextElem[0]
                filteredHighRes.append(tuple(elem))
                """print elem
                print nextElem
                print newHighRes"""
        if len(filteredHighRes) != len(highResGenome) + len(listOfFakes):
            print "WOW"
            print len(filteredHighRes)
            print len(highResGenome), len(listOfFakes)
        
        highRes = formFullMapAdj(filteredHighRes, maximum)
        length2 = len(highRes)

                    

    aDict = defaultdict(list)
    aAdjDict = defaultdict(list)
    
    i = 0 # low res counter
    j = 0 # high res counter
    currentLow = lowRes[i]
    currentHigh = highRes[j]
    #print currentLow, currentHigh

    overlapStats = defaultdict(list)
    
    # always add first telomere end of high res in low res
    appending(aDict, aAdjDict, currentLow, currentLow)
    j += 1

    while i < length1 and j < length2:
        if len(currentLow) == 3 and currentLow not in aDict:
            aDict[currentLow]
        if len(currentLow) == 2 and currentLow not in aAdjDict:
            aAdjDict[currentLow]
        currentLow = lowRes[i]
        if len(currentLow) == 2 and currentLow[1] - currentLow[0] == 0\
           and currentLow not in aAdjDict:
            aAdjDict[currentLow]
            i += 1
            continue

        currentHigh = highRes[j]

        med = (currentHigh[0] + currentHigh[1]) / 2
        if len(currentHigh) == 2 and currentHigh[1] - currentHigh[0] == 0\
           and isInside(currentLow, med):
            appending(aDict, aAdjDict, currentLow, currentHigh)
            j += 1
            continue
        #print currentLow, currentHigh

        if currentHigh[0] >= currentLow[1]:
            i += 1
            continue

        # if the block of high resolution lays within the block of low resolution
        elif currentLow[0] <= currentHigh[0] and currentLow[1] >= currentHigh[1]:
            appending(aDict, aAdjDict, currentLow, currentHigh)
            j += 1

        # if the block of low resolution lays within the block of high resolution    
        elif currentLow[0] >= currentHigh[0] and currentLow[1] < currentHigh[1]:
            i += 1
            if isInside(currentLow, med):
                appending(aDict, aAdjDict, currentLow, currentHigh)
            elif len(currentHigh) == 3 and len(currentLow) == 3:
                aDict[currentLow]
            elif len(currentHigh) == 2 and len(currentLow) == 3:
                aDict[currentLow]
            elif len(currentHigh) == 2 and len(currentLow) == 2:
                aAdjDict[currentLow]
            elif len(currentHigh) == 3 and len(currentLow) == 2:
                appending(aDict, aAdjDict, currentLow, currentHigh)
            if currentLow[1] == currentHigh[1]:
                j += 1
                i += 1
                


        # if the block of low resolution overlaps with block of high resolution
        elif currentLow[0] > currentHigh[0] and currentLow[0] < currentHigh[1]:
            if len(currentLow) == 3 and len(currentHigh) == 3:
                overlapStats[currentLow].append(currentHigh)
            if isInside(currentLow, med):
                appending(aDict, aAdjDict, currentLow, currentHigh)
                j += 1
            else:
                appending(aDict, aAdjDict, lowRes[i-1], currentHigh)
                j += 1
        elif currentLow[1] > currentHigh[0] and currentLow[1] < currentHigh[1]:
            if len(currentLow) == 3 and len(currentHigh) == 3:
                overlapStats[currentLow].append(currentHigh)
            if isInside(currentLow, med):
                appending(aDict, aAdjDict, currentLow, currentHigh)
                i += 1
            else:
                if len(currentLow) == 3 and currentLow not in aDict:
                    aDict[currentLow]
                elif isInside(lowRes[i+1], med):
                    appending(aDict, aAdjDict, lowRes[i+1], currentHigh)
                i += 1
        elif currentLow[1] == currentHigh[0]:
            i += 1
        elif currentLow[0] == currentHigh[1]:
            j += 1
        else:
            print currentLow, currentHigh
            sys.exit("Error message")

    if lowResGenome[-1] not in aDict:
        print lowResGenome[-1]
        print aDict
        

    #while currentLow[1] == maximum and currentHigh[1] != maximum:
        
    #print aDict, aAdjDict
    #print len(aDict), len(lowResGenome)

        
    

    return (aDict, aAdjDict, overlapStats, len(listOfFakes))














def formMapSmooth12(lowResGenome, highResGenome, sixDimVect):
    """returns statistics for M1, M2, M3, N1, N2, N3 for different genomes/chromosomes"""
    maximum = 30000000000
    lowRes = formFullMapAdj(lowResGenome, maximum)
    highRes = formFullMapAdj(highResGenome, maximum)

    length1 = len(lowRes)
    length2 = len(highRes)

    aDict = defaultdict(list)
    aAdjDict = defaultdict(list)
    
    i = 0
    j = 0
    currentLow = lowRes[i]
    currentHigh = highRes[j]

    overlapStats = 0

    appending(aDict, aAdjDict, currentLow, currentHigh)
    i += 1
    if currentHigh[1] <= currentLow[1]:
        j += 1
    
    while (currentLow[1] != maximum and currentHigh[1] != maximum):
        currentLow = lowRes[i]
        currentHigh = highRes[j]
        currHighMed = float(currentHigh[1] + currentHigh[0]) / 2
        if ((currentLow[0] <= currentHigh[0] and currentLow[1] >= currentHigh[1])):
            appending(aDict, aAdjDict, currentLow, currentHigh)
            j += 1
        elif (currentHigh[0] <= currentLow[0] and currentHigh[1] >= currentLow[1]):
            appending(aDict, aAdjDict, currentLow, currentHigh)
            j += 1
            i += 1

        elif currentHigh[0] < currentLow[0] and currentHigh[1] > currentLow[0]:
            overlapStats += 1
            j += 1
            if currHighMed >= currentLow[0] and currHighMed < currentLow[1]:
                appending(aDict, aAdjDict, currentLow, currentHigh)
            else:
                prevLow = lowRes[i-1]
                appending(aDict, aAdjDict, prevLow, currentHigh)
                
        elif currentHigh[0] < currentLow[1] and currentHigh[1] > currentLow[1]:
            overlapStats += 1
            j += 1
            if currHighMed < currentLow[1] and currentHigh[1] != maximum:
                appending(aDict, aAdjDict, currentLow, currentHigh)
            elif currentHigh[1] != maximum:
                nextLow = lowRes[i + 1]
                appending(aDict, aAdjDict, nextLow, currentHigh)

        elif currentHigh[0] >= currentLow[1]:
            i += 1
        flagOfBlocksWithin = False
        if currentLow in aDict:
            for elem in aDict[currentLow]:
                if len(elem) == 3:
                    flagOfBlocksWithin = True
        if ((len(currentLow) == 3 and currentLow not in aDict)
            or (flagOfBlocksWithin == False and len(currentLow) == 3  and currentLow != lowRes[i])):
            aDict[currentLow]


    
    while not (currentLow[1] == maximum and currentHigh[1] == maximum):
        if (currentLow[1] != maximum):
            currentLow = lowRes[i]
        elif (currentHigh[1] != maximum):
            currentHigh = highRes[j]
        if currentHigh[1] == maximum:
            currentLow = lowRes[i]
            if currentLow[1] > currentHigh[0] and currentLow[1] == maximum:
                appending(aDict, aAdjDict, currentLow, currentHigh)
            elif currentLow not in aDict and currentLow not in aAdjDict and len(currentLow) == 3:
                aDict[currentLow]
            i += 1
        elif currentLow[1] == maximum:
            currentHigh = highRes[j]
            if (currentHigh[0] >= currentLow[0]):
                appending(aDict, aAdjDict, currentLow, currentHigh)
            else:
                overlapStats += 1
                currHighMed = float(currentHigh[1] + currentHigh[0]) / 2
                if currHighMed > currentLow[0]:
                    appending(aDict, aAdjDict, currentLow, currentHigh)
                else:
                    prevLow = lowRes[i-1]
                    appending(aDict, aAdjDict, prevLow, currentHigh)
            j += 1

    return (aDict, aAdjDict, overlapStats)



            

def formMapStrict(lowResGenome, highResGenome):
    """Returns elements of class that lays between the borders of elements of class b. Class b:
    elements of syntheny blocks and adj list"""
    flag = False
    maximum = "infty" #3 * (10 ** 9) # 3 billions of bp length - no one chromosome can reach this limit! =)
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
    overlapStats = 0
    anotherOverlap = 0
    # TODO: rewrite this part of code with while loops
    # don't forget to process 1) first part 2) main part 3) end
    while (currLow[1] is not done) or (currHigh[1] is not done):
        if currLow[1] == currHigh[1] == maximum and currHigh[0] >= currLow[0]:
            aAdjDict[currLow].append(currHigh)
            break
        elif currLow[1] == currHigh[1] == maximum and currHigh[0] < currLow[0]:
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
        elif currHigh[1] == maximum and currLow[1] != maximum:
            if len(currLow) == 3:
                currLow = next(aAdj, done)
            elif len(currLow) == 2:
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
                currHigh = next(b, done)
                anotherOverlap += 1

            elif len(currHigh) == 3:
                if len(currLow) == 3:
                    flag = True
                    overlapStats += 1
                    # print currHigh, currLow
                else:
                    flag = True
                    overlapStats += 1
                    # print currHigh, currLow
                currHigh = next(bAdj, done)

        elif (currHigh[0] < currLow[0] and currHigh[1] > currLow[0]):
            if len(currHigh) == 2:
                currHigh = next(b, done)
                anotherOverlap += 1
            elif len(currHigh) == 3:
                if len(currLow) == 3:
                    flag = True
                    overlapStats += 1
                    # print currHigh, currLow

                else:
                    anotherOverlap += 1
                currHigh = next(bAdj, done)

        elif currHigh[0] >= currLow[1]:
            if len(currLow) == 2:
                currLow = next(a, done)
            elif len(currLow) == 3:
                currLow = next(aAdj, done)

    return (aDict, aAdjDict, overlapStats, anotherOverlap, flag)
