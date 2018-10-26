# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:36:13 2018

@author: DIC SKG
"""
import numpy as np
import math

# Number of mismatches between strings is called hamming distance 
def HammingDistance (p, q):
    count = 0
    l = len(p)
    for i in range (l):
        if p[i] != q[i]:
            count += 1
    return count

def NumberToSymbol (k):
    if k == 0:
        return 'A'
    elif k == 1:
        return 'C'
    elif k == 2:
        return 'G'
    elif k == 3:
        return 'T'


def NumberToPattern (num, k):
    if k == 1:
        return NumberToSymbol (num)
    else:
        r = num % 4
        prefixnum = num // 4
        # print(r)
        prefix = NumberToPattern (prefixnum, k - 1)
        symbol = NumberToSymbol (r)
        return prefix + symbol

def motif(dna,pat,d):
    k = len(pat)
    l = len(dna)
    for i in range(0,l-k+1):
        pat1 = dna[i:i+k]
        if(HammingDistance(pat,pat1)) <= d:
            return True
    return False

# Brute force (also known as exhaustive search) is a general problem-solving technique that explores all possible solution candidates and checks whether each candidate solves the problem. 
def MotifEnumeration(Dna, k, d):
    Patterns = [ ]
    l1 = len(Dna)
    #print(l1)
    
    for i in range(4**k):
        pat1 = NumberToPattern(i,k)
        #print(pat1)
        a = 0
        for j in range(0,l1):
            if motif(Dna[j],pat1,d):
                a = a + 1
                #print(pat1)
                #print(Dna[j])
        if a == l1 :
                Patterns.append(pat1)
    Patterns = set(Patterns)
        
    return Patterns

def Motifs(Profile, Dna):
    motifs = []
    t = len(Dna)
    k = len(Profile['A'])
    for i in range(t):
        motifs.append(ProfilemostProbable(Dna[i], k, Profile))
    return motifs

# The most popular letters in each column of the motif matrix
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

# Counting the number of occurrences of each nucleotide in each column of the motif matrix; the (i, j)-th element of Count(Motifs) stores the number of times that nucleotide i appears in column j of Motifs
def Count(Motifs):
    l = len(Motifs)
    m = len(Motifs[0])
    a = {}
    for symbol in "ACGT":
        a[symbol] = []
        for j in range(m):
            a[symbol].append(0)

    for i in range(l):
        for j in range(m):
            symbol = Motifs[i][j]
            a[symbol][j] += 1

    return a

# A profile matrix P = Profile(Motifs) for which Pi,j is the frequency of the i-th nucleotide in the j-th column of the motif matrix. Note that the elements of any column of the profile matrix sum to 1
def Profile(Motifs):
    profile = {}
    t = len(Motifs)
    CountMotifs = Count(Motifs)

    for symbol in "ACGT":
        profile[symbol] = []

    for x in CountMotifs:
        for y in CountMotifs[x]:
            z = y/float(t)
            profile[x].append(z)
    return profile


# Entropy is a measure of the uncertainty of a probability distribution (p1, …, pN).
def Entropy(Motifs):
    a = Profile(Motifs)
    m = len(Motifs[0])
    #e = [0] * 1
    #for i in range(1):
    #    e[i] = [0] * m
    e = 0
    for i in range(m):
        for j in range(4):
            if a[j][i] != 0:
                l = math.log2(a[j][i])
                #e[0][i] = e[0][i]*(-1) - (a[j][i]*l)
                e = e - a[j][i]*l
    return e

# Number of unpopular (lower case) letters in the motif matrix Motifs
def Score(Motifs):
    count = 0
    k = len(Motifs[0])
    t = len(Motifs)
     ConsensusMotif = Consensus(Motifs)
    for i in range(t):
        for j in range(k):
            if Motifs[i][j] != ConsensusMotif[j]:
                count += 1
    return count

   
# Sum of Hamming distances between Pattern and each Motifi
def d(pattern, motifs):
    a = 0
    for m in motifs:
        a = a + HammingDistance(m,pattern)
    return a


# Minimum Hamming distance between Pattern and any k-mer in Text
def dtext(pattern, text):
    l = len(text)
    k = len(pattern)
    m = [ ]
    for j in range(l-k+1):
            pat = text[j:j+k]
            m.append(HammingDistance(pattern, pat))
    minm = m.index(min(m))
    motimi = text[minm:minm+k]
    return motimi    

# Sum of distances between Pattern and all strings in Dna
def ddna(pattern, dna):
    n = len(dna)
    l = len(dna[0])
    k = len(pattern)
    motimi = 0
    for i in range(n):
        m = [ ]
        for j in range(l-k+1):
                pat = dna[i][j:j+k]
                m.append(HammingDistance(pattern, pat))
        minm = m.index(min(m))
        motimi = motimi + HammingDistance(dna[i][minm:minm+k], pattern)
    return motimi

                
# Finding a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern    
def MedianString(Dna, k):
    m = []
    for i in range(4**k):
        pat = NumberToPattern(i,k)
        m.append(ddna(pat, Dna))
    minm = m.index(min(m))
    Median = NumberToPattern(minm,k)        
    return Median
  
# Fast heuristics that trade accuracy for speed in order to find an approximate solution    
def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p
   
# k-mer that was most likely to have been generated by Profile among all k-mers in Text    
def ProfilemostProbable(text,k,profile):
    p_dict = {}
    for i in range (len (text) - k + 1):
        p = Pr (text[i: i + k], Profile)
        p_dict[i] = p
    m = max (p_dict.values ())
    keys = [k for k, v in p_dict.items () if v == m]
    ind = keys[0]
    return text[ind: ind + k]


# Algorithm to find motifs  
def GreedyMotifSearch(dna,k,t):
    bestmotifs = []
    l = len(dna[0])
    for i in range(t):
        bestmotifs.append(dna[i][0:k])
    base = dna[0]
    motifs = []
    for i in range(l-k+1):
        motifs.append(base[i:i+k])
        for j in range(1,t):
            profilematrix = Profile(motifs[0:j])
            nextmotif = ProfilemostProbable(dna[j],k,profilematrix)
            motifs.append(nextmotif)           
        if Score(motifs) < Score(bestmotifs):
                bestmotifs = motifs
    index = ' '.join(str(e) for e in bestmotifs)            
    return index




