# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# Where in the genome does replication begin ?


# An error during replication can lead to various diseases, including cancer. To understand how replication initiation works and what causes it to malfunction, we must first know where to look for replication origins. For this reason, we must accurately locate ori sites in the genome to study their replication initiation. Things are made even more difficult when we move from bacteria to more complex organisms; the human genome has thousands of origins of replication.
# Legrand's method to look for frequent words in the ori because for various biological processes, certain nucleotide strings often appear surprisingly often in small regions of the genome. This is often because certain proteins can only bind to DNA if a specific string of nucleotides is present, and if there are more occurrences of the string, then it is more likely that binding will successfully occur.

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count + 1
            
    return count


# Algorithm for finding the most frequent k-mers (DnaA boxes are 9-mers) in a string Text checks all k-mers appearing in this string and then computes how many times each k-mer appears in Text. 

def FrequentWords(Text, k):
    Count = []
    Freq = []
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count.append(PatternCount(Text,Pattern))
    maxCount = max(Count)
    for i in range(len(Count)):
        if Count[i] == maxCount:
            Freq.append(Text[i:i+k])
    Freq = list(set(Freq))
    return Freq
    #string = ' '.join(str(e) for e in Freqarray)


#  Having one strand and a supply of “free floating” nucleotides, one can imagine the synthesis of a complementary strand on a template strand. This model of replication was confirmed rigorously by Meselson and Stahl in 1958. The reverse complement of a string Pattern = p1 … pn is the string Patternrc = pn* … p1* formed by taking the complement of each nucleotide in Pattern, then reversing the resulting string. 

def ReverseComplement(Pattern):
    rcomplement = ""
    for char in Pattern:
        if char == "A":
            rcomplement = "T" + rcomplement
        elif char == "T":
            rcomplement = "A" + rcomplement 
        elif char == "G":
            rcomplement = "C" + rcomplement 
        elif char == "C":
            rcomplement = "G" + rcomplement 
            
    return rcomplement


# DnaA protein that binds to DnaA boxes and initiates replication does not care which of the two strands it binds to, also reverse compliment of the same k-mer are considered for DnaA box determination algorithm but before concluding about the DnaA boxes will find it's occurence in the ori
# Algorithm for finding occurence of a perticular nucleotide sequence

def PatternMatching(Pattern, Genome):
    Ck = []
    for i in range(len(Genome)-len(Pattern)+1):
        if Pattern == Genome[i:i+len(Pattern)]:
            Ck.append(i)
    #index = ' '.join(str(e) for e in Ck)
    print(Ck)
    

# All lexicographically ordered k-mers, the resulting list is still ordered lexicographically 

def SymbolToNumber(Pattern):
    if Pattern == "A":
        return 0
    elif Pattern == "C":
        return 1
    elif Pattern == "G":
        return 2
    elif Pattern == "T":
        return 3
def PatternToNumber(Pattern):
    k = len(Pattern)
    if Pattern == "A" or Pattern == "C" or Pattern == "G" or Pattern == "T":
        return SymbolToNumber(Pattern)
    else:
        Prefix = Pattern[:k-1]
        Symbol = Pattern[k-1]
        return 4 * PatternToNumber(Prefix) + SymbolToNumber(Symbol)

# Vice versa
        
def NumberToSymbol(k):
    if k == 0:
        return 'A'
    elif k == 1:
        return 'C'
    elif k == 2:
        return 'G'
    elif k == 3:
        return 'T'
    
    
def NumberToPattern(num, k):
    if k == 1:
        return NumberToSymbol(num)
    else:
        r = num % 4
        prefixnum = num//4
        #print(r)
        prefix = NumberToPattern(prefixnum, k - 1)
        symbol = NumberToSymbol(r)
        return prefix + symbol
        

# Algorithm to compute frequencies of all the k-mers in lexicographic order    

def computing_frequencies(Text, k):
    Freqarray = []
    m = []
    for i in range(4**k):
        Freqarray.append(0)
    for i in range(len(Text)-k+1):
        m.append(PatternToNumber(Text[i:i+k]))
    for i in range(len(Text)-k+1):
        j = PatternToNumber(Text[i:i+k])
        Freqarray[j] = m.count(j)
    #string = ' '.join(str(e) for e in Freqarray)
    return Freqarray  


# Clumps, i.e., part of genome where nucleotide sequence appear close to each other in a small region of the genome which may represent the hidden message to DnaA to start replication
# A k-mer Pattern forms an (L, t)-clump inside a (longer) string Genome if there is an interval of Genome of length L in which this k-mer appears at least t times. (This definition assumes that the k-mer completely fits within the interval

def ClumpFinding( Genome, k, L, t ):
    out = []
    for start in range(len(Genome)-L+1):
        window = Genome[start:start+L]
        counts = {}
        for i in range(len(window)-k+1):
            if window[i:i+k] not in counts:
                counts[window[i:i+k]] = 0
            counts[window[i:i+k]] += 1
        for kmer in counts:
            if counts[kmer] >= t and kmer not in out:
                out.append(kmer)
    out = " ".join(str(e) for e in out)
    print(out)

# Improved
    
def ClumpFindingImp(Genome, k, L, t):
    Freqpat = []
    clump = []
    for i in range (4**k):
        clump.append(0)
    for i in range(len(Genome)-L+1):
        text = Genome[i:i+L]
        Freqarray = ComputingFrequencies(text,k)
        for j in range(4**k):
            if int (Freqarray[j]) >= t:
                clump[j] = 1
    for i in range (4**k):
        if clump[i] == 1:
            Pattern = NumberToPattern(i,k)
            Freqpat.append(Pattern)
    return Freqpat



















