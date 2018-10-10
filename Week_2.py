# coding=utf-8

# Where in the genome does the replication begin 


# Skew diagram from traversing the genome, keeping a running total of the difference between the counts of G and C. If this difference starts increasing, then we guess that we are on the forward half-strand; on the other hand, if this difference starts decreasing, then the we guess that we are on the reverse half-strand.

def skew(text):
    string = " "
    list = []
    list.append(0)
    for i in range(1,len(text)+1):
        if text[i-1] == 'G':
            k = list[i-1] + 1
            list.append(k)
        elif text[i-1] == 'C':
            j = list[i-1] - 1
            list.append(j)
        else:
             l = list[i-1]
             list.append(l)
    string = " ".join(str(e) for e in list)
    return list


# The skew is decreasing along the reverse half-strand and increasing along the forward half-strand. Thus, the skew should achieve a minimum at the position where the reverse half-strand ends and the forward half-strand begins, which is exactly the location of ori.

def MinimumSkew(text):
    l = []
    s = []
    l = skew(text)
    k = min(skew(text))
    for i in range(len(l)):
        if l[i] == k:
            s.append(i)
    #s = " ".join(str(e) for e in s)
    return s


# Number of mismatches between strings is called hamming distance 
    
def HammingDistance (p, q):
    count = 0
    for i in range (len (p)):
        if not p[i] == q[i]:
            count += 1

    return count


# Algorithm to find starting positions of k-mer Pattern which appears as a substring of Text with at most d mismatches if there is some k-mer substring Pattern of Text having d or fewer mismatches with Pattern

def ApproxPatternmatching(p,q,d):
    t = [ ]
    for i in range(len(q)-len(p)+1):
        if HammingDistance(q[i:i+len(p)],p) <= d:
            t.append(i)
    #t = " ".join(str(e) for e in t)
    return(t)


# Algorithm to find total number of occurrences of Pattern in Text with at most d mismatches to find DnaA boxes by identifying frequent k-mers, possibly with mismatches
    
def ApproximatePatternCount(p,q,t):
    count = 0
    for i in range(len(q)-len(p)+1):
        pn = q[i:i+len(p)]
        if HammingDistance(p,pn) <= t:
            count += 1
    return(count)


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


def SymbolToNumber (Pattern):
    if Pattern == "A":
        return 0
    elif Pattern == "C":
        return 1
    elif Pattern == "G":
        return 2
    elif Pattern == "T":
        return 3


def PatternToNumber (Pattern):
    k = len (Pattern)
    if Pattern == "A" or Pattern == "C" or Pattern == "G" or Pattern == "T":
        return SymbolToNumber (Pattern)
    else:
        Prefix = Pattern[:k - 1]
        Symbol = Pattern[k - 1]
        return 4 * PatternToNumber (Prefix) + SymbolToNumber (Symbol)


# Algorithm to find most frequent k-mer with up to d mismatches
        
def Freqmermismatch (p, k, d):
    freq = [ ]
    mk = [ ]
    for i in range(4**k):
        freq.append(0)
    for i in range(4**k):
        for j in range(len(p)-k+1):
            q = p[j:j+k]
            freq[i] += ApproximatePatternCount(NumberToPattern(i,k),q,d)
    m = max(freq)
    for i in range(len(freq)):
        if m == freq[i]:
            mk.append(NumberToPattern(i,k))
    #mk = " ".join(e for e in mk)
    return mk


def ReverseComplement (Pattern):
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


# Algorithm to find most frequent k-mer including respective reverse compliment with up to d mismatches

def Freqmismatchrc(p,k,d):
    freq1 = [ ]
    mk = [ ]
    freq = [ ]
    for i in range(4**k):
        freq1.append(0)
    for i in range(4**k):
        for j in range(len(p)-k+1):
            q = p[j:j+k]
            freq1[i] += ApproximatePatternCount(NumberToPattern(i,k),q,d)
    freq2 = [ ]
    for i in range(4**k):
        freq2.append(0)
    for i in range(4**k):
        for j in range(len(p)-k+1):
            q2 = p[j:j+k]
            freq2[i] += ApproximatePatternCount(ReverseComplement(NumberToPattern(i,k)),q2,d)

    freq = map (sum, zip (freq1, freq2))
    freq = list(map(int,freq))
    
    m = max(freq)
    for i in range(len(freq)):
        if freq[i] == m:
            mk.append(NumberToPattern(i,k))
    mk = " ".join(e for e in mk)
    return(mk)



