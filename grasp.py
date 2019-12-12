from functools import reduce
from operator import mul
import random
from profile import Profile
from common import *

def main():
    for s in alignNW('A','-')[1]:
        print(s)
    MSA(['GAAC', 'CAAGAC', 'GTCCA'])

def MSA(seqs):
    alignedSeqs = [None]*len(seqs)
    i,j,res = weightedRandomChoice(alignAllPairs(seqs))
    profile = Profile(*res[1])
    for idx,x in enumerate([i,j]):
        seqs.pop(x-idx)
        alignedSeqs[idx]=res[1][idx]
    print(profSeqAlign(profile,seqs[0]))

def alignAllPairs(seqs):
    return sorted([ (i,j,alignNW(s,t)) for i,j,s,t in allPairs(seqs)], key=lambda t:t[2][0],reverse=True)

def weightedRandomChoice(rank):
    rankWeights = [ a[2][0]+abs(rank[-1][2][0])+1 for a in rank ] #transforma los scores en pesos positivos
    return random.choices(rank,weights=rankWeights)[0]

def allPairs(seqs):
    return [(a,b,seqs[a],seqs[b]) for a in range(len(seqs)) for b in range(a+1,len(seqs))]

def alignNW(s, t):
    m, iMax, jMax = initialize(s, t)
    alignMap = createAlignMap(s, t)
    for i in range(1, iMax):
        for j in range(1, jMax):
            m[i][j] = transform(s[i-1], t[j-1]) + m[i-1][j-1]
            alignMap[i][j] = (i-1, j-1, s[i-1], t[j-1])

            newScore = transform(s[i-1],'-') + m[i][j-1]
            if newScore > m[i][j]:
                m[i][j] = newScore
                alignMap[i][j] = (i, j-1, s[i-1], '-')

            newScore = transform('-',t[j-1]) + m[i-1][j]
            if newScore > m[i][j]:
                m[i][j] = newScore
                alignMap[i][j] = (i-1, j, '-', t[j-1])
    return (m[-1][-1],buildResult(alignMap))

def profSeqAlign(p,s):
    m, iMax, jMax = initializePSA(p,s)
    for i in range(1, iMax):
        for j in range(1, jMax):
            m[i][j] = calcScore(s[i-1], p.getColumn(j-1)) + m[i-1][j-1]
            
            newScore = transform(s[i-1],'-') + m[i][j-1]
            if newScore > m[i][j]:
                m[i][j] = newScore

            newScore = calcScore('-',p.getColumn(j-1)) + m[i-1][j]
            if newScore > m[i][j]:
                m[i][j] = newScore
    return m[-1][-1]

def calcScore(ch,pCol):
    return sum([ p*transform(char[ch], i) for i,p in enumerate(pCol) if p>0])

def initializePSA(p,s):
    iMax, jMax = len(s)+1, p.len()+1
    return ([[initValPSA(i, j, s[i-1], p.get(j-1)) for j in range(jMax)] for i in range(iMax)], iMax, jMax)

def initValPSA(i, j, sChar, pCol):
    if i == j == 0:
        return 0
    elif i == 0:
        return j*calcScore('-', pCol)
    elif j == 0:
        return i*transform(sChar,'-')
    else:
        return None


def createAlignMap(s, t):
    iMax, jMax = len(s)+1, len(t)+1
    alignMap = [[None]*jMax for i in range(iMax)]
    for i in range(iMax):
        for j in range(jMax):
            if (i == 0 and j != 0):
                alignMap[i][j] = (i, j-1, '-', t[j-1]) 
            elif (i != 0 and j==0):
                alignMap[i][j] = (i-1, j, s[i-1], '-')
    return alignMap


def buildResult(m):
    path = []
    i = (-1, -1)
    res = ['']*2
    while i != (0, 0):
        aux = m[i[0]][i[1]]
        res[0] = aux[2]+res[0]
        res[1] = aux[3]+res[1]
        i = (aux[0], aux[1])
        path.insert(0,i)
    print('Path:',path)
    return res

def initialize(s, t):
    iMax, jMax = len(s)+1, len(t)+1
    return ([[initVal(i, j, s[i-1], t[j-1]) for j in range(jMax)] for i in range(iMax)], iMax, jMax)


def initVal(i, j, sChar, tChar):
    if i == j == 0:
        return 0
    elif i == 0:
        return j*transform('-', tChar)
    elif j == 0:
        return i*transform(sChar,'-')
    else:
        return None

if __name__ == "__main__":
    main()
