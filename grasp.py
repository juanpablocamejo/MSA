from functools import reduce
from operator import mul
import random
from profile import Profile
from copy import deepcopy
from common import *


def main():
    print( MSA(['GAAC', 'CAAGAC', 'GTCCA']))


def MSA(seqs):
    #GREEDY RANDOM PONDERADO
    profile, remSeqs = initial_random_greedy_alignment(seqs)
    #ALINEAMIENTOS PROFILE-SECUENCIA
    profile = profile_sequence_alignments(profile,remSeqs)
    #BUSQUEDA LOCAL
    return localSearch(profile.seqs)


def initial_random_greedy_alignment(seqs):
    remSeqs = seqs.copy()  # secuencias sin alinear
    alignedSeqs = [None]*len(seqs)
    i, j, res = weightedRandomChoice(alignAllPairs(seqs))
    score, alignment = res
    profile = Profile(*alignment, score)
    for idx, x in enumerate([i, j]):
        remSeqs.pop(x-idx)  # elimino las que ya alineé
        alignedSeqs[idx] = alignment[idx]
    return (profile,remSeqs)

def profile_sequence_alignments(profile,remSeqs):
    while remSeqs:
        i, p = sorted([(i, profSeqAlign(deepcopy(profile), s)) for i, s in enumerate(
            remSeqs)], key=lambda r: r[1].score, reverse=True)[0]
        remSeqs.pop(i)
    return p

def localSearch(solution):
    colScores = [columnScore([s[i] for s in solution]) for i in range(len(solution[0]))]
    alternatives = allSeqAlternatives(solution)
    while alternatives:
        i, swp= alternatives.pop()
        a,b = swp
        alt = swap(solution,i,a,b)
        if colScores[a]+colScores[b] < columnScoreAt(alt,a,b):
            solution = alt
            alternatives = allSeqAlternatives(solution)
    return solution
    
def allSeqAlternatives(seqs):
    return [(i, a) for i,s in enumerate(seqs) for a in alternativeSeqs(s) ]

def alternativeSeqs(seq):
    res = []
    for i in range(len(seq)-1):
        a = seq[i]
        b = seq[i+1]
        if (a == '-') ^ (b == '-'):
            res.append((i, i+1))
    return res

def columnScoreAt(alt,a,b):
    return columnScore([x[a] for x in alt]) + columnScore([x[b] for x in alt])

def swap(seqs, i, a, b):
    res = seqs.copy()
    tmp = res[i][a]
    res[i] = insertChar(res[i],a,res[i][b])
    res[i] = insertChar(res[i],b,tmp) 
    return res

def insertChar(seq,i,r):
    return seq[:i] + r + seq[i + 1:]
    


def alignAllPairs(seqs):
    return sorted([(i, j, alignNW(s, t)) for i, j, s, t in allPairs(seqs)], key=lambda t: t[2][0], reverse=True)


def weightedRandomChoice(rank):
    # transforma los scores en pesos positivos
    rankWeights = [a[2][0]+abs(rank[-1][2][0])+1 for a in rank]
    return random.choices(rank, weights=rankWeights)[0]


def allPairs(seqs):
    return [(a, b, seqs[a], seqs[b]) for a in range(len(seqs)) for b in range(a+1, len(seqs))]


def alignNW(s, t):
    m, iMax, jMax = initialize(s, t)
    alignMap = createAlignMap(s, t)
    for i in range(1, iMax):
        for j in range(1, jMax):
            m[i][j] = transform(s[i-1], t[j-1]) + m[i-1][j-1]
            alignMap[i][j] = (i-1, j-1, s[i-1], t[j-1])

            newScore = transform(s[i-1], '-') + m[i][j-1]
            if newScore > m[i][j]:
                m[i][j] = newScore
                alignMap[i][j] = (i, j-1, s[i-1], '-')

            newScore = transform('-', t[j-1]) + m[i-1][j]
            if newScore > m[i][j]:
                m[i][j] = newScore
                alignMap[i][j] = (i-1, j, '-', t[j-1])
    return (m[-1][-1], buildResult(alignMap))


def profSeqAlign(p, s):
    p.initializeAlignMap(s)
    m, iMax, jMax = initializePSA(p, s)
    for i in range(1, iMax):
        for j in range(1, jMax):
            # match/missmatch
            m[i][j] = calcScore(s[i-1], p.getColumn(j-1)) + m[i-1][j-1]
            p.alignMap[i][j] = ((i-1, j-1), p.getSeqsColumn(j-1) + [s[i-1]])
            # inserta gap en el profile
            newScore = transform(s[i-1], '-') + m[i][j-1]
            if newScore > m[i][j]:
                m[i][j] = newScore
                p.alignMap[i][j] = ((i, j-1), ['-']*len(p.seqs) + [s[i-1]])
            # inserta gap en la secuencia
            newScore = calcScore('-', p.getColumn(j-1)) + m[i-1][j]
            if newScore > m[i][j]:
                m[i][j] = newScore
                p.alignMap[i][j] = ((i-1, j), p.getSeqsColumn(j-1) + ['-'])
    p.buildResult()
    p.score = m[-1][-1]
    return p


def calcScore(ch, pCol):
    return sum([p*transform(char[ch], i) for i, p in enumerate(pCol) if p > 0])


def initializePSA(p, s):
    iMax, jMax = len(s)+1, p.len()+1
    return ([[initValPSA(i, j, s[i-1], p.getColumn(j-1)) for j in range(jMax)] for i in range(iMax)], iMax, jMax)


def initValPSA(i, j, sChar, pCol):
    if i == j == 0:
        return 0
    elif i == 0:
        return j*calcScore('-', pCol)
    elif j == 0:
        return i*transform(sChar, '-')
    else:
        return None


def createAlignMap(s, t):
    iMax, jMax = len(s)+1, len(t)+1
    alignMap = [[None]*jMax for i in range(iMax)]
    for i in range(1, iMax):
        alignMap[i][0] = (i-1, 0, s[i-1], '-')
    for j in range(1, jMax):
        alignMap[0][j] = (0, j-1, '-', t[j-1])
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
        path.insert(0, i)
    # print('Path:',path)
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
        return i*transform(sChar, '-')
    else:
        return None


if __name__ == "__main__":
    main()
