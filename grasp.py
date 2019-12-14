from functools import reduce
from operator import mul
import random
from profile import Profile
from copy import deepcopy
from common import *
from dynamic_programming import msa_dp

def msa_grasp(seqs):
    iterMax = len(seqs)*20
    solution, score = grasp_iteration(seqs)
    for _ in range(1, iterMax):
        newSol, newScore = grasp_iteration(seqs)
        if newScore > score:
            solution = newSol
            score = newScore
            print(_,newScore)
    return (solution, score)


def grasp_iteration(seqs):
    # GREEDY RANDOM PONDERADO
    profile, remSeqs = initial_random_greedy_alignment(seqs)
    # ALINEAMIENTOS PROFILE-SECUENCIA
    profile = profile_sequence_alignments(profile, remSeqs)
    # BUSQUEDA LOCAL
    s = local_search(profile.seqs)
    return (s, total_score(s))


# obtiene un primer alineamiento de entre todos los pares de secuencias
# eligiendolo a través de un random ponderado por el score de dichos alineamientos
# devuelve la tupla (profile, secuenciasRestantes)
def initial_random_greedy_alignment(seqs):
    remSeqs = seqs.copy()  # secuencias sin alinear
    i, j, res = weighted_random_choice(alignAllPairs(seqs))
    alignment,score = res
    profile = Profile(*alignment, score)
    for idx, x in enumerate([i, j]):
        remSeqs.pop(x-idx)  # elimino las que ya alineé
    return (profile, remSeqs)

# Agrega al profile todas las secuencia restantes
# tomando en cada paso la secuencia de mayor similitud con el profile
def profile_sequence_alignments(profile, remSeqs):
    p = profile
    while remSeqs:
        allProfSeqAlignments = [(i, prof_seq_align(deepcopy(p), s)) for i, s in enumerate(remSeqs)]
        i, p = sorted(allProfSeqAlignments, key=lambda r: r[1].score, reverse=True)[0]
        remSeqs.pop(i)
    return p

# Encuentra el óptimo local en el vecindario de la solución
def local_search(solution):
    s_scores = solution_scores(solution)
    N = neighborhood(solution)
    while N:
        i, swapIdx = N.pop()
        a, b = swapIdx
        neighbor_sol = neighbor_solution(solution, i, a, b)
        newScores = [columnScoreAt(neighbor_sol, i) for i in (a, b)]
        if s_scores[a]+s_scores[b] < sum(newScores):
            solution = neighbor_sol
            s_scores[a] = newScores[0]
            s_scores[b] = newScores[1]
            N = neighborhood(solution)
    return solution

# devuelve todas las soluciones vecinas con la forma (i,(a,b))
# donde 'i' es el indice de la secuencia a modificar
# 'a' y 'b' son las posiciones de la secuencia que podrian swapearse
def neighborhood(seqs):
    return [(i, a) for i, s in enumerate(seqs) for a in allGapSwaps(s)]

    
# devuelve una la solucion vecina en la que la secuencia 'i'
# tiene los caracteres de indices 'a' y 'b' swapeados
def neighbor_solution(seqs, i, a, b):
    res = seqs.copy()
    tmp = res[i][a]
    res[i] = insertChar(res[i], a, res[i][b])
    res[i] = insertChar(res[i], b, tmp)
    return res

# devuelve todos los pares de indices consecutivos de una secuencia
# que contengan un swap y un caracter
def allGapSwaps(seq):
    res = []
    for i in range(len(seq)-1):
        a = seq[i]
        b = seq[i+1]
        if (a == '-') ^ (b == '-'):
            res.append((i, i+1))
    return res

# devuelve los scores por columna de un alineamiento
# pre: 'seqs' no está vacía
def solution_scores(seqs):
    return [columnScore([s[i] for s in seqs]) for i in range(len(seqs[0]))]

#devuelve el score total de un alineamiento
def total_score(seqs):
    return sum(solution_scores(seqs))

# obtiene el score de la columna 'i' en la lista de secuencias alineadas
def columnScoreAt(seqs, i):
    return columnScore([x[i] for x in seqs])

def insertChar(seq, i, r):
    return seq[:i] + r + seq[i + 1:]

def alignAllPairs(seqs):
    res = []
    for i, j, s, t in allPairs(seqs):
        res.append((i, j, needleman_wunsch(s, t)))
    return sorted(res, key=lambda t: t[2][0], reverse=True)

def weighted_random_choice(rank):
    # transforma los scores en pesos positivos
    rankWeights = [a[2][1]+abs(rank[-1][2][1])+1 for a in rank]
    return random.choices(rank, weights=rankWeights)[0]

def allPairs(seqs):
    return [(a, b, seqs[a], seqs[b]) for a in range(len(seqs)) for b in range(a+1, len(seqs))]

# alinea 2 secuencias aplicando el algoritmo de Needlman-Wunsch
# devuelve la tupla (score, alineamiento)
def needleman_wunsch(s, t):
    m, iMax, jMax = initialize(s, t)
    alignMap = createAlignMap(s, t)
    for i in range(1, iMax):
        for j in range(1, jMax):
            # match/missmatch
            m[i][j] = transform(s[i-1], t[j-1]) + m[i-1][j-1]
            alignMap[i][j] = (i-1, j-1, s[i-1], t[j-1])
            # gap en 't'
            newScore = transform(s[i-1], '-') + m[i-1][j]
            if newScore > m[i][j]:
                m[i][j] = newScore
                alignMap[i][j] = (i-1, j, s[i-1], '-')
            # gap en 's'
            newScore = transform('-', t[j-1]) + m[i][j-1]
            if newScore > m[i][j]:
                m[i][j] = newScore
                alignMap[i][j] = (i, j-1, '-', t[j-1])
    return ( buildResult(alignMap), m[-1][-1])

#alinea un profile con una secuencia
#devuelve una instancia de Profile
def prof_seq_align(p, s):
    p.initializeAlignMap(s)
    m, iMax, jMax = initializePSA(p, s)
    for i in range(1, iMax):
        for j in range(1, jMax):
            # match/missmatch
            m[i][j] = calcScore(s[i-1], p.getColumn(j-1)) + m[i-1][j-1]
            p.alignMap[i][j] = ((i-1, j-1), p.getSeqsColumn(j-1) + [s[i-1]])
            # gap en el profile
            newScore = calcScore(s[i-1], [1]) + m[i-1][j]
            if newScore > m[i][j]:
                m[i][j] = newScore
                p.alignMap[i][j] = ((i-1, j), ['-']*len(p.seqs) + [s[i-1]])
            # gap en la secuencia
            newScore = calcScore('-', p.getColumn(j-1)) + m[i][j-1]
            if newScore > m[i][j]:
                m[i][j] = newScore
                p.alignMap[i][j] = ((i, j-1), p.getSeqsColumn(j-1) + ['-'])
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
        return i*transform(sChar, '-')
    else:
        return None


seqs2 =[
'CCGGTTTATACCGACTTCATCATGTTCAA',
'ACATTATCTTTAAGGATGAGGAGAACACG'#,'TTTAAGGAACTCACCCCCAATTACAACCC'
]

def main():
    for i in seqs2:
        print(i)
    r,s =msa_grasp(seqs2)
    for j in r:
        print(j)

    
def printAlignment(title, seqs,score):
    w = len(seqs[0])+2
    print('='*w)
    print(title, 'Score:', score)
    print('='*w)
    for x in seqs:
        print(x)
    print('='*w)


if __name__ == '__main__':
    main()