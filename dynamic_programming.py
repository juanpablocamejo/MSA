from common import *
import numpy as np


# m: matriz n-dimensional para almacenar las computaciones
# n: Cant de dimensiones / secuencias
# l: Lista de longitudes de las dimensiones / len(s)+1 para cada secuencia
def msa_dp(seqs):
    m = emptyMatN(seqs)
    n = len(seqs)
    l = mapAdd(1, map(len, seqs))    
    alignMap = initialize(seqs, m, n, l) #Casos base
    i = startIdx = [0]*n
    while True:
        if getM(m, i) is None:
            bestScore = prevIdx = col = None
            for j in allPrevIndex(i):
                currPrevIdx = mapSub(i, j)
                currCol = calcChoiceCol(seqs, i, j)
                currScore = columnScore(currCol) + getM(m, currPrevIdx)
                if bestScore is None or bestScore < currScore:
                    bestScore = currScore
                    prevIdx = currPrevIdx
                    col = currCol
            alignMap[listToStr(i)] = (listToStr(prevIdx), col)
            setM(m, i, bestScore)
        i = nextIndex(i, l)
        if i == startIdx:
            break
    return (buildResult(alignMap, l),getM(m,[-1]*n))


def allPrevIndex(idx):
    return [ i for i in mkAllBinIndex(len(idx), 1) if validPrevIdx(idx,i)]

def validPrevIdx(idx,prev):
    for i in range(len(idx)):
        if idx[i]==0 and prev[i]==1:
            return False
    return True


# Arma el resultado de la alineación de acuerdo al map que contiene
# el camino almacenado y la columna de caracteres asociada a la decisión tomada en cada paso.
def buildResult(alignMap, l):
    k = listToStr(mapAdd(-1, l))
    path = [k]
    res = ['']*len(l)
    while k != '0'*len(l):
        for i, c in enumerate(alignMap[k][1]):
            res[i] = c + res[i]
        k = alignMap[k][0]
        path.insert(0, k)
    print('path_DP:',path)
    return res

# devuelve la lista de caracteres/gaps de acuerdo a la
# opción explorada (idx - j)
def calcChoiceCol(seqs, idx, j):
    return [seqs[x][idx[x]-1] if idx[x]>0 and j[x] == 1 else '-' for x in range(len(idx))]

def initialize(seqs, m, n, l):
    alignMap = {}
    setM(m, [0]*n, 0)
    for i in range(0, n):
        for j in range(1, l[i]):
            idx = [0]*n
            idx[i] = j
            prevIdx = idx.copy()
            prevIdx[i] = j-1
            col = ['-' if x == 0 else seqs[i][x-1] for x in idx]
            score = getM(m, prevIdx) + columnScore(col)
            setM(m, idx, score)
            alignMap[listToStr(idx)] = (listToStr(prevIdx), col)
    return alignMap