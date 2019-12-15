from common import *
from nmatrix import NMatrix
from align_map import AlignMap

# m: matriz n-dimensional para almacenar las computaciones
# n: Cant de dimensiones / secuencias
# l: Lista de longitudes de las dimensiones / len(s)+1 para cada secuencia
def msa_dp(seqs):
    n = len(seqs)
    l = map_add(1, map(len, seqs))    
    m = NMatrix(l)
    alignMap = initialize(seqs, m, n, l) #Casos base
    i = startIdx = [0]*n
    while True:
        if m[i] is None:
            bestScore = prevIdx = col = None
            for j in all_prev_index(i):
                currPrevIdx = map_sub(i, j)
                currCol = calc_choice_column(seqs, i, j)
                currScore = column_score(currCol) + m[currPrevIdx]
                if bestScore is None or bestScore < currScore:
                    bestScore = currScore
                    prevIdx = currPrevIdx
                    col = currCol
            alignMap[i] = (prevIdx, col)
            m[i] = bestScore
        i = next_index(i, l)
        if i == startIdx:break
    return (alignMap.build_result(), m[[-1]*n])

#Genera todos los posibles indices de la matriz
#que podrían ser la posición anterior en el camino óptimo
def all_prev_index(idx):
    return [ i for i in all_bin_index(len(idx), 1) if valid_prev_idx(idx,i)]

# Verifica que el numero binario a restar
# no resulte en indices negativos
def valid_prev_idx(idx,prev):
    for i in range(len(idx)):
        if idx[i]==0 and prev[i]==1:
            return False
    return True

# devuelve la lista de caracteres/gaps de acuerdo a la
# opción explorada (idx - j)
def calc_choice_column(seqs, idx, j):
    return [seqs[x][idx[x]-1] if idx[x]>0 and j[x] == 1 else '-' for x in range(len(idx))]

#inicializa la matriz n-dimensional cargando los casos base
def initialize(seqs, m, n, l):
    alignMap = AlignMap(l)
    m[[0]*n]= 0
    for i in range(0, n):
        for j in range(1, l[i]):
            idx = [0]*n
            idx[i] = j
            prevIdx = idx.copy()
            prevIdx[i] = j-1
            col = ['-' if x == 0 else seqs[i][x-1] for x in idx]
            m[idx] = m[prevIdx] + column_score(col)
            alignMap[idx] = (prevIdx, col)
    return alignMap