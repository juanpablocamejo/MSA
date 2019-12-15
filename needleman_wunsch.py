from common import *

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

#inicializa una matriz para almacenar el camino
#del alineamiento entre 2 secuencias
def createAlignMap(s, t):
    iMax, jMax = len(s)+1, len(t)+1
    alignMap = [[None]*jMax for i in range(iMax)]
    for i in range(1, iMax):
        alignMap[i][0] = (i-1, 0, s[i-1], '-')
    for j in range(1, jMax):
        alignMap[0][j] = (0, j-1, '-', t[j-1])
    return alignMap

#arma el resultado del alineamiento entre 2 secuencias a partir de la matriz (alignMap)
#que almacena el camino y la columna de caracteres en cada paso.
def buildResult(m):
    i = (-1, -1)
    res = ['']*2
    while i != (0, 0):
        aux = m[i[0]][i[1]]
        res[0] = aux[2]+res[0]
        res[1] = aux[3]+res[1]
        i = (aux[0], aux[1])
    return res

#inicializa la matriz con los casos base para el alineamiento entre 2 secuencias
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