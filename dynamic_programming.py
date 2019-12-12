from common import *

def main():
    seqs = ['A','A','A','A']
    for s in align(seqs):
        print(' '+s)
    print('-------------')

# m: matriz n-dimensional para almacenar las computaciones
# n: Cant de dimensiones / secuencias
# l: Lista de longitudes de las dimensiones / len(s)+1 para cada secuencia
def align(seqs):
    m = emptyMatN(seqs)
    n = len(seqs)
    l = mapAdd(1, map(len, seqs))
    
    alignMap = initialize(seqs, m, n, l) #Casos base
    print('-------------')
    print('Casos base')
    print('-------------')
    printM(m,l,True)
    print('-------------')
    i = startIdx = [0]*n
    while True:
        if getM(m, i) is None:
            bestScore = prevIdx = col = None
            for j in mkAllBinIndex(n, 1):
                currPrevIdx = mapSub(i, j)
                if currPrevIdx != i:
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
    print('Matriz:')
    print('-------------')
    printM(m,l)
    print('-------------')
    print('Score:', getM(m, [-1]*n))
    print('-------------')
    return buildResult(alignMap, l)




# Arma el resultado de la alineación de acuerdo al map que contiene
# el camino almacenado y las columnas de caracteres asociada a la decisión tomada en cada paso.
def buildResult(alignMap, l):
    res = ['']*len(l)
    k = listToStr(mapAdd(-1, l))
    path =[k]
    while k != '0'*len(l):
        for i, c in enumerate(alignMap[k][1]):
            res[i] = c + res[i]
        k = alignMap[k][0]
        path.insert(0,k)
    print('Path:',path)
    return res

# devuelve la lista de caracteres/gaps de acuerdo a la
# opción explorada (i - j)
def calcChoiceCol(seqs, i, j):
    return [seqs[x][i[x]-1] if j[x] == 1 else '-' for x in range(len(i))]


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


def columnScore(col):
    return sum([transform(col[a], col[b]) for a in range(len(col)) for b in range(a+1, len(col))])




if __name__ == "__main__":
    main()