from functools import reduce
from operator import mul, add
from math import log2

char = { '-':0, 'A':1, 'C':2, 'G':3, 'T':4 }
scoring = [[1 if i==j else -1 for i in range(5)] for j in range(5)] # matriz de scoring
gapPenalty = [-1]*5


def main():
    seqs = ['ATGCATGC','GCATGAAA']
    for s in align(seqs):
        print(' '+s)
    print('-------------')

# m: matriz n-dimensional para almacenar las computaciones
# n: Cant de dimensiones / secuencias
# l: Lista de longitudes de las dimensiones / len(s)+1 para cada secuencia
def align(seqs, scoreMatrix=None):
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
                    currScore = columnScore(currCol, scoreMatrix) + getM(m, currPrevIdx)
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


def printMap(m):
    for k in m.keys():
        print(k, m[k])

def writeGap(c):
    return '-' if c=='_' else c

# Arma el resultado de la alineación de acuerdo al map que contiene
# el camino almacenado y las columnas de caracteres asociada a la decisión tomada en cada paso.
def buildResult(alignMap, l):
    res = ['']*len(l)
    k = listToStr(mapAdd(-1, l))
    path =[k]
    while k != '0'*len(l):
        for i, c in enumerate(alignMap[k][1]):
            res[i] = writeGap(c) + res[i]
        k = alignMap[k][0]
        path.insert(0,k)
    print('Path:',path)
    return res

# devuelve la lista de caracteres/gaps de acuerdo a la
# opción explorada (i - j)
def calcChoiceCol(seq, i, j):
    return [seq[x][i[x]-1] if j[x] == 1 else '_' for x in range(len(i))]


def initialize(seqs, m, n, l):
    alignMap = {}
    setM(m, [0]*n, 0)
    for i in range(0, n):
        for j in range(1, l[i]):
            idx = [0]*n
            idx[i] = j
            prevIdx = idx.copy()
            prevIdx[i] = j-1
            col = ['_' if x == 0 else seqs[i][x-1] for x in idx]
            score = getM(m, prevIdx) + columnScore(col)
            setM(m, idx, score)
            alignMap[listToStr(idx)] = (listToStr(prevIdx), col)
    return alignMap


def columnScore(col, m=None):
    return sum([charPairScore(col[a], col[b], m) for a in range(len(col)) for b in range(a+1, len(col))])

def transform(a, b):
    return scoring[char[a]][char[b]]

def gap(c):    
    return gapPenalty[char[writeGap(c)]]

#NOTA: uso '_' para diferenciar la accion de insertar un gap de la de comparar un gap que ya contengan las secuencias
def charPairScore(a, b, m=None):
    if a=='_':
        return gap(b)
    elif b=='_':
        return gap(a)
    else:
        return transform(a,b)

# crea una matriz n-dimensional para almacenar las computaciones
def emptyMatN(seqs, i=0):
    dimLen = len(seqs[i])+1
    if i == len(seqs)-1:
        return [None]*(dimLen)
    return [emptyMatN(seqs, i+1) for _ in range(dimLen)]

# suma n a cada elemento de la lista xs
def mapAdd(n, xs):
    return list(map(lambda x: x+n, xs))

# convierte una lista en un string de todos los elementos concatenados
def listToStr(idx):
    return ''.join(map(str, idx))

# resta 2 listas elemento a elemento, omitiendo los resultados negativos
# pre: las listas tienen el mismo tamaño
def mapSub(a, b):
    res = []
    for i in range(len(a)):
        res.append(0 if a[i] == 0 else a[i]-b[i])
    return res

# devuelve todas las listas de naturales posibles
# que tengan el mismo tamaño que 'lengths' y cuyos digitos
# para cada posición 'i' esten en el intervalo [0, lengths[i])
def mkAllIndex(lenghts):
    res = [[0]*len(lenghts)]
    for _ in range(1, reduce(mul, lenghts, 1)):
        res.append(nextIndex(res[-1], lenghts))
    return res


def nextIndex(prev, lenghts):
    n = len(lenghts)
    aux = prev.copy()
    for j in range(n-1, -1, -1):
        if aux[j] == lenghts[j]-1:
            aux[j] = 0
        else:
            aux[j] += 1
            break
    return aux


def mkAllBinIndex(length, skip=0):
    return mkAllIndex([2]*length)[skip:]

# Obtiene un valor de la matriz n-dimensional 'm'
# usando los valores de 'idx' como indices para cada dimensión
def getM(m, idx):
    d = m
    for i in range(len(idx)):
        d = d[idx[i]]
    return d

# Guarda un valor en la matriz n-dimensional 'm'
# en la posición determinada por los indices de 'idx'
def setM(m, idx, val):
    res = m
    for i in range(len(idx)-1):
        res = res[idx[i]]
    res[idx[-1]] = val


cons = lambda v: lambda:v
# imprime el contenido de la matriz n-dimensional
# como una lista de indices sequidos de sus valores ej: [0][1]...[6] = None
def printM(m, lengths, skipEmpty=False):
    for i in mkAllIndex(lengths):
        v = getM(m, i)
        if not skipEmpty or v is not None:
            print(reduce(lambda acc, c: acc + '[' + str(c) + ']', i, ''), '=', v)


if __name__ == "__main__":
    main()
