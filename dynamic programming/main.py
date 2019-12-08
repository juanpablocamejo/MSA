from functools import reduce
from operator import mul, add
from math import log2

def main():
    seqs = [
     'ACGT'
    ,'AAAACGT'
    ,'TTACGT'
    ,'ACGTGGGG']
    for s in align(seqs):
        print(s)

def align(seqs):
    m = emptyMatN(seqs)
    n = len(seqs)
    l = mapAdd(1, map(len,seqs))
    alignMap = {}
    initialize(seqs, m, n, l, alignMap)
    for i in mkAllIndex(l):
        if getM(m, i) is None:
            scores = []
            prevChoices = []
            cols=[]
            for j in mkAllBinIndex(n, 1):
                idx = mapSub(i, j)
                if idx != i:
                    prevChoices.append(idx)
                    cols.append(getChoiceCol(seqs, i, j))
                    scores.append(columnScore(cols[-1]) + getM(m, idx))
            maxS = scores[0]
            prevIdx = prevChoices[0]
            col = cols[0]
            for c in range(1, len(scores)):
                if scores[c] > maxS:
                    maxS = scores[c]
                    prevIdx = prevChoices[c]
                    col = cols[c]
            alignMap[listToStr(i)] = (listToStr(prevIdx), col)
            setM(m, i, maxS)
    print(getM(m,mapAdd(-1,l)))
    return buildResult(alignMap,l)

def buildResult(alignMap,l):
    res = ['']*len(l)
    k = listToStr(mapAdd(-1,l))
    while k != '0'*len(l):
        for i,c in enumerate(alignMap[k][1]):
            res[i]=c + res[i]
        k = alignMap[k][0]
    return res

def getChoiceCol(seq, i, j):
    return [seq[x][i[x]-1] if j[x] == 1 else '-' for x in range(len(i))]


def initialize(seqs, m, n, l,am):
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
            am[listToStr(idx)]=(listToStr(prevIdx),col)

def columnScore(col):
    return sum([charCost(col[a], col[b]) for a in range(len(col)) for b in range(a+1, len(col))])


# TODO: Reemplazar con la funcion de costos real
def charCost(a, b):
    if a == b:
        return 1
    else:
        return -1
            
### FUNCIONES AUXILIARES ###

# crea una matriz n-dimensional para almacenar las computaciones
def emptyMatN(seqs, i=0):
    dimLen = len(seqs[i])+1
    if i == len(seqs)-1:
        return [None]*(dimLen)
    return [emptyMatN(seqs, i+1) for _ in range(dimLen)]
    
#suma n a cada elemento de la lista xs
def mapAdd(n,xs):
    return list(map(lambda x: x+n, xs))

#convierte una lista en un string de todos los elementos concatenados 
def listToStr(idx):
    return ''.join(map(str, idx))

#resta 2 listas elemento a elemento, omitiendo los resultados negativos
#pre: las listas tienen el mismo tamaño
def mapSub(a, b):
    res = []
    for i in range(len(a)):
        res.append(0 if a[i] == 0 else a[i]-b[i])
    return res

#devuelve todas las listas de naturales posibles
#que tengan el mismo tamaño que 'lengths' y cuyos digitos 
#para cada posición 'i' esten en el intervalo [0, lengths[i])
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

#Obtiene un valor de la matriz n-dimensional 'm'
#usando los valores de 'idx' como indices para cada dimensión
def getM(m, idx):
    d = m
    for i in range(len(idx)):
        d = d[idx[i]]
    return d

#Guarda un valor en la matriz n-dimensional 'm'
#en la posición determinada por los indices de 'idx'
def setM(m, idx, val):
    res = m
    for i in range(len(idx)-1):
        res = res[idx[i]]
    res[idx[-1]] = val

#imprime el contenido de la matriz n-dimensional
#como una lista de indices sequidos de sus valores ej: [0][1]...[6] = None
def printM(m, lengths):
    for i in mkAllIndex(lengths):
        print(reduce(lambda acc, c: acc +
                     '[' + str(c) + ']', i, ''), '=', getM(m, i))


if __name__ == "__main__":
    main()
