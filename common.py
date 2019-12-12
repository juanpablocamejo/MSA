from functools import reduce
from operator import mul, add


char = { '-':0, 'A':1, 'C':2, 'G':3, 'T':4 }
scoring = [[-1 if i!=j or 0 in(i,j) else 1 for i in range(5)] for j in range(5)] # matriz de scoring

def transform(a, b):
    return scoring[char[a] if isinstance(a,str) else a][char[b] if isinstance(b,str) else b]


def printMap(m):
    for k in m.keys():
        print(k, m[k])

# imprime el contenido de una matriz n-dimensional
# como una lista de indices seguidos de sus valores ej: [0][1]...[6] = None
def printM(m, lengths, skipEmpty=False):
    for i in mkAllIndex(lengths):
        v = getM(m, i)
        if not skipEmpty or v is not None:
            print(reduce(lambda acc, c: acc + '[' + str(c) + ']', i, ''), '=', v)

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

# crea una matriz n-dimensional para almacenar las computaciones
def emptyMatN(seqs, i=0):
    dimLen = len(seqs[i])+1
    if i == len(seqs)-1:
        return [None]*(dimLen)
    return [emptyMatN(seqs, i+1) for _ in range(dimLen)]

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