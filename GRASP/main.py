from functools import reduce
from operator import mul

char = { '-':0, 'A':1, 'C':2, 'G':3, 'T':4 }
scoring = [[1 if i==j else -1 for i in range(5)] for j in range(5)] # matriz de scoring
gapPenalty = [-1]*5

def allPairs(seqs):
    return [(seqs[a],seqs[b]) for a in range(len(seqs)) for b in range(a+1,len(seqs))]

def gap(c):
    return gapPenalty[char[c]]

def main():
    print(alignNW('ATGCATGC','GCATGAAA'))
    MSA(['GAAC', 'CAAGAC', 'GTCCA'])
    
def profile(alignment):
    return None

def MSA(seqs):
    rank = sorted([ alignNW(*p) for p in allPairs(seqs)], key=lambda t:t[0],reverse=True)
    print(rank)


def alignNW(s, t):
    m, iMax, jMax = initialize(s, t)
    alignMap = createAlignMap(s, t)
    for i in range(1, iMax):
        for j in range(1, jMax):
            m[i][j] = transform(s[i-1], t[j-1]) + m[i-1][j-1]
            alignMap[i][j] = (i-1, j-1, s[i-1], t[j-1])

            newScore = gap(s[i-1]) + m[i][j-1]
            if newScore > m[i][j]:
                m[i][j] = newScore
                alignMap[i][j] = (i, j-1, s[i-1], '-')

            newScore = gap(t[j-1]) + m[i-1][j]
            if newScore > m[i][j]:
                m[i][j] = newScore
                alignMap[i][j] = (i-1, j, '-', t[j-1])
    return (m[-1][-1],buildResult(alignMap))


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


# resta 2 listas elemento a elemento, omitiendo los resultados negativos
# pre: las listas tienen el mismo tamaño
def mapSub(a, b):
    res = []
    for i in range(len(a)):
        res.append(0 if a[i] == 0 else a[i]-b[i])
    return res

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

def transform(a, b):
    return scoring[char[a]][char[b]]


def initialize(s, t):
    iMax, jMax = len(s)+1, len(t)+1
    return ([[initVal(i, j, s[i-1], t[j-1]) for j in range(jMax)] for i in range(iMax)], iMax, jMax)


def initVal(i, j, sChar, tChar):
    if i == j == 0:
        return 0
    elif i == 0:
        return j*gap(tChar)
    elif j == 0:
        return i*gap(sChar)
    else:
        return None

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

# imprime el contenido de la matriz n-dimensional
# como una lista de indices sequidos de sus valores ej: [0][1]...[6] = None
def printM(m, lengths):
    for i in mkAllIndex(lengths):
        print(reduce(lambda acc, c: acc +
                     '[' + str(c) + ']', i, ''), '=', getM(m, i))

# Obtiene un valor de la matriz n-dimensional 'm'
# usando los valores de 'idx' como indices para cada dimensión
def getM(m, idx):
    d = m
    for i in range(len(idx)):
        d = d[idx[i]]
    return d
if __name__ == "__main__":
    main()
