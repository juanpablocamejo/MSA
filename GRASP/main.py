from functools import reduce
from operator import mul

def main():
    s, t = 'GAAC', 'CAAGAC'
    for x in align(s, t):
        print(x)

def align(s, t, scoring=None):
    m, iMax, jMax = initialize(s, t, scoring)
    alignMap = createAlignMap(s, t)
    for i in range(1, iMax):
        for j in range(1, jMax):
            m[i][j] = transform(s[i-1], t[j-1], scoring) + m[i-1][j-1]
            alignMap[i][j] = (i-1, j-1, s[i-1], t[j-1])

            newScore = transform(s[i-1], '-', scoring) + m[i][j-1]
            if newScore > m[i][j]:
                m[i][j] = newScore
                alignMap[i][j] = (i, j-1, s[i-1], '-')

            newScore = transform('-', t[j-1], scoring) + m[i-1][j]
            if newScore > m[i][j]:
                m[i][j] = newScore
                alignMap[i][j] = (i-1, j, '-', t[j-1])
    print('Score:',m[-1][-1])
    printM(m,[iMax,jMax])
    return buildResult(alignMap)


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

def buildResult(m):
    i = (-1, -1)
    res = ['']*2
    while i != (0, 0):
        aux = m[i[0]][i[1]]
        res[0] = aux[2]+res[0]
        res[1] = aux[3]+res[1]
        i = (aux[0], aux[1])
    return res


def transform(a, b, m):
    if m is None:  # default
        if a == b:
            return 1
        else:
            return -1
    else:
        return m[a][b]


def initialize(s, t, scoring):
    iMax, jMax = len(s)+1, len(t)+1
    return ([[initVal(i, j, s[i-1], t[j-1], scoring) for j in range(jMax)] for i in range(iMax)], iMax, jMax)


def initVal(i, j, sChar, tChar, scoring):
    if i == j == 0:
        return 0
    elif i == 0:
        return j*transform('-', tChar, scoring)
    elif j == 0:
        return i*transform(sChar, '-', scoring)
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
