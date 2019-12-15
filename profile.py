from common import char

# Modela el resultado de una alineación de multiples secuencias
# Se inicializa con las 2 primeras secuencias alineadas y el score de dicho alineamiento
class Profile:
    def __init__(self, s, t, score):
        self.total = 2
        self.seqs = [s, t]
        self.matrix = self.__calcMatrix()
        self.alignMap = None
        self.score = score
    
    # calcula la matriz de ocurrencias
    def __calcMatrix(self):
        return [[[s[i] for s in self.seqs].count(v) for v, k in char.items()] for i in range(len(self.seqs[0]))]

    def len(self):
        return len(self.matrix)

    def getColumnProbs(self, i):
        return [c/self.total for c in self.matrix[i]]

    # inicializa la matriz que almacenará el camino del alineamiento
    def initializeAlignMap(self, newSeq):
        iMax, jMax = len(newSeq)+1, self.len()+1
        self.alignMap = [[None]*jMax for i in range(iMax)]
        for i in range(1, iMax):
            self.alignMap[i][0] = (
                (i-1, 0), ['-']*(len(self.seqs)) + [newSeq[i-1]])
        for j in range(1, jMax):
            self.alignMap[0][j] = ((0, j-1), self.getColumnChars(j-1) + ['-'])

    def getColumnChars(self, i):
        return list(map(lambda x: x[i], self.seqs))

    def buildResult(self):
        i, j = (-1, -1)
        res = ['']*(len(self.seqs)+1)
        while (i, j) != (0, 0):
            aux = self.alignMap[i][j]
            for idx, c in enumerate(aux[1]):
                res[idx] = c + res[idx]
            i, j = aux[0]
        self.seqs = res
        self.matrix = self.__calcMatrix()
        return res
