char = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}

class Profile:
    def __init__(self, s, t):
        self.total = 2
        self.seqs = [s, t]
        self.calcMatrix()
        self.alignMap = self.createAlignMap(s, t)

    def calcMatrix(self):
        self.matrix = [[[s[i] for s in self.seqs].count(v) for v, k in char.items()] for i in range(len(self.seqs[0]))]

    def len(self):
        return len(self.matrix)

    
      

    def getColumn(self, i):
        return [c/self.total for c in self.matrix[i]]

    

    def createAlignMap(self, s, t):
        iMax, jMax = len(s)+1, len(t)+1
        alignMap = [[None]*jMax for i in range(iMax)]
        for i in range(iMax):
            for j in range(jMax):
                if (i == 0 and j != 0):
                    alignMap[i][j] = (i, j-1, ['-', t[j-1]])
                elif (i != 0 and j == 0):
                    alignMap[i][j] = (i-1, j, [s[i-1], '-'])
        return alignMap


def buildResult(self):
    path = []
    i = (-1, -1)
    res = ['']*self.seqs
    while i != (0, 0):
        aux = self.alignMap[i[0]][i[1]]
        for i, c in enumerate(aux[1]):
            res[i] = c + res[i]
        i = aux[0]
        path.insert(0, i)
    print('Path:', path)
    return res
