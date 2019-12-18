from common import char


# Modela el resultado de una alineación de multiples secuencias
# Se inicializa con las 2 primeras secuencias alineadas y el score de dicho alineamiento
class Profile:
    def __init__(self, s, t, score):
        self.total = 2
        self.seqs = [s, t]
        self.matrix = self.__calc_matrix()
        self.alignMap = None
        self.score = score

    # calcula la matriz de ocurrencias
    def __calc_matrix(self):
        return [[[s[i] for s in self.seqs].count(v) for v, k in char.items()] for i in range(len(self.seqs[0]))]

    def len(self):
        return len(self.matrix)

    def get_column_probs(self, i):
        return [c / self.total for c in self.matrix[i]]

    # inicializa la matriz que almacenará el camino del alineamiento
    def initialize_align_map(self, new_seq):
        i_max, j_max = len(new_seq) + 1, self.len() + 1
        self.alignMap = [[None] * j_max for i in range(i_max)]
        for i in range(1, i_max):
            self.alignMap[i][0] = (
                (i - 1, 0), ['-'] * (len(self.seqs)) + [new_seq[i - 1]])
        for j in range(1, j_max):
            self.alignMap[0][j] = ((0, j - 1), self.get_column_chars(j - 1) + ['-'])

    def get_column_chars(self, i):
        return list(map(lambda x: x[i], self.seqs))

    def build_result(self):
        i, j = (-1, -1)
        res = [''] * (len(self.seqs) + 1)
        while (i, j) != (0, 0):
            aux = self.alignMap[i][j]
            for idx, c in enumerate(aux[1]):
                res[idx] = c + res[idx]
            i, j = aux[0]
        self.seqs = res
        self.matrix = self.__calc_matrix()
        return res
