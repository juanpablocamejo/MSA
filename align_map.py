from nmatrix import NMatrix
from common import map_add


class AlignMap(NMatrix):
    def __init__(self, lenghts):
        super().__init__(lenghts)
        self.path = []
        self.result = []

    def build_result(self):
        idx = tuple(map_add(-1, self.lenghts))
        self.path = [idx]
        res = [''] * len(self.lenghts)
        zero = (0,) * self.n
        while tuple(idx) != zero:
            for i, c in enumerate(self[idx][1]):
                res[i] = c + res[i]
            idx = self[idx][0]
            self.path.insert(0, idx)
        return res
