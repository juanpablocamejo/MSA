from functools import reduce
from common import next_index

class NMatrix:
    def __init__(self,lenghts):
        self.lenghts = lenghts
        self.n = len(lenghts)
        self.m = self.__emptyMatN(lenghts)

    #crea una matriz n-dimensional
    def __emptyMatN(self, lenghts, i=0):
        dimLen = lenghts[i]
        if i == len(lenghts)-1:
            return [None]*(dimLen)            
        return [self.__emptyMatN(lenghts, i+1) for _ in range(dimLen)]

    def __getitem__(self,idx):
        d = self.m
        for i in range(len(idx)):
            d = d[idx[i]]
        return d        

    def __setitem__(self, idx, val):
        res = self.m
        for i in range(len(idx)-1):
            res = res[idx[i]]
        res[idx[-1]] = val

    def __repr__(self):
        start = idx = [0]*self.n
        res=''
        while True:
            res += reduce(lambda acc, c: acc + '[' + str(c) + ']', idx, '') + ' = ' + str(self.__getitem__(idx)) + '\n'
            idx = next_index(idx,self.lenghts)
            if idx == start: break
        return res