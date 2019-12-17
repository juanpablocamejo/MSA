import logging
from functools import reduce
from operator import mul

char = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}
scoring = [[-1 if i != j or 0 in (i, j) else 1 for i in range(5)]
           for j in range(5)]  # matriz de scoring


def log(*args):
    logging.info(reduce(lambda acc, x: acc + str(x), args, ''))


# devuelve el score de una lista de caracteres
# calculado como la sumatoria de los scores de todos los posibles pares
def column_score(col):
    return sum([transform(col[a], col[b]) for a in range(len(col)) for b in range(a + 1, len(col))])


def transform(a, b):
    return scoring[char[a] if isinstance(a, str) else a][char[b] if isinstance(b, str) else b]


# devuelve todas las listas de naturales posibles
# que tengan el mismo tamaño que 'lengths' y cuyos digitos
# para cada posición 'i' esten en el intervalo [0, lengths[i])
def mk_all_index(lenghts):
    res = [[0] * len(lenghts)]
    for _ in range(1, reduce(mul, lenghts, 1)):
        res.append(next_index(res[-1], lenghts))
    return res


def next_index(prev, lenghts):
    n = len(lenghts)
    aux = prev.copy()
    for j in range(n - 1, -1, -1):
        if aux[j] == lenghts[j] - 1:
            aux[j] = 0
        else:
            aux[j] += 1
            break
    return aux


def all_bin_index(length, skip=0):
    return mk_all_index([2] * length)[skip:]


# suma n a cada elemento de la lista xs
def map_add(n, xs):
    return list(map(lambda x: x + n, xs))


# resta 2 listas elemento a elemento
# pre: las listas tienen el mismo tamaño
def map_sub(a, b):
    return [a[i] - b[i] for i in range(len(a))]


def replace_char(seq, i, r):
    return seq[:i] + r + seq[i + 1:]
