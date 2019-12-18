from common import *


# alinea 2 secuencias aplicando el algoritmo de Needlman-Wunsch
# devuelve la tupla (score, alineamiento)
def needleman_wunsch(s, t):
    m, i_max, j_max = initialize(s, t)
    align_map = create_align_map(s, t)
    for i in range(1, i_max):
        for j in range(1, j_max):
            # match/missmatch
            m[i][j] = transform(s[i - 1], t[j - 1]) + m[i - 1][j - 1]
            align_map[i][j] = (i - 1, j - 1, s[i - 1], t[j - 1])
            # gap en 't'
            new_score = transform(s[i - 1], '-') + m[i - 1][j]
            if new_score > m[i][j]:
                m[i][j] = new_score
                align_map[i][j] = (i - 1, j, s[i - 1], '-')
            # gap en 's'
            new_score = transform('-', t[j - 1]) + m[i][j - 1]
            if new_score > m[i][j]:
                m[i][j] = new_score
                align_map[i][j] = (i, j - 1, '-', t[j - 1])
    return build_result(align_map), m[-1][-1]


# inicializa una matriz para almacenar el camino
# del alineamiento entre 2 secuencias
def create_align_map(s, t):
    i_max, j_max = len(s) + 1, len(t) + 1
    align_map = [[None] * j_max for i in range(i_max)]
    for i in range(1, i_max):
        align_map[i][0] = (i - 1, 0, s[i - 1], '-')
    for j in range(1, j_max):
        align_map[0][j] = (0, j - 1, '-', t[j - 1])
    return align_map


# arma el resultado del alineamiento entre 2 secuencias a partir de la matriz (alignMap)
# que almacena el camino y la columna de caracteres en cada paso.
def build_result(m):
    i = (-1, -1)
    res = [''] * 2
    while i != (0, 0):
        aux = m[i[0]][i[1]]
        res[0] = aux[2] + res[0]
        res[1] = aux[3] + res[1]
        i = (aux[0], aux[1])
    return res


# inicializa la matriz con los casos base para el alineamiento entre 2 secuencias
def initialize(s, t):
    i_max, j_max = len(s) + 1, len(t) + 1
    return [[init_val(i, j, s[i - 1], t[j - 1]) for j in range(j_max)] for i in range(i_max)], i_max, j_max


def init_val(i, j, s_char, t_char):
    if i == j == 0:
        return 0
    elif i == 0:
        return j * transform('-', t_char)
    elif j == 0:
        return i * transform(s_char, '-')
    else:
        return None
