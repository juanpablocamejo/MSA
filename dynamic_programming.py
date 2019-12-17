from common import *
from nmatrix import NMatrix
from align_map import AlignMap


# m: matriz n-dimensional para almacenar las computaciones
# n: Cant de dimensiones / secuencias
# l: Lista de longitudes de las dimensiones / len(s)+1 para cada secuencia
# noinspection PyPep8
def msa_dp(seqs):
    n = len(seqs)
    l = map_add(1, map(len, seqs))
    m = NMatrix(l)
    align_map = initialize(seqs, m, n, l)  # Casos base
    i = start_idx = [0] * n
    while True:
        if m[i] is None:
            best_score = prev_idx = col = None
            for j in all_prev_index(i):
                curr_prev_idx = map_sub(i, j)
                curr_col = calc_choice_column(seqs, i, j)
                curr_score = column_score(curr_col) + m[curr_prev_idx]
                if best_score is None or best_score < curr_score:
                    best_score = curr_score
                    prev_idx = curr_prev_idx
                    col = curr_col
            align_map[i] = (prev_idx, col)
            m[i] = best_score
        i = next_index(i, l)
        if i == start_idx: break
    return align_map.build_result(), m[[-1] * n]


# Genera todos los posibles indices de la matriz
# que podrían ser la posición anterior en el camino óptimo
def all_prev_index(idx):
    return [i for i in all_bin_index(len(idx), 1) if valid_prev_idx(idx, i)]


# Verifica que el numero binario a restar
# no resulte en indices negativos
def valid_prev_idx(idx, prev):
    for i in range(len(idx)):
        if idx[i] == 0 and prev[i] == 1:
            return False
    return True


# devuelve la lista de caracteres/gaps de acuerdo a la
# opción explorada (idx - j)
def calc_choice_column(seqs, idx, j):
    return [seqs[x][idx[x] - 1] if idx[x] > 0 and j[x] == 1 else '-' for x in range(len(idx))]


# inicializa la matriz n-dimensional cargando los casos base
# noinspection PyPep8
def initialize(seqs, m, n, l):
    align_map = AlignMap(l)
    m[[0] * n] = 0
    for i in range(0, n):
        for j in range(1, l[i]):
            idx = [0] * n
            idx[i] = j
            prev_idx = idx.copy()
            prev_idx[i] = j - 1
            col = ['-' if x == 0 else seqs[i][x - 1] for x in idx]
            m[idx] = m[prev_idx] + column_score(col)
            align_map[idx] = (prev_idx, col)
    return align_map
