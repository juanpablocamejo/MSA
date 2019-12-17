import random
from copy import deepcopy

from common import *
from needleman_wunsch import needleman_wunsch
from profile import Profile


def msa_grasp(seqs):
    iter_max = 2 * len(seqs)
    scores = []
    solution, score = grasp_iteration(seqs, 0)
    scores.append(score)
    for _ in range(1, iter_max):
        new_sol, new_score = grasp_iteration(seqs, _)
        if new_score > score:
            solution = new_sol
            score = new_score
        scores.append(score)
    print(reduce(lambda acc, x: acc + str(x) + '|', scores, ''))
    return solution, score


def grasp_iteration(seqs, i):
    log('GRASP Iteration:', i)
    # GREEDY RANDOM PONDERADO
    log('  greedy random')
    profile, rem_seqs = initial_random_greedy_alignment(seqs)
    # ALINEAMIENTOS PROFILE-SECUENCIA
    log('  profile to sequence')
    profile = profile_sequence_alignments(profile, rem_seqs)
    # BUSQUEDA LOCAL
    s = local_search(profile.seqs)
    score = total_score(s)
    log('  score:', score)

    return s, total_score(s)


# obtiene un primer alineamiento de entre todos los pares de secuencias
# eligiendolo a través de un random ponderado por el score de dichos alineamientos
# devuelve la tupla (profile, secuenciasRestantes)
def initial_random_greedy_alignment(seqs):
    rem_seqs = seqs.copy()  # secuencias sin alinear
    i, j, res = weighted_random_choice(align_all_pairs(seqs))
    alignment, score = res
    profile = Profile(*alignment, score)
    for idx, x in enumerate([i, j]):
        rem_seqs.pop(x - idx)  # elimino las que ya alineé
    return profile, rem_seqs


# Agrega al profile todas las secuencia restantes
# tomando en cada paso la secuencia de mayor similitud con el profile
def profile_sequence_alignments(profile, rem_seqs):
    p = profile
    while rem_seqs:
        all_prof_seq_alignments = [(i, prof_seq_align(deepcopy(p), s)) for i, s in enumerate(rem_seqs)]
        i, p = sorted(all_prof_seq_alignments, key=lambda r: r[1].score, reverse=True)[0]
        rem_seqs.pop(i)
    return p


# Encuentra el óptimo local en el vecindario de la solución
def local_search(solution):
    s_scores = solution_scores(solution)
    neib = neighborhood(solution)
    cnt = 0
    while neib:
        cnt += 1
        i, swap_idx = neib.pop()
        a, b = swap_idx
        neighbor_sol = neighbor_solution(solution, i, a, b)
        new_scores = [column_score_at(neighbor_sol, i) for i in (a, b)]
        if s_scores[a] + s_scores[b] < sum(new_scores):
            solution = neighbor_sol
            s_scores[a] = new_scores[0]
            s_scores[b] = new_scores[1]
            neib = neighborhood(solution)
    log('  local search iterations:', cnt)
    return solution


# devuelve todas las soluciones vecinas con la forma (i,(a,b))
# donde 'i' es el indice de la secuencia a modificar
# 'a' y 'b' son las posiciones de la secuencia que podrian swapearse
def neighborhood(seqs):
    return [(i, a) for i, s in enumerate(seqs) for a in all_gap_swaps(s)]


# devuelve una la solucion vecina en la que la secuencia 'i'
# tiene los caracteres de indices 'a' y 'b' swapeados
def neighbor_solution(seqs, i, a, b):
    res = seqs.copy()
    tmp = res[i][a]
    res[i] = replace_char(res[i], a, res[i][b])
    res[i] = replace_char(res[i], b, tmp)
    return res


# devuelve todos los pares de indices consecutivos de una secuencia
# que contengan un swap y un caracter
def all_gap_swaps(seq):
    res = []
    for i in range(len(seq) - 1):
        a = seq[i]
        b = seq[i + 1]
        if (a == '-') ^ (b == '-'):
            res.append((i, i + 1))
    return res


# devuelve los scores por columna de un alineamiento
# pre: 'seqs' no está vacía
def solution_scores(seqs):
    return [column_score([s[i] for s in seqs]) for i in range(len(seqs[0]))]


# devuelve el score total de un alineamiento
def total_score(seqs):
    return sum(solution_scores(seqs))


# obtiene el score de la columna 'i' en la lista de secuencias alineadas
def column_score_at(seqs, i):
    return column_score([x[i] for x in seqs])


# alinea todos los posibles pares de seqs
# aplicando needleman-wunsch y los devuelve ordenados por mayor score
# el resultado es la tupla (i,j,alignment,score) donde i y j son
# los indices de las secuencias en la lista original
def align_all_pairs(seqs):
    res = []
    for i, j, s, t in all_pairs(seqs):
        res.append((i, j, needleman_wunsch(s, t)))
    return sorted(res, key=lambda x: x[2][0], reverse=True)


def weighted_random_choice(rank):
    # transforma los scores en pesos positivos
    rank_weights = [a[2][1] + abs(rank[-1][2][1]) + 1 for a in rank]
    return random.choices(rank, weights=rank_weights)[0]


# devuelve todos los pares posibles de secuencias
# el resultado es una lista de tuplas (i,j,s,t)
# donde i y j son los indices, y  s y t las secuencias
def all_pairs(seqs):
    return [(a, b, seqs[a], seqs[b]) for a in range(len(seqs)) for b in range(a + 1, len(seqs))]


# alinea un profile con una secuencia
# devuelve una instancia de Profile
def prof_seq_align(p, s):
    p.initialize_align_map(s)
    m, len_i, len_j = initialize_psa(p, s)
    for i in range(1, len_i):
        for j in range(1, len_j):
            # match/missmatch
            m[i][j] = calc_score(s[i - 1], p.get_column_probs(j - 1)) + m[i - 1][j - 1]
            p.alignMap[i][j] = ((i - 1, j - 1), p.get_column_chars(j - 1) + [s[i - 1]])
            # gap en el profile
            new_score = calc_score(s[i - 1], [1]) + m[i - 1][j]
            if new_score > m[i][j]:
                m[i][j] = new_score
                p.alignMap[i][j] = ((i - 1, j), ['-'] * len(p.seqs) + [s[i - 1]])
            # gap en la secuencia
            new_score = calc_score('-', p.get_column_probs(j - 1)) + m[i][j - 1]
            if new_score > m[i][j]:
                m[i][j] = new_score
                p.alignMap[i][j] = ((i, j - 1), p.get_column_chars(j - 1) + ['-'])
    p.build_result()
    p.score = m[-1][-1]
    return p


# calcula el score para un paso del alineamiento profile-secuencia
# ch es el caracter de la secuencia y pCol la columna de la matriz del profile
def calc_score(ch, p_col):
    return sum([p * transform(char[ch], i) for i, p in enumerate(p_col) if p > 0])


# inicializa la matriz y carga los casos base, para la alineación profile-secuencia
def initialize_psa(p, s):
    len_i, len_j = len(s) + 1, p.len() + 1
    return (
        [[init_val_psa(i, j, s[i - 1], p.get_column_probs(j - 1)) for j in range(len_j)] for i in range(len_i)], len_i, len_j)


def init_val_psa(i, j, s_char, p_col):
    if i == j == 0:
        return 0
    elif i == 0:
        return j * calc_score('-', p_col)
    elif j == 0:
        return i * transform(s_char, '-')
    else:
        return None
