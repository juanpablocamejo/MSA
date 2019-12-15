from grasp import *
from dynamic_programming import msa_dp
import pytest

seqs3 = [
    'AGAGTACTGGCCGCGGATTACGACGGAAACCCGGTTTATACCGACTTCATCATGTTCAA',
    'TTAAAGAACTTACTCCAAATTACAATCCCGACATTATCTTTAAGGATGAGGAGAACACG',
    'GAGGGGAAGATCACCAGAAACTCAGAGCGATTTAAGGAACTCACCCCCAATTACAACCC'
]

seqs2 = [
    'AGAGTACTGGCCGCGGATTACGACGGAAACCCGGTTTATACCGACTTCATCATGTTCAA',
    'TTAAAGAACTTACTCCAAATTACAATCCCGACATTATCTTTAAGGATGAGGAGAACACG'
    ]

#verifico que no se pierdan caracteres de las secuencias originales
#y que las secuencias del resultado tenga el mismo length
def test_msa_dp_valid_alignment():
    res,_ = msa_dp(seqs3)
    for s in res:
        assert len(s) == len(res[0])
    for s in res:
        assert s.replace('-', '') in seqs3

def test_msa_grasp_valid_alignment():
    res,_ = msa_grasp(seqs3)
    for s in res:
        assert len(s) == len(res[0])
    for s in res:
        assert s.replace('-', '') in seqs3

def test_needlman_wunsch_valid_alignment():
    res,_ = needleman_wunsch(*seqs2)
    for s in res:
        assert len(s) == len(res[0])
    for s in res:
        assert s.replace('-', '') in seqs2

def test_local_search():
    sol = ['AG-','A-G']
    neig = ['AG-','AG-']
    newSol = local_search(sol)
    assert newSol == neig

def test_neighborhood():
    seqs = ['AG-','A-G']
    res = [(0,(1,2)),(1,(0,1)),(1,(1,2))]
    N =  neighborhood(seqs)
    assert len(res)==len(N)
    for n in N:
        assert n in res

def test_neighborhood_solution():
    seqs = ['AG-','A-G']
    n = neighbor_solution(seqs,0,1,2)
    assert n==['A-G','A-G']

def test_all_gap_swaps():
    s = '-A-CG-T-'
    expected = [(0,1),(1,2),(2,3),(4,5),(5,6),(6,7)]
    res = all_gap_swaps(s)
    assert len(res)==len(expected)
    for r in res:
        assert r in expected
