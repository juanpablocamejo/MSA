from dynamic_programming import msa_dp
from grasp import msa_grasp, needleman_wunsch
from common import *
import random
import logging

logging.basicConfig(filename='app.log', format='%(asctime)s - %(message)s', level=logging.INFO)


seqs =['GCGGGTCACTGAGGGCTGGGATGAGGACGGCCACCACTTCGAGGAGTCCCTTCACTACGAGGGCAGGGCCGTGGACATCACCACGTCAGACAGGGACAAGAGCAAGTACGGCACCCTGTCCAGACTGGCGGTGGAAGCTGGGTTCGACTGGGTCTACTATGAGTCCAAAGCGCACATCCACTGCTCTGTGAAAGCAGAAAGCTCAGTCGCTGCAAAGTCGGGCGGTTGCTTCCCAGGATCCTCCACGGTCACCCTGGAAAATGGCACCCAGAGGCCCGTCAAAGATCTCCAACCCGGGGACAGAGTACTGGCCGCGGATTACGACGGAAACCCGGTTTATACCGACTTCATCATGTTCAA',
        'CTACGGCAGAAGAAGACATCCGAAAAAGCTGACACCTCTCGCCTACAAGCAGTTCATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATAACGCGCAATTCGGAGAGATTTAAAGAACTTACTCCAAATTACAATCCCGACATTATCTTTAAGGATGAGGAGAACACG',
        'CTACGGCAGAAGAAGACATCCCAAGAAGCTGACACCTCTCGCCTACAAGCAGTTTATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATCACGCGCAATTCGGAGAGATTTAAAGAACTTACTCCAAATTACAATCCCGACATTATCTTTAAGGATGAGGAGAACACT',
        'TGCTGCTGCTGGCGAGATGTCTGCTGGTGCTGCTTGTCTCCTCGCTGTTGATGTGCTCGGGGCTGGCGTGCGGACCCGGCAGGGGATTTGGCAAGAGGCGGAACCCCAAAAAGCTGACCCCTTTAGCCTACAAGCAGTTTATCCCCAACGTGGCGGAGAAGACCCTAGGGGCCAGTGGAAGATATGAGGGGAAGATCACCAGAAACTCAGAGCGATTTAAGGAACTCACCCCCAATTACAACCC']


def generateInstance(seqs, cnt, lenght):
    res = []
    for n in range(cnt):
        src = seqs[n%len(seqs)]
        minStart = len(src)-lenght-1
        ini = random.randint(0,minStart)
        fin = ini+lenght
        res.append(src[ini:fin])
    return res      
instances = [generateInstance(seqs,4,60) for _ in range(10)]


def main():
    alpha = 'abcdefghijklmnñopqrstuvwxyz'
    for j,i in enumerate(instances):
        log('='*60)
        log('Instance:',alpha[j])
        log('='*60)
        for s in i: log(' ',s)
        log('='*60)
        sol, score = msa_grasp(i)
        log('='*60)
        log('Score:',score)
        for r in sol: log(' ',r)
        log('='*60)


if __name__ == '__main__':
    main()