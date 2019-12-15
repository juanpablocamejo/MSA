from dynamic_programming import msa_dp
from grasp import msa_grasp, needleman_wunsch
from common import *
import random
import logging

logging.basicConfig(filename='app.log', filemode='w', format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)


seqs =['GCGGGTCACTGAGGGCTGGGATGAGGACGGCCACCACTTCGAGGAGTCCCTTCACTACGAGGGCAGGGCCGTGGACATCACCACGTCAGACAGGGACAAGAGCAAGTACGGCACCCTGTCCAGACTGGCGGTGGAAGCTGGGTTCGACTGGGTCTACTATGAGTCCAAAGCGCACATCCACTGCTCTGTGAAAGCAGAAAGCTCAGTCGCTGCAAAGTCGGGCGGTTGCTTCCCAGGATCCTCCACGGTCACCCTGGAAAATGGCACCCAGAGGCCCGTCAAAGATCTCCAACCCGGGGACAGAGTACTGGCCGCGGATTACGACGGAAACCCGGTTTATACCGACTTCATCATGTTCAA',
        'CTACGGCAGAAGAAGACATCCGAAAAAGCTGACACCTCTCGCCTACAAGCAGTTCATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATAACGCGCAATTCGGAGAGATTTAAAGAACTTACTCCAAATTACAATCCCGACATTATCTTTAAGGATGAGGAGAACACG',
        'CTACGGCAGAAGAAGACATCCCAAGAAGCTGACACCTCTCGCCTACAAGCAGTTTATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATCACGCGCAATTCGGAGAGATTTAAAGAACTTACTCCAAATTACAATCCCGACATTATCTTTAAGGATGAGGAGAACACT',
        'TGCTGCTGCTGGCGAGATGTCTGCTGGTGCTGCTTGTCTCCTCGCTGTTGATGTGCTCGGGGCTGGCGTGCGGACCCGGCAGGGGATTTGGCAAGAGGCGGAACCCCAAAAAGCTGACCCCTTTAGCCTACAAGCAGTTTATCCCCAACGTGGCGGAGAAGACCCTAGGGGCCAGTGGAAGATATGAGGGGAAGATCACCAGAAACTCAGAGCGATTTAAGGAACTCACCCCCAATTACAACCC']


def generateInstance(seqs, cnt, lenght):
    res = []
    for n in range(cnt):
        src = seqs[n%4]
        minStart = len(src)-lenght-1
        ini = random.randint(0,minStart)
        fin = ini+lenght
        res.append(src[ini:fin])
    return res      
instances = [generateInstance(seqs,3,30) for _ in range(2)]


def main():
    for i in instances:
        log('='*40)
        for s in i: log(' ',s)
        log('='*40)
        sol, score = msa_dp(i)
        log('='*40)
        log('Score:',score)
        for r in sol: log(' ',r)
        log('='*40)
        print('score:',score)
        for seq in sol:
            print(seq)

    
def printAlignment(title, seqs,score):
    w = len(seqs[0])+2
    print('='*w)
    print(title, 'Score:', score)
    print('='*w)
    for x in seqs:
        print(x)
    print('='*w)


if __name__ == '__main__':
    main()