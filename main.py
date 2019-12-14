from dynamic_programming import msa_dp
from grasp import msa_grasp, needleman_wunsch


seqs2 =[
'AAACGGCCGGGTTT',
'ACGTACGGGTACGT']

def main():
    printAlignment('GRASP',seqs2,*msa_grasp(seqs2))
    printAlignment('Dynamic Programming', seqs2,*msa_dp(seqs2))
    #printAlignment('Needleman-Wunsch',seqs2, *needleman_wunsch(*seqs2))
    
def printAlignment(title,src, seqs,score):
    w = len(seqs[0])+2
    print('='*w)
    print(title, 'Score:', score)
    print('='*w)
    for x in seqs:
        print(x)
    print('='*w)
    for i in range(len(seqs)):
        b =seqs[i].replace('-','') in src
        if not b:
             print(src[i])
             print(seqs[i])


if __name__ == '__main__':
    main()