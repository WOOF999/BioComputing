#2021245044 소프트웨어학부 서은하
import sys
import time

# Custom error class
class NoCorrectFormatErr(Exception): 
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg
    
# read fasta format file
def read_fasta_file(filename):
    sequences = []
    comment_cnt=0
    with open(filename, 'r') as f:
        seq = ''
        for line in f:
            if line.startswith('>'):
                comment_cnt+=1
                if seq != '':
                    whitespace_tmp=seq.replace(' ',"").upper()
                    is_DNA_seq(whitespace_tmp)
                    sequences.append(whitespace_tmp)
                    seq = ''
                if len(sequences) == 2:
                    break
            # raise no correct format error
            elif comment_cnt==0 and not line.startswith('>'):
                raise NoCorrectFormatErr("No Correct Format")
            else:
                seq += line.strip()
        if seq != '':
            whitespace_tmp=seq.replace(' ',"").upper()
            is_DNA_seq(whitespace_tmp)
            sequences.append(whitespace_tmp)
    
    # check exception
    if len(sequences) == 0:
        print("No DNA sequence")
        sys.exit()
    elif len(sequences) == 1:
        print("Need one more sequence")
        sys.exit()
    elif len(sequences) > 2:
        sequences = sequences[:2]
    return sequences

# DNA sequence validation
def is_DNA_seq(seq):
    if not set(seq).issubset({'A', 'C', 'G', 'T'}):
        print("No DNA sequence")
        sys.exit()
    else:
        return 
    
#find longest common subsequence
def lcs(seq1, seq2):
    # return lcs!
    m, n = len(seq1), len(seq2)
    # initialize first column, row with 0
    grid = [[0 for j in range(n+1)] for i in range(m+1)]
    for i in range(1, m+1):
        for j in range(1, n+1):
            if seq1[i-1] == seq2[j-1]:
                # diagonal move (weight 1)
                grid[i][j] = grid[i-1][j-1] + 1
            else:
                # choose the maximum score
                max_score = 0 # initialize max_score to 0
                if grid[i-1][j] >= grid[i][j-1]:
                    max_score = grid[i-1][j] # vertical move
                else:
                    max_score = grid[i][j-1] # horizontal move
                if max_score == grid[i-1][j-1]:
                    # diagonal move 
                    if seq1[i-1] == seq2[j-1]:
                        grid[i][j] = grid[i-1][j-1] + 1
                    else:
                        grid[i][j] = grid[i-1][j-1]
                else:
                    grid[i][j] = max_score
    # backtrack to find the longest common subsequence
    i, j = m, n
    subseq = ''
    while i > 0 and j > 0:
        if seq1[i-1] == seq2[j-1]:
            subseq = seq1[i-1] + subseq
            i -= 1
            j -= 1
        elif grid[i-1][j] >= grid[i][j-1]:
            i -= 1
        else:
            j -= 1
    return subseq

if __name__ == '__main__':
    input_file = sys.argv[1]
    #input_file = "dna.txt"
    output_file = 'lcs_output.txt'
    try:
        sequences = read_fasta_file(input_file)
    except FileNotFoundError:
        print("No DNA sequence")
        sys.exit()
    except NoCorrectFormatErr:
        print("No correct format")
        sys.exit()

    # check elapsed time
    start_time = time.time()
    lcs_seq = lcs(sequences[0], sequences[1])

    with open(output_file, 'w') as f:
        f.write(lcs_seq)

    # print elapsed time
    print("Elapsed time (microseconds): %f" % (time.time() - start_time))
