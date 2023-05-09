#2021245044 소프트웨어학부 서은하
import sys
import time

# custom error class
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
                    is_Protein_seq(whitespace_tmp)
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
            is_Protein_seq(whitespace_tmp)
            sequences.append(whitespace_tmp)
    
    # check exception
    if len(sequences) == 0:
        print("No Protein sequence")
        sys.exit()
    elif len(sequences) == 1:
        print("Need one more sequence")
        sys.exit()
    elif len(sequences) > 2:
        sequences = sequences[:2]
    return sequences

# BLOSUM62 txt file to dictionary
def blosum_matrix_to_dict(file_name):
    with open(file_name, "r") as f:
        lines = f.readlines()
    # Get protein symbols
    amino_acids = lines[0].split("\t")
    # Initialize BLOSUM dictionary
    blosum_matrix = {}

    # Parse BLOSUM matrix
    for line in lines[1:]:
        tokens = line.split('\t')
        aa1 = tokens[0]
        for i, score in enumerate(tokens[1:]):
            aa2 = amino_acids[i+1]
            blosum_matrix[(aa1, aa2)] = int(score)
    return blosum_matrix

# protein sequence validation
def is_Protein_seq(seq):
    if not set(seq).issubset({'C','S','T','P','A','G','N','D','E','Q','H','R','K','M','I','L','V','F','Y','W'}):
        print("No Protein sequence")
        sys.exit()
    else:
        return 

def local_align(sequence1, sequence2, blosum, gap_penalty):
    # Initialize variables
    m, n = len(sequence1), len(sequence2)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0
    max_i, max_j = 0, 0
    
    # Fill in the score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i - 1][j - 1] + blosum[(sequence1[i - 1], sequence2[j - 1])]
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(0, match, delete, insert)
            
            # Update maximum score
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_i, max_j = i, j
    
    # Traceback alignment
    aligned_sequence1, aligned_sequence2 = '', ''
    i, j = max_i, max_j
    while score_matrix[i][j] != 0:
        current_score = score_matrix[i][j]
        diagonal_score = score_matrix[i - 1][j - 1]
        up_score = score_matrix[i - 1][j]
        left_score = score_matrix[i][j - 1]
        
        if current_score == diagonal_score + blosum[(sequence1[i - 1], sequence2[j - 1])]:
            aligned_sequence1 = sequence1[i - 1] + aligned_sequence1
            aligned_sequence2 = sequence2[j - 1] + aligned_sequence2
            i -= 1
            j -= 1
        elif current_score == up_score + gap_penalty:
            aligned_sequence1 = sequence1[i - 1] + aligned_sequence1
            aligned_sequence2 = '-' + aligned_sequence2
            i -= 1
        elif current_score == left_score + gap_penalty:
            aligned_sequence1 = '-' + aligned_sequence1
            aligned_sequence2 = sequence2[j - 1] + aligned_sequence2
            j -= 1    
    return max_score, aligned_sequence1, aligned_sequence2

if __name__ == '__main__':
    input_file = sys.argv[1]
    #input_file = "dna.txt"
    output_file = 'Assignment5_output.txt'

    #read protein sequence file
    try:
        sequences = read_fasta_file(input_file)
    except FileNotFoundError:
        print("No input file")
        sys.exit()
    except NoCorrectFormatErr:
        print("No correct format")
        sys.exit()

    # read blosum txt file
    try:
        #input_BLOSUM=sys.argv[2]
        input_BLOSUM="BLOSUM62.txt"
        blosum=blosum_matrix_to_dict(input_BLOSUM)     
    except FileNotFoundError:
        print("There is no such file : ",end='')
        print(input_BLOSUM)
        time.sleep(5)
        sys.exit()

    # start timer
    start_time = time.time()
    score,seq1,seq2= local_align(sequences[0], sequences[1], blosum, -5)

    # check elapsed time
    print("Elapsed time (microseconds): %f" % (time.time() - start_time))

    # make output file
    if len(seq1)>60:
        newline=len(seq1)//60
        i=0
        with open(output_file, 'w') as f:
            score_str="Score : %d\n" % score
            f.write(score_str)
            while newline!=0:
                f.write(seq1[i:i+60]+"\n")
                f.write(seq2[i:i+60]+"\n")
                f.write("------------------------------------------------------------\n")
                i+=60
                newline-=1
                if newline==0:
                    f.write(seq1[i:len(seq1)]+"\n")
                    f.write(seq2[i:len(seq2)])        
    else:
        with open(output_file, 'w') as f:
            score_str="Score : %d\n" % score
            f.write(score_str)
            f.write(seq1+"\n")
            f.write(seq2)
        f.close()



