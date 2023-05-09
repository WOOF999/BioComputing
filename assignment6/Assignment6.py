#2021245044 소프트웨어학부 서은하
import sys
import time
import re
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
    pattern=re.compile('[A-Z]+')
    match=pattern.match(seq)
    if not match:
        print("No Protein sequence")
        sys.exit()
    else:
        return 

def global_align(sequence1, sequence2, blosum, gap_penalty):
    # make dp score matrix
    m = len(sequence1)
    n = len(sequence2)

    dp_matrix = [[0 for j in range(n+1)] for i in range(m+1)]
    
    # initialize
    for i in range(1, m+1):
        dp_matrix[i][0] = i * gap_penalty
    for j in range(1, n+1):
        dp_matrix[0][j] = j * gap_penalty
        
    # score dp matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            diagonal = dp_matrix[i-1][j-1] + blosum[(sequence1[i - 1], sequence2[j - 1])]
            left = dp_matrix[i-1][j] + gap_penalty
            up = dp_matrix[i][j-1] + gap_penalty
            dp_matrix[i][j] = max(diagonal, left, up)
    # trace back
    align1 = ''
    align2 = ''
    i = m
    j = n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp_matrix[i][j] == dp_matrix[i-1][j-1] + blosum[(sequence1[i - 1], sequence2[j - 1])]:
            align1 = sequence1[i-1] + align1
            align2 = sequence2[j-1] + align2
            i -= 1
            j -= 1
        elif i > 0 and dp_matrix[i][j] == dp_matrix[i-1][j] + gap_penalty:
            align1 = sequence1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = sequence2[j-1] + align2
            j -= 1   
    return dp_matrix[m][n],align1, align2

def multiple_seq_align(sequences):
    msl=[]
    msl.append(sequences[0][0])
    msl.append(sequences[0][1])
    for i in range(1,len(sequences)):
        psa_seq_c=sequences[i][0]
        psa_seq=sequences[i][1]
        msl=star_align(msl,psa_seq_c,psa_seq)
    return msl
            

def star_align(msl,psa_seq_c,psa_seq):
    center_seq=msl[0]
    num=len(center_seq)
    i=0
    while num!=0:
        if psa_seq_c[i]!='-' and center_seq[i]=='-' :
            psa_seq=insert_gap(psa_seq,i)
            psa_seq_c=insert_gap(psa_seq_c,i)
            i+=1
            num-=1
        elif psa_seq_c[i]=='-' and center_seq[i]!='-':
            num+=1
            center_seq=insert_gap(center_seq,i)
            for j in range(1,len(msl)):
                msl[j]=insert_gap(msl[j],i)
        else:
            num-=1
            i+=1
    msl[0]=center_seq
    msl.append(psa_seq)
    return msl

def insert_gap(sequence, index):
    return sequence[:index] + '-' + sequence[index:]

def make_star_sign(msl):
    center=msl[0]
    star=''
    for i in range(len(center)):
        check=center[i]
        is_star=True
        for seq in msl:
            if seq[i]!=check or seq[i]=='-':
                star=star+' '
                is_star=False
                break
        if is_star:
            star=star+'*'
    msl.append(star)
    return msl

if __name__ == '__main__':
    input_file = sys.argv[1]
    #input_file = "test8.txt"
    output_file = 'Assignment6_output.txt'

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
        input_BLOSUM="BLOSUM62.txt"
        blosum=blosum_matrix_to_dict(input_BLOSUM)     
    except FileNotFoundError:
        print("There is no such file : ",end='')
        print(input_BLOSUM)
        time.sleep(5)
        sys.exit()

    # start timer
    start_time = time.time()
    gsl=[]
    gsl_score=[]
    sequences_len=len(sequences)
    for i in range(sequences_len):
        tmp_score=0
        for j in range(sequences_len):
            temp=[]
            if i!=j:
                score,seq1,seq2= global_align(sequences[i], sequences[j], blosum, -5)
                tmp_score+=score
                temp.append(seq1)
                temp.append(seq2)
                gsl.append(temp)
        gsl_score.append(tmp_score)

    max_idx=gsl_score.index(max(gsl_score))
    if max_idx==0:
        center_gsl=gsl[0:sequences_len-1]
    else:
        center_gsl=gsl[max_idx*(sequences_len-1):(max_idx+1)*(sequences_len-1)]

    msl=multiple_seq_align(center_gsl)
    msl=make_star_sign(msl)
    # check elapsed time
    print("Elapsed time (microseconds): %f" % (time.time() - start_time))

    # make output file
    if len(msl[0])>60:
        newline=len(msl[0])//60
        i=0
        with open(output_file, 'w') as f:
            while newline!=0:
                for seq in msl:
                    f.write(seq[i:i+60]+"\n")
                f.write("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
                i+=60
                newline-=1
                if newline==0:
                    for seq in msl:
                        f.write(seq[i:len(seq)]+"\n")       
    else:
        with open(output_file, 'w') as f:
            for seq in msl:
                f.write(seq+"\n")
        f.close()



