#2021245044 소프트웨어학부 서은하
import re
import sys

check_comment=False
DNASeq=[]

file=open(sys.argv[1],'r')
#file=open('dna.txt','r')

def text_preprocess(line):
    strip_tmp=line.strip()
    whitespace_tmp=strip_tmp.replace(' ',"")
    upper_tmp=whitespace_tmp.upper()
    return upper_tmp 

def concat(DNASeq_lst):
    DNASeq="".join(DNASeq_lst)
    return DNASeq

def check_DNA_sequence(DNASeq):
    DNASeq=re.sub("[^ATCG]","*",DNASeq) 
    if '*' not in DNASeq: #DNA sequence exist
        DNASeq=DNASeq.replace("T","U")
        return DNASeq
    else:
        print("No DNA sequence")
        return False

while True:
    line=file.readline()
    #check a single line of comment
    if check_comment is False:
        check_comment=True
    else:
        if not line:
            break
        DNASeq.append(text_preprocess(line))       
file.close()  

DNAStr=concat(DNASeq)

#make output file 
file=open("Assignment1_output.txt",'w')    
if check_DNA_sequence(DNAStr):
    file.writelines(check_DNA_sequence(DNAStr))
file.close()