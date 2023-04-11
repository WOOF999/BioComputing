#2021245044 소프트웨어학부 서은하
import re
import sys
import time

# read input file
try:
    with open(sys.argv[1], "r") as file:
    #with open("dna.txt", "r") as file:
        lines = file.readlines()
except:
    print("No input file")
    sys.exit()

# read first dna sequence
check_seq_num=0
seq = ""
for line in lines:
    if line[0] == ">":
        check_seq_num+=1
        continue
    elif line[0] !=">"and check_seq_num==0:
        print("No correct format")
        sys.exit()
    if check_seq_num==1 :
        strip_tmp= line.strip().upper()
        whitespace_tmp=strip_tmp.replace(' ',"")
        seq += whitespace_tmp

# check empty file
if len(seq)==0:
    print("No DNA sequence")
    sys.exit()

# check DNA sequence
if not set(seq).issubset({'A', 'C', 'G', 'T'}):
    print("No DNA sequence")
    sys.exit()

# check elapsed time 
start_time = time.time()

# Find all matches in the sequence
pattern=re.compile(r"((\w{2,5})\2{2,})")
matches = pattern.finditer(seq)
low_complexity_regions_index=[]

# Find starting index
for match in matches:
    low_complexity_regions_index.append(match.start())

# make output file
if len(low_complexity_regions_index) == 0:
    print("No low-complexity region found") 
else:
    with open("Assignment3-1_output.txt", 'w') as f:
        for region in low_complexity_regions_index:
            f.write(str(region))
            f.write("\n")
    print("Elapsed time (microseconds): %f" % (time.time() - start_time))