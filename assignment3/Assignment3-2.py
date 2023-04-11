#2021245044 소프트웨어학부 서은하
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
    if check_seq_num==1:
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

# finding low-complexity regions without regex
def find_low_complexity_regions(seq):
    low_complexity_regions = []
    for starting_idx in range(len(seq)-6):
        maximum_length=[]
        for i in range(2,6):
            repeat_cnt=0
            while(True):
                if starting_idx==0:
                    repeat_cnt=finding_repeat_seq(seq,starting_idx,i,repeat_cnt)
                elif low_complexity_regions and starting_idx<=low_complexity_regions[-1][-1]:
                    break
                else:                    
                    repeat_cnt=finding_repeat_seq(seq,starting_idx,i,repeat_cnt)
                if repeat_cnt>1:
                    temp=[]
                    end_idx=starting_idx+(repeat_cnt+1)*i-1
                    temp.append(starting_idx)
                    temp.append(end_idx)
                    maximum_length.append(temp)
                    break
                else:
                    break
        if maximum_length:
            low_complexity_regions.append(maximum_length[-1])  
    return low_complexity_regions

# finding repeat count
def finding_repeat_seq(seq,s_idx,s_len,repeat_cnt):
    if seq[s_idx:s_idx+s_len]==seq[s_idx+s_len:s_idx+2*s_len]:
        repeat_cnt+=1
        return finding_repeat_seq(seq,s_idx+s_len,s_len,repeat_cnt)
    else:
        return repeat_cnt

matches=find_low_complexity_regions(seq)

# writing results to output file
if matches:
    with open("Assignment3-2_output.txt", "w") as file:
        for match in matches:
            file.write(str(match[0]) + "\n")
else:
    print("No low-complexity region found")
print("Elapsed time (microseconds): %f" % (time.time() - start_time))