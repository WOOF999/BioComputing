import sys
import time

def kmp_match(text, pattern):
    n = len(text)
    m = len(pattern)
    pi = compute_prefix(pattern)
    matches = []
    j = 0
    for i in range(n):
        while j > 0 and pattern[j] != text[i]:
            j = pi[j-1]
        if pattern[j] == text[i]:
            j += 1
        if j == m:
            matches.append(i-m+1)
            j = pi[j-1]
    return matches

def compute_prefix(pattern):
    m = len(pattern)
    pi = [0] * m
    j = 0
    for i in range(1, m):
        while j > 0 and pattern[j] != pattern[i]:
            j = pi[j-1]
        if pattern[j] == pattern[i]:
            j += 1
        pi[i] = j
    return pi

if __name__ == '__main__':
    pattern_file = sys.argv[1]
    text_file = sys.argv[2]
    output_file = sys.argv[3]

    with open(pattern_file, 'r') as f:
        pattern = f.read().strip()

    with open(text_file, 'r') as f:
        text = f.read().strip()

    start_time = time.perf_counter()
    matches = kmp_match(text, pattern)
    elapsed_time = (time.perf_counter() - start_time) * 1e6

    with open(output_file, 'w') as f:
        for match in matches:
            f.write(str(match) + '\n')

    print(f'Time elapsed: {elapsed_time:.2f} microseconds')
