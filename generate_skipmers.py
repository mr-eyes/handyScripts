import sys

def getSkipmers(m, n, seq):
    seq_len = len(seq)
    orfs = []
    for start in range(n):
        new_seq = ""
        for i in range(start, seq_len, n):
            new_seq += seq[i:i+m]
        
        orfs.append(new_seq)
    
    return orfs


seq = sys.stdin.read()

reading_frames = getSkipmers(2,3,seq)

kSize = 21

for i in range(3):
    for j in range(0, len(reading_frames[i]) - kSize + 1, 2):
        print(reading_frames[i][j:j+kSize])