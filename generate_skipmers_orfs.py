from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os

if len(sys.argv) < 2:
    exit("python generate.py <fasta> <names>")

else:
    input_fasta = sys.argv[1]
    old_names = sys.argv[2]

output_fasta = "skipmers_" + os.path.basename(input_fasta)
names_file = output_fasta + ".names"
m = 2
n = 3

def getSkipmers(m, n, seq):
    seq_len = len(seq)
    orfs = []
    for start in range(n):
        new_seq = ""
        for i in range(start, seq_len, n):
            new_seq += seq[i:i+m]
        
        orfs.append(new_seq)
    
    return orfs


sequences = SeqIO.parse(open(input_fasta), 'fasta')
output_handle = open(output_fasta, 'w')
namesFile =  open(names_file, 'w')
oldNames = open(old_names, 'r')
        
for seq in sequences:
    original_header = seq.description + "|"
    skipmers_frames = getSkipmers(m, n, str(seq.seq))
    oldGeneName = next(oldNames).strip().split("\t")[1]
    
    print("Seq: %s Length: %s bp" % (seq.description,len(seq.seq)))
    
    for _i in range(len(skipmers_frames)):
        frame = skipmers_frames[_i]
        new_header = original_header + "RF-" + str(_i)
        _temp_seq = SeqRecord(Seq(frame))
        _temp_seq.id = _temp_seq.name = _temp_seq.description = new_header
        SeqIO.write(_temp_seq, output_handle, 'fasta')
        print("Seq: %s Length: %s bp" % (_temp_seq.name,len(_temp_seq.seq)))
        namesFile.write(new_header + "\t" + oldGeneName + "\n")


output_handle.close()
namesFile.close()
oldNames.close()