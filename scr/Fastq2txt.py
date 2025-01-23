#Convert a fastq file to a txt file (seq_name\tseq)
#We only keep the name of the sequence and the sequence itself, we do not keep the quality score.
import sys

Arg = sys.argv[:]

if len(Arg) not in [3]:
    print("Use : " + Arg[0] + " input.fastq output.txt")
    exit()

with open(Arg[1], 'r') as f:
    with open(Arg[2], 'w') as f_out:
        for line in f:
            if line[0] == '@':
                f_out.write(">"+line[1:])
                #now getting the complet sequence before adding it to the output file
                seq = ''
            elif line[0] == '+':
                f_out.write(seq + '\n')
            else:
                seq += line[:-1]

