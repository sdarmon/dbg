import sys

Arg = sys.argv[:]

def homomorphic_compression(read,w):
    #Removing consecutive duplicates window by window
    S = read[0:w]
    current_w = read[0:w]
    for i in range(1, len(read)-w):
        next_w = read[i:i+w]
        if current_w != next_w:
            S += read[i+w-1]
            current_w = next_w
    return S

if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " TE.fa name_TE.txt output")
    exit()

#Reading name_TE.txt into a list
TE = []
with open(Arg[2], 'r') as f:
    for line in f:
        TE.append(line[:-1])

#Reading TE.fa into a list

on_a_TE = False

with open(Arg[1], 'r') as f:
    with open(Arg[3], 'w') as f_out:
        for line in f:
            if line[0] == '>' and line[0:-1] in TE:
                on_a_TE = True
                f_out.write(line)
                continue
            if line[0] == '>' :
                on_a_TE = False
                continue
            if on_a_TE:
                f_out.write(line)