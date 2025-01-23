# Ce programme permet de savoir quelles bulles (d'un fichier bulle.fa) traversent les composantes (du dossier dirComp)
# On indique si les chemins traversent les composantes par le chemin du haut, du bas ou des deux.

import sys

Arg = sys.argv[:]


# Fonction renvoyant le reverse complement d'une séquence de nucléotides
def reverseC(seq):
    s = ""
    for el in seq:
        if el == 'A':
            s = 'T' + s
        elif el == 'T':
            s = 'A' + s
        elif el == 'C':
            s = 'G' + s
        elif el == 'G':
            s = 'C' + s
    return s


# Renvoie une séquence de nucléotides en majuscule / En fait on s'en fiche, les minuscules correspondent au bout de
# l'unitig.
def maj(sequ):
    s = ""
    for el in sequ:
        if el == 'a':
            s = s + 'A'
        elif el == 'c':
            s = s + 'C'
        elif el == 'g':
            s = s + 'G'
        elif el == 't':
            s = s + 'T'
        else:
            s = s + el
    return s


if len(Arg) not in [5, 6, 7, 8, 9]:
    print("Use : " + Arg[
        0] + " dirComp bulle.fa nbComp threshold -k kmer [-rapport/-label event.tsv/-compare type.fa]\n\t k = 41 by default.")
    exit()
if len(Arg) in [5, 6, 7, 8, 9]:
    titre = ""  # Variable contenant le titre de la séquence dans le fichier .fa
    seq = ""  # Sequence de nucléotides
    event = dict()
    sequences = dict()
    seq_vu=set()
    bubble_vu = set()
    if Arg[-2] == "-label":
        with open(Arg[-1], 'r') as f:
            for line in f:
                L = line.split("\t")
                type = L[4]
                bubble = L[15]
                event[bubble] = type

    if Arg[-2] == "-compare":
        with open(Arg[-1], 'r') as f:
            odd = 1
            for line in f:
                if odd:
                    L = line.split("|")
                    bubble = L[0][1:] + "|" + L[1]
                    odd = 0
                else:
                    key = maj(line[:-1])
                    if sequences.get(key,-1) != -1:
                        sequences[key].append(bubble)
                    else:
                        sequences[key]= [bubble]
                    odd = 1

    if (len(Arg) >= 7 and Arg[5][1] == 'k'):
        k = int(Arg[6])
    else:
        k = 41
    threshold = int(Arg[4])
    kmerFrom = dict()  # Dictionnaire associant à un kmer donné la composante dont il est issu
    # Récupération des séquences des composantes
    for i in range(int(Arg[3])):
        with open(Arg[1] + "comp" + str(i) + ".txt", 'r') as f:
            for line in f:
                if len(line) < 2:
                    continue
                L = line.split("\t")
                # On l'ajoute au dictionnaire ainsi que son complémentaire
                compl = reverseC(L[1])
                for pos in range(len(L[1]) - k + 1):
                    kmerFrom[L[1][pos:k + pos]] = i
                    kmerFrom[compl[pos:k + pos]] = i
    compVu = [0 for i in range(int(Arg[3]))]
    # On parcourt maintenant les bulles
    with open(Arg[2], 'r') as f:
        upper = False  # Variable booléenne indiquant si on est un upper path ou un under path
        for line in f:
            if len(line) < 2:
                continue
            if titre == "":  # Cas où la ligne est un titre
                titre = line[:-1]
                upper = not upper  # On change la valeur de upper
                if upper:
                    L = line.split("|")
                    bubbles = L[0][1:] + "|" + L[1]
                    type = event.get(bubbles, "NA")
                    titreUpper = titre[:-1]  # On stocke le titre
                else:
                    titreUnder = titre[:-1]
                continue
            seq = line[:-1]
            seq = maj(seq)  # Les séquences ayant des nucléotides écrits en miniscule, on les passe en majuscule
            # On récupère les numéros des composantes contenant des kmers de la séquence
            comp_possible = set()
            trouve = 0
            for pos in range(len(seq) - k + 1):
                mer = seq[pos: pos + k]
                ind = kmerFrom.get(mer, -1)
                if ind != -1:
                    comp_possible.add(ind)
                    trouve = 1
            if upper:  # Cas chemin du haut
                seqUpper = line[:-1]
                upperComp = comp_possible
                debut = kmerFrom.get(seq[0:k], -1)
                end = kmerFrom.get(seq[-k:], -1)
                trouveUpper = trouve
                intersect_connu_upper = sequences.get(seq,"UpperNotFound")
                if intersect_connu_upper != "UpperNotFound":
                    seq_vu.add(seq)
                    for el in intersect_connu_upper:
                        bubble_vu.add(el)
                    intersect_connu_upper = "||".join(intersect_connu_upper)
            else:  # Cas chemin du bas
                seqUnder = line[:-1]
                underComp = comp_possible
                trouveUnder = trouve
                intersect_connu_under = sequences.get(seq,"UnderNotFound")
                if intersect_connu_under != "UnderNotFound":
                    seq_vu.add(seq)
                    for el in intersect_connu_under:
                        bubble_vu.add(el)
                    intersect_connu_under = "||".join(intersect_connu_under)
                # On peut écrit la bulle et son rapport si c'est intéressant
                if True: #if trouveUnder or trouveUpper:
                    text = ""
                    printing = True #False initialement
                    if len(upperComp) != 0 or len(underComp) != 0:
                        A = upperComp & underComp
                        B = upperComp - underComp
                        C = underComp - upperComp
                        if debut != -1:
                            printing = True
                            text += "Start in " + str(debut) + "\t"
                        if end != -1:
                            printing = True
                            text += "End in " + str(end) + "\t"
                        # if debut != end : #debut>= 0 and fin >= 0 and
                        #     printing = False
                        #     titre = ""
                        #     continue
                        if len(A) != 0:
                            t = ""
                            # printing = False
                            # titre = ""
                            # continue
                            for el in A:
                                t += " " + str(el)
                            text += "in both path :" + t + "\t"
                        if len(B) != 0:
                            t = ""
                            # vu = False
                            # for el in B:
                            #     if compVu[el] == 0:
                            #         compVu[el] = 1
                            #         vu = True
                            #         break
                            # if vu == False:
                            #     printing = False
                            #     titre = ""
                            #     continue
                            printing = True
                            for el in B:
                                t += " " + str(el)
                            text += "only in upper :" + t + "\t"
                        if len(C) != 0:
                            t = ""
                            # vu = False
                            # for el in B:
                            #     if compVu[el] == 0:
                            #         compVu[el] = 1
                            #         vu = True
                            #         break
                            # if vu == False:
                            #     printing = False
                            #     titre = ""
                            #     continue
                            printing = True
                            for el in C:
                                t += " " + str(el)
                            text += "only in under :" + t + "\t"
                    if printing :
                        if Arg[-1] == "-rapport":
                            print(titreUpper)
                            print(seqUpper)
                            print(titreUnder)
                            print(seqUnder)

                        text += bubbles + "\t" + type + "\t" + intersect_connu_upper + "\t" + intersect_connu_under
                        print(text)
            titre = ""  # On part pour la ligne suivante qui sera un titre
    """
    for key,value in sequences.items():
        if key not in seq_vu:
            for bubble in value:
                if bubble not in bubble_vu:
                    print("missing",bubble)
    
    """