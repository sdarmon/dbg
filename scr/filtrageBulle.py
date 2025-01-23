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


# Renvoie une séquence de nucléotides en masjucule / En fait on s'en fiche, les miniscules correspondent au bout de
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
        0] + " dirComp bulle.fa nbComp threshold -k kmer [-rapport/-label event.tsv]\n\t k = 41 by default.")
    exit()
if len(Arg) in [5, 6, 7, 8, 9]:
    titre = ""  # Variable contenant le titre de la séquence dans le fichier .fa
    seq = ""  # Sequence de nucléotides
    event = dict()
    if Arg[-2] == "-label":
        with open(Arg[-1], 'r') as f:
            for line in f:
                L = line.split("\t")
                type = L[4]
                bubble = L[15]
                event[bubble] = type

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
                    bubble = L[0][1:] + "|" + L[1]
                    type = event.get(bubble, "NA")
                    titreUpper = titre[:-1]  # On stocke la titre
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
            else:  # Cas chemin du bas
                seqUnder = line[:-1]
                underComp = comp_possible
                trouveUnder = trouve
                # On peut écrit la bulle et son rapport si c'est intéressant
                if trouveUnder or trouveUpper:
                    text = ""
                    printing = False
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
                            printing = False
                            titre = ""
                            continue
                            for el in A:
                                t += " " + str(el)
                            text += "in both path :" + t + "\t"
                        if len(B) != 0:
                            t = ""
                            vu = False
                            for el in B:
                                if compVu[el] == 0:
                                    compVu[el] = 1
                                    vu = True
                                    break
                            if vu == False:
                                printing = False
                                titre = ""
                                continue
                            printing = True
                            for el in B:
                                t += " " + str(el)
                            text += "only in upper :" + t + "\t"
                        if len(C) != 0:
                            t = ""
                            vu = False
                            for el in B:
                                if compVu[el] == 0:
                                    compVu[el] = 1
                                    vu = True
                                    break
                            if vu == False:
                                printing = False
                                titre = ""
                                continue
                            printing = True
                            for el in C:
                                t += " " + str(el)
                            text += "only in under :" + t + "\t"
                        if printing:
                            if Arg[-1] == "-rapport":
                                print(titreUpper)
                                print(seqUpper)
                                print(titreUnder)
                                print(seqUnder)
                            text += bubble + "\t" + type
                            print(text)
            titre = ""  # On part pour la ligne suivante qui sera un titre
