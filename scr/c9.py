#Le but de ce code est dénumérer tous les graphes planaires à 9 sommets,
#en partant d'un C9, et en ajoutant un nombre maximal d'arêtes.



#Les sommets sont numérotés de 0 à 8.
#Les arêtes (i,i+1) pour i de 0 à 8 sont déjà présentes.
#Un graphe est représenté par l'ensemble des arêtes (i,j,bool) avec i < j ajoutées
# et bool = True si l'arête est à l'intérieur du C9, False sinon.
#Un graphe est valide s'il n'y a pas de pair d'arêtes croisées :
#(i1,j1,b1) et (i2,j2,b2) tels que b1 == b2 et i1 < i2 < j1 < j2.
#Un graphe est maximal s'il est valide et qu'on ne peut pas ajouter d'arête
def verification(graphe):
    #On verifie que toutes les pairs d'arêtes du C9 sont bien distantes de 1
    for a1 in range(8):
        for a2 in range(a1+3,9):
            i1 = a1
            j1 = (a1 +1) %9
            i2 = a2
            j2 = (a2 +1 )%9
            aretes_possibles = [(i1,i2,True),(i1,j2,True),(j1,i2,True),(j1,j2,True),(i1,i2,False),(i1,j2,False),(j1,i2,False),(j1,j2,False)]
            trouve = False
            for arete in aretes_possibles:
                if arete in graphe:
                    trouve = True
                    break
            if not trouve:
                return False
    return True

def verification_simple(graphe,perm):
    #On verifie que toutes les pairs d'arêtes du C9 sont bien distantes de 1
    for a1 in range(8):
        for a2 in range(a1+3,9):
            i1 = perm[a1]
            j1 = perm[(a1 +1) %9]
            i2 = perm[a2]
            j2 = perm[(a2 +1)%9]
            aretes_possibles = [(i1,i2),(i1,j2),(j1,i2),(j1,j2)]
            trouve = False
            for arete in aretes_possibles:
                if arete[1] in graphe[arete[0]]:
                    trouve = True
                    break
            if not trouve:
                return False
    return True



def verification_simple8(graphe,perm):
    #On trouve le 9eme sommets qui n'est pas dans la permutation
    reste = [x for x in range(9) if x not in perm]
    s = reste[0]
    aretes_s_possibles = [(s,el) for el in graphe[s]]
    voisins_s = graphe[s]
    aretes_s_fonctionnelles = []
    #On verifie que toutes les pairs d'arêtes du C9 sont bien distantes de 1
    for a1 in range(7):
        for a2 in range(a1+3,7):
            i1 = perm[a1]
            j1 = perm[(a1 +1) %8]
            i2 = perm[a2]
            j2 = perm[(a2 +1)%8]
            aretes_possibles = [(i1,i2),(i1,j2),(j1,i2),(j1,j2)]
            trouve = False
            for arete in aretes_possibles:
                if arete[1] in graphe[arete[0]]:
                    trouve = True
                    break
            if not trouve:
                return False,[]
    #On verifie que l'arête (s,el) est bien fonctionnelle
    for arete in aretes_s_possibles:
        fonctionne = True
        for a1 in range(8):
            i1 = perm[a1]
            j1 = perm[(a1 +1) %8]
            if i1 in voisins_s or j1 in voisins_s:
                continue
            if arete[1] == perm[(a1+7)%8] or arete[1] == perm[(a1+2)%8]:
                continue
            fonctionne = False
            break
        if fonctionne :
            aretes_s_fonctionnelles.append(arete)

    return True,aretes_s_fonctionnelles

def verification_simple8_10(graphe,perm):
    #On trouve le 9eme sommets qui n'est pas dans la permutation
    reste = [x for x in range(10) if x not in perm]
    s1 = reste[0]
    s2 = reste[1]
    aretes_s1_possibles = [(s1,el) for el in graphe[s1]]
    aretes_s2_possibles = [(s2,el) for el in graphe[s2] if el != s1]
    aretes_s_possibles =[]
    if s2 in graphe[s1]:
        aretes_s_possibles.append((s1,s2))
    voisins_s1 = [x for x in graphe[s1] if x != s2]
    voisins_s2 = [x for x in graphe[s2] if x != s1]
    aretes_s_fonctionnelles = []
    #On verifie que toutes les pairs d'arêtes du C9 sont bien distantes de 1
    for a1 in range(7):
        for a2 in range(a1+3,7):
            i1 = perm[a1]
            j1 = perm[(a1 +1) %8]
            i2 = perm[a2]
            j2 = perm[(a2 +1)%8]
            aretes_possibles = [(i1,i2),(i1,j2),(j1,i2),(j1,j2)]
            trouve = False
            for arete in aretes_possibles:
                if arete[1] in graphe[arete[0]]:
                    trouve = True
                    break
            if not trouve:
                return False,[]
    #On verifie que l'arête (s,el) est bien fonctionnelle
    for arete in aretes_s1_possibles:
        fonctionne = True
        for a1 in range(8):
            i1 = perm[a1]
            j1 = perm[(a1 +1) %8]
            if i1 in voisins_s1 or j1 in voisins_s1:
                continue
            if arete[1] == perm[(a1+7)%8] or arete[1] == perm[(a1+2)%8]:
                continue
            fonctionne = False
            break
        if fonctionne :
            aretes_s_fonctionnelles.append(arete)
    for arete in aretes_s2_possibles:
        fonctionne = True
        for a1 in range(8):
            i1 = perm[a1]
            j1 = perm[(a1 +1) %8]
            if i1 in voisins_s2 or j1 in voisins_s2:
                continue
            if arete[1] == perm[(a1+7)%8] or arete[1] == perm[(a1+2)%8]:
                continue
            fonctionne = False
            break
        if fonctionne :
            aretes_s_fonctionnelles.append(arete)
    for arete in aretes_s_possibles:
        fonctionne = True
        for a1 in range(8):
            i1 = perm[a1]
            j1 = perm[(a1 +1) %8]
            if i1 in voisins_s1 or j1 in voisins_s2 or i1 in voisins_s2 or j1 in voisins_s1:
                continue
            fonctionne = False
            break
        if fonctionne :
            aretes_s_fonctionnelles.append(arete)

    return True,aretes_s_fonctionnelles




def ajoutable(graphe,arete):
    #renvoie True si l'arête est ajoutable, False sinon
    i1,j1,b1= arete[0],arete[1],arete[2]
    for arete2 in graphe:
        i2,j2,b2 = arete2[0],arete2[1],arete2[2]
        if b1 == b2 and ( i1 < i2 < j1 < j2 or i2 < i1 < j2 < j1):
            return False
    return True


lst_graphes_max = []
print("Début de l'énumération")
nb_trouve = 0

#Lire le fichier des graphes contenu dans le fichier ~/Documents/data/planar_graphes_9_good.txt
#le fichier ressemble à :
# Graph 23436, order 9.
# 0 : 8;
# 1 : 4 8;
# 2 : 5 6 7;
# 3 : 6 7 8;
# 4 : 1 5 6 8;
# 5 : 2 4 7 8;
# 6 : 2 3 4 8;
# 7 : 2 3 5 8;
# 8 : 0 1 3 4 5 6 7;
#
# Graph 23437, order 9.
# 0 : 8;
# 1 : 4 8;
# 2 : 5 6 7;
# 3 : 6 7 8;
# 4 : 1 5 6 8;
# 5 : 2 4 7 8;
# 6 : 2 3 4 7;
# 7 : 2 3 5 6;
# 8 : 0 1 3 4 5;

lst_graphes_max = []
graphe = [[] for _ in range(10)]
with open("/home/sdarmon/Documents/data/planar_conn.10.txt","r") as f:
    for line in f:
        if line[0] == "G":
            graphe = [[] for _ in range(10)]
        elif line[0] == " ":
            s = line.split(":")[1].split(";")[0]
            s = [int(x) for x in s if x != " "]
            for x in s:
                graphe[int(line[2])].append(x)
        else:
            lst_graphes_max.append(graphe)

print("Fin de la lecture")


def test_ajout(graphe,perm,restant):
    s = perm[-1]
    o = perm[0]
    #Condition d'arret
    if len(perm) == 8:
        if o in graphe[s]:
            return perm
    if len(perm) == 9:
        if o in graphe[s]:
            return perm
        else:
            return None
    #Induction
    voisin = [x for x in graphe[s] if x in restant]
    if len(voisin) == 0:
        return None
    for v in voisin:
        p = test_ajout(graphe,perm + [v], [x for x in restant if x != v])
        if p != None:
            return p
    return None

print("Nombre de graphes à tester",len(lst_graphes_max))
for graphe in lst_graphes_max:
    #Trouver la permutation des sommets qui correspond à la permutation du C9
    # s'il y a un sommet de degree 1, on skip
    skip = False
    deg2 = -1
    deg3 = -1
    deg4 = -1
    for i in range(10):
        d = len(graphe[i])
        if d == 1:
            skip = True
            break
        if d == 2:
            deg2 = i
        if d == 3:
            deg3 = i
        if d == 4:
            deg4 = i

    if skip :
        continue
    if deg2 == -1 and deg3 == -1 and deg4 == -1:
        print("Erreur, pas de sommets de degré < 5", graphe)
        continue

    if deg2 != -1:
        #On commence par ce sommet et ses voisins
        a = graphe[deg2][0]
        b = graphe[deg2][1]
        lst_perm = [[a,deg2,b]]
        for perm in lst_perm :
            restant = [x for x in range(10) if x not in perm]
            #On ajoute les sommets restants en testant toutes les permutations
            p = test_ajout(graphe,perm,restant)
            if p == None:
                continue
            # if len(p) == 9 and verification_simple(graphe,p):
            #     print(graphe,p)
            if len(p) == 8 :
                b,aretes = verification_simple8_10(graphe,p)
                if b and len(aretes) > 2:
                    print(graphe,p,aretes)

    elif deg3 != -1:
        #On commence par ce sommet et ses voisins
        a = graphe[deg3][0]
        b = graphe[deg3][1]
        c = graphe[deg3][2]
        lst_perm = [[a,deg3,b],[a,deg3,c],[c,deg3,b]]
        for perm in lst_perm :
            restant = [x for x in range(10) if x not in perm]
            #On ajoute les sommets restants en testant toutes les permutations
            p = test_ajout(graphe,perm,restant)
            if p == None:
                continue
            #On permute les arêtes pour obtenir le graphes
            # if len(p) == 9 and verification_simple(graphe,p):
            #     print(graphe,p)
            if len(p) == 8 :
                b,aretes = verification_simple8_10(graphe,p)
                if b and len(aretes) > 2:
                    print(graphe,p,aretes)
    elif deg4 != -1:
        #On commence par ce sommet et ses voisins
        a = graphe[deg4][0]
        b = graphe[deg4][1]
        c = graphe[deg4][2]
        d = graphe[deg4][3]
        lst_perm = [[a,deg4,b],[a,deg4,c],[c,deg4,b],[a,deg4,d],[d,deg4,b],[c,deg4,d]]
        for perm in lst_perm :
            restant = [x for x in range(10) if x not in perm]
            #On ajoute les sommets restants en testant toutes les permutations
            p = test_ajout(graphe,perm,restant)
            if p == None:
                continue
            #On permute les arêtes pour obtenir le graphes
            # if len(p) == 9 and verification_simple(graphe,p):
            #     print(graphe,p)
            #     break
            if len(p) == 8 :
                b,aretes = verification_simple8_10(graphe,p)
                if b and len(aretes) > 2:
                    print(graphe,p,aretes)
