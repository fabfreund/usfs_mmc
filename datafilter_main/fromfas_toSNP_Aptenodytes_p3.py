#!/usr/bin/python
# -*- coding: utf-8 -*-

prim=">Aptenodytes_patagonicus" #population name in fastafile
out=">Aptenodytes_forsteri" #outgroup name in fasta file

n=20 #Number of individual in the population
lab="Aptenodytes" #label of the population

with open("../example_fasta/Aptenodytes_concat_alignment_clean.fas","r") as file1:
    nb_seq=0
    seq=''
    seqs=[]
    anc=0
    seq_ancestrale=[]
    for ligne in file1:
            
            if ligne[0]=='>' and nb_seq==0:
                nb_seq=1
            elif ligne[0]=='>' and nb_seq>0:
                if anc==1:
                    seq_ancestrale.append(seq)
                    anc=0
                else:
                    seqs.append(seq)
                seq=''
                nb_seq+=1
            else:
                seq+=ligne[:-1]
            
            if ligne[:len(out)]==out:
                anc=1
            
    seqs.append(seq)

print("Number of sequences",len(seqs))
print("Lengt",len(seqs[-1]))

liste_dico=[]

for i in range(len(seqs[0])):
    dico={}
    liste_dico.append(dico)

for j in range(0,len(seqs)):
    for i in range(len(seqs[j])):
        if seqs[j][i]=='A' or seqs[j][i]=='T' or seqs[j][i]=='C' or seqs[j][i]=='G':
            if seqs[j][i] in liste_dico[i] :
                liste_dico[i][seqs[j][i]]+=1
            else:
                liste_dico[i][seqs[j][i]]=1
        else:
            if 'N' in liste_dico[i]:
                liste_dico[i]['N']+=1
            else:
                liste_dico[i]['N']=1
SNP=[] ##diallelic SNP
nb_un=0 #conversed positions
nbSNP=0 #Number of diallelic SNP
nb_trois=0 #Number of triallelic positions
nb_quatre=0 #Number of positions with 4 alleles
nb_N=0 #Number of unkown positions

for dico in liste_dico:
    if 'N' in dico:
        nb_N+=1
    elif len(dico)==1:
        nb_un+=1
    elif len(dico)==2:
        nbSNP+=1
    elif len(dico)==3:
        nb_trois+=1
    elif len(dico)==4:
        nb_quatre+=1
    else:
        print(dico)
        
for i in range(len(liste_dico)):
    if 'N' not in liste_dico[i] and len(liste_dico[i])==2:
        SNP.append([list(liste_dico[i].keys())[0],str(liste_dico[i][list(liste_dico[i].keys())[0]]),list(liste_dico[i].keys())[1],str(liste_dico[i][list(liste_dico[i].keys())[1]]),seq_ancestrale[0][i],seq_ancestrale[1][i]])
        

			
print("Conserved",nb_un)
print("2 variants",nbSNP)
print("3 variants",nb_trois)
print("4 variants",nb_quatre)
print("unknown",nb_N)

print("Number of SNP in file",len(SNP))



with open("SNPfile_"+lab+".txt","w") as file2:
	for snp in SNP:
		file2.write(str(snp[0])+'\t'+str(snp[1])+'\t'+str(snp[2])+'\t'+str(snp[3])+'\t'+str(snp[4])+'\t'+str(snp[5])+'\n')
	
