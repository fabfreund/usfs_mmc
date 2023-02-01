#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import os


n=20 #Number of individual in the population
lab="Aptenodytes" #label of the population
sfs=[0]*(n-1)
dif_anc=0
anc_diff_N=0
ancestral_diff=0
anc_N=0
nb_SNP=0

##SNP file has to have format REF,#REF, ALT, #ALT, Allele_Outgroup1, Allele_Outgroup2
with open("SNPfile_"+lab+".txt","r") as file1:
    for ligne in file1:
        liste=re.split('\t|\n',ligne)
        if liste[4]!=liste[5]: #If two outgroup alleles different
            if liste[4]=="N" or liste[4]=="-" or liste[5]=="N" or liste[5]=="-": #check if outgroup allele known
                anc_diff_N+=1
            ancestral_diff+=1
        else : #if outgroup allele similar: check if similar to REF or ALT
            if liste[0]==liste[4]:
                sfs[int(liste[3])-1]+=1
                nb_SNP+=1
            elif liste[2]==liste[4]:
                sfs[int(liste[1])-1]+=1
                nb_SNP+=1
            else:
                if liste[4]=="N" or liste[4]=="-":
                    anc_N+=1
                dif_anc+=1


print(sfs)
print("Number of Polarized SNP",nb_SNP)
print("Number of Diallelic outgroup",ancestral_diff-anc_diff_N)
print("Number of site with different allele in the outgroup (not diallelic)",dif_anc-anc_N)

##Write SFS in file
with open(lab+"_SFS.txt","w") as file2:
    for i in range(len(sfs)):
        file2.write(str(i+1)+'\t'+str(sfs[i])+'\n')
