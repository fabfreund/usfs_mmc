#!/usr/bin/python
# -*- coding: utf-8 -*-

prim="Aptenodytes_patagonicus" #population name in fastafile
out="Aptenodytes_forsteri" #outgroup name in fasta file

n=20 #Number of individual in the population
lab="Aptenodytes" #label of the population

prim_contig={}
out_contig={}

contig=''
seq=''
espece=''

with open("Aptenodytes_patagonicus+Aptenodytes_forsteri.fas","r") as file1:
	for ligne in file1:
		if ligne[0]==">":
			if espece==prim:
				prim_contig[contig].append(seq)
			elif espece==out:
				out_contig[contig].append(seq)
			liste=ligne.split('|')
			contig=liste[0]
			if prim_contig.has_key(contig)==False:
				prim_contig[contig]=[]
				out_contig[contig]=[]
			espece=liste[1]
			seq=''
			
		else:
			seq+=ligne[:-1]
	
if espece==prim:
	prim_contig[contig].append(seq)
elif espece==out:
	out_contig[contig].append(seq)
	
print "Number of contigs before filtering :",len(prim_contig), len(out_contig)


#Only keep contigs if obtained for all individuals in the population (n indiv)
		
for k, v in prim_contig.items():
	if len(v) != n:
		del prim_contig[k]
		del out_contig[k]
	
print "Number of contigs after filtering",len(prim_contig), len(out_contig)

ltot=0
for key in prim_contig:
	ltot+=len(prim_contig[key][0])

print "Total length",ltot


SNP=[] ##diallelic SNP
nb_un=0 #conversed positions
nbSNP=0 #Number of diallelic SNP
nb_trois=0 #Number of triallelic positions
nb_quatre=0 #Number of positions with 4 alleles
nb_N=0 #Number of unkown positions

for contig in prim_contig:
	liste_dico=[]

	for i in xrange(len(prim_contig[contig][0])):
		dico={}
		liste_dico.append(dico)

	for j in xrange(0,len(prim_contig[contig])):
		for i in xrange(len(prim_contig[contig][j])):
			if prim_contig[contig][j][i]=='A' or prim_contig[contig][j][i]=='T' or prim_contig[contig][j][i]=='C' or prim_contig[contig][j][i]=='G':
				if liste_dico[i].has_key(prim_contig[contig][j][i]):
					liste_dico[i][prim_contig[contig][j][i]]+=1
				else:
					liste_dico[i][prim_contig[contig][j][i]]=1
			else:
				if liste_dico[i].has_key('N'):
					liste_dico[i]['N']+=1
				else:
					liste_dico[i]['N']=1

	for dico in liste_dico:
		if dico.has_key('N'):
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
			print "/!\ ",dico

	for i in xrange(len(liste_dico)):
		if liste_dico[i].has_key('N')==False and len(liste_dico[i])==2:
			SNP.append([liste_dico[i].keys()[0],str(liste_dico[i][liste_dico[i].keys()[0]]),liste_dico[i].keys()[1],str(liste_dico[i][liste_dico[i].keys()[1]]),out_contig[contig][0][i],out_contig[contig][1][i]])
	
		
print "Conserved",nb_un
print "2 variants",nbSNP
print "3 variants",nb_trois
print "4 variants",nb_quatre
print "unknown",nb_N

print "Number of SNP in file",len(SNP)



with open("SNPfile_"+lab+".txt","w") as file2:
	for snp in SNP:
		file2.write(str(snp[0])+'\t'+str(snp[1])+'\t'+str(snp[2])+'\t'+str(snp[3])+'\t'+str(snp[4])+'\t'+str(snp[5])+'\n')
	
