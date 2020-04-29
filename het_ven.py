#/usr/bin/env python

import sys
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles
# from matplotlib.pyplot as plt

coarse_all = ['Pvag_2_107.sort', 'Pvag_2_108.sort','Pvag_2_117.sort','Pvag_3_114.sort','Pvag_2_120.sort','Pvag_2_207.sort','Pvag_2_212.sort', 'Pvag3_179.sort','Pvag_2_76.sort','Pvag_3_123.sort','Pvag_3_179.sort', 'Pvag_3_183.sort','Pvag_3_202.sort']
coarse_dip = []
coarse_poly = []
turf = []
dist = ['Pvag_1_284500.sort','Pvag_1_364977.sort','Pvag_1_647916.sort','Pvag_1_ur5.sort','Pvag_3_139.sort','Pvag_3_64.sort']
admix = ['Pvag_3_208.sort','612771.sort']
marker_dict ={}
individual_order =[]
vendia = {"dis" : 0, "tur" : 0, "cor": 0, 'distur': 0, 'discor' : 0, 'turcor' : 0, "disturcor" : 0,'' : 0}
with open(sys.argv[2],'r') as f:
	for line in f:
		line=line.split()
		individual_order.append (line[1])
		if line[1] not in coarse_all and line[1] not in dist and line[1] not in admix:
			turf.append(line[1])




geno_dict = {'11':'A','22':'B','12':'H','21':'H'}
with open(sys.argv[1],'r') as f:
	for line in f:
		line=line.split()
		marker_dict[line[0]+'_'+line[3]] = []
		for i in range(4,len(line),2):
			if line[i] +line[i+1] in geno_dict:

				marker_dict[line[0]+'_'+line[3]].append( geno_dict[line[i] +line[i+1]])

			else:
					# print '!' + line[i] +line[i+1] 
				marker_dict[line[0]+'_'+line[3]].append('_')
for x in marker_dict:
	stuff = ['dis','tur', 'cor']
	things = [False,False,False]
	count = 0
	vennstuff = ""
	print marker_dict[x].count('H') / float(len(marker_dict[x]))
	if (marker_dict[x].count('H') / float(len(marker_dict[x]))) > float(sys.argv[3]):
		for i in marker_dict[x]:
			if i == 'H':
				if individual_order[count] in dist:
					things[0] = True
				elif individual_order[count] in turf:
					things[1] = True
				elif individual_order[count] in coarse_all:
					things[2] = True
			
			count += 1
	for z in [0,1,2]:
		if things[z]:
			vennstuff += stuff[z]
	vendia[vennstuff] += 1

print vendia

venn3(subsets = (vendia['tur'], vendia['cor'], vendia['turcor'], vendia["dis"], vendia['distur'], vendia['discor'], vendia['disturcor']), set_labels = ('Fine-textured', 'Coarse-textured', 'P. distichum'), alpha = 0.5);

plt.savefig("foo.pdf", bbox_inches='tight')
plt.show()


