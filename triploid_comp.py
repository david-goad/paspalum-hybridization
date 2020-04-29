#!/usr/bin/env python
import time
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import operator
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms



count_dict = {}
ratio_dict = {}
hist_dict = {}
order = []
ratio_dict2 = {}
site_dict = {}

with open(sys.argv[1], 'r') as f:
	for line in f:
		line = line.strip()
		line = line.split()
		count_dict[line[1]] = [0,0]
		hist_dict[line[1]] = [0,0,0,0,0,0,0,0,0,0,0]
		ratio_dict[line[1]] = []
		order.append(line[1])
		ratio_dict2[line[1]] = []

print order
with open(sys.argv[2], 'r') as f:
	for line in f:
		line = line.strip()
		line = line.split()
		LG = line[0]
		pos = line[3]


		all_het = True

		site = LG + '_' + pos

		if line[0]  :#== 'LG10' and int(pos) <= 5000000:



			for i in range(len(order)):
				
				if line[4+(i*2)] == line[5+(i*2)]:

					all_het = False
					


			if all_het == True:
				site_dict[site] = []


				for i in range(len(order)):
					site_dict[site].append(0)

					# print order[i]

					if line[4+(i*2)] != line[5+(i*2)]:
						command = "samtools view " + order[i] + '.bam "' + LG +':' + pos + '-' + pos + '" '
						p = subprocess.Popen(command, shell =True,stdout=subprocess.PIPE)
						

						x = p.communicate()
						


						print LG + ' ' + pos

						print count_dict
						print hist_dict
						read_count = [0,0]
						read_count_dict = {}
					
					
						# print x[0]
						

						for il_read in x[0].split('\n'):
							if len (il_read) > 10:
						

									# print order[i]
								subs = il_read.split()[12].split(':')[2]
								if subs in read_count_dict:
									read_count_dict[subs] += 1
								else:
									read_count_dict[subs] = 1

						sorted_x = sorted(read_count_dict.items(), key=operator.itemgetter(1))

						if len(sorted_x) > 1:

							highest = [sorted_x[-1],sorted_x[-2]]

							has_ref = False
							has_alt = False
							count_total = 0
							for y in highest:
								print y

								count_total += y[1]
								if len(y[0]) == 2:
									read_count[0] = y[1]
									has_ref = True

								else:
									read_count[1] = y[1]
									has_alt = True

							print count_total




							if has_ref == True and count_total > 20 and has_alt == True:
								count_dict[order[i]][0] += read_count[0]
								count_dict[order[i]][1] += read_count[1]
								ratio_dict[order[i]].append(read_count)

								ratio = float(read_count[0]) / float(read_count[0]+read_count[1])
								hist_dict[order[i]][int(ratio * 10)] += 1
								ratio_dict2[order[i]] .append( ratio)


								site_dict[site][i] = ratio

								if ratio > .75:
									print order[i] +' ' + LG + ' ' + pos






					

print count_dict
print len (site_dict)

for i in ratio_dict:


	data = np.array(ratio_dict[i])
	x, y = data.T
	plt.scatter(x,y,color = 'k', alpha=0.1)

	plt.ylabel('X genome')
	plt.xlabel('D genome')
	plt.title(i)
	plt.axis('square')

	y_lim = plt.ylim()
	x_lim = plt.xlim()
	plt.plot([0,x_lim[1]], [0,y_lim[1]], color='r', linestyle='--', linewidth=1)
	plt.plot([0,x_lim[1]], [0,y_lim[1]*2], color='r', linestyle='--', linewidth=1)
	plt.plot([0,x_lim[1]], [0,y_lim[1]/2], color='r', linestyle='--', linewidth=1)

	

	file_name = i + '_scatter.pdf'

	plt.savefig(file_name) 

	plt.show()


for i in ratio_dict2:



	
	plt.hist(ratio_dict2[i], color = 'black', bins = 50, range=(0,1))

	plt.title(i)
	plt.xlabel('Percent of reads identical to reference genome')
	plt.ylabel('number of loci')


	file_name = i + '_hist.pdf'

	plt.savefig(file_name) 

	plt.show()




# x = []
# y = []
# for i in site_dict:
	
# 	if site_dict[i][1] not in [0,1] and site_dict[i][2] not in [0,1]:
# 		x.append(site_dict[i][1])
# 		y.append(site_dict[i][2])

# plt.scatter(x,y, color = 'k', alpha=0.1)
# plt.gca().set_aspect('equal', adjustable='box')
# plt.ylabel(order[2])
# plt.xlabel(order[1])
# plt.show()



# x = []
# y = []
# for i in site_dict:
# 	if site_dict[i][1] not in [0,1] and site_dict[i][0] not in [0,1]:
# 		x.append(site_dict[i][1])
# 		y.append(site_dict[i][0])

# plt.scatter(x,y, color = 'k', alpha=0.1)
# plt.gca().set_aspect('equal', adjustable='box')
# plt.ylabel(order[0])
# plt.xlabel(order[1])
# plt.show()



# x = []
# y = []
# for i in site_dict:
# 	if site_dict[i][0] not in [0,1] and site_dict[i][2] not in [0,1]:
# 		x.append(site_dict[i][0])
# 		y.append(site_dict[i][2])

# plt.scatter(x,y, color = 'k', alpha=0.1)
# plt.gca().set_aspect('equal', adjustable='box')
# plt.ylabel(order[2])
# plt.xlabel(order[0])
# plt.show()


fig, axs = plt.subplots(3, 2)

#117
i = 'Pvag_2_117.sort'
data = np.array(ratio_dict[i])
x, y = data.T

axs[0,0].scatter(x,y,color = 'k', alpha=0.1)

axs[0,0].ylabel('X genome')
axs[0,0].xlabel('D genome')
axs[0,0].title(i)
axs[0,0].gca().set_aspect('equal','box')
x_lim = axs[0,0].xlim()
axs[0,0].plot([0,x_lim[1]], [0,y_lim[1]], color='r', linestyle='--', linewidth=1)
axs[0,0].plot([0,x_lim[1]], [0,y_lim[1]*2], color='r', linestyle='--', linewidth=1)
axs[0,0].plot([0,x_lim[1]], [0,y_lim[1]/2], color='r', linestyle='--', linewidth=1)





axs[0,1].hist(ratio_dict2[i], color = 'black', bins = 100)

axs[0,1].title(i)
axs[0,1].xlabel('Percent of reads identical to reference genome')
axs[0,1].ylabel('number of loci')









#076
i = 'Pvag_2_76.sort'
data = np.array(ratio_dict[i])
x, y = data.T

axs[1,0].scatter(x,y,color = 'k', alpha=0.1)

axs[1,0].ylabel('X genome')
axs[1,0].xlabel('D genome')
axs[1,0].title(i)
axs[1,0].gca().set_aspect('equal','box')
x_lim = axs[1,0].xlim()
axs[1,0].plot([0,x_lim[1]], [0,y_lim[1]], color='r', linestyle='--', linewidth=1)
axs[1,0].plot([0,x_lim[1]], [0,y_lim[1]*2], color='r', linestyle='--', linewidth=1)
axs[1,0].plot([0,x_lim[1]], [0,y_lim[1]/2], color='r', linestyle='--', linewidth=1)



axs[1,1].hist(ratio_dict2[i], color = 'black', bins = 100)

axs[1,1].title(i)
axs[1,1].xlabel('Percent of reads identical to reference genome')
axs[1,1].ylabel('number of loci')






#107

i = 'Pvag_2_107.sort'
data = np.array(ratio_dict[i])
x, y = data.T

axs[2,0].scatter(x,y,color = 'k', alpha=0.1)

axs[2,0].ylabel('X genome')
axs[2,0].xlabel('D genome')
axs[2,0].title(i)
axs[2,0].gca().set_aspect('equal','box')
x_lim = axs[2,0].xlim()
axs[2,0].plot([0,x_lim[1]], [0,y_lim[1]], color='r', linestyle='--', linewidth=1)
axs[2,0].plot([0,x_lim[1]], [0,y_lim[1]*2], color='r', linestyle='--', linewidth=1)
axs[2,0].plot([0,x_lim[1]], [0,y_lim[1]/2], color='r', linestyle='--', linewidth=1)



axs[2,1].hist(ratio_dict2[i], color = 'black', bins = 100)

axs[2,1].title(i)
axs[2,1].xlabel('Percent of reads identical to reference genome')
axs[2,1].ylabel('number of loci')




fig.tight_layout()
plt.show()		

