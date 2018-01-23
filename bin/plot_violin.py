import os
import numpy as np
from subprocess import call
import matplotlib.pyplot as plt
matplotlib.use('Agg')
################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return data0

################################################################################################
### plot violin
def plot_violin(input_file_list, outputname, log2, small_num):
	input_file_list = open(input_file_list, 'r')
	signal_track_list = []
	filename_list = []
	for records in input_file_list:
		filename = records.split()[0]
		filename_list.append(records.split()[1])
		### read file
		signal_track = read2d_array(filename, float)
		if log2=='T':
			signal_track = np.log2(signal_track + small_num)
		signal_track_list.append(signal_track[:,0])

	print((signal_track_list[0]))
	print((signal_track_list[1]))
	### plot violin plot
	pos = range(1,len(filename_list)+1)
	print('plot violinplot of index:' + outputname)
	plt.figure()
	plt.violinplot(signal_track_list, pos, points=20, widths=0.5, showmeans=True, showextrema=False, showmedians=False)
	plt.savefig(outputname + '.violin.png')

############################################################################
#time python plot_violin.py -i  -a 4 -s all5cell.binary_matrix.png.binary.counts.thresh.bed -c 2 -o B_SPL.meansig.signif.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:o:l:s:")
	except getopt.GetoptError:
		print 'time python index_label2meansig.py -i input_file_list -o outputname -l log2 -s small_num'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python index_label2meansig.py -i input_file_list -o outputname -l log2 -s small_num'
			sys.exit()
		elif opt=="-i":
			input_file_list=str(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())
		elif opt=="-l":
			log2=str(arg.strip())
		elif opt=="-s":
			small_num=float(arg.strip())		


	plot_violin(input_file_list, outputname, log2, small_num)

if __name__=="__main__":
	main(sys.argv[1:])

