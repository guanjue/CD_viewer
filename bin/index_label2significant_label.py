import os
import numpy as np
################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	### get colnames
	colnames = data.readline()
	colnames = [x.strip() for x in colnames.split(' ')]
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split(' ')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return {'data0': data0, 'colnames': colnames}

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

################################################################################################
### ideas_matrix2ideas_bed
def index_label2significant_label(od_index_label_table, significant_index_label, significant_index_label_col, outputname):
	### read od_index_label matrix & significant label
	od_index_label = read2d_array(od_index_label_table, str)
	significant_index_label = read2d_array(significant_index_label, str)
	### signif index_label to dict
	significant_index_label_dict = {}
	for records in significant_index_label:
		label_tmp = records[significant_index_label_col-1]
		significant_index_label_dict[label_tmp] = ''
	### get insignif label
	insignif_label = ''
	print(od_index_label[0])
	eg = od_index_label[0][0].split('_')
	for i in range(0,len(eg)-1):
		insignif_label = insignif_label + 'X_'
	insignif_label = insignif_label + 'X'
	### od to signif
	od2signif_index_label = []
	for records in od_index_label:
		if records[0] in significant_index_label_dict:
			od2signif_index_label.append(records[0])
		else:
			od2signif_index_label.append(insignif_label)
	### write output 
	r1=open(outputname,'w')
	for records in od2signif_index_label:
		r1.write(str(records)+'\n')
	r1.close()

############################################################################
#time python index_label2significant_label.py -i B_SPL.binary_matrix.txt -s all5cell.binary_matrix.png.binary.counts.thresh.bed -c 2 -o B_SPL.binary_matrix.signif.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:s:c:o:")
	except getopt.GetoptError:
		print 'time python index_label2significant_label.py -i od_index_label_table -s significant_index_label -c significant_index_label_col -o outputname'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python index_label2significant_label.py -i od_index_label_table -s significant_index_label -c significant_index_label_col -o outputname'
			sys.exit()
		elif opt=="-i":
			od_index_label_table=str(arg.strip())
		elif opt=="-s":
			significant_index_label=str(arg.strip())		
		elif opt=="-c":
			significant_index_label_col=int(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())

	index_label2significant_label(od_index_label_table, significant_index_label, significant_index_label_col, outputname)

if __name__=="__main__":
	main(sys.argv[1:])