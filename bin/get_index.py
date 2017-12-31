import os
import numpy as np
from subprocess import call
from scipy import stats
from collections import Counter
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
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

################################################################################################
### get convert bedtools window output to matrix of pk and intersect ideas label info
def ideas_label_info(input_bedtools_window, id_col, lb_col, pk_col, ideas_col):
	data_info0=open(input_bedtools_window, 'r')
	### read DNA region orders
	data_info=[]
	for records in data_info0:
		tmp=[x.strip() for x in records.split('\t')]
		### get intersect region; midpoint dist; TF peak length
		if ((int(tmp[ideas_col-1]) - int(tmp[pk_col-1]))>=0) and ((int(tmp[ideas_col]) - int(tmp[pk_col]))<=0) :
			### IDEAS Bin >= pk region
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[ideas_col])-int(tmp[ideas_col-1]), (float(tmp[ideas_col])+float(tmp[ideas_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[ideas_col])-int(tmp[ideas_col-1]) ]
		elif ((int(tmp[ideas_col-1]) - int(tmp[pk_col-1]))<0) and ((int(tmp[ideas_col]) - int(tmp[pk_col]))>0) :
			### IDEAS Bin < pk region
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[pk_col])-int(tmp[pk_col-1]), (float(tmp[ideas_col])+float(tmp[ideas_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[ideas_col])-int(tmp[ideas_col-1]) ]
		elif ((int(tmp[ideas_col-1]) - int(tmp[pk_col-1]))<0) and ((int(tmp[ideas_col]) - int(tmp[pk_col]))<=0) :
			### IDEAS Bin upstream < pk region upstream & IDEAS Bin downstream <= pk region downstream
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[ideas_col])-int(tmp[pk_col-1]), (float(tmp[ideas_col])+float(tmp[ideas_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[ideas_col])-int(tmp[ideas_col-1]) ]
		elif ((int(tmp[ideas_col-1]) - int(tmp[pk_col-1]))>=0) and ((int(tmp[ideas_col]) - int(tmp[pk_col]))>0) :
			### IDEAS Bin upstream >= pk region upstream & IDEAS Bin downstream > pk region downstream
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[pk_col])-int(tmp[ideas_col-1]), (float(tmp[ideas_col])+float(tmp[ideas_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[ideas_col])-int(tmp[ideas_col-1]) ]
		data_info.append(tmp_vec)
	data_info0.close()
	return (data_info)

################################################################################################
### get peak's ideas labels
def get_cRE_ideas_state(data_info_matrix, id_col, lb_col, cover_col, middist_col, ideaslen, bed_od_file, bed_od_idcol, outputname):
	### read DNA region IDEAS state info matrix
	pk_id_list = []
	data_ideas1={}
	data_ideas1_maxcover={} ### coverage size
	data_ideas1_middist={} ### midpoint dist
	data_ideas1_statelen={} ### ideas state len
	### initialize problem counter
	k=0
	for info in data_info_matrix:
		pk_id = info[id_col-1]
		### creat pk_id_list for keeping the id order in output
		pk_id_list.append(pk_id)
		if not (pk_id in data_ideas1):
			data_ideas1[pk_id] = info[lb_col-1]
			data_ideas1_maxcover[pk_id] = info[cover_col-1]
			data_ideas1_middist[pk_id] = info[middist_col-1]
			data_ideas1_statelen[pk_id] = info[ideaslen-1]
		elif info[cover_col-1] > data_ideas1_maxcover[pk_id]:
			### if multiple cover; select the highest covering state
			data_ideas1[pk_id] = info[lb_col-1]
			data_ideas1_maxcover[pk_id] = info[cover_col-1]
			data_ideas1_middist[pk_id] = info[middist_col-1]
			data_ideas1_statelen[pk_id] = info[ideaslen-1]
		elif info[cover_col-1] == data_ideas1_maxcover[pk_id]: ### if 2 states cover the same region with same length
			if info[middist_col-1] < data_ideas1_middist[pk_id]: 
				### if cover the same; check mid point distance
				data_ideas1[pk_id] = info[lb_col-1]
				data_ideas1_maxcover[pk_id] = info[cover_col-1]
				data_ideas1_middist[pk_id] = info[middist_col-1]
				data_ideas1_statelen[pk_id] = info[ideaslen-1]
			elif info[middist_col-1] == data_ideas1_middist[pk_id]: ### if 2 states cover the same region with same length; with same midpoint dist
				if info[ideaslen-1] < data_ideas1_statelen[pk_id]:
					### if cover same & mid point distance same; check state len 
					data_ideas1[pk_id] = info[lb_col-1]
					data_ideas1_maxcover[pk_id] = info[cover_col-1]
					data_ideas1_middist[pk_id] = info[middist_col-1]
					data_ideas1_statelen[pk_id] = info[ideaslen-1]
				else: ### if 2 states cover the same region with same length; with same midpoint dist; with same state length ...give attention!
					k=k+1
					print('problem!')
					print(k)
	### read original bed file to get the pk id list
	bed_od_file=open(bed_od_file,'r')
	bed_od_id_list = []
	for records in bed_od_file:
		bed_od_id_list.append(records.split()[bed_od_idcol-1])
	bed_od_file.close()
	### write ideas label output
	result=open(outputname,'w')
	for pkid in bed_od_id_list:
		if pkid in data_ideas1:
			tmp=data_ideas1[pkid]
			result.write(tmp+'\n')
		else:
			tmp=records
			result.write('NA'+'\n')
	result.close()

################################################################################################
### get index/signal matrix
def get_mark_matrix(peak_bed, peak_bed_colnum, mark_list, output_file, signal_col, method, sort):
	### sort input bed files
	sort_bed_file = peak_bed + '.sort.bed'
	call('sort -k1,1 -k2,2n ' + peak_bed + ' > ' + sort_bed_file, shell=True)
	call('cp ' + sort_bed_file + ' ' + output_file, shell=True)
	### generate index mark matrix
	mark_list_vec = open(mark_list, 'r')
	celltype_list = []
	for mark_bed in mark_list_vec:	
		tmp = [x.strip() for x in mark_bed.split('\t')]
		### read bianry label file list
		mark_bed_file = tmp[0]
		print(mark_bed_file)
		### add cell type name to cell type list
		celltype_list.append(tmp[1])
		### sort bianry label bed files
		if sort == 'T':
			call('sort -k1,1 -k2,2n ' + mark_bed_file + ' > ' + mark_bed_file+'.sort.bed', shell=True)
		else:
			call('cp ' + mark_bed_file + ' ' + mark_bed_file+'.sort.bed', shell=True)
		### use bedtools to generate the index/signal matrix
		if method == 'intersect':
			### used bedtools intersect to get the binary label of each peak
			call('bedtools intersect -c -a ' + sort_bed_file + ' -b ' + mark_bed_file+'.sort.bed' + ' > ' + mark_bed_file+'.tmp01.txt', shell=True)
		elif method == 'map':
			### used bedtools map to get the average signal of each peak
			call('bedtools map -c ' + str(signal_col) + ' -null 0 -F 0.5 -o mean -a ' + sort_bed_file + ' -b ' + mark_bed_file+'.sort.bed' + ' > ' + mark_bed_file+'.tmp01.txt', shell=True)
		elif method == 'window':
			### used bedtools map to get the average signal of each peak
			call('bedtools window -a ' + sort_bed_file + ' -b ' + mark_bed_file+'.sort.bed' + ' -w 0 > ' + mark_bed_file+'.tmp01.txt', shell=True)
			### convert bedtools window output to matrix of pk and intersect ideas label info (intersect region; midpoint dist; TF peak length)
			data_info_matrix = ideas_label_info(mark_bed_file+'.tmp01.txt', 4, 8, 2, 6)
			### get peak's ideas labels based on intersect region; midpoint dist; TF peak length
			get_cRE_ideas_state(data_info_matrix, 1, 2, 3, 4, 5, sort_bed_file, 4, mark_bed_file+'.tmp01.txt')
		### cut the map number column
		call('cut -f'+ str(peak_bed_colnum+1) +" -d$'\t' " + mark_bed_file+'.tmp01.txt' + ' > ' + mark_bed_file+'.tmp02.txt', shell=True)
		### cbind to matrix
		call('paste ' + output_file + ' ' + mark_bed_file+'.tmp02.txt' + ' > ' + output_file+'.tmp.txt' + ' && mv ' + output_file+'.tmp.txt ' + output_file, shell=True)
		### remove tmp files
		call('rm ' + mark_bed_file+'.tmp01.txt' + ' ' + mark_bed_file+'.tmp02.txt' + ' ' + mark_bed_file+'.sort.bed', shell=True)
	mark_list_vec.close()

################################################################################################
### convert index matrix to index counts dict
def index_matrix2index_count_dict(index_matrix, index_matrix_start_col):
	### binary matrix to index X_X_X
	index_vector = []
	for vec in index_matrix:
		index_tmp = ''
		for i in range(index_matrix_start_col+1, len(vec)-1):
			index_tmp = index_tmp + vec[i] + '_'
		index_tmp = index_tmp + vec[len(vec)-1]
		index_vector.append( index_tmp )
	### index_vector 2 np array
	index_vector = np.array(index_vector)
	### index peak counts (dict)
	index_uniq_count_dict = Counter(index_vector)
	### index peak counts dict 2 index peak counts np array
	index_vector_count_vec = []
	for index in index_uniq_count_dict:
		index_vector_count_vec.append(index_uniq_count_dict[ index ])
	index_vector_count_vec = np.array(index_vector_count_vec)
	return {'index_vector': index_vector.reshape(index_vector.shape[0], 1), 'index_uniq_count_dict': index_uniq_count_dict, 'index_vector_count_vec': index_vector_count_vec}

################################################################################################
### get index set counts threshold
def index_count_thresh(count_vec, thesh):
	mean = np.mean(count_vec)
	std = np.std(count_vec)
	n = len(count_vec)
	t_stat = stats.t.ppf(thesh, n)
	index_count_thresh = mean + t_stat * std / np.sqrt(n)
	return(index_count_thresh)

################################################################################################
### use while loop to select the threshold of index set counts
def select_index_set_counts_thresh(index_matrix, index_matrix_start_col, siglevel_counts):
	index_count_dict = index_matrix2index_count_dict(index_matrix, index_matrix_start_col)
	index_vector = index_count_dict['index_vector']
	index_uniq_count_dict = index_count_dict['index_uniq_count_dict']
	index_vector_count_vec = index_count_dict['index_vector_count_vec']
	print('calculating index peak counts thresh...')
	### initalize thresh 
	index_vector_count_vec_select = index_vector_count_vec
	index_count_thresh_1 = max(index_vector_count_vec)
	index_count_thresh_2 = 0
	i = 1
	while index_count_thresh_1 != index_count_thresh_2:
		print('select index peak counts thresh - round: ' + str(i))
		i = i+1
		### use  significant level to find insignificant index counts 
		index_count_thresh_1 = index_count_thresh(index_vector_count_vec_select, siglevel_counts)
		index_vector_count_vec_select = index_vector_count_vec_select[index_vector_count_vec_select<index_count_thresh_1]
		index_count_thresh_2 = index_count_thresh(index_vector_count_vec_select, siglevel_counts)
	print('select index peak counts thresh: ' + str(index_count_thresh_2))
	### extract insignificant index names
	insig_index = []
	for index in index_uniq_count_dict:
		if index_uniq_count_dict[ index ] < index_count_thresh_2:
			insig_index.append( index )
	return { 'index_vector': index_vector, 'insig_index': insig_index, 'index_count_thresh': index_count_thresh_2, 'index_vector_count_vec': index_vector_count_vec }

################################################################################################
### calculating multiple variable norm density score
def mvn_density_score(signal_matrix_od, signal_matrix_start_col, log_signal, small_value, qda_round, index_vector, insig_index):
	print('calculating multiple variable norm density score...')
	### initalize signal matrix for each index (dict)
	index_signal_matrix_dict = {}
	### initalize uniq index vector
	uniq_index = []
	### extract bed files
	signal_matrix_bed = signal_matrix_od[:,range(0, signal_matrix_start_col-1)]
	### extract signal matrix
	signal_matrix = signal_matrix_od[:,range(signal_matrix_start_col-1,signal_matrix_od.shape[1])]
	### convert string matrix to float matrix
	signal_matrix = signal_matrix.astype(float)
	### log transform the data
	if log_signal == 'T':
		signal_matrix = np.log2(signal_matrix+small_value)
	###############
	### QDA start
	for l in range(0, qda_round):
		print('Round: '+ str(l))
		### initialize the index name vector for all peaks
		if l==0:
			### in the 1st round, we use the orginal index to label the peaks
			index_vector_loop = index_vector
		else:
			### after 1st round, we use previous round generated index to label the peaks
			index_vector_loop = index_name_vec
		### convert index to index and insignificant index set
		print('filter insig_index...')
		if l == qda_round-1:
			index_vector_filter = []	
		for i in range(0, len(index_vector_loop)):
			index = index_vector_loop[i, 0]
			signal_vector = signal_matrix[i,:]
			### if not in the insig_index vector, we use the index as one index set
			if not (index in insig_index):
				if index in index_signal_matrix_dict:
					index_signal_matrix_dict[ index ].append( signal_vector )
				else:
					uniq_index.append( index )
					index_signal_matrix_dict[ index ] = [ signal_vector ]
			### if in the insig_index vector, we use the insig_index ('X_X_X...') as one index set
			else:
				index = 'X_X_X...'
				if index in index_signal_matrix_dict:
					index_signal_matrix_dict[ index ].append( signal_vector )
				else:
					uniq_index.append( index )
					index_signal_matrix_dict[ index ] = [ signal_vector ]
			### in last round, we append to index_vector_filter for original index vector relabeling
			if l == qda_round-1:
				index_vector_filter.append(index)
	print('calculating mean vector and cov matrix...')
	### initialize cov dict & mean dict
	index_signal_cov_dict = {}
	index_signal_mean_dict = {}
	### initialize mean matrix & counts matrix for index set output
	index_set_signal_mean_matrix = []
	index_set_peak_counts_matrix = []
	### loop uniq index
	for index in uniq_index:
		print(index)
		### extract index signal matrix
		one_index_matrix = np.array(index_signal_matrix_dict[ index ], dtype = float)
		### calculating index matrix covariance matrix & mean vector
		one_index_matrix_cov = np.cov(one_index_matrix, rowvar = False)
		one_index_matrix_mean = np.mean(one_index_matrix, axis = 0)
		### append to covariance matrix dict & mean vector dict & index set matrix
		index_signal_cov_dict[ index ] = one_index_matrix_cov
		index_signal_mean_dict[ index ] = one_index_matrix_mean
		index_set_signal_mean_matrix.append(one_index_matrix_mean)
		index_set_peak_counts_matrix.append(one_index_matrix.shape[0])

	print('calculating Quadratic Scores...')
	score_i_exp_matrix = np.empty((signal_matrix.shape[0], 0), float)
	for index in uniq_index:
		### extract index signal matrix of one index set
		cov_i = index_signal_cov_dict[ index ]
		cov_i_inverse = np.linalg.inv(cov_i)
		mean_i = index_signal_mean_dict[ index ]
		### calculate log scale score
		d = np.sum(- 0.5 * np.log( abs(cov_i) ))
		score_i = d - 0.5 * np.sum( np.dot((signal_matrix-mean_i), cov_i_inverse) * (signal_matrix-mean_i), axis = 1 )
		### convert to exp scale
		score_i_exp = np.exp(score_i).reshape((score_i.shape[0],1))
		### cbind exp_scale score to score_i_exp_matrix
		score_i_exp_matrix = np.concatenate((score_i_exp_matrix, score_i_exp), axis=1)
	print('calculating Quadratic Scores...DONE')
	print('calculating max p and index...') 
	index_name_vec = []
	index_p_vec = []
	for p_vec in score_i_exp_matrix:
		index_i = uniq_index[np.argmax(p_vec)]
		index_name_vec.append(index_i)
		p_max = np.max(p_vec) / np.sum(p_vec)
		index_p_vec.append(p_max)
	print('calculating max p and index...DONE') 
	index_name_vec = np.array(index_name_vec).reshape(len(index_name_vec),1)

	return { 'signal_matrix_bed': signal_matrix_bed, 'index_name_vec': index_name_vec, 'index_p_vec': index_p_vec, 'signal_matrix': signal_matrix, 'uniq_index': uniq_index, 'index_set_peak_counts_matrix': index_set_peak_counts_matrix, 'index_set_signal_mean_matrix': index_set_signal_mean_matrix, 'index_vector_filter': index_vector_filter }

################################################################################################


################################################################################################
### get Binary index matrix
peak_bed = 'homerTable3.peaks.filtered.interval.bed'
peak_bed_colnum = 4
mark_list_index = 'peak_list.txt'
output_file = 'homerTable3.peaks.filtered.interval.bed.index.matrix.txt'
signal_col = 'N/A'
method = 'intersect'
sort_sigbed = 'T'
print('get binary matrix...')
get_mark_matrix(peak_bed, peak_bed_colnum, mark_list_index, output_file, signal_col, method, sort_sigbed)

### get TF ChIP-seq matrix
peak_bed = 'homerTable3.peaks.filtered.interval.bed'
peak_bed_colnum = 4
mark_list_chip = 'chip_list.txt'
output_file = 'homerTable3.peaks.filtered.interval.bed.chip.matrix.txt'
signal_col = 'N/A'
method = 'intersect'
sort_sigbed = 'T'
print('get chip matrix...')
get_mark_matrix(peak_bed, peak_bed_colnum, mark_list_chip, output_file, signal_col, method, sort_sigbed)

### get ideas label matrix
peak_bed = 'homerTable3.peaks.filtered.interval.bed'
peak_bed_colnum = 0
mark_list_chip = 'ideas_list.txt'
output_file = 'homerTable3.peaks.filtered.interval.bed.ideas.matrix.txt'
signal_col = 'N/A'
method = 'window'
sort_sigbed = 'F'
print('get ideas matrix...')
get_mark_matrix(peak_bed, peak_bed_colnum, mark_list_chip, output_file, signal_col, method, sort_sigbed)

### get signal matrix
peak_bed = 'homerTable3.peaks.filtered.interval.bed'
peak_bed_colnum = 4
mark_list_signal = 'signal_list.txt'
output_file = 'homerTable3.peaks.filtered.interval.bed.signal.matrix.txt'
signal_col = 5
method = 'map'
sort_sigbed = 'T'
print('get signal matrix...')
get_mark_matrix(peak_bed, peak_bed_colnum, mark_list_signal, output_file, signal_col, method, sort_sigbed)

### Multi-variable norm p-value (QDA)
index_matrix_start_col = 5
signal_matrix_start_col = 5
siglevel_counts = 0.95
small_value = 0.001
log_signal = 'T'
qda_round = 2
bins_folder = '/Volumes/MAC_Data/data/labs/zhang_lab/01projects/CD_viewer/bin/'
index_matrix = read2d_array('homerTable3.peaks.filtered.interval.bed.index.matrix.txt', 'str')
signal_matrix_od = read2d_array('homerTable3.peaks.filtered.interval.bed.signal.matrix.txt', 'str')

### use while loop to select the threshold of index set counts
index_matrix = read2d_array('homerTable3.peaks.filtered.interval.bed.index.matrix.txt', 'str')
insig_index_dict = select_index_set_counts_thresh(index_matrix, index_matrix_start_col)
index_vector = insig_index_dict['index_vector']
insig_index = insig_index_dict['insig_index']
index_count_thresh_2 = insig_index_dict['index_count_thresh']
index_vector_count_vec = insig_index_dict['index_vector_count_vec']

### calculating multiple variable norm density score
mvn_density_score_dict = mvn_density_score(signal_matrix_od, signal_matrix_start_col, log_signal, small_value, qda_round, index_vector, insig_index)
signal_matrix_bed = mvn_density_score_dict['signal_matrix_bed']
index_name_vec = mvn_density_score_dict['index_name_vec']
index_p_vec = mvn_density_score_dict['index_p_vec']
signal_matrix = mvn_density_score_dict['signal_matrix']
uniq_index = mvn_density_score_dict['uniq_index']
index_set_peak_counts_matrix = mvn_density_score_dict['index_set_peak_counts_matrix']
index_set_signal_mean_matrix = mvn_density_score_dict['index_set_signal_mean_matrix']
index_vector_filter = mvn_density_score_dict['index_vector_filter']


print('write output matrix files')
### index sort peak signal matrix
### reshape vector to N X 1 matrix
index_p_vec = np.array(index_p_vec).reshape(signal_matrix.shape[0],1)
index_name_vec = np.array(index_name_vec).reshape(signal_matrix.shape[0],1)
### cbind
signal_matrix_indexed = np.concatenate((signal_matrix_bed, index_name_vec, index_p_vec, signal_matrix), axis = 1)
#signal_matrix_indexed = np.concatenate((index_name_vec, index_p_vec), axis = 1)
#signal_matrix_indexed = np.concatenate((signal_matrix_indexed, signal_matrix), axis = 1)
### write output
write2d_array(signal_matrix_indexed, 'signal_matrix_indexed.txt')
## sort based on index and p-value
call('sort -k5,5 -k6,6rn signal_matrix_indexed.txt > signal_matrix_indexed_sort.txt', shell=True)
call('rm signal_matrix_indexed.txt', shell=True)
call('time Rscript ' + bins_folder + 'plot_pk_p_sig.R signal_matrix_indexed_sort.txt 6 orange 7 red '+ mark_list_signal + ' 2 signal_matrix_indexed_sort_hm', shell=True)	

### index set mean matrix
### reshape vector to N X 1 matrix
uniq_index = np.array(uniq_index).reshape(len(uniq_index),1)
index_set_peak_counts_matrix = np.array(index_set_peak_counts_matrix).reshape(len(index_set_peak_counts_matrix),1)
index_set_signal_mean_matrix = np.array(index_set_signal_mean_matrix)
### cbind
index_set_signal_mean_matrix = np.concatenate((uniq_index, index_set_peak_counts_matrix, index_set_signal_mean_matrix), axis = 1)
### write output
write2d_array(index_set_signal_mean_matrix, 'index_set_signal_mean_matrix.txt')
## sort based on index 
call('sort -k1,1 index_set_signal_mean_matrix.txt > index_set_signal_mean_matrix_sort.txt', shell=True)
call('rm index_set_signal_mean_matrix.txt', shell=True)

### orignal index with insig_index relabeled
### reshape vector to N X 1 matrix
index_vector_filter = np.array(index_vector_filter).reshape(len(index_vector_filter),1)
### cbind
signal_matrix_indexed_od = np.concatenate((signal_matrix_bed, index_vector_filter, signal_matrix), axis = 1)
### write output
write2d_array(signal_matrix_indexed_od, 'signal_matrix_indexed_od.txt')
## sort based on index 
call('sort -k5,5 signal_matrix_indexed_od.txt > signal_matrix_indexed_od_sort.txt', shell=True)
call('rm signal_matrix_indexed_od.txt', shell=True)	

### index_vector_count_vec and counts threshold
index_vector_count_vec = np.append(index_vector_count_vec, index_count_thresh_2).reshape(index_vector_count_vec.shape[0]+1, 1)
write2d_array(index_vector_count_vec, 'index_vector_count_vec.txt')
### plot histogram
call('time Rscript ' + bins_folder + 'plot_index_set_counts_hist.R ' + 'index_vector_count_vec.txt' + ' ' + 'index_vector_count_vec.png', shell=True)	









