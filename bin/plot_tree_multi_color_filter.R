library(networkD3, lib.loc='/storage/home/gzx103/work/r_package/')
#library(data.tree)
library(igraph, lib.loc='/storage/home/gzx103/work/r_package/')

####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_matrix_file = args[1]
filter_matrix_file = args[2]
cd_tree = args[3]
signal_range_color_file = args[4]
signal_input_list = args[5]
signal_matrix_start_col = args[6]
filter_matrix_start_col = args[7]
log2 = args[8]
smallnum = as.numeric(args[9])

#################################################### 
############ read input files
####################################################
######
### read filter matrix file
filter_matrix_od = as.matrix(read.table(filter_matrix_file, header=FALSE))
### extract filter matrix without info
filter_matrix = filter_matrix_od[ , filter_matrix_start_col:dim(filter_matrix_od)[2] ]
### convert to numeric matrix
class(filter_matrix) = 'numeric'
### filter range
filter_range = max(filter_matrix) - min(filter_matrix)
### subtract min(filter matrix) (Entropy smaller -> less white filter)
filter_percent = (filter_matrix - min(filter_matrix) ) / filter_range

### read signal matrix file
signal_matrix_od = as.matrix(read.table(signal_matrix_file, header=FALSE))
### extract signal matrix without info
signal_matrix = signal_matrix_od[ , signal_matrix_start_col:dim(signal_matrix_od)[2] ]
### get index_set name
index_set_name = signal_matrix_od[,1]
### convert to numeric matrix
class(signal_matrix) = 'numeric'

###### read colnames file
colname_file = read.table(signal_input_list, header=F)
colname = colname_file[,2]

### read cell development tree file
tree = read.table(cd_tree, header = F, sep=',')
tree.df = as.data.frame(tree)
colnames(tree.df) = c('Node.1', 'Node.2')

### get color list

signal_matrix_color = signal_matrix
###### read color file
signal_range_color = as.matrix(read.table(signal_range_color_file, header=F))

### replace NA by background signal in background color table
signal_matrix[is.na(signal_matrix)] = as.numeric(signal_range_color[dim(signal_range_color)[1],1])

### log2 transform
if (log2=='T'){
	signal_matrix = log2(signal_matrix+smallnum)
}

### add colnames
for (i in seq(1,dim(signal_range_color)[1])){
	### read signal range and color range
	sig_range_high = signal_range_color[i,1]
	class(sig_range_high) = 'numeric'
	sig_range_low = signal_range_color[i,2]
	class(sig_range_low) = 'numeric'
	sig_range_high_col = unlist(strsplit(signal_range_color[i,3], split = ','))
	class(sig_range_high_col) = 'numeric'
	sig_range_high_col = sig_range_high_col / 255
	sig_range_low_col = unlist(strsplit(signal_range_color[i,4], split = ','))
	class(sig_range_low_col) = 'numeric'
	sig_range_low_col = sig_range_low_col / 255
	### bg color
	sig_range_bg_col = unlist(strsplit(signal_range_color[dim(signal_range_color)[1],4], split = ','))
	class(sig_range_bg_col) = 'numeric'
	sig_range_bg_col = sig_range_bg_col / 255
	#########
	#########
	### identify signal within this signal range of the color
	signal_color_cell = as.logical((signal_matrix<=sig_range_high) * (signal_matrix>sig_range_low))
	#########
	#########
	col_range = sig_range_high - sig_range_low
	###### get r channel score
	### signal to color
	r = sig_range_high_col[1] - ((sig_range_high-signal_matrix) / col_range) * (sig_range_high_col[1] - sig_range_low_col[1])
	### add filter: od_color * (1-filter) + bg_color * filter
	r = r * (1-filter_percent)+ sig_range_bg_col * filter_percent
	### set upper & lower limit for r
	r[r>1] = 1
	r[r<0] = 0
	###### get g channel score
	### signal to color
	g = sig_range_high_col[2] - ((sig_range_high-signal_matrix) / col_range) * (sig_range_high_col[2] - sig_range_low_col[2])
	### add filter: od_color * (1-filter) + bg_color * filter
	g = g * (1-filter_percent)+ sig_range_bg_col * filter_percent
	### set upper & lower limit for g
	g[g>1] = 1
	g[g<0] = 0
	###### get b channel score
	### signal to color
	b = sig_range_high_col[3] - ((sig_range_high-signal_matrix) / col_range) * (sig_range_high_col[3] - sig_range_low_col[3])
	### add filter: od_color * (1-filter) + bg_color * filter
	b = b * (1-filter_percent)+ sig_range_bg_col * filter_percent
	### set upper & lower limit for b
	b[b>1] = 1
	b[b<0] = 0
	###### creat rgb color matrix
	color_matrix_tmp = NULL
	for (i in seq(1,dim(r)[2])){
		#print(i)
		### convert each r,g,b vector to a rgb matrix column
		cor_col_tmp = rgb( r[,i], g[,i], b[,i] )
		### cbind column to a rgb matrix
		color_matrix_tmp = cbind(color_matrix_tmp, cor_col_tmp)
	}
	###### replace color matrix cell values
	signal_matrix_color[ signal_color_cell ] = color_matrix_tmp[signal_color_cell]
}

### plot trees
for (i in seq(1,dim(signal_matrix_color)[1])){
	### get color vector from color matrix
	value_col = signal_matrix_color[i,]

	### get tree
	tree.igraph = graph.data.frame(tree.df, directed=TRUE)
	tree_names = V(tree.igraph)$name
	### sort colnames by tree nodes id
	match_id = match(tree_names, colname)
	V(tree.igraph)$color = value_col[match_id]
	V(tree.igraph)$size = 25

	png(paste(toString(i-1), '.', signal_input_list, index_set_name[i], '.tree.sh.png', sep = ''), width = 1200, height = 1200)
	plot(tree.igraph, layout = layout_as_tree(tree.igraph, root=c(1)))
	dev.off()
}

