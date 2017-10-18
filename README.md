# CD_viewer

##################################
### required library

#### python
###### numpy

#### R
###### pheatmap


##################################
###### set script folder; input folder; output folder
###### in the "runall.sh" file

#### CD-veiwer bin folder
```
script_folder='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/cell_develop_visual/bin/'
```
#### input data folder
```
input_folder='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/cell_develop_visual/test_data/'
```
#### ouput data folder
```
output_folder='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/cell_develop_visual/test_data_output/'
```

##################################
###### set script folder; input folder; output folder

#### input matrix: "celltype.signal.txt"
###### The first column is the unique peak id
###### The rest columns have the target epigenetic mark signal is the corresponding cell type
###### The first line is the header line:
```
name	cell1	cell2	cell3
1	0.0	1.2	0.5
2	0.1	0.0	0.6
3	0.5	1.6	1.6
...
```

#### cell type order table: "celltype.order.txt"
###### The first column is the column number in the input matrix
###### The second the column is the cell type name in the ouput matrix
###### The row order is the output cell type order
```
1	LSK
2	CMP
4	MEP
3	GMP
...
```

#### signal level table: "signal_level_range.txt"
#### This method needs to put the epigenetic marker signal into different clusters.
#### The "signal_level_range.txt" includes the user difined signal cluster and the corresponding signal range:
```
cluster	low_lim	high_lim
0.0	0.0	0.6
1.0	0.6	1.0
2.0	1.0	1.6
...
```

###### scripts in "runall.sh"
#### get index sets (CORE!!!)
###### input signal matrix: $input_folder'celltype.signal.txt'
###### user defined cell type order table: $input_folder'celltype.order.txt'
###### user defined signal cluster and the corresponding signal range: $input_folder'signal_level_range.txt'
###### sorted index matrix output file name : $output_folder'celltype.index.sorted.txt'
###### index set matrix output file name: $output_folder'celltype.index_set.sorted.txt'
```
time python $script_folder'get_index_set.py' -i $input_folder'celltype.signal.txt' -r $input_folder'celltype.order.txt' -l $input_folder'signal_level_range.txt' -f $output_folder'celltype.index.sorted.txt' -s $output_folder'celltype.index_set.sorted.txt'
```

#### plot histogram of number of DNA regions in each index set
###### all sorted index set matrix (ready for plot): $output_folder'celltype.index_set.sorted.txt'

###### histogram file name: $output_folder'index_hist.pdf'
```
time Rscript $script_folder'plot_index_set_region_hist.R' $output_folder'celltype.index_set.sorted.txt' $output_folder'index_hist.pdf' 
```

#### plot index set
###### all sorted index set matrix (ready for plot): $output_folder'celltype.index_set.sorted.txt'
###### sorted index matrix (ready for plot): $output_folder'celltype.index.sorted.txt'
###### index set matrix that passed the threshold (ready for plot): $output_folder'celltype.index_set_filtered.sorted.txt'
###### The heatmap file name of all index sets: $output_folder'index_set_all.pdf'
###### The heatmap file name of the index sets that pass the threshold: $output_folder'index_set_thresh.pdf'
###### The heatmap file name of all DNA regions in the input matrix: $output_folder'index.png'
###### The threshold for number of DNA regions in index set: 200 (the plotted index sets require more than >200 DNA regions in each of them)
###### Output heatmap color: black

```
time Rscript $script_folder'plot_index_set_module.R' $output_folder'celltype.index_set.sorted.txt' $output_folder'celltype.index.sorted.txt' $output_folder'celltype.index_set_filtered.sorted.txt' $output_folder'index_set_all.pdf' $output_folder'index_set_thresh.pdf' $output_folder'index.png' 200 black
```

### run CD-viewer
```
time bash runall.sh
```