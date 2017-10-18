##################################
script_folder='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/cell_develop_visual/bin/'
input_folder='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/cell_develop_visual/test_data/'
output_folder='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/cell_develop_visual/test_data_output/'

### mkdir output folder
if [ -d "$output_folder" ]; then  
rm -r $output_folder
fi
mkdir $output_folder

### get index sets (CORE!!!)
echo 'get index sets (CORE!!!)'
time python $script_folder'get_index_set.py' -i $input_folder'celltype.signal.txt' -r $input_folder'celltype.order.txt' -l $input_folder'signal_level_range.txt' -f $output_folder'celltype.index.sorted.txt' -s $output_folder'celltype.index_set.sorted.txt'

### plot index set
echo 'plot index set peak number distribution'
time Rscript $script_folder'plot_index_set_region_hist.R' $output_folder'celltype.index_set.sorted.txt' $output_folder'index_hist.pdf' $output_folder'index_hist_noxlim.pdf'
echo 'plot index set'
time Rscript $script_folder'plot_index_set_module.R' $output_folder'celltype.index_set.sorted.txt' $output_folder'celltype.index.sorted.txt' $output_folder'celltype.index_set_filtered.sorted.txt' $output_folder'index_set_all.pdf' $output_folder'index_set_thresh.pdf' $output_folder'index.png' 200 black
