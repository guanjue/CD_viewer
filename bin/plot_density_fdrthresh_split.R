####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
bed_file = args[1]
output_file = args[2]
uselog2 = args[3]
smallnum = as.numeric(args[4])
fdr_pval_thresh = as.numeric(args[5])
line_num = as.numeric(args[6])
file_list = args[7]
### read file
bed = as.data.frame(read.table(bed_file, header=FALSE))
sig = as.numeric(bed[,5])
print(bed_file)
if (uselog2=='T'){
	sig = log2(sig+smallnum)
}

### get mean & std
sig_mean = mean(sig)
sig_sd = sd(sig)
print(paste('mean: ', toString(sig_mean)))
print(paste('std: ', toString(sig_sd)))
### get z score
sig_z = (sig - sig_mean) / sig_sd

### get p-value
sig_pval = pnorm(sig_z, lower.tail = FALSE)

### get fdr p-value
sig_fdr_pval = p.adjust(sig_pval, 'fdr')

### one-side test
fdr_pval_thresh = fdr_pval_thresh

sig_0 = sig[sig_pval >= fdr_pval_thresh]
sig_1 = sig[sig_pval < fdr_pval_thresh]

print(sum(sig_pval >= fdr_pval_thresh) / length(sig_pval))
### plot histogram
png(output_file)
plot(density(sig))
abline(v=(max(sig_0)), col='red')
print((max(sig_0)))

### write binary bed file
sig_binary = sig
sig_binary[sig_fdr_pval >= fdr_pval_thresh] = 0
sig_binary[sig_fdr_pval < fdr_pval_thresh] = 1

bed_sig_binary = as.data.frame(cbind(bed[,c(1,2,3,4)], sig_binary))

head(bed_sig_binary)

### read file list
file_list = read.table(file_list, header=FALSE)
n = 0
for (file in file_list){
	write.table(bed_sig_binary[1:line_num,], paste(file_list[n+1], '.binary.allbins.bed', sep=''), sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
	n = n+1
}

