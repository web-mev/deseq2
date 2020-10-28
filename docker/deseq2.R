if(!require("DESeq2", character.only=T)) stop("Please install the DESeq2 package first.")

# args from command line:
args<-commandArgs(TRUE)
RAW_COUNT_MATRIX<-args[1]
BASE_CONDITION_SAMPLES <- args[2]
EXPERIMENTAL_CONDITION_SAMPLES <- args[3]
CONDITION_A<-args[4]
CONDITION_B<-args[5]
OUTPUT_DESEQ_FILE_BASE <- 'deseq2_results'

# change the working directory to co-locate with the counts file:
working_dir <- dirname(RAW_COUNT_MATRIX)
setwd(working_dir)

# convert names in case they are not "proper" for R
CONDITION_A <- make.names(CONDITION_A)
CONDITION_B <- make.names(CONDITION_B)

# create a string to identify the contrast:
contrast_str = paste0(CONDITION_B, '_vs_', CONDITION_A)

# the sample names are given as a comma-delimited string. Split them
base_samples <- make.names(strsplit(BASE_CONDITION_SAMPLES, ',')[[1]])
exp_samples <- make.names(strsplit(EXPERIMENTAL_CONDITION_SAMPLES, ',')[[1]])
all_samples <- c(base_samples, exp_samples)

condition_a_list <- rep(CONDITION_A, length(base_samples))
condition_b_list <- rep(CONDITION_B, length(exp_samples))
all_conditions <- c(condition_a_list, condition_b_list)

# full annotation matrix:
annotations <- as.data.frame(cbind(all_samples, all_conditions), stringsAsFactors = F)
colnames(annotations) <- c('sample', 'condition')

# read the raw count matrix, genes as row names:
count_data <- read.table(RAW_COUNT_MATRIX, sep='\t', header = T, row.names = 1, stringsAsFactors = F)

# subset to keep only the samples in the count table.  This is important if the annotation
# file has more samples than the count matrix:
count_mtx_cols = colnames(count_data)
annotations <- annotations[annotations$sample %in% count_mtx_cols,]

# subset to only keep samples corresponding to the current groups in the count_data dataframe
count_data <- count_data[,annotations$sample]

# DESeq2 expects that the rownames of the annotation data frame are the sample names.  Set the rownames and drop that col
rownames(annotations) <- annotations$sample
annotations <- annotations[-1]

# Need to set the condition as a factor since it's going to be used as a design matrix
annotations$condition <- as.factor(annotations$condition)

# run the actual differential expression:
dds <- DESeqDataSetFromMatrix(countData = count_data,
							  colData = annotations,
							  design = ~condition)

dds <- DESeq(dds)
res <- results(dds, cooksCutoff=F, contrast=c("condition", CONDITION_B, CONDITION_A))
original_colnames = colnames(res)
n = length(original_colnames)
baseMeanA = rowMeans(counts(dds,normalized=TRUE)[,dds$condition == CONDITION_A]) 
baseMeanB = rowMeans(counts(dds,normalized=TRUE)[,dds$condition == CONDITION_B]) 
res = cbind(rownames(res), res[,1],baseMeanA, baseMeanB, as.data.frame(res[,2:n])) 
colnames(res) = c('Gene', 'overall_mean', CONDITION_A, CONDITION_B, original_colnames[2:n])
resOrdered <- res[order(res$padj),]

# extract and output the normalized counts:
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
nc <- cbind(gene=rownames(nc), nc)

# merge to create a single table
m <- merge(resOrdered, nc, by="row.names")
rownames(m) <- m[,'Row.names']
drops <- c("Row.names", "Gene", "gene")
m <- m[, !(names(m) %in% drops)]
output_filename <- paste(OUTPUT_DESEQ_FILE_BASE, contrast_str, 'tsv', sep='.')
output_filename <- paste(working_dir, output_filename, sep='/')
write.table(m, output_filename, sep='\t', quote=F)

json_str = paste0(
       '{"dge_results":"', output_filename, '"}'
)
write(json_str, 'outputs.json')

