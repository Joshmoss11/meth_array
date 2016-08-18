group_by_cancer_type <- function(parent_dir, data_source)
{
	# Read sample sheet csv
	parent.dir <- file.path(parent_dir,data_source)
	ss <- read.csv(file.path(parent.dir,'sample_sheet_all.csv'),stringsAsFactors=F)

	# Create directory
	betas.by_cancer.dir <- file.path(parent.dir,'by_cancer')
	dir.create(betas.by_cancer.dir, showWarnings=FALSE)

	# Replace sample type column to one of 3 values
	ss <- ss[!ss$sample_type =="no_info",]
	ss$sample_type <- gsub("Additional Metastatic", "Metastatic",ss$sample_type)
	ss$sample_type <- gsub("Recurrent Tumor", "Primary Tumor",ss$sample_type)
	ss$sample_type <- gsub("Additional - New Primary", "Primary Tumor",ss$sample_type)

	# Get array numbers
	array.nums <- substr(ss$Array.Data.File,1,10)

	cancerT <- unique(ss$cancer_type)

	for (i in 1:length(cancerT)){
		# Create cancer directory
		cancer.dir <- file.path(betas.by_cancer.dir, cancerT[i])
		dir.create(cancer.dir, showWarnings=FALSE)

		# Get all sample types of current cancer type	
		sample.type <- unique(ss$sample_type[ss$cancer_type==cancerT[i]])

		for (j in 1:length(sample.type)) {
	
			cancer.arrays <- unique(array.nums[ss$cancer_type==cancerT[i] & ss$sample_type == sample.type[j]])	
			for ( k in 1:length(cancer.arrays)){
				# Read current array
				curr.array <- read.csv(file.path(parent.dir,"betas","by_array", paste0(cancer.arrays[k],".csv")),stringsAsFactors=F, row.names=1, check.names=F)
	
				curr.sample.type <- ss$Array.Data.File[array.nums==cancer.arrays[k] & ss$cancer_type==cancerT[i] & ss$sample_type==sample.type[j]]
				if (k==1){
					data <- curr.array[,curr.sample.type]
				}else {
					if (length(curr.sample.type) == 1) {
						temp.data <- as.data.frame(curr.array[,curr.sample.type])
						colnames(temp.data) <- curr.sample.type
						data <- cbind(data, temp.data)	
					}else{
						data <- cbind(data, curr.array[,curr.sample.type])
					}
				}	
	
			}	

			sample.type.file.dir <- file.path(cancer.dir, paste0(sample.type[j],".csv"))
			write.csv(data,sample.type.file.dir,quote=FALSE)
		}
	}
}

args<-commandArgs(trailingOnly=TRUE)
group_by_cancer_type(args[1], args[2])



