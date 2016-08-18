group_by_cancer_type <- function(parent_dir, cancer_name)
{
	# Read sample sheet csv
	print(paste0("group_by_cancer_type script - READING SAMPLE SHEET - ",Sys.time()))
	parent.dir <- file.path(parent_dir)
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

	# Create cancer directory
	print(paste0("group_by_cancer_type script - CREATE CANCER DIR - ",Sys.time()))
	cancer.dir <- file.path(betas.by_cancer.dir, cancer_name)
	dir.create(cancer.dir, showWarnings=FALSE)

	# Get all sample types of current cancer type	
	print(paste0("group_by_cancer_type script - LOOP THROUGH SAMPLE TYPES - ",Sys.time()))
	sample.type <- unique(ss$sample_type[ss$cancer_type==cancer_name])
	
	if (length(sample.type) > 0) {
		for (j in 1:length(sample.type)) {
			
			print(paste0("group_by_cancer_type script - CREATING DATA FOR ",sample.type[j]," - ",Sys.time()))
			cancer.arrays <- unique(array.nums[ss$cancer_type==cancer_name & ss$sample_type == sample.type[j]])	
			for ( k in 1:length(cancer.arrays)){
				# Read current array
				curr.array <- read.csv(file.path(parent.dir,"betas","by_array", paste0(cancer.arrays[k],".csv")),stringsAsFactors=F, row.names=1, check.names=F)
		
				curr.sample.type <- ss$Array.Data.File[array.nums==cancer.arrays[k] & ss$cancer_type==cancer_name & ss$sample_type==sample.type[j]]
				sample.type.filt <- intersect(curr.sample.type,colnames(curr.array))
	
				if (length(sample.type.filt) > 0) {
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
	
			}	
			
			print(paste0("group_by_cancer_type script - CREATING FILE FOR ",sample.type[j]," - ",Sys.time()))
			sample.type.file.dir <- file.path(cancer.dir, paste0(sample.type[j],".csv"))
			write.csv(data,sample.type.file.dir,quote=FALSE)
		}
	}
	print(paste0("group_by_cancer_type script - FINISHED - ",Sys.time()))
}

args<-commandArgs(trailingOnly=TRUE)
group_by_cancer_type(args[1], args[2])



