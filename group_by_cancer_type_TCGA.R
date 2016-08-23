group_by_cancer_type <- function(parent.dir, cancer_name)
{
	# Read sample sheet csv
	print(paste0("group_by_cancer_type script - READING SAMPLE SHEET - ",Sys.time()))
	ss <- read.csv(file.path(parent.dir,'sample_sheet_all.csv'),stringsAsFactors=F)
	parent.dir <- file.path(parent.dir,'betas')

	# Create directory path
        betas.by_cancer.dir <- file.path(parent.dir,'by_cancer')

	# Replace sample type column to one of 3 values
	ss <- ss[!ss$sample_type =="no_info",]
	ss$sample_type <- gsub("Additional Metastatic", "Metastatic",ss$sample_type)
	ss$sample_type <- gsub("Recurrent Tumor", "Primary Tumor",ss$sample_type)
	ss$sample_type <- gsub("Additional - New Primary", "Primary Tumor",ss$sample_type)
	ss$sample_type <- gsub(" ", "_", ss$sample_type)

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
				curr.array <- read.csv(file.path(parent.dir,"by_array", paste0(cancer.arrays[k],".csv")),stringsAsFactors=F, row.names=1, check.names=F)
		
				curr.sample.type <- ss$Array.Data.File[array.nums==cancer.arrays[k] & ss$cancer_type==cancer_name & ss$sample_type==sample.type[j]]
				sample.type.filt <- intersect(curr.sample.type,colnames(curr.array))
	
				if (length(sample.type.filt) > 0) {
					if (k==1){
						if (length(sample.type.filt) == 1) {
                                                        data <- as.data.frame(curr.array[,sample.type.filt])
                                                        colnames(data) <- sample.type.filt
							rownames(data) <- rownames(curr.array)
                                                }else{
	                                                data <- curr.array[,sample.type.filt]
						}
					}else {
						if (length(sample.type.filt) == 1) {
							temp.data <- as.data.frame(curr.array[,sample.type.filt])
							colnames(temp.data) <- sample.type.filt 
							data <- cbind(data, temp.data)	
						}else{
							data <- cbind(data, curr.array[,sample.type.filt])
						}
					}	
				}
	
			}	
			
			print(paste0("group_by_cancer_type script - CREATING FILE FOR ",sample.type[j]," - ",Sys.time()))
			sample.type.file.dir <- file.path(cancer.dir, paste0(cancer_name,"_",sample.type[j],".csv"))
			write.csv(data,sample.type.file.dir,quote=FALSE)
		}
	}
	print(paste0("group_by_cancer_type script - FINISHED - ",Sys.time()))
}

args<-commandArgs(trailingOnly=TRUE)
group_by_cancer_type(args[1], args[2])



