group_by_tissue_type <- function(parent.dir, tissue_name)
{
	# Read sample sheet csv
	print(paste0("group_by_tissue_type script - READING SAMPLE SHEET - ",Sys.time()))
	ss <- read.csv(file.path(parent.dir,'sample_sheet.csv'),stringsAsFactors=F)
	parent.dir <- file.path(parent.dir,'betas')
	ss$Type <- gsub(" ", "_", ss$Type)
	ss$Chip.ID[is.na(ss$Chip.ID)] <- "ENCSR"

	# Create directory
        betas.by_tissue.dir <- file.path(parent.dir,'by_tissue')
	
	# Get all arrays of current tissue type	
	print(paste0("group_by_tissue_type script - LOOP THROUGH ARRAYS - ",Sys.time()))
	tissue.arrays <- unique(ss$Chip.ID[ss$Type==tissue_name])
	for ( k in 1:length(tissue.arrays)){
		# Read current array
		curr.array <- read.csv(file.path(parent.dir,"by_array", paste0(tissue.arrays[k],".csv")),stringsAsFactors=F, row.names=1, check.names=F)
		curr.sample.type <- ss$Basename[ss$Chip.ID==tissue.arrays[k] & ss$Type==tissue_name]
		sample.type.filt <- intersect(curr.sample.type,colnames(curr.array))
	

		if (length(sample.type.filt) > 0) {
			if (k==1){
				 if (length(sample.type.filt) == 1) {
                                        data <- as.data.frame(curr.array[,sample.type.filt])
                                        colnames(data) <- sample.type.filt
					rownames(data) <- rownames(curr.array)
                                } else {
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
	
	print(paste0("group_by_tissue_type script - CREATING FILE - ",Sys.time()))
	sample.type.file.dir <- file.path(betas.by_tissue.dir, paste0(tissue_name,".csv"))
	write.csv(data,sample.type.file.dir,quote=FALSE)
	print(paste0("group_by_tissue_type script - FINISHED - ",Sys.time()))
}

args<-commandArgs(trailingOnly=TRUE)
group_by_tissue_type(args[1], args[2])
