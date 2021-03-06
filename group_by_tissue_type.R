group_by_tissue_type <- function(parent.dir, tissue_name)
{
	# Read sample sheet csv
	print(paste0("group_by_tissue_type script - READING SAMPLE SHEET - ",Sys.time()))
	ss <- read.csv(file.path(parent.dir,'sample_sheet.csv'),stringsAsFactors=F)
	parent.dir <- file.path(parent.dir,'betas')
	ss$Type <- gsub(" ", "_", ss$Type)
	ss$Chip.ID[is.na(ss$Chip.ID)] <- "ENCSR"

	# Create path
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
				data <- as.matrix(curr.array[,sample.type.filt])
				 if (length(sample.type.filt) == 1) {
                                        colnames(data) <- sample.type.filt
					rownames(data) <- rownames(curr.array)
                                }
				# In case there is only one array in data
				full.data <- data
			}else {
				# Union row names for arrays with 450 or 850 rows
				full.row.names <- union(rownames(curr.array),rownames(data))
				full.data <- matrix(NA, nrow=length(full.row.names), ncol=ncol(data)+length(sample.type.filt))
		                rownames(full.data) <- full.row.names
                		colnames(full.data) <- c(colnames(data), sample.type.filt)
				
				# Match previous data
				match_idx <- match(rownames(data), rownames(full.data))
				full.data[match_idx,1:ncol(data)] <- data

				# Match current data
				match_idx <- match(rownames(curr.array), rownames(full.data))
                                full.data[match_idx,ncol(data)+1:length(sample.type.filt)] <- as.matrix(curr.array[,sample.type.filt])

				# Copy current data for next loop union
				data <- full.data
			}
		}
	}	
	
	print(paste0("group_by_tissue_type script - CREATING FILE - ",Sys.time()))
	sample.type.file.dir <- file.path(betas.by_tissue.dir, paste0(tissue_name,".csv"))
	write.csv(full.data,sample.type.file.dir,quote=FALSE)
	print(paste0("group_by_tissue_type script - FINAL FILE ROWS, COLS - ", nrow(full.data), ",", ncol(full.data), " - ", Sys.time()))
	print(paste0("group_by_tissue_type script - FINISHED - ",Sys.time()))
}

args<-commandArgs(trailingOnly=TRUE)
group_by_tissue_type(args[1], args[2])
