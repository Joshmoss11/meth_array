sum_by_tissue <- function(betas_file, sum_betas_file) {

	# Calculate statistics per sample_type file
	print(paste0("sum_by_tissue script - READING BETAS - ",Sys.time()))
	betas <- read.csv(betas_file,stringsAsFactors=F, row.names=1, check.names=F)
	
	# Calculate statistics per sample_type file
	print(paste0("sum_by_tissue script - CALC STATISTICS - ",Sys.time()))
	samples_num <- ncol(betas)
	means <- rowMeans(betas, na.rm = TRUE)
	std <- apply(betas, 1, sd, na.rm = TRUE)
	quantile_0.1 <- apply(betas, 1, quantile, probs=c(0.1), na.rm = TRUE)
	quantile_0.9 <- apply(betas, 1, quantile, probs=c(0.9), na.rm = TRUE)

	# Merge and write to file
	print(paste0("sum_by_tissue script - MERGE AND WRITE FILE - ",Sys.time()))
	summary_data <- data.frame(means, std, quantile_0.1, quantile_0.9, samples_num)
	write.csv(summary_data, sum_betas_file, quote=FALSE)
		
	print(paste0("sum_by_tissue script - FINISHED - ",Sys.time()))
}


args<-commandArgs(trailingOnly=TRUE)
sum_by_tissue(args[1], args[2])
