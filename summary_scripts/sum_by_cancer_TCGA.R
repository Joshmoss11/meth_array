sum_by_cancer_TCGA <- function(betas_folder, sum_betas_output_folder) {

# Get all sample_type files in current cancer folder
print(paste0("sum_by_cancer_TCGA script - LISTING FILES IN CURRENT CANCER FOLDER - ",Sys.time()))
beta_files <- list.files(betas_folder, pattern="*.csv")

if(length(beta_files) > 0) {
	for (i in 1:length(beta_files)){
		# Read beta file
		print(paste0("sum_by_cancer_TCGA script - READING BETAS - ",beta_files[i], "-",Sys.time()))
		betas <- read.csv(file.path(betas_folder,beta_files[i]),stringsAsFactors=F, row.names=1, check.names=F)
	
		# Calculate statistics per sample_type file
		print(paste0("sum_by_cancer_TCGA script - CALC SAMPLES NUM - ",Sys.time()))
		samples_num <- ncol(betas)
		print(paste0("sum_by_cancer_TCGA script - CALC MEANS - ",Sys.time()))
		means <- rowMeans(betas, na.rm = TRUE)
		print(paste0("sum_by_cancer_TCGA script - CALC STD - ",Sys.time()))
		std <- apply(betas, 1, sd, na.rm = TRUE)
		print(paste0("sum_by_cancer_TCGA script - CALC QUANTILE 0.1 - ",Sys.time()))
		quantile_0.1 <- apply(betas, 1, quantile, probs=c(0.1), na.rm = TRUE)
		print(paste0("sum_by_cancer_TCGA script - CALC QUANTILE 0.9 - ",Sys.time()))
		quantile_0.9 <- apply(betas, 1, quantile, probs=c(0.9), na.rm = TRUE)

		# Merge and write to file
		print(paste0("sum_by_cancer_TCGA script - MERGE AND WRITE FILE - ",Sys.time()))
		summary_data <- data.frame(means, std, quantile_0.1, quantile_0.9, samples_num)
		summary_file_name <- file.path(sum_betas_output_folder, beta_files[i])
		write.csv(summary_data, summary_file_name, quote=FALSE)
	}
  }

print(paste0("sum_by_cancer_TCGA script - FINISHED - ",Sys.time()))
}


args<-commandArgs(trailingOnly=TRUE)
sum_by_cancer_TCGA(args[1], args[2])
