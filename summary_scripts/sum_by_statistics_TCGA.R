sum_by_statistics_TCGA <- function(sample_type, parent_dir, statistics_dir) 
{

	# Get all tissue summary files by sample_type
	print(paste0("sum_by_statistics_TCGA script - LISTING SUMMARY FILES - ",Sys.time()))
	sum_files <- list.files(parent_dir, pattern="*.csv", recursive=TRUE)
	sum_files <- sum_files[grepl(sample_type,sum_files)]

	if (length(sum_files) > 0) {
        	for (i in 1:length(sum_files)){
		        # Read summary file
		        print(paste0("sum_by_statistics_TCGA script - READING SUMMARY FILE - ", sum_files[i], " - ", Sys.time()))
			sums <- read.csv(file.path(parent_dir,sum_files[i]),stringsAsFactors=F, row.names=1, check.names=F)	
			cancer_name <- strsplit(sum_files[i], "/")[[1]][1]
			
			if (i==1){
				print(paste0("sum_by_statistics_TCGA script - CALC FIRST MEANS - ",Sys.time()))
				means <- as.data.frame(sums$means)
                                colnames(means) <- cancer_name
				rownames(means) <- rownames(sums)
				
				print(paste0("sum_by_statistics_TCGA script - CALC FIRST STD - ",Sys.time()))
				std <- as.data.frame(sums$std)
                                colnames(std) <- cancer_name
                                rownames(std) <- rownames(sums)

				print(paste0("sum_by_statistics_TCGA script - CALC FIRST QUANTILE 0.1 - ",Sys.time()))
				quantile_0.1 <- as.data.frame(sums$quantile_0.1)
                                colnames(quantile_0.1) <- cancer_name
                                rownames(quantile_0.1) <- rownames(sums)
			
				print(paste0("sum_by_statistics_TCGA script - CALC FIRST QUANTILE 0.9 - ",Sys.time()))
				quantile_0.9 <- as.data.frame(sums$quantile_0.9)
                                colnames(quantile_0.9) <- cancer_name
                                rownames(quantile_0.9) <- rownames(sums)

				print(paste0("sum_by_statistics_TCGA script - CALC FIRST SAMPLES_NUM - ",Sys.time()))
				samples_num <- as.data.frame(sums$samples_num)
                                colnames(samples_num) <- cancer_name
                                rownames(samples_num) <- rownames(sums)
			} else {
				print(paste0("sum_by_statistics_TCGA script - CALC MEANS - ",Sys.time()))
				temp_means <- as.data.frame(sums$means)
                                colnames(temp_means) <- cancer_name

				print(paste0("sum_by_statistics_TCGA script - CALC STD - ",Sys.time()))
				temp_std <- as.data.frame(sums$std)
                                colnames(temp_std) <- cancer_name
				
				print(paste0("sum_by_statistics_TCGA script - CALC QUANTILE 0.1 - ",Sys.time()))
				temp_quantile_0.1 <- as.data.frame(sums$quantile_0.1)
                                colnames(temp_quantile_0.1) <- cancer_name

				print(paste0("sum_by_statistics_TCGA script - CALC QUANTILE 0.9 - ",Sys.time()))
				temp_quantile_0.9 <- as.data.frame(sums$quantile_0.9)
                                colnames(temp_quantile_0.9) <- cancer_name

				print(paste0("sum_by_statistics_TCGA script - CALC SAMPLES_NUM - ",Sys.time()))
				temp_samples_num <- as.data.frame(sums$samples_num)
                                colnames(temp_samples_num) <- cancer_name

				print(paste0("sum_by_statistics_TCGA script - BIND MATRICES - ",Sys.time()))
				means <- cbind(means, temp_means)
                                std <- cbind(std, temp_std)
                                quantile_0.1 <- cbind(quantile_0.1, temp_quantile_0.1)
                                quantile_0.9 <- cbind(quantile_0.9, temp_quantile_0.9)
                               	samples_num <- cbind(samples_num, temp_samples_num)
			}
		}
	}
	
	# Write files
	print(paste0("sum_by_statistics_TCGA script - WRITE FILES FOR SAMPLE_TYPE - ", sample_type, " - ", Sys.time()))
	write.csv(means, file.path(statistics_dir,paste0("means_",sample_type, ".csv")), quote=FALSE)
	write.csv(std, file.path(statistics_dir,paste0("std_", sample_type, ".csv")), quote=FALSE)
	write.csv(quantile_0.1, file.path(statistics_dir,paste0("quantile_0.1_", sample_type, ".csv")), quote=FALSE)
	write.csv(quantile_0.9, file.path(statistics_dir,paste0("quantile_0.9_", sample_type, ".csv")), quote=FALSE)
	write.csv(samples_num, file.path(statistics_dir,paste0("samples_num_", sample_type, ".csv")), quote=FALSE)
		
	print(paste0("sum_by_statistics_TCGA script - FINISHED - ",Sys.time()))
}


args<-commandArgs(trailingOnly=TRUE)
sum_by_statistics_TCGA(args[1], args[2], args[3])
