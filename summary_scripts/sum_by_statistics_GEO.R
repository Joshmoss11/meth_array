sum_by_statistics <- function(sum_dir, statistics_dir) {

	# Get all tissue beta files
	print(paste0("sum_by_statistics script - LISTING SUMMARY FILES - ",Sys.time()))
	sum_files <- list.files(sum_dir, pattern="*.csv")

	if (length(sum_files) > 0) {
        	for (i in 1:length(sum_files)){
		        # Read summary file
		        print(paste0("sum_by_statistics script - READING SUMMARY FILE - ", sum_files[i], " - ", Sys.time()))
			sums <- read.csv(file.path(sum_dir,sum_files[i]),stringsAsFactors=F, row.names=1, check.names=F)	
			tissue_name <- substr(sum_files[i],1,nchar( sum_files[i])-nchar(".csv"))

			if (i==1){
				print(paste0("sum_by_statistics script - CALC FIRST MEANS - ",Sys.time()))
				means <- as.data.frame(sums$means)
                                colnames(means) <- tissue_name
				rownames(means) <- rownames(sums)
				
				print(paste0("sum_by_statistics script - CALC FIRST STD - ",Sys.time()))
				std <- as.data.frame(sums$std)
                                colnames(std) <- tissue_name
                                rownames(std) <- rownames(sums)

				print(paste0("sum_by_statistics script - CALC FIRST QUANTILE 0.1 - ",Sys.time()))
				quantile_0.1 <- as.data.frame(sums$quantile_0.1)
                                colnames(quantile_0.1) <- tissue_name
                                rownames(quantile_0.1) <- rownames(sums)
			
				print(paste0("sum_by_statistics script - CALC FIRST QUANTILE 0.9 - ",Sys.time()))
				quantile_0.9 <- as.data.frame(sums$quantile_0.9)
                                colnames(quantile_0.9) <- tissue_name
                                rownames(quantile_0.9) <- rownames(sums)

				print(paste0("sum_by_statistics script - CALC FIRST SAMPLES_NUM - ",Sys.time()))
				samples_num <- as.data.frame(sums$samples_num)
                                colnames(samples_num) <- tissue_name
                                rownames(samples_num) <- rownames(sums)
			} else {
				print(paste0("sum_by_statistics script - CALC MEANS - ",Sys.time()))
				temp_means <- as.data.frame(sums$means)
                                colnames(temp_means) <- tissue_name
                                rownames(temp_means) <- rownames(sums)

				print(paste0("sum_by_statistics script - CALC STD - ",Sys.time()))
				temp_std <- as.data.frame(sums$std)
                                colnames(temp_std) <- tissue_name
				rownames(temp_std) <- rownames(sums)
				
				print(paste0("sum_by_statistics script - CALC QUANTILE 0.1 - ",Sys.time()))
				temp_quantile_0.1 <- as.data.frame(sums$quantile_0.1)
                                colnames(temp_quantile_0.1) <- tissue_name
				rownames(temp_quantile_0.1) <- rownames(sums)

				print(paste0("sum_by_statistics script - CALC QUANTILE 0.9 - ",Sys.time()))
				temp_quantile_0.9 <- as.data.frame(sums$quantile_0.9)
                                colnames(temp_quantile_0.9) <- tissue_name
                                rownames(temp_quantile_0.9) <- rownames(sums)

				print(paste0("sum_by_statistics script - CALC SAMPLES_NUM - ",Sys.time()))
				temp_samples_num <- as.data.frame(sums$samples_num)
                                colnames(temp_samples_num) <- tissue_name
				rownames(temp_samples_num) <- rownames(sums)

				if (nrow(sums) != nrow(means)){
					print(paste0("sum_by_statistics script - MERGE MATRICES - ",Sys.time()))
                                        means <- merge(temp_means, means, by=0, all=TRUE)
                                        std <- merge(temp_std, std, by=0, all=TRUE)
                                        quantile_0.1 <- merge(temp_quantile_0.1, quantile_0.1, by=0, all=TRUE)
                                        quantile_0.9 <- merge(temp_quantile_0.9, quantile_0.9, by=0, all=TRUE)
                                        samples_num <- merge(temp_samples_num, samples_num, by=0, all=TRUE)
				} else {			
					print(paste0("sum_by_statistics script - BIND MATRICES - ",Sys.time()))
					means <- cbind(means, temp_means)
        	                        std <- cbind(std, temp_std)
        	                        quantile_0.1 <- cbind(quantile_0.1, temp_quantile_0.1)
        	                        quantile_0.9 <- cbind(quantile_0.9, temp_quantile_0.9)
                                	samples_num <- cbind(samples_num, temp_samples_num)
				}
			}
		}
	}
	
	# Write files
	print(paste0("sum_by_statistics script - WRITE FILE - ",Sys.time()))
	write.csv(means, file.path(statistics_dir,paste0("means.csv")), quote=FALSE)
	write.csv(std, file.path(statistics_dir,paste0("std.csv")), quote=FALSE)
	write.csv(quantile_0.1, file.path(statistics_dir,paste0("quantile_0.1.csv")), quote=FALSE)
	write.csv(quantile_0.9, file.path(statistics_dir,paste0("quantile_0.9.csv")), quote=FALSE)
	write.csv(samples_num, file.path(statistics_dir,paste0("samples_num.csv")), quote=FALSE)
		
	print(paste0("sum_by_statistics script - FINISHED - ",Sys.time()))
}


args<-commandArgs(trailingOnly=TRUE)
sum_by_statistics(args[1], args[2])
