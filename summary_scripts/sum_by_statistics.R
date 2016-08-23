sum_by_statistics <- function(sum_dir, statistics_dir) {

	# Get all tissue beta files
	print(paste0("sum_by_statistics script - LISTING SUMMARY FILES - ",Sys.time()))
	sum_files <- list.files(sum_dir, pattern="*.csv")
	
	if (length(sum_files) > 0) {
		# Flag that indicates whether match is needed. wil be TRUE after the first time that rownum is changed
		is_match_needed <- FALSE

        	for (i in 1:length(sum_files)){
		        # Read summary file
		        print(paste0("sum_by_statistics script - READING SUMMARY FILE - ", sum_files[i], " - ", Sys.time()))
			sums <- read.csv(file.path(sum_dir,sum_files[i]),stringsAsFactors=F, row.names=1, check.names=F)	
			tissue_name <- substr(sum_files[i],1,nchar( sum_files[i])-nchar(".csv"))

			if (i==1){
				print(paste0("sum_by_statistics script - CALC FIRST MEANS - ",Sys.time()))
				means <- matrix(sums$means)
                                colnames(means) <- tissue_name
				rownames(means) <- rownames(sums)
				
				print(paste0("sum_by_statistics script - CALC FIRST STD - ",Sys.time()))
				std <- matrix(sums$std)
                                colnames(std) <- tissue_name
                                rownames(std) <- rownames(sums)

				print(paste0("sum_by_statistics script - CALC FIRST QUANTILE 0.1 - ",Sys.time()))
				quantile_0.1 <- matrix(sums$quantile_0.1)
                                colnames(quantile_0.1) <- tissue_name
                                rownames(quantile_0.1) <- rownames(sums)
			
				print(paste0("sum_by_statistics script - CALC FIRST QUANTILE 0.9 - ",Sys.time()))
				quantile_0.9 <- matrix(sums$quantile_0.9)
                                colnames(quantile_0.9) <- tissue_name
                                rownames(quantile_0.9) <- rownames(sums)

				print(paste0("sum_by_statistics script - CALC FIRST SAMPLES_NUM - ",Sys.time()))
				samples_num <- matrix(sums$samples_num)
                                colnames(samples_num) <- tissue_name
                                rownames(samples_num) <- rownames(sums)
			} else {
				print(paste0("sum_by_statistics script - CALC MEANS - ",Sys.time()))
				temp_means <- matrix(sums$means)
                                colnames(temp_means) <- tissue_name
                                rownames(temp_means) <- rownames(sums)

				print(paste0("sum_by_statistics script - CALC STD - ",Sys.time()))
				temp_std <- matrix(sums$std)
                                colnames(temp_std) <- tissue_name
				rownames(temp_std) <- rownames(sums)
				
				print(paste0("sum_by_statistics script - CALC QUANTILE 0.1 - ",Sys.time()))
				temp_quantile_0.1 <- matrix(sums$quantile_0.1)
                                colnames(temp_quantile_0.1) <- tissue_name
				rownames(temp_quantile_0.1) <- rownames(sums)

				print(paste0("sum_by_statistics script - CALC QUANTILE 0.9 - ",Sys.time()))
				temp_quantile_0.9 <- matrix(sums$quantile_0.9)
                                colnames(temp_quantile_0.9) <- tissue_name
                                rownames(temp_quantile_0.9) <- rownames(sums)

				print(paste0("sum_by_statistics script - CALC SAMPLES_NUM - ",Sys.time()))
				temp_samples_num <- matrix(sums$samples_num)
                                colnames(temp_samples_num) <- tissue_name
				rownames(temp_samples_num) <- rownames(sums)

				if ((nrow(sums) != nrow(means)) && (is_match_needed == FALSE)){
					is_match_needed <- TRUE

					# Create final rownames for all data
					full.row.names <- union(rownames(temp_means), rownames(means))

					# Create matrix of final data
					full.means <- matrix(NA, nrow=length(full.row.names), ncol=1)
					rownames(full.means) <- full.row.names
					colnames(full.means) <- tissue_name
					full.std <- matrix(NA, nrow=length(full.row.names), ncol=1)
                                        rownames(full.std) <- full.row.names
					colnames(full.std) <- tissue_name
					full.quantile_0.1 <- matrix(NA, nrow=length(full.row.names), ncol=1)
                                        rownames(full.quantile_0.1) <- full.row.names
					colnames(full.quantile_0.1) <- tissue_name
					full.quantile_0.9 <- matrix(NA, nrow=length(full.row.names), ncol=1)
                                        rownames(full.quantile_0.9) <- full.row.names
					colnames(full.quantile_0.9) <- tissue_name
					full.samples_num <- matrix(NA, nrow=length(full.row.names), ncol=1)
                                        rownames(full.samples_num) <- full.row.names
					colnames(full.samples_num) <- tissue_name

					match_idx <- match(rownames(temp_means), rownames(full.means))
					full.means[match_idx,1] <- temp_means[,1]				
					full.std[match_idx,1] <- temp_std[,1]
					full.quantile_0.1[match_idx,1] <- temp_quantile_0.1[,1]
					full.quantile_0.9[match_idx,1] <- temp_quantile_0.9[,1]
					full.samples_num[match_idx,1] <- temp_samples_num[,1]
					
						

					match_idx <- match(rownames(means), rownames(full.means))
				} else if (is_match_needed == TRUE) {
					
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
