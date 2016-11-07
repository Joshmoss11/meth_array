sum_by_statistics <- function(sum_dir, statistics_dir) {

	# Get all tissue beta files
	print(paste0("sum_by_statistics script - LISTING SUMMARY FILES - ",Sys.time()))
	sum_files <- list.files(sum_dir, pattern="*.csv")
	
	if (length(sum_files) > 0) {
		i <- 1
		found_450 <- FALSE
		found_850 <- FALSE
		print(paste0("sum_by_statistics script - SEARCHING CG LIST (450+850)", Sys.time()))
		while ((i <= length(sum_files)) && (!found_450 || !found_850)) {
			row.count <- length(count.fields(file.path(sum_dir,sum_files[i])))
			if ((!found_450) && (row.count > 400000) && (row.count < 500000)){
				sums <- read.csv(file.path(sum_dir,sum_files[i]),stringsAsFactors=F, row.names=1, check.names=F)
				cg_450 <- rownames(sums)
				found_450 <- TRUE
			} else if ((!found_850) && (row.count > 800000) && (row.count < 900000)){
				sums <- read.csv(file.path(sum_dir,sum_files[i]),stringsAsFactors=F, row.names=1, check.names=F)
				cg_850 <- rownames(sums)
                                found_850 <- TRUE
			}
			i <- i+1
		}

		print(paste0("sum_by_statistics script - CREATING CG LIST (450+850)", Sys.time()))
		full.row.names <- union(cg_450, cg_850)
		tissue_names <- substr(sum_files,1,nchar( sum_files)-nchar(".csv"))
                full.means <- matrix(NA, nrow=length(full.row.names), ncol=length(sum_files))
                rownames(full.means) <- full.row.names
		colnames(full.means) <- tissue_names
                full.std <- matrix(NA, nrow=length(full.row.names), ncol=length(sum_files))
                rownames(full.std) <- full.row.names
		colnames(full.std) <- tissue_names
                full.quantile_0.1 <- matrix(NA, nrow=length(full.row.names), ncol=length(sum_files))
                rownames(full.quantile_0.1) <- full.row.names
		colnames(full.quantile_0.1) <- tissue_names
                full.quantile_0.9 <- matrix(NA, nrow=length(full.row.names), ncol=length(sum_files))
                rownames(full.quantile_0.9) <- full.row.names
		colnames(full.quantile_0.9) <- tissue_names
                full.samples_num <- matrix(NA, nrow=length(full.row.names), ncol=length(sum_files))
                rownames(full.samples_num) <- full.row.names
		colnames(full.samples_num) <- tissue_names

		for (i in 1:length(sum_files)) {
			print(paste0("sum_by_statistics script - READING TISSUE FILE - ", sum_files[i], " - ", Sys.time()))
			sums <- read.csv(file.path(sum_dir,sum_files[i]),stringsAsFactors=F, row.names=1, check.names=F)
			match_idx <- match(rownames(sums), rownames(full.means))
			full.means[match_idx,i] <- sums[,"means"]
			full.std[match_idx,i] <- sums[,"std"]	
			full.quantile_0.1[match_idx,i] <- sums[,"quantile_0.1"]
			full.quantile_0.9[match_idx,i] <- sums[,"quantile_0.9"]
			full.samples_num[match_idx,i] <- sums[,"samples_num"]
		}
		
		# Write files
        	print(paste0("sum_by_statistics script - WRITE FILES - ",Sys.time()))
	        write.csv(full.means, file.path(statistics_dir,paste0("means.csv")), quote=FALSE)
        	write.csv(full.std, file.path(statistics_dir,paste0("std.csv")), quote=FALSE)
	        write.csv(full.quantile_0.1, file.path(statistics_dir,paste0("quantile_0.1.csv")), quote=FALSE)
	        write.csv(full.quantile_0.9, file.path(statistics_dir,paste0("quantile_0.9.csv")), quote=FALSE)
	        write.csv(full.samples_num, file.path(statistics_dir,paste0("samples_num.csv")), quote=FALSE)	
	}
	
	print(paste0("sum_by_statistics script - FINISHED - ",Sys.time()))
}


args<-commandArgs(trailingOnly=TRUE)
sum_by_statistics(args[1], args[2])
