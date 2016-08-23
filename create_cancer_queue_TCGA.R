create_cancer_queue <- function(parent_dir)
{
	# Read sample sheet csv
	parent.dir <- file.path(parent_dir)
	ss <- read.csv(file.path(parent.dir,'sample_sheet_all.csv'),stringsAsFactors=F)

        cancerT <- unique(ss$cancer_type)
	fileConn <- file(file.path(parent_dir,"cancer_queue"))
	writeLines(cancerT,fileConn)

}

args<-commandArgs(trailingOnly=TRUE)
create_cancer_queue(args[1])
