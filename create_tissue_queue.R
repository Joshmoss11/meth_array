create_tissue_queue <- function(parent_dir)
{
	# Read sample sheet csv
	parent.dir <- file.path(parent_dir)
	ss <- read.csv(file.path(parent.dir,'sample_sheet.csv'),stringsAsFactors=F)

        tissueT <- unique(ss$Type)
	tissueT <- gsub(" ", "_" ,tissueT)
	fileConn <- file(file.path(parent_dir,"tissue_queue"))
	writeLines(tissueT,fileConn)

}

args<-commandArgs(trailingOnly=TRUE)
create_tissue_queue(args[1])
