library(minfi)
library(wateRmelon)

process_array <- function(parent_dir,array_id)
{
	array_path <- file.path(parent_dir, "idats/by_array",array_id)
	out_path <- file.path(parent_dir, "betas/by_array",paste0(array_id,".csv"))
	
	print(paste0("1. process_array script - READING ARRAY - ",Sys.time()))
	rgSet <- read.metharray.exp(array_path,extended=T)
	
	# Get p-values	
	print(paste0("2. process_array script - P_VALUES CALCULATION - ",Sys.time()))
	pVal <- detectionP(rgSet, type = "m+u")

	# ssNoob
	print(paste0("3. process_array script - PREPROCESSNOOB - ",Sys.time()))
	MethSet <- preprocessNoob(rgSet,dyeMethod="single")
	
	# Context removal
	print(paste0("4. process_array script - CONTEXT_REMOVAL - ",Sys.time()))
	MethDrop <- dropMethylationLoci(MethSet, dropRS = TRUE, dropCH = TRUE)

	# Get RatioSet
	print(paste0("5. process_array script - CONVERT TO RATIOSET - ",Sys.time()))
	RSet <- ratioConvert(MethDrop, what = "beta", keepCN = FALSE)

	# Get GenomicRatioSet
	print(paste0("6. process_array script - CONVERT TO GENOMICRATIOSET - ",Sys.time()))
        grSet <- mapToGenome(RSet)
	
	# Remove SNPs
	print(paste0("7. process_array script - REMOVE SNPs - ",Sys.time()))
        grSet <- dropLociWithSnps(grSet, snps=c("SBE","CpG"), maf=0)
	
	# Get beta values
	print(paste0("8. process_array script - BETA VALUES CALCULATION - ",Sys.time()))
	betas<-getBeta(grSet)

	# Remove sex chromosomes
	print(paste0("9. process_array script - REMOVE SEX CHROMOSOMES - ",Sys.time()))
	chrom <- as.character(seqnames(granges(grSet)))
	chrom.sex <- chrom=='chrX' | chrom=='chrY'
	betas.aut <- betas[!chrom.sex,]
	
	# Convert probe types from 'I' and 'II' to 1,2
	print(paste0("10. process_array script - GETTING PROBE TYPES - ",Sys.time()))
	ProbeType <- getProbeType(grSet)
	ProbeType <- ProbeType[!chrom.sex]
    	ProbeTypeNum <- integer(length(ProbeType))
	ProbeTypeNum[ProbeType=='I'] <- 1
    	ProbeTypeNum[ProbeType=='II'] <- 2

	# Filter p-values
	print(paste0("12. process_array script - FIX ACCORDING TO P_VALUES - ",Sys.time()))
	pVal.Filt <- pVal[match(rownames(betas.aut),rownames(pVal))]
	betas.aut[pVal.Filt>0.05] <- NA
	
	print(paste0("13. process_array script - WRITING ARRAY TO CSV - ",Sys.time()))
	if (length(sampleNames(rgSet)) == 1){
		betas.aut <- as.data.frame(betas.aut)
		colnames(betas.aut)<-sampleNames(rgSet)
	}
	
	write.csv(betas.aut,out_path,quote=FALSE)

	print(paste0("14. process_array script - FINISHED - ",Sys.time()))	
}

args<-commandArgs(trailingOnly=TRUE)
process_array(args[1], args[2])

