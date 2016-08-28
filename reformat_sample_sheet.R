reformat <- function(ss.old.f,ss.new.f){
	ss.old <- read.table(ss.old.f,header=T,sep='\t',row.names=1,stringsAsFactors=F,check.names=F,skip=8)
	ss.new <- data.frame(Sample = ss.old[,1], Type = ss.old[,2], Chip.id = ss.old[,3], 
		Chip.CO.position = substr(ss.old[,4],4,6), Chip.Row.Pos = substr(ss.old[,4],1,3), 
		Ref = rep(0,nrow(ss.old)),barcode = paste(ss.old[,3],ss.old[,4],sep='_'), 
		Basename = paste(ss.old[,3],ss.old[,4],sep='_'))
	write.csv(ss.new,ss.new.f,row.names=F,quote=F)
}
args<-commandArgs(trailingOnly=TRUE)
reformat(args[1], args[2])
