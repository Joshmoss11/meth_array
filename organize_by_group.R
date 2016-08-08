organize_by_group <- function(sample_sheet,parent_dir,sample_id_col,array_id_col,group_col){
	ss<- read.csv(sample_sheet,header=T,stringsAsFactors=F)
	betas.by_array.dir <- file.path(parent_dir,'beta_vals','by_array')
	betas.by_group.dir <- file.path(parent_dir,'beta_vals','by_group')
	dir.create(betas.by_group.dir)
	groups <- ss[,group_col]
	groups.unique <- unique(groups)

}
