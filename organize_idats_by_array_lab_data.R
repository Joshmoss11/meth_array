parent.dir <- file.path('/mnt','lustre','hms-01','fs01','joshua.moss','dor','meth_array_data','lab_data')
idats.dir <- file.path('/cs','icore','joshua.moss','dor','probe_selection','all_idats')
idats.by_array.dir <- file.path(parent.dir,'idats','by_array')
ss.dir <- file.path(parent.dir,'sample_sheets')
ss.by_array.dir <- file.path(ss.dir,'by_array')
dir.create(idats.by_array.dir)
dir.create(ss.dir)
dir.create(ss.by_array.dir)
ss <- read.csv(file.path(parent.dir,'sample_sheet.csv'),stringsAsFactors=F)
arrays <- ss$Chip.ID
arrays.unique <- unique(arrays)
for (i in 1:length(arrays.unique)){
  array.cur.dir <- file.path(idats.by_array.dir,arrays.unique[i])
  dir.create(array.cur.dir)
  system(paste0('cp ',idats.dir,'/*',arrays.unique[i],'* ',array.cur.dir))
  ss.cur <- ss[arrays %in% arrays.unique[i],]
  write.csv(ss.cur,file.path(ss.by_array.dir,paste0(arrays.unique[i],'.csv')),row.names=F,quote=F)
}
