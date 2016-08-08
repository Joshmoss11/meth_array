parent.dir <- file.path('/cs','icore','joshua.moss','dor','atlas','TCGA')
out.dir <- file.path('/mnt','lustre','hms-01','fs01','joshua.moss','dor','meth_array_data','TCGA')

tcga.dirs <- c('ACC','BLCA','BRCA','CESC','COAD','DLBC','ESCA','GBM','HNSC',
               'KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO',
               'OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT',
               'THCA','THYM','UCEC','UCS','UVM')

i <- 1
ss.file <- file.path(parent.dir,tcga.dirs[i],'DNA_Methylation','JHU_USC__HumanMethylation450',paste(tcga.dirs[i],'_sample_annotation.txt',sep=''))
ss <- read.table(ss.file,sep='\t',header=T)
ss$cancer_type <- rep(tcga.dirs[i])

for (i in 2:length(tcga.dirs)){
  ss.file <- file.path(parent.dir,tcga.dirs[i],'DNA_Methylation','JHU_USC__HumanMethylation450',paste(tcga.dirs[i],'_sample_annotation.txt',sep=''))
  ss.temp <- read.table(ss.file,sep='\t',header=T)
  if (!("tumor_tissue_site" %in% colnames(ss.temp))){
    ss.temp$tumor_tissue_site <- rep('blood')
  }
  ss.temp$cancer_type <- rep(tcga.dirs[i])
  ss <- rbind(ss,ss.temp)
}

out.file <- file.path(out.dir,'sample_sheet_all.csv')
write.csv(ss,out.file,row.names=F,quote=F)
