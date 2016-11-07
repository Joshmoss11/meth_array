BLOOD_METH_THRESHOLD <- 0.8
BLOOD_UNMETH_THRESHOLD <- 0.2
OTHER_TISSUES_METH_THRESHOLD <- 0.8
OTHER_TISSUES_UNMETH_THRESHOLD <- 0.2
TISSUE_UNMETH_THRESHOLD <- 0.5
TISSUE_METH_THRESHOLD <- 0.5
NUM_ALLOWED_BLOOD <- 0
nrows <- -1

merge_matrix <- function(...){
        input_objects <- list(...)

        # Create list of rownames and colnames for merged data
        rownames.full <- rownames(input_objects[[1]])
        colnames.full <- colnames(input_objects[[1]])
        for (i in 2:length(input_objects)){
                rownames.full <- union(rownames.full, rownames(input_objects[[i]]))
                colnames.full <- c(colnames.full, colnames(input_objects[[i]]))
        }

        # Create final object
        merged_data <- matrix(NA, nrow=length(rownames.full), ncol=length(colnames.full))
        rownames(merged_data) <- rownames.full
        colnames(merged_data) <- colnames.full

        # Merge each object to final object
        start.border <- 1
        end.border <- 0
        for (i in 1:length(input_objects)){
                curr_object <- input_objects[[i]]
                start.border <- 1 + end.border
                end.border <- ncol(curr_object) + end.border

                match_idx <- match(rownames(curr_object), rownames(merged_data))
                merged_data[match_idx, start.border:end.border] <- curr_object
        }

        return(merged_data)
}


get_rows_methylated_in_blood <- function(blood.10,BLOOD_METH_THRESHOLD){
	num.blood <- ncol(blood.10)
	blood.methylated <- blood.10 > BLOOD_METH_THRESHOLD
	blood.methylated.count <- apply(blood.methylated,1,sum)
	rows.methylated.in.blood <- which(blood.methylated.count >= (num.blood-NUM_ALLOWED_BLOOD))
	return(rows.methylated.in.blood)
}

get_rows_unmethylated_in_blood <- function(blood.90,BLOOD_UNMETH_THRESHOLD){
	num.blood <- ncol(blood.90)
	blood.unmethylated <- blood.90 < BLOOD_UNMETH_THRESHOLD
	blood.unmethylated.count <- apply(blood.unmethylated,1,sum)
	rows.unmethylated.in.blood <- which(blood.unmethylated.count >= (num.blood-NUM_ALLOWED_BLOOD))
	return(rows.unmethylated.in.blood)
}

get_unmeth_markers <- function(tissue,blood.10.meth_blood,tissues.10.meth_blood,tissues.90.meth_blood,TISSUE_UNMETH_THRESHOLD,OTHER_TISSUES_METH_THRESHOLD){
	other.tissues <- tissues.10.meth_blood[,setdiff(colnames(tissues.10.meth_blood),tissue)]
	num.tissues <- ncol(other.tissues)
	target.tissue <- tissues.90.meth_blood[,tissue,drop=FALSE]
	target.tissue.unmeth <- target.tissue < TISSUE_UNMETH_THRESHOLD
	target.tissue.unmeth[is.na(target.tissue.unmeth)]<-F
	other.tissues.meth <- tissues.10.meth_blood > OTHER_TISSUES_METH_THRESHOLD
	other.tissues.meth.count <- data.frame(apply(other.tissues.meth,1,sum,na.rm=T,drop=FALSE))
	other.tissues.meth.many <- other.tissues.meth.count > (num.tissues/2)

	markers <- data.frame(target.tissue,blood.10.meth_blood,tissues.10.meth_blood,other.tissues.meth.count, num.tissues)
	markers <- markers[target.tissue.unmeth & other.tissues.meth.many,]
	names(markers)[ncol(markers)-1] <- 'num_tissues_threshold'
	markers <- markers[order(markers[,ncol(markers)-1],decreasing=T),]
	return(markers)
}

get_meth_markers <- function(tissue,blood.table.90.unmeth_blood,tissues.table.10.unmeth_blood,tissues.table.90.unmeth_blood,TISSUE_METH_THRESHOLD,OTHER_TISSUES_UNMETH_THRESHOLD){
	other.tissues <- tissues.90.unmeth_blood[,setdiff(colnames(tissues.90.unmeth_blood),tissue)]
	num.tissues <- ncol(other.tissues)
	target.tissue <- tissues.10.unmeth_blood[,tissue,drop=FALSE]
	target.tissue.meth <- target.tissue > TISSUE_METH_THRESHOLD
	target.tissue.meth[is.na(target.tissue.meth)]<-F
	other.tissues.unmeth <- tissues.90.unmeth_blood < OTHER_TISSUES_UNMETH_THRESHOLD
	other.tissues.unmeth.count <- data.frame(apply(other.tissues.unmeth,1,sum,na.rm=T,drop=FALSE))
	other.tissues.unmeth.many <- other.tissues.unmeth.count > (num.tissues/2)

	markers <- data.frame(target.tissue,blood.90.unmeth_blood,tissues.90.unmeth_blood,other.tissues.unmeth.count, num.tissues)
	markers <- markers[target.tissue.meth & other.tissues.unmeth.many,]
	names(markers)[ncol(markers)-1] <- 'num_tissues_threshold'
	markers <- markers[order(markers[,ncol(markers)-1],decreasing=T),]
	return(markers)
}

write_markers_file <- function(markers, annotation, markers_dir, cancer, marker_type, tissue, annotation.id, cancer.data){
	markers <- merge(annotation,markers,by.x=annotation.id,by.y="row.names",sort=F)
	markers <- merge(markers, cancer.data, by.x=annotation.id, by.y="row.names",sort=F)
	markers <- markers[order(markers$num_tissues,decreasing=T),]
	markers.file <- file.path(markers_dir, cancer, paste0(tissue,"_", marker_type,".txt"))
	write.table(markers,markers.file,sep='\t',col.names=NA,quote=F)
	num.markers <- nrow(markers)
	print(paste('Found',num.markers, marker_type,'markers for',tissue))
}

barplot_eb <- function(means,sds,names=1:length(means),ylim=c(0,1),ylab='y_label',main='tile'){
  bar.plot <- barplot(means,las=3,names.arg=names,ylim=ylim,ylab=ylab,main=main)
  segments(bar.plot, means - sds, bar.plot, means + sds, lwd = 1.5)
  arrows(bar.plot, means - sds, bar.plot, means + sds, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
}


TCGA_parent_dir <- "/mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/TCGA"
GEO_parent_dir <- "/mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/GEO"
lab_data_parent_dir <- "/mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/lab_data"

summary_dir <- "summary/by_statistics"
annot.file <- "/mnt/lustre/hms-01/fs01/hadargolan/markers/Illumina_450k_annotation.txt"
markers_dir <-"/mnt/lustre/hms-01/fs01/hadargolan/markers/"
cancer <- "LC"

# Blood columns definition
geo.blood.cols <- c("Whole_blood", "PBMC", "Granulocytes", "CD4._T.cells", "CD8._T.cells", "CD19._B.cells", "CD14._Monocytes", "CD56._NK.cells", "Neutrophils", "Eosinophils", "Whole_maternal_blood")
lab_data.blood.cols <- c("Leukocyte")
blood.cols <- c(geo.blood.cols,lab_data.blood.cols)

# Read annotation file
print(paste("Reading annotation file -", Sys.time()))
annotation <- read.table(annot.file,sep='\t',header=T)
annotation.id <- 'IlmnID'

# Read cpg's data
print(paste("Reading GEO quantiles -", Sys.time()))
geos.10 <- as.matrix(read.csv(file.path(GEO_parent_dir,summary_dir,"quantile_0.1.csv"), header=T, row.names=1, nrows=nrows))
geos.90 <- as.matrix(read.csv(file.path(GEO_parent_dir,summary_dir,"quantile_0.9.csv"), header=T, row.names=1, nrows=nrows))

print(paste("Reading lab_data quantiles -", Sys.time()))
lab_data.10 <-  as.matrix(read.csv(file.path(lab_data_parent_dir,summary_dir,"quantile_0.1.csv"), header=T, row.names=1, nrows=nrows))
lab_data.90 <-  as.matrix(read.csv(file.path(lab_data_parent_dir,summary_dir,"quantile_0.9.csv"), header=T, row.names=1, nrows=nrows))

print(paste("Reading TCGA quantiles -", Sys.time()))
tcga.10.Metastatic <- as.matrix(read.csv(file.path(TCGA_parent_dir,summary_dir, "Metastatic", "quantile_0.1_Metastatic.csv"), header=T, row.names=1, nrows=nrows))
tcga.90.Metastatic <- as.matrix(read.csv(file.path(TCGA_parent_dir,summary_dir, "Metastatic", "quantile_0.9_Metastatic.csv"), header=T, row.names=1, nrows=nrows))
tcga.10.PT <- as.matrix(read.csv(file.path(TCGA_parent_dir,summary_dir, "Primary_Tumor", "quantile_0.1_Primary_Tumor.csv"), header=T, row.names=1, nrows=nrows))
tcga.90.PT <- as.matrix(read.csv(file.path(TCGA_parent_dir,summary_dir, "Primary_Tumor", "quantile_0.9_Primary_Tumor.csv"), header=T, row.names=1, nrows=nrows))
tcga.10.Normal <- as.matrix(read.csv(file.path(TCGA_parent_dir,summary_dir, "Solid_Tissue_Normal", "quantile_0.1_Solid_Tissue_Normal.csv"), header=T, row.names=1, nrows=nrows))
tcga.90.Normal <- as.matrix(read.csv(file.path(TCGA_parent_dir,summary_dir, "Solid_Tissue_Normal", "quantile_0.9_Solid_Tissue_Normal.csv"), header=T, row.names=1, nrows=nrows))

# Update colnames to contain tissue type
colnames(tcga.10.Metastatic) <- paste0(colnames(tcga.10.Metastatic),"_Metastatic")
colnames(tcga.90.Metastatic) <- paste0(colnames(tcga.90.Metastatic),"_Metastatic")
colnames(tcga.10.PT) <- paste0(colnames(tcga.10.PT),"_Primary_Tumor")
colnames(tcga.90.PT) <- paste0(colnames(tcga.90.PT),"_Primary_Tumor")
colnames(tcga.90.Normal) <- paste0(colnames(tcga.90.Normal),"_Solid_Tissue_Normal")
colnames(tcga.10.Normal) <- paste0(colnames(tcga.10.Normal),"_Solid_Tissue_Normal")

# Merge all data
print(paste("Merging data -", Sys.time()))
all.10 <- merge_matrix(geos.10, lab_data.10, tcga.10.Normal, tcga.10.PT, tcga.10.Metastatic)
all.90 <- merge_matrix(geos.90, lab_data.90, tcga.90.Normal, tcga.90.PT, tcga.90.Metastatic)

# Remove irrelevant tissues, including healthy tissue (solid tissue normal) of target tissue
all.10 <- all.10[,!colnames(all.10) %in% c(paste0(cancer,"_Solid_Tissue_Normal"),"Mix","Negative","cfDNA","Cord_blood","Leukocytes_cord_blood")]
all.90 <- all.90[,!colnames(all.90) %in% c(paste0(cancer,"_Solid_Tissue_Normal"),"Mix","Negative","cfDNA","Cord_blood","Leukocytes_cord_blood")]

# Separating blood tissues
blood.cols <- intersect(blood.cols, colnames(all.10))

blood.10 <- all.10[,blood.cols]
blood.90 <- all.90[,blood.cols]

# Separating healthy and cancer tissues
tissues.cols <- setdiff(colnames(all.10),blood.cols)
healthy.cols <- tissues.cols[!(!grepl(paste0(cancer,'_'), tissues.cols) & grepl("_Primary_Tumor", tissues.cols))]
healthy.cols <- healthy.cols[!(!grepl(cancer, healthy.cols) & grepl("_Metastatic", healthy.cols))]
cancer.cols <- setdiff(tissues.cols, healthy.cols)

tissues.10 <- all.10[,healthy.cols]
tissues.90 <- all.90[,healthy.cols]
cancer.10 <- all.10[,cancer.cols]
cancer.90 <- all.90[,cancer.cols]

# Filter markers from blood
print(paste("Getting blood markers -", Sys.time()))
rows.meth.blood <- get_rows_methylated_in_blood(blood.10,BLOOD_METH_THRESHOLD)
rows.unmeth.blood <- get_rows_unmethylated_in_blood(blood.90,BLOOD_UNMETH_THRESHOLD)
blood.10.meth_blood <- blood.10[rows.meth.blood,]
blood.90.unmeth_blood <- blood.90[rows.unmeth.blood,]

# Filter tissues according to markers from blood
print(paste("Filtering tissues according to blood markers -", Sys.time()))
tissues.10.meth_blood <- tissues.10[rows.meth.blood,]
tissues.90.meth_blood <- tissues.90[rows.meth.blood,]
tissues.10.unmeth_blood <- tissues.10[rows.unmeth.blood,]
tissues.90.unmeth_blood <- tissues.90[rows.unmeth.blood,]

# Find unmeth markers for Primary_Tumor
tissue <- paste0(cancer,"_Primary_Tumor")
print(paste("Find unmeth markers for", tissue, "-", Sys.time()))
unmeth.markers <- get_unmeth_markers(tissue,blood.10.meth_blood,tissues.10.meth_blood,tissues.90.meth_blood,TISSUE_UNMETH_THRESHOLD,OTHER_TISSUES_METH_THRESHOLD)
write_markers_file(unmeth.markers, annotation, markers_dir, cancer, "unmeth", tissue, annotation.id, cancer.10)

print(paste("Find meth markers for", tissue, "-", Sys.time()))
meth.markers <- get_meth_markers(tissue,blood.table.90.unmeth_blood,tissues.table.10.unmeth_blood,tissues.table.90.unmeth_blood,TISSUE_METH_THRESHOLD,OTHER_TISSUES_UNMETH_THRESHOLD)
write_markers_file(meth.markers, annotation, markers_dir, cancer, "meth", tissue, annotation.id, cancer.90)

if (FALSE){
read_matrix <- function(path){
	return(as.matrix(read.csv(path, header=T, row.names=1, nrows=nrows)))
}


# Read markers means
print(paste("Reading markers means -", Sys.time()))
geos.means <- read_matrix(file.path(GEO_parent_dir,summary_dir,"means.csv"))
lab_data.means <- read_matrix(file.path(lab_data_parent_dir,summary_dir,"means.csv"))
tcga.means.PT <- read_matrix(file.path(TCGA_parent_dir,summary_dir, "Primary_Tumor", "means_Primary_Tumor.csv"))
tcga.means.Normal <- read_matrix(file.path(TCGA_parent_dir,summary_dir, "Solid_Tissue_Normal", "means_Solid_Tissue_Normal.csv"))

colnames(tcga.means.PT) <- paste0(colnames(tcga.means.PT),"_Primary_Tumor")
colnames(tcga.means.Normal) <- paste0(colnames(tcga.means.Normal),"_Solid_Tissue_Normal")

all.means <- merge_matrix(geos.means, lab_data.means, tcga.means.PT, tcga.means.Normal)
unmeth.means <- all.means[rownames(unmeth.markers), healthy.cols]
meth.means <- all.means[rownames(meth.markers), healthy.cols]

# Read markers std
print(paste("Reading markers std -", Sys.time()))
geos.std <- read_matrix(file.path(GEO_parent_dir,summary_dir,"std.csv"))
lab_data.std <- read_matrix(file.path(lab_data_parent_dir,summary_dir,"std.csv"))
tcga.std.PT <- read_matrix(file.path(TCGA_parent_dir,summary_dir, "Primary_Tumor", "std_Primary_Tumor.csv"))
tcga.std.Normal <- read_matrix(file.path(TCGA_parent_dir,summary_dir, "Solid_Tissue_Normal", "std_Solid_Tissue_Normal.csv"))

colnames(tcga.std.PT) <- paste0(colnames(tcga.std.PT),"_Primary_Tumor")
colnames(tcga.std.Normal) <- paste0(colnames(tcga.std.Normal),"_Solid_Tissue_Normal")

all.std <- merge_matrix(geos.std, lab_data.std, tcga.std.PT, tcga.std.Normal)
unmeth.std <- all.std[rownames(unmeth.markers), healthy.cols]
meth.std <- all.std[rownames(meth.markers), healthy.cols]

pdf('unmeth_marker_barplot.pdf')
barplot_eb(unmeth.means,unmeth.std,names=colnames(unmeth.means))
dev.off()
}
print(paste("FINISHED -", Sys.time()))
