library("gplots")
BLOOD_METH_THRESHOLD <- 0.8
BLOOD_UNMETH_THRESHOLD <- 0.2
OTHER_TISSUES_METH_THRESHOLD <- 0.8
OTHER_TISSUES_UNMETH_THRESHOLD <- 0.2
TISSUE_UNMETH_THRESHOLD <- 0.6
TISSUE_METH_THRESHOLD <- 0.4
#TISSUE_UNMETH_THRESHOLD <- 0.5
#TISSUE_METH_THRESHOLD <- 0.5
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

write_markers_files <- function(markers, annotation, markers_dir, marker_type, tissue, annotation.id, cancer.data, markers.means){
	# Write means files -------------TODO:FINISH-------
	col.order <- c(blood.cols,setdiff(colnames(unmeth.means),blood.cols))
	col.order <- c(blood.cols,healthy.cols,setdiff(setdiff(colnames(unmeth.means),blood.cols),healthy.cols))	
	means <- merge(annotation,markers.means,by.x=annotation.id,by.y="row.names",sort=F)
	means.file <- file.path(markers_dir, tissue, paste0(tissue,"_", marker_type, "_means.csv"))
	write.csv(means,means.file, row.names=FALSE)

	# Write quantile files
	markers$id <- nrow(markers):1
	markers <- merge(annotation,markers,by.x=annotation.id,by.y="row.names",sort=F)
	markers <- merge(markers, cancer.data, by.x=annotation.id, by.y="row.names",sort=F)
	markers <- markers[order(markers$id,decreasing=T),]
	markers <- markers[,!names(markers) %in% "id"]
	markers.file <- file.path(markers_dir, tissue, paste0(tissue,"_", marker_type,"quantiles.csv"))
	write.csv(markers,markers.file, row.names=FALSE)
	num.markers <- nrow(markers)
	print(paste('Found',num.markers, marker_type,'markers for',tissue))
}

barplot_eb <- function(means,sds,names=1:length(means),ylim=c(0,1),ylab='y_label',main='title'){
  bar.plot <- barplot(means,las=3,names.arg=names,ylim=ylim,ylab=ylab,main=main)
  segments(bar.plot, means - sds, bar.plot, means + sds, lwd = 1.5)
  arrows(bar.plot, means - sds, bar.plot, means + sds, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
}

histogram <- function(num_tissues,breaks,main='title'){
	hist(num_tissues, breaks=breaks, xaxt='n', xlab="num tissues", ylab="markers", main=main)
	axis(side=1, at=breaks)
}

read_matrix <- function(path){
        return(as.matrix(read.csv(path, header=T, row.names=1, nrows=nrows)))
}

convert_tcga_names <- function(tcga_ids){
	tcga.id <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP",
		"LIHC", "LUAD", "LUSC", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
		"THCA", "THYM", "UCEC", "LC", "ACC", "DLBC", "KICH", "LGG", "MESO", "OV", "TGCT", "UCS", "UVM")
	tcga.names <- c('Bladder','Breast','Cervix','Colon','Esophagus','Brain','Head_and_neck','Kidney(KIRC)','Kidney(KIRP)','Liver', 'Lung(LUAD)','Lung(LUSC)','Pancreas','Testes','Prostate','Rectum','Mesenchyme','Skin','Stomach','Thyroid','Thymus','Uterus', 'Larynx', 'Adrenocortical', 'Diffuse_Large_B_cell', 'Kidney_Chromophobe', 'Low_Grade_Glioma', 'Mesothelioma', 'Ovary', 'Testicular_Germ_Cell', 'Uterine_Carcinosarcoma', 'Uveal_Melanoma')
tcga_names <- tcga.names[match(tcga_ids,tcga.id)]
return(tcga_names)
}


TCGA_parent_dir <- "/mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/TCGA"
GEO_parent_dir <- "/mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/GEO"
lab_data_parent_dir <- "/mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/lab_data"

summary_dir <- "summary/by_statistics"
annot.file <- "/mnt/lustre/hms-01/fs01/hadargolan/markers/Illumina_450k_annotation.txt"
markers_dir <-"/mnt/lustre/hms-01/fs01/hadargolan/markers/"

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
colnames(tcga.10.Metastatic) <- paste0(convert_tcga_names(colnames(tcga.10.Metastatic)), "_Metastatic")
colnames(tcga.90.Metastatic) <- paste0(convert_tcga_names(colnames(tcga.90.Metastatic)), "_Metastatic")
colnames(tcga.10.PT) <- paste0(convert_tcga_names(colnames(tcga.10.PT)), "_Primary_Tumor")
colnames(tcga.90.PT) <- paste0(convert_tcga_names(colnames(tcga.90.PT)), "_Primary_Tumor")
colnames(tcga.10.Normal) <- paste0(convert_tcga_names(colnames(tcga.10.Normal)), "_Solid_Tissue_Normal")
colnames(tcga.90.Normal) <- paste0(convert_tcga_names(colnames(tcga.90.Normal)), "_Solid_Tissue_Normal")

# Merge all data
print(paste("Merging data -", Sys.time()))
all.10 <- merge_matrix(geos.10, lab_data.10, tcga.10.Normal, tcga.10.PT, tcga.10.Metastatic)
all.90 <- merge_matrix(geos.90, lab_data.90, tcga.90.Normal, tcga.90.PT, tcga.90.Metastatic)

# Remove irrelevant tissues, including healthy tissue (solid tissue normal) of target tissue
filtered_tissues <- c("Mix","Negative","cfDNA","Cord_blood","Leukocytes_cord_blood","Whole_maternal_blood","Leukocyte","Amniotic_epithelial_cell","Amniotic_membrane","Aortic_smooth_muscle_cell","Astrocyte","Bronchial_epithelial_cell","Cardiac_fibroblast","Cardiac_muscle_cell","Chorionic_membrane","Choroid_plexus_epithelial_cell","Epithelial_cell_of_alveolus_of_lung","Epithelial_cell_of_esophagus","Epithelial_cell_of_prostate","Epithelial_cell_of_proximal_tubule","Fibroblast","Iris_pigment_epithelial_cell","Kidney_epithelial_cell","Mammary_epithelial_cell","Non.pigmented_ciliary_epithelial_cell","Omentum","Placental_chorionic_villi","Placental_decidua","Placental_mesenchyme","Placental_trophoblast","Renal_cortical_epithelial_cell","Skeletal_muscle_cell","Umbilical_cord",)
all.10 <- all.10[,!colnames(all.10) %in% filtered_tissues]
all.90 <- all.90[,!colnames(all.90) %in% filtered_tissues]

# Blood columns definition
geo.blood.cols <- paste0(c("Whole_blood", "PBMC", "Granulocytes", "CD4._T.cells", "CD8._T.cells", "CD19._B.cells", "CD14._Monocytes", "CD56._NK.cells", "Neutrophils", "Eosinophils"),"_GEO")
blood.cols <- c(geo.blood.cols)

# Separating healthy and cancer tissues
tissues.cols <- setdiff(colnames(all.10),blood.cols)
healthy.cols <- tissues.cols[!grepl("_Primary_Tumor", tissues.cols) & !grepl("_Metastatic", tissues.cols)]
cancer.cols <- setdiff(tissues.cols, healthy.cols)

# Separating blood tissues
blood.cols <- intersect(blood.cols, colnames(all.10))
blood.10 <- all.10[,blood.cols]
blood.90 <- all.90[,blood.cols]

# Filter markers from blood
print(paste("Getting blood markers -", Sys.time()))
rows.meth.blood <- get_rows_methylated_in_blood(blood.10,BLOOD_METH_THRESHOLD)
rows.unmeth.blood <- get_rows_unmethylated_in_blood(blood.90,BLOOD_UNMETH_THRESHOLD)
blood.10.meth_blood <- blood.10[rows.meth.blood,]
blood.90.unmeth_blood <- blood.90[rows.unmeth.blood,]

tissues.10 <- all.10[,healthy.cols]
tissues.90 <- all.90[,healthy.cols]
cancer.10 <- all.10[,cancer.cols]
cancer.90 <- all.90[,cancer.cols]

# Filter tissues according to markers from blood
print(paste("Filtering tissues according to blood markers -", Sys.time()))
tissues.10.meth_blood <- tissues.10[rows.meth.blood,]
tissues.90.meth_blood <- tissues.90[rows.meth.blood,]
tissues.10.unmeth_blood <- tissues.10[rows.unmeth.blood,]
tissues.90.unmeth_blood <- tissues.90[rows.unmeth.blood,]

# Create matrix to save amount of markers
markers_count <- matrix(NA, nrow=length(healthy.cols), ncol=2)
rownames(markers_count) <- healthy.cols
colnames(markers_count) <- c("unmeth", "meth")


# Read markers means - to draw barplots
print(paste("Reading markers means -", Sys.time()))
geos.means <- read_matrix(file.path(GEO_parent_dir,summary_dir,"means.csv"))
lab_data.means <- read_matrix(file.path(lab_data_parent_dir,summary_dir,"means.csv"))
tcga.means.PT <- read_matrix(file.path(TCGA_parent_dir,summary_dir, "Primary_Tumor", "means_Primary_Tumor.csv"))
tcga.means.Normal <- read_matrix(file.path(TCGA_parent_dir,summary_dir, "Solid_Tissue_Normal", "means_Solid_Tissue_Normal.csv"))

colnames(tcga.means.PT) <- paste0(convert_tcga_names(colnames(tcga.means.PT)),"_Primary_Tumor")
colnames(tcga.means.Normal) <- paste0(convert_tcga_names(colnames(tcga.means.Normal)),"_Solid_Tissue_Normal")
colnames(geos.means) <- paste0(colnames(geos.means), '_GEO')
colnames(lab_data.means) <- paste0(colnames(lab_data.means), '_LAB')

# Read markers std
print(paste("Reading markers std -", Sys.time()))
geos.std <- read_matrix(file.path(GEO_parent_dir,summary_dir,"std.csv"))
lab_data.std <- read_matrix(file.path(lab_data_parent_dir,summary_dir,"std.csv"))
tcga.std.PT <- read_matrix(file.path(TCGA_parent_dir,summary_dir, "Primary_Tumor", "std_Primary_Tumor.csv"))
tcga.std.Normal <- read_matrix(file.path(TCGA_parent_dir,summary_dir, "Solid_Tissue_Normal", "std_Solid_Tissue_Normal.csv"))

colnames(tcga.std.PT) <- paste0(convert_tcga_names(colnames(tcga.std.PT)),"_Primary_Tumor")
colnames(tcga.std.Normal) <- paste0(convert_tcga_names(colnames(tcga.std.Normal)),"_Solid_Tissue_Normal")
colnames(geos.std) <- paste0(colnames(geos.std), '_GEO')
colnames(lab_data.std) <- paste0(colnames(lab_data.std), '_LAB')

all.means <- merge_matrix(geos.means, lab_data.means, tcga.means.PT, tcga.means.Normal)
all.std <- merge_matrix(geos.std, lab_data.std, tcga.std.PT, tcga.std.Normal)

for (i in 1:length(healthy.cols)) {

	curr_tissue <- healthy.cols[i]
	dir.create(file.path(markers_dir, curr_tissue))
	print(paste("Find markers for", curr_tissue, "-", Sys.time()))

	unmeth.markers <- get_unmeth_markers(curr_tissue,blood.10.meth_blood,tissues.10.meth_blood,tissues.90.meth_blood,TISSUE_UNMETH_THRESHOLD,OTHER_TISSUES_METH_THRESHOLD)
	meth.markers <- get_meth_markers(curr_tissue,blood.90.unmeth_blood,tissues.10.unmeth_blood,tissues.90.unmeth_blood,TISSUE_METH_THRESHOLD,OTHER_TISSUES_UNMETH_THRESHOLD)
	markers_count[i,"unmeth"] <- nrow(unmeth.markers)
	markers_count[i,"meth"] <- nrow(meth.markers)

	unmeth.means <- all.means[rownames(unmeth.markers), , drop=FALSE]
	meth.means <- all.means[rownames(meth.markers), ,drop=FALSE]
	unmeth.std <- all.std[rownames(unmeth.markers), ,drop=FALSE]
	meth.std <- all.std[rownames(meth.markers), ,drop=FALSE]

	unmeth_amount <- ifelse(nrow(unmeth.means) > 10, 10, nrow(unmeth.means))
	meth_amount <- ifelse(nrow(meth.means) > 10, 10, nrow(meth.means))	
	
	# Create unmeth markers and barplots
	if (nrow(unmeth.markers) > 0) {
		write_markers_files(unmeth.markers, annotation, markers_dir, "unmeth", curr_tissue, annotation.id, cancer.10, unmeth.means)
		pdf(file.path(markers_dir, curr_tissue, paste0(curr_tissue,"_unmeth_markers_plots.pdf")), width = 20)
		par(cex.main=2,cex.lab=2)

		# Create histogram according to num_tissues
		breaks <-seq(floor(unmeth.markers$num.tissues[1]/2),unmeth.markers$num.tissues[1], 1)
		histogram(unmeth.markers$num_tissues_threshold, breaks=breaks, main=paste0(curr_tissue, " unmeth markers num tissues"))
		axis(side=1, at=breaks)
		par(mar = c(20,5,2,2), cex.main=1.8,cex.lab=1.8)
		for (j in 1:unmeth_amount) {
			unmeth.curr <- data.frame(unmeth.means[j,], unmeth.std[j,])
			rownames(unmeth.curr) <- colnames(unmeth.means)	
			unmeth.curr <- unmeth.curr[order(unmeth.curr[,1], decreasing=T),]

			barplot_eb(unmeth.curr[,1],unmeth.curr[,2],names=rownames(unmeth.curr), ylab="methylation", main=paste0("unmeth marker for ",curr_tissue,' ', rownames(unmeth.means)[j]))
		}
		dev.off()
	}

	# Create meth barplots
	if (nrow(meth.markers) > 0) {
		write_markers_quantile_file(meth.markers, annotation, markers_dir, "meth", curr_tissue, annotation.id, cancer.90)
		pdf(file.path(markers_dir, curr_tissue,paste0(curr_tissue,"_meth_markers_plots.pdf")), width = 20)
		par(cex.main=2,cex.lab=2)

		# Create histogram according to num_tissues
                breaks <-seq(floor(meth.markers$num.tissues[1]/2),meth.markers$num.tissues[1], 1)
                histogram(meth.markers$num_tissues_threshold, breaks=breaks, main=paste0(curr_tissue, " meth markers num tissues"))
                axis(side=1, at=breaks)
		par(mar = c(20,5,2,2), cex.main=1.8,cex.lab=1.8)
		for (j in 1:meth_amount) {
			meth.curr <- data.frame(meth.means[j,], meth.std[j,])
	                rownames(meth.curr) <- colnames(meth.means)
	                meth.curr <- meth.curr[order(meth.curr[,1], decreasing=T),]
		
			barplot_eb(meth.curr[,1],meth.curr[,2],names=rownames(meth.curr), ylab="methylation", main=paste0("meth marker for ",curr_tissue,' ', rownames(meth.means)[j]))
		}
		dev.off()
	}

		
	# Create heatmap
	print(paste("Create heatmaps ", curr_tissue," ",unmeth_amount," ",meth_amount, "-", Sys.time()))
	create.heatmap <- FALSE
	if ((meth_amount > 0) && (unmeth_amount > 0)){
		markers.heat <- rbind(unmeth.means[1:unmeth_amount,],meth.means[1:meth_amount,])
		create.heatmap <- TRUE
	} else if (meth_amount > 1){
		markers.heat <- meth.means[1:meth_amount,]
		create.heatmap <- TRUE
	} else if (unmeth_amount > 1){
		markers.heat <- unmeth.means[1:unmeth_amount,]
		create.heatmap <- TRUE
	}
	
	if (isTRUE(create.heatmap)){
		pdf(file.path(markers_dir, curr_tissue,paste0(curr_tissue,"_markers_heatmap.pdf")), width = 20)
		heatmap.2(markers.heat, trace="none", keysize=1, density.info="none", margin=c(13, 7), main=paste0(curr_tissue,"_markers"))
		dev.off()

		# Create heatmap - NO CLUSTERING
	        pdf(file.path(markers_dir, curr_tissue,paste0(curr_tissue,"_markers_heatmap_no_cluster.pdf")), width = 20)
        	heatmap.2(markers.heat, trace="none", keysize=1, density.info="none", margin=c(13, 7), Rowv=FALSE, dendrogram="none", Colv="NA", main=paste0(curr_tissue,"_markers"))
	        dev.off()
	}
}

write.csv(markers_count,file.path(markers_dir, "markers_count.csv"),quote=FALSE)

print(paste("FINISHED -", Sys.time()))
