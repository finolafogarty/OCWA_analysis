# Introduction ----

# Started by Darren Irwin on 12 November 2017.
# Based on R code written for Greenish Warbler analysis, and then NA warbler analysis.
# The functions are in a file called genomics_R_functions.R
# This script will give general info and use those functions.

##### Orange-crowned Warbler GBS analysis
### Finola Fogarty
### April 2020
# Based on code proveded by Kenny Askelsen and using functions written by Darren Irwin

#Initial Histogram to determine samples to remove

OCWA_FF_March2020_MQ20_het06 <- read.table("/Users/finolafogarty/Desktop/OCWA 2020/Analysis/OCWA_FF_March2020_MQ20_het06.imiss", header = TRUE)

print(OCWA_FF_March2020_MQ20_het06)

hist(OCWA_FF_March2020_MQ20_het06$F_MISS)


# Initial setup ----

#load the following files to computer (command is in Atom):
#"OCWA_FF_March2020_MQ20_het06_missing_80_sample_names_biallelic_GQ10_maf005_miss70_noIndel.012"
#"OCWA_FF_March2020_MQ20_het06_missing_80_sample_names_biallelic_GQ10_maf005_miss70_noIndel.012.indv"
#"OCWA_FF_March2020_MQ20_het06_missing_80_sample_names_biallelic_GQ10_maf005_miss70_noIndel.012.pos”

# set directory where files stored
setwd("~/Desktop/OCWA_2020/Analysis/PCA_April2020")

# Load functions
source("~/Desktop/OCWA_2020/Analysis/PCA_April2020/genomics_R_functions.R")

Analysis_set <- 1  # 1: 137 individuals with Europe;   #2: 131 inds, no Europe;   #3 126 inds, filtered pruned

# PCA whole-genome  ----
# Load 012NA file containing only variable sites throughout genome;
# construct a PCA based on all sites passing an Fst threshold between the "groups" below;
# and all individuals in "groups.to.plot.PCA" according to colors in "group.colors.PCA"

if (Analysis_set == 1) {      # all populations
  groups <- c("Lut", "Ore", "Cel")  # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)
  group.colors.PCA <- c("yellow", "red", "blue")
  groups.to.plot.PCA <- c("Lut", "Ore", "Cel")  
  base.file.name <- "OCWA_FF_March2020_MQ20_het06_missing_80_sample_names_biallelic_GQ10_maf005_miss70_noIndel"
  pos <- read.table(paste0(base.file.name, ".012.pos"))
  column_names <- c("null", paste("c", pos$V1, pos$V2, sep="."))
  geno <- read.table(paste0(base.file.name, ".012"), colClasses = "integer", col.names = column_names)
  SNPnum <- length(geno[1,]) -1   # because the first column is not a SNP (just a count from zero)
  ind <- read.table(paste0(base.file.name, ".012.indv"))
  locations <- read.table("OCWA.Fst_groups.txt", header=TRUE)
  num_loc_cols <- length(locations[1,])
  ind_with_locations <- cbind(ind,locations)
  combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])
  combo[combo == -1] <- NA #Replace -1 with NA
}

# determine number of missing SNPs per bird, and filter out those with more than X% missing SNPs
X <- 100   # this is the percentage threshold
threshold_NA <- SNPnum * X/100
numNAs <- rowSums(is.na(combo))
selection <- (numNAs <= threshold_NA)
combo.NApass.all <- combo[selection,]
ind[which(selection==F),]  # prints out the individuals being left out.

#individuals being left out: 
#factor(0)
#69 Levels: sp_ocwa_plate1__OCWA055 sp_ocwa_plate1__OCWA056 sp_ocwa_plate1__OCWA057 sp_ocwa_plate1__OCWA058 ... sp_ocwa_plate1_OCWA053

# filter out SNPs with more than SNP_missing_max_percent % missing genotypes
SNP_missing_max_percent <- 100
threshold_SNP_missing_max <- length(combo.NApass.all[,1]) * SNP_missing_max_percent/100
numNAs <- colSums(is.na(combo.NApass.all[,(num_loc_cols+1):length(combo.NApass.all[1,])]))
selection <- (numNAs <= threshold_SNP_missing_max)
selection <- c(rep(TRUE, num_loc_cols), selection)
combo.NApass.all <- combo.NApass.all[,selection]

# filter out all but selected chromosome (or set of them):
choose.chrom <- FALSE
if (choose.chrom == TRUE) {
  chrom <- "CM018259.1"
  loci.selection <- (pos$V1 == chrom)
  loci.selection <- c(TRUE, TRUE, TRUE, loci.selection)  # add placeholders for info columns
  combo.NApass <- combo.NApass.all[,loci.selection]
  region.text <- paste0("chr", chrom)	
}	else {
  region.text <- "whole_genome"
  combo.NApass <- combo.NApass.all
}

# Calculate allele freqs and sample sizes (use column Fst_group)
temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
freqs <- temp.list$freqs
sample_size <- temp.list$sample_size
rm(temp.list)

# calculate WC84_Fst for each SNP
temp.list <- getWC84Fst(freqs, sample_size, groups, among=FALSE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
WC84_Fst <- temp.list$WC84_Fst
WC84_Fst_numerator <- temp.list$WC84_Fst_numerator
WC84_Fst_denominator <- temp.list$WC84_Fst_denominator
rm(temp.list)

Fst.filter <- FALSE
Fst.cutoff <- 0.15
# if filtering (Fst.filter == T), choose whether to filter by Fst between pair of populations, or by Fst_among (as defined above)
groups.to.compare <- "Lut_Ore_Cel"
axes <- 3
PCA_results <- plotPCA(Fst.filter, Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass, num_loc_cols, region.text,
                       groups.to.plot.PCA, group.colors.PCA, axes) 
legend("topleft", legend = groups.to.plot.PCA, fill = group.colors.PCA, cex = 0.9)



######################
#Lut_Ore comparison

fst_data<-((as.data.frame(t(WC84_Fst))))
clr <- ifelse(seq(-0.07,1,by=0.01) < 0.07, "red", "grey")[-length(seq(-0.07,1,by=0.01))]
hist(fst_data$Lut_Ore, xlim = c(0,1), breaks = seq(-0.07,1,by=0.01), axes = FALSE, 
     col = clr, main = expression(paste("Distribution of ", italic("F"["ST"]))), xlab =expression(paste(italic("F"["ST"]))))
axis(1, at=seq(0,1,by=0.1))
axis(2, at=seq(0,1000,by=200))
abline(v=0.07,col="black",lty=5)

hist(fst_data$Lut_Ore, xlim = c(0,1), axes = FALSE, 
     col = "grey", main = expression(paste("Distribution of ", italic("F"["ST"]))), xlab =expression(paste(italic("F"["ST"]))))


#Lut_Cel comparison

fst_data<-((as.data.frame(t(WC84_Fst))))
clr <- ifelse(seq(-0.07,1,by=0.01) < 0.07, "red", "grey")[-length(seq(-0.07,1,by=0.01))]
hist(fst_data$Lut_Cel, xlim = c(0,1), breaks = seq(-0.07,1,by=0.01), axes = FALSE, 
     col = clr, main = expression(paste("Distribution of ", italic("F"["ST"]))), xlab =expression(paste(italic("F"["ST"]))))
axis(1, at=seq(0,1,by=0.1))
axis(2, at=seq(0,1000,by=200))
abline(v=0.07,col="black",lty=5)

hist(fst_data$Lut_Cel, xlim = c(0,1), axes = FALSE, 
     col = "grey", main = expression(paste("Distribution of ", italic("F"["ST"]))), xlab =expression(paste(italic("F"["ST"]))))


#Ore_Cel comparison

fst_data<-((as.data.frame(t(WC84_Fst))))
clr <- ifelse(seq(-0.07,1,by=0.01) < 0.07, "red", "grey")[-length(seq(-0.07,1,by=0.01))]
hist(fst_data$Ore_Cel, xlim = c(0,1), breaks = seq(-0.07,1,by=0.01), axes = FALSE, 
     col = clr, main = expression(paste("Distribution of ", italic("F"["ST"]))), xlab =expression(paste(italic("F"["ST"]))))
axis(1, at=seq(0,1,by=0.1))
axis(2, at=seq(0,1000,by=200))
abline(v=0.07,col="black",lty=5)

hist(fst_data$Ore_Cel, xlim = c(0,1), axes = FALSE, 
     col = "grey", main = expression(paste("Distribution of ", italic("F"["ST"]))), xlab =expression(paste(italic("F"["ST"]))))



# PCA_results is a list containing var_explained, scores, and data 
PCA_results$var_explained
View(as.data.frame(PCA_results$scores$PC2))

dat <- as.data.frame(PCA_results$scores$PC2)
dat$PC2 <- PCA_results$scores$PC2
dat$names <- PCA_results$data$ID

write.csv(dat, "dat.csv")

#scores for PC1 in new cvs file (don't know how to combine)
PCA_results$var_explained
View(as.data.frame(PCA_results$scores$PC1))

dat <- as.data.frame(PCA_results$scores$PC1)
dat$PC1 <- PCA_results$scores$PC1
dat$names <- PCA_results$data$ID

write.csv(dat, "datPC1.csv")

#var_explained: [1] 0.03652755 0.06423642 0.08844150

###############
##############
#########
#PCA by SAMPLING LOCATION:

if (Analysis_set == 1) {      # all populations
  groups <- c("Van", "EastBC", "NWBC", "SK", "Rockies", "CentAB", "Kluane", "HG", "Unk")  # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)
  group.colors.PCA <- c("blue", "cyan", "red", "gray", "pink", "yellow", "orange", "green", "brown")
  groups.to.plot.PCA <- c("Van", "EastBC", "NWBC", "SK", "Rockies", "CentAB", "Kluane", "HG", "Unk")  
  base.file.name <- "OCWA_FF_March2020_MQ20_het06_missing_80_sample_names_biallelic_GQ10_maf005_miss70_noIndel"
  pos <- read.table(paste0(base.file.name, ".012.pos"))
  column_names <- c("null", paste("c", pos$V1, pos$V2, sep="."))
  geno <- read.table(paste0(base.file.name, ".012"), colClasses = "integer", col.names = column_names)
  SNPnum <- length(geno[1,]) -1   # because the first column is not a SNP (just a count from zero)
  ind <- read.table(paste0(base.file.name, ".012.indv"))
  locations <- read.table("OCWA.Fst_groups_locations.txt", header=TRUE)
  num_loc_cols <- length(locations[1,])
  ind_with_locations <- cbind(ind,locations)
  combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])
  combo[combo == -1] <- NA #Replace -1 with NA
}

# determine number of missing SNPs per bird, and filter out those with more than X% missing SNPs
X <- 100   # this is the percentage threshold
threshold_NA <- SNPnum * X/100
numNAs <- rowSums(is.na(combo))
selection <- (numNAs <= threshold_NA)
combo.NApass.all <- combo[selection,]
ind[which(selection==F),]  # prints out the individuals being left out.

#individuals being left out: 
#factor(0)
#69 Levels: sp_ocwa_plate1__OCWA055 sp_ocwa_plate1__OCWA056 sp_ocwa_plate1__OCWA057 sp_ocwa_plate1__OCWA058 ... sp_ocwa_plate1_OCWA053

# filter out SNPs with more than SNP_missing_max_percent % missing genotypes
SNP_missing_max_percent <- 100
threshold_SNP_missing_max <- length(combo.NApass.all[,1]) * SNP_missing_max_percent/100
numNAs <- colSums(is.na(combo.NApass.all[,(num_loc_cols+1):length(combo.NApass.all[1,])]))
selection <- (numNAs <= threshold_SNP_missing_max)
selection <- c(rep(TRUE, num_loc_cols), selection)
combo.NApass.all <- combo.NApass.all[,selection]

# filter out all but selected chromosome (or set of them):
choose.chrom <- FALSE
if (choose.chrom == TRUE) {
  chrom <- "CM018259.1"
  loci.selection <- (pos$V1 == chrom)
  loci.selection <- c(TRUE, TRUE, TRUE, loci.selection)  # add placeholders for info columns
  combo.NApass <- combo.NApass.all[,loci.selection]
  region.text <- paste0("chr", chrom)	
}	else {
  region.text <- "whole_genome"
  combo.NApass <- combo.NApass.all
}

# Calculate allele freqs and sample sizes (use column Fst_group)
temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
freqs <- temp.list$freqs
sample_size <- temp.list$sample_size
rm(temp.list)

# calculate WC84_Fst for each SNP
temp.list <- getWC84Fst(freqs, sample_size, groups, among=FALSE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
WC84_Fst <- temp.list$WC84_Fst
WC84_Fst_numerator <- temp.list$WC84_Fst_numerator
WC84_Fst_denominator <- temp.list$WC84_Fst_denominator
rm(temp.list)

Fst.filter <- FALSE
Fst.cutoff <- 0.15
# if filtering (Fst.filter == T), choose whether to filter by Fst between pair of populations, or by Fst_among (as defined above)
groups.to.compare <- "Van_EastBC_NWBC_SK_Rockies_CentAB_Kluane_HG_Unk"
axes <- 3
PCA_results <- plotPCA(Fst.filter, Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass, num_loc_cols, region.text,
                       groups.to.plot.PCA, group.colors.PCA, axes)
legend("topleft", legend = groups.to.plot.PCA, fill = group.colors.PCA, cex = 0.9)

# PCA_results is a list containing var_explained, scores, and data 
PCA_results$var_explained
View(as.data.frame(PCA_results$scores$PC2))

dat <- as.data.frame(PCA_results$scores)
dat$PC2 <- PCA_results$scores$PC2
dat$names <- PCA_results$data$ID

write.csv(dat, "dat.csv")






### Emily
dat <- as.data.frame(PCA_results$scores)
dat$names <- PCA_results$data$ID

write.csv(dat, "dat.csv")

PCA_results$var_explained
View(as.data.frame(PCA_results$scores$PC2))


dat$PC2 <- PCA_results$scores$PC2






