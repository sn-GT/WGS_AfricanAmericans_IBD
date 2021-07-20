setwd("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/Genotypes_cases")

library(data.table)
controls <- fread("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/Genotypes_cases/Genotypes_150kControls_IDS.txt", header = T)
controls$PHENOTYPE <- 0
names(controls)[1] <- "FID"

for(i in 1:5){
case <- fread(paste0("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/Genotypes_cases/Test",i,"_case.raw"), header = T)
names(case) <- gsub('_.*(.*)','\\1', names(case))
case <- case[,-(2:5)]

train_ids <- read.table(paste0("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/AA_5Sets/Test",i,".txt"))
train_ids$V1 <- as.character(as.vector(train_ids$V1))

full <- rbind(case,controls)
df <- full[(full$FID %in% train_ids$V1),]
df$PHENOTYPE[df$PHENOTYPE == 2] <- 1

dim(df)
table(df$PHENOTYPE)

df$PHENOTYPE <- as.factor(as.vector(df$PHENOTYPE))
df <- data.frame(df)
temp <- df

profile <- read.table(paste0("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/liability/Dave_Correction/Dave_Corrected_liability_AA_Set",i,".txt"), header = T)
#profile <- read.table(paste0("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/liability/Dave_Correction/Dave_Corrected_liability_UKB_Set",i,".txt"), header = T)


comm <- intersect(colnames(df), profile$SNP)

df <- df[,comm]

for(j in 1:ncol(df)){
	print(j)
	id <- colnames(df)[j]
	alpha1 <- profile[profile$SNP == id, "alpha1"]
	alpha2 <- profile[profile$SNP == id, "alpha2"]
	OR <- profile[profile$SNP == id, "OR_p_AA"]

	df[,id][df[,id] == 2] <- 2*alpha1
	df[,id][df[,id] == 1] <- alpha1+alpha2
	df[,id][df[,id] == 0] <- 2*alpha2
	df[,id][is.na(df[,id])] <- 0
}

df$PRS <- rowSums(df)
f <- cbind(temp, df)

write.table(f[,c("FID", "PHENOTYPE", "PRS")], paste0("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/liability/PRS_AA_Set",i,"_liability_AA.txt"), row.names=F, col.names=T, quote = F)
#write.table(f[,c("FID", "PHENOTYPE", "PRS")], paste0("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/liability/Dave_Correction/PRS/Dave_Corrected_PRS_AA_Set",i,"_liability_UKB.txt"), row.names=F, col.names=T, quote = F)
}