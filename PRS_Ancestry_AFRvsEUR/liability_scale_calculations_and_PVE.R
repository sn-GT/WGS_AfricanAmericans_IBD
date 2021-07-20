# Liability scale calulations and PVE
# On training set with 1242 cases and 105000 controls

setwd("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/Genotypes_cases")

library(data.table)
controls <- fread("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/Genotypes_cases/Genotypes_150kControls_IDS.txt", header = T)
controls$PHENOTYPE <- 0
names(controls)[1] <- "FID"

for(i in 1:5){
	case <- fread(paste0("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/Genotypes_cases/Train",i,"_case.raw"), header = T)
	names(case) <- gsub('_.*(.*)','\\1', names(case))
	case <- case[,-(2:5)]

	train_ids <- read.table(paste0("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/AA_5Sets/Train",i,".txt"))

	full <- rbind(case,controls)
	df <- full[(full$FID %in% train_ids$V1),]
	df$PHENOTYPE[df$PHENOTYPE == 2] <- 1
	df$PHENOTYPE <- as.factor(as.vector(df$PHENOTYPE))
	df <- data.frame(df)

	table(df$PHENOTYPE)

	profile <- read.table(paste0("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/PROFILES/Profile_AA_Set",i,".txt"), header = F)
	profile$OR <- exp(profile$V3)

	disease_prev <- 0.01
	thresh_disease_prev <- qnorm(disease_prev, mean = 0, sd = 1, lower.tail = F)

	final <- data.frame()
	for(id in profile$V1){
		snp <- as.character(id)
		print(snp)

		cont <- df[df$PHENOTYPE == 0,]  # using the observed allele frequency in controls corresponding to odds ratio of effect allele
		allele_freq <- sqrt(table(cont[,snp])/nrow(cont))

		# p is risk allele [coded as 2], q is protective allele coded as 0
		# OR is the odds ratio
		freq_p <- allele_freq[3]

		OR_p <- profile[profile$V1 == snp, "OR"]

		T1 <- (OR_p * freq_p) / (1 - freq_p)

		Pcase <- T1/(1+T1)

		Case_pp <- Pcase^2
		Case_pq <- 2*Pcase*(1 - Pcase)
		Case_qq <- (1 - Pcase)^2

		Control_pp <- freq_p^2
		Control_pq <- 2*freq_p*(1 - freq_p)
		Control_qq <- (1 - freq_p)^2

		Poverall <- (disease_prev * Pcase) + ((1 - disease_prev)*freq_p)

		Overall_pp <-Poverall^2
		Overall_pq <- 2*Poverall*(1 - Poverall)
		Overall_qq <- (1 - Poverall)^2

		penetrance_pp <- (Case_pp * disease_prev)/Overall_pp
		penetrance_pq <-  (Case_pq * disease_prev)/Overall_pq
		penetrance_qq <-  (Case_qq * disease_prev)/Overall_qq

		displacement_pp <- thresh_disease_prev - qnorm(penetrance_pp, lower.tail = F)
		displacement_pq <- thresh_disease_prev - qnorm(penetrance_pq, lower.tail = F)
		displacement_qq <- thresh_disease_prev - qnorm(penetrance_qq, lower.tail = F)

		mean <- (displacement_pp*Overall_pp) + (displacement_pq*Overall_pq) + (displacement_qq*Overall_qq)

		central_displacement_pp <- displacement_pp - mean
		central_displacement_pq <- displacement_pq - mean
		central_displacement_qq <- displacement_qq - mean

		alpha1 <- (Poverall*central_displacement_pp) + ((1 - Poverall)*central_displacement_pq)
		alpha2 <- (Poverall*central_displacement_pq) + ((1 - Poverall)*central_displacement_qq)

		# PVE
		total_variance <- 2*((alpha1^2 * Poverall) + (alpha2^2*(1 - Poverall)))
		 
		dx <- data.frame("SNP" = snp, "OR_p_AA" = OR_p, "P_cont" = freq_p, "P_case" =  Pcase, "Poverall" = Poverall, "alpha1" = alpha1, "alpha2" = alpha2, "total_variance" = total_variance)
		final <- rbind(final,dx)

	}

	write.table(final, paste0("/gpfs/scratch1/3/snagpal3/WGS_AA/Analysis/UPSAMPLING/liability_PVE_AA_Set",i,".txt"), row.names = F, col.names = T, quote = F)
}