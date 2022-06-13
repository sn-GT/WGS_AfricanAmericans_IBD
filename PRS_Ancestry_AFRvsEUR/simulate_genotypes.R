
# Data with SNP,  allele, allele frequency for control individuals
mydata = read.csv("allele_freq_controls.csv")
t1 <- (1-mydata$Freq)^2
t2 <- 1-mydata$Freq^2
set.seed(7251963)

# Number of individuals to simulate genotypes for
num_individuals <- 150000
genolist = list()
for(i in 1:num_individuals) {
  n1000 <- runif(1000,0,1)
  geno <- ifelse(n1000<t1, 0, ifelse(n1000<t2,1,2))
  genolist[[i]] <- geno
}

# dataframe with simulated genotypes for num_individuals and SNPs provided in input data
sim_geno_df <- data.frame(genolist)
colnames(sim_geno_df) <- c(paste0("ID", c(1:num_individuals)))
sim_geno_df$SNP <- mydata$SNP

write.table(sim_geno_df, "Simulatedgenotypes.txt", row.names = F, col.names = T, quote = F, sep = "\t")

#############################################################################################

# Check simulated genotypes
# Check if allele frequency of simulated genotypes match with 
# allele frequency provided in input real data

sim_geno_df$AF <- (rowSums(sim_geno_df))/(2*num_individuals)
sim_geno_df$SNP <- rownames(sim_geno_df)

match_sim_actual <- merge(mydata, sim_geno_df, by = "SNP")

# Compare allele freqeuncies of simulated and real input data
plot(match_sim_actual$Freq, match_sim_actual$AF, xlab = "Allele freq simulated data", ylab = "Allele freq real data")
cor.test(match_sim_actual$Freq, match_sim_actual$AF)

#############################################################################################
