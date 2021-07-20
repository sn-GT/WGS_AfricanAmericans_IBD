setwd("C:\\Users\\snagpal3\\Dropbox (GaTech)\\IBD\\WGS_AA\\Upsampling_AA")
mydata = read.csv("allele_freq_controls.csv")
t1 <- (1-mydata$Freq)^2
t2 <- 1-mydata$Freq^2
set.seed(7251963)

genolist = list()
weightlist = list()
for(i in 1:100000) {
  n1000 <- runif(1000,0,1)
  geno <- ifelse(n1000<t1, 0, ifelse(n1000<t2,1,2))
  wgeno <- ifelse(n1000<t1, 0, ifelse(n1000<t2,mydata$LogOdds,2*mydata$LogOdds))
  genolist[[i]] <- geno
  weightlist[[i]] <- wgeno }
pop_geno = do.call(cbind,datalist)
pop_wgeno = do.call(cbind,weightlist)
write.csv(pop_geno, file="Population Genotypes G0")				#	//OPTIONAL 

prs <- colSums(pop_wgeno)
quantile(prs,c(0.5))
rm(pop_wgeno)
rm(weightlist)
pop_geno <- rbind(pop_geno,prs)
Sub <- subset(pop_geno, , pop_geno[1001,] >= quantile(prs,c(0.2)))			#//SELECTION STEP
Sub <- Sub[-1001,]

new_freq <- rowSums(Sub)/(2*80000)
hist(prs, xlim=c(-3,5))
t1 <- (1-new_freq)^2
t2 <- 1-new_freq^2
prob_dist <- exp(prs)/135.9
plot(prs, prob_dist, xlim=c(-3,5), ylim=c(0,1))

sort_pd <- sort(prob_dist)
rr <- runif(100000,0,1)
case <- ifelse(rr<sort_pd,1,0)
sum(case[1:100000])
dlist = list()
for(i in 1:50) {
  sumcase <- sum(case[(2000*(i-1)+1):(2000*(i-1)+2000)])
  dlist[[i]] <- sumcase }
prevalences = do.call(cbind,dlist)/20
percentile <- seq(1,50,by=1)
plot(percentile,prevalences, ylim=c(0,30))

Exponential model
prob_dist <- exp(prs*0.02)*0.000000004
prob_dist <- exp(prs*0.03)*0.0000000000002			steeper
prob_dist <- 0.1 + ifelse(prs<950, 0, exp(prs*0.03)*0.0000000000004)		higher prevalence

Importing the Frequencies and Weights

mydata = read.csv("CADtest.csv")
t1 <- (1-mydata$Freq)^2
t2 <- 1-mydata$Freq^2
set.seed(7251963)

genolist = list()
weightlist = list()
for(i in 1:100000) {
  n1000 <- runif(1000,0,1)
  geno <- ifelse(n1000<t1, 0, ifelse(n1000<t2,1,2))
  wgeno <- ifelse(n1000<t1, 0, ifelse(n1000<t2,mydata$LogOdds,2*mydata$LogOdds))
  genolist[[i]] <- geno
  weightlist[[i]] <- wgeno }
pop_geno = do.call(cbind,genolist)
pop_wgeno = do.call(cbind,weightlist)
write.csv(pop_geno, file="Population Genotypes G0")

prs <- colSums(pop_wgeno)
quantile(prs,c(0.5))
pop_geno <- rbind(pop_geno,prs)
Sub <- subset(pop_geno, , pop_geno[1001,] <= quantile(prs,c(0.8)))
Sub <- Sub[-1001,]

new_freq <- rowSums(Sub)/(2*50000)
hist(prs, xlim=c(-3,5))
t1 <- (1-new_freq)^2
t2 <- 1-new_freq^2
prob_dist <- exp(prs)/135.9
plot(prs,prob_dist, xlim=c(-3,5), ylim=c(0,1))

sort_pd <- sort(prob_dist)
rr <- runif(100000,0,1)
case <- ifelse(rr<sort_pd,1,0)
sum(case[1:100000])
dlist = list()
for(i in 1:50) {
  sumcase <- sum(case[(2000*(i-1)+1):(2000*(i-1)+2000)])
  dlist[[i]] <- sumcase }
prevalences = do.call(cbind,dlist)/20
percentile <- seq(1,50,by=1)
plot(percentile,prevalences, ylim=c(0,30))

