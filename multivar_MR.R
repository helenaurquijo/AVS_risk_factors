################################################################################
           
              # Multivariable analyses of AVS risk factors #

################################################################################

library(dplyr)
library(stringr)
library(tidyverse)
library(TwoSampleMR)
library(vroom)
library(data.table)
library(readr)
library(tidyr)
library(tibble)
library(MVMR)

setwd("~/wd/")

load("AVS_GWAS.RData") 

u <- 4807/(4807+458120)
v <- u*(1-u)

################################################################################
                  # Analyses with local exposure datasets
################################################################################

############################# 3-lipid model  ###################################
#each lipid gnetically adjusted for the two others 

exp1 <- fread('~/wd/triglycerides_multivar.txt', header=T)
exp2 <- fread('~/wd/LDL_multivar.txt', header=T)
exp3 <- fread('~/wd/apolipoprotein_B_multivar.txt', header=T)


exp1$Phenotype <- 'Triglycerides'
exp2$Phenotype <- 'LDL cholesterol'
exp3$Phenotype <- 'Apolipoprotein B'

exp_dat1 <- format_data(exp1, type="exposure")
exp_dat2 <- format_data(exp2, type="exposure")
exp_dat3 <- format_data(exp3, type="exposure")
exp_dat <- rbind(exp_dat1, exp_dat2, exp_dat3)

mvdat <- mv_harmonise_data(exp_dat, out_dat, harmonise_strictness = 1)

mvres <- mv_multiple(mvdat)

#BOLT LMM correction and convert to OR
mv_res <- mvres$result %>%
  select(-id.exposure, -id.outcome) %>% 
  mutate(outcome = "AVS") %>% 
  mutate(b = b/v,
         se = se/v) %>% 
  mutate(OR = exp(b),
         UCI = exp(b + 1.96 * se),
         LCI = exp(b - 1.96 * se))


write.table(all_lipids, file = "lipids_mvmr.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE)

#Conditional F stats 
colnames(exp1)<-paste(colnames(exp1),"1",sep = "_")
colnames(exp2)<-paste(colnames(exp2),"2",sep = "_")
colnames(exp3)<-paste(colnames(exp3),"3",sep = "_")

out<-out_dat[,c(1,6,7)]

dat1<- inner_join(exp1,exp2, by=c("SNP_1"="SNP_2")) 
dat1<- inner_join(dat1,exp3, by=c("SNP_1"="SNP_3"))
dat1<- inner_join(dat1,out, by=c("SNP_1"="SNP"))

bx <- as.matrix(dat1[,c("beta_1","beta_2","beta_3")])
bxse <- as.matrix(dat1[,c("se_1","se_2","se_3")])
dat <- MendelianRandomization::mr_mvinput(
  bx = bx, bxse = bxse, by = dat1$beta.outcome, byse = dat1$se.outcome,
  snps = dat1$SNP_1
)
mvdat <- mrmvinput_to_mvmr_format(dat)

cov <- fread("~/wd/correlation_matrix.txt")
Pcov <- as.matrix(cov)

phe_cov <- phenocov_mvmr(Pcov, bxse)

sres <- strength_mvmr(r_input = mvdat, gencov = phe_cov)
pres <- pleiotropy_mvmr(r_input = mvdat, gencov = phe_cov)

res <- as.data.frame(ivw_mvmr(r_input = mvdat)) %>% 
  mutate(`Estimate` = `Estimate`/v)

res1 <- qhet_mvmr(mvdat, Pcov, CI = F, iterations = 100) 

################################################################################
# ApoB and Lp(a)

exp1 <- fread('~/wd/apo_B_multivar_lipoprotein_A.txt', header=T)
exp2 <- fread('~/wd/lipoprotein_A_multivar_apo_B.txt', header=T)

exp1$Phenotype <- 'Apolipoprotein B'
exp2$Phenotype <- 'Lipoprotein(a)'

exp_dat1 <- format_data(exp1, type="exposure")
exp_dat2 <- format_data(exp2, type="exposure")

exp_dat <- rbind(exp_dat1, exp_dat2)

mvdat <- mv_harmonise_data(exp_dat, out_dat, harmonise_strictness = 1)

mvres <- mv_multiple(mvdat)

mv_res <- mvres$result %>%
  select(-id.exposure, -id.outcome) %>% 
  mutate(outcome = "AVS") %>% 
  mutate(b = b/v,
         se = se/v) %>% 
  mutate(OR = exp(b),
         UCI = exp(b + 1.96 * se),
         LCI = exp(b - 1.96 * se))


################################################################################
# Multivariable paiwrise analyses according to pearson's correlation of >0.1

phe_corr <- fread("~/wd/correlation_matrix.txt")

get_pairs <- function(phe_cor,thres){
  pairs <- as.data.frame(which((phe_cor<1) & (phe_cor>thres), arr.ind = TRUE))
  key<-as.data.frame(colnames(phe_cor)) 
  key <- key %>% mutate(n = seq.int(nrow(key)))
  paired<-left_join(pairs, key, by=c("row"="n"), keep = FALSE) %>% 
    select(-row) %>% rename(exp_1=`colnames(phe_cor)`)
  paired<-left_join(paired, key, by=c("col"="n"), keep = FALSE) %>% 
    select(-col) %>% rename(exp_2=`colnames(phe_cor)`)
  unique_pairs <- paired[!duplicated(t(apply(paired[c("exp_1", "exp_2")], 1, sort))), ] 
  return(unique_pairs)
}

pairwise <- get_pairs(phe_cor = phe_corr, thres = 0.1)

exps <- list.files('~/wd/pairwise_mvmr/')

exps_path<- "~/wd/pairwise_mvmr/"

pairwise_res <- data.frame()

key<-as.data.frame(colnames(phe_corr)) 

for (i in 1:nrow(pairwise)) {
  exp_dat <- NULL
  mvdat<- NULL
  exp1 <- pairwise$exp_1[i]
  exp2 <- pairwise$exp_2[i]
  exposure1 <- fread(paste(exps_path, exp1,"_multivar_", exp2,".txt", sep=""))
  exposure2 <- fread(paste(exps_path, exp2,"_multivar_", exp1,".txt", sep=""))
  
  
  col_names <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval', 'tmp', 'tmp1')
  col_n<-ncol(exposure1)
  names(exposure1)<-col_names[1:col_n]
  col_n<-ncol(exposure2)
  names(exposure2)<-col_names[1:col_n]
  
  exposure1$Phenotype <- exp1
  exposure2$Phenotype <- exp2
  exp_dat1 <- format_data(exposure1, type='exposure')
  exp_dat2 <- format_data(exposure2, type='exposure')
  
  print(paste(exp1, "and",exp2), sep = " ")
  exp_dat <- rbind(exp_dat1,exp_dat2)
  
  mvdat <- mv_harmonise_data(exp_dat, out_dat)
  
  mvres <- mv_multiple(mvdat)
  
  mv_res <- mvres$result %>%
    select(-id.exposure, -id.outcome) %>% 
    mutate(outcome = "AVS") %>% 
    mutate(b = b/v,
           se = se/v) %>% 
    mutate(OR = exp(b),
           UCI = exp(b + 1.96 * se),
           LCI = exp(b - 1.96 * se)) %>% 
    mutate(`Adjusted for` = 
             case_when(
               str_detect(exposure, fixed(exp1)) ~ paste(exp2),
               str_detect(exposure, fixed(exp2)) ~ paste(exp1)
             ))
  
  mv_res <- mv_res[,c(1,10,2:9)]
  
  pairwise_res <- rbind(pairwise_res, mv_res)
  
  #Conditional F stats
  colnames(exposure1)<-paste(colnames(exposure1),"1",sep = "_")
  colnames(exposure2)<-paste(colnames(exposure2),"2",sep = "_")
  
  out<-out_dat[,c(1,6,7)]
  
  dat1<- inner_join(exposure1,exposure2, by=c("SNP_1"="SNP_2")) 
  dat1<- inner_join(dat1,out, by=c("SNP_1"="SNP"))
  
  bx <- as.matrix(dat1[,c("beta_1","beta_2")])
  bxse <- as.matrix(dat1[,c("se_1","se_2")])
  dat <- MendelianRandomization::mr_mvinput(
    bx = bx, bxse = bxse, by = dat1$beta.outcome, byse = dat1$se.outcome,
    snps = dat1$SNP_1
  )
  
  mvdat <- mrmvinput_to_mvmr_format(dat)
  
  cov <- as.matrix(phe_corr)
  cor <- cov[as.double(which(str_detect(key$`colnames(phe_corr)`,paste(exp2)))),]
  cor <- cor[as.double(which(str_detect(key$`colnames(phe_corr)`,paste(exp1))))]
  Pcov <- matrix(c(1,cor,cor,1), nrow = 2)
  
  phe_cov <- phenocov_mvmr(Pcov, bxse)
  
  sres <- strength_mvmr(r_input = mvdat, gencov = phe_cov)
  pres <- pleiotropy_mvmr(r_input = mvdat, gencov = phe_cov)
  
  res <- as.data.frame(ivw_mvmr(r_input = mvdat)) %>% 
    mutate(`Estimate` = `Estimate`/v)
  
  res1 <- qhet_mvmr(mvdat, Pcov, CI = F, iterations = 100) 
}

write.table(pairwise_res, file = "pairwise.mvmr.txt", row.names = F, col.names = T)
