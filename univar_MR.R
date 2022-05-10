################################################################################

            # Univariable MR analyses of AVS risk factors #

################################################################################
library(dplyr)
library(stringr)
library(tidyverse)
library(TwoSampleMR)
library(vroom)
library(data.table)

setwd("~/wd/")

load("AVS_GWAS.RData") 

################################################################################
                  # Analyses with local exposure datasets
################################################################################
#List of exposures - LDL, trigs, ApoB, Lp(a), SBP, BMI, smoking, 
# alcohol, hba1c, calcium, phosphate

var_names <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
               'chisq_lm','p','beta','se','chisq_bolt','pval', 'tmp', 'tmp1', 'Phenotype')

exp1 <- fread('~/wd/LDL_direct_clumped.txt', header=T)
exp2 <- fread('~/wd/apolipoprotein_B_clumped.txt', header=T)
exp3 <- fread('~/wd/triglycerides_clumped.txt', header=T)
exp4 <- fread('~/wd/lipoprotein_A_clumped1.txt', header = T)
exp5 <- fread('~/wd/HbA1c_clumped.txt', header=T)
exp6 <- fread('~/wd/calcium_clumped.txt', header=T)
exp7 <- fread('~/wd/phosphate_clumped.txt', header=T)
exp8 <- fread('~/wd/SBP_clumped.txt', header=T)
exp9 <- fread('~/wd/BMI_clumped.txt', header=T)
exp10 <- fread('~/wd/alcohol_frequency_clumped.txt', header=T)
exp11 <- fread('~/wd/lifetime_smoking_clumped.txt', header=T) %>% 
  mutate(BETA = BETA/0.6940093,
         SE = SE/0.6940093) #GWAS was not standardised



exp1$Phenotype <- 'LDL cholesterol'
exp2$Phenotype <- 'Apolipoprotein B'
exp3$Phenotype <- 'Triglycerides'
exp4$Phenotype <- 'Lipoprotein(a)'
exp5$Phenotype <- 'HbA1c'
exp6$Phenotype <- 'Serum calcium'
exp7$Phenotype <- 'Serum phosphate'
exp8$Phenotype <- 'Systolic blood pressure'
exp9$Phenotype <- 'Body mass index'
exp10$Phenotype <- 'Alcohol intake frequency'
exp11$Phenotype <- 'Lifetime smoking index'


exp <- rbind(exp1, exp2, exp3, exp4, exp5, exp6, exp7, exp8, exp9, exp10, exp11, fill = TRUE)
names(exp) <- var_names
exp_dat <- format_data(exp, type = "exposure")


univar_results <- data.frame()
all_sens <- data.frame()

dat <- NULL 
mr_results <- NULL
out_dat <- gwas

dat <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat,
  action = 1
  )
  
mr_results <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median",
                                        "mr_egger_regression"))  # main MR analysis

#BOLT LMM correction and convert to OR
u <- 4807/(4807+458120)
v <- u*(1-u)
mr_results <- mr_results %>% 
  dplyr::mutate(b = b/v,
                se = se/v) %>% 
  dplyr::mutate(OR = exp(b),
                UCI = exp(b + 1.96 * se),
                LCI = exp(b - 1.96 * se))

#Sensitivity analyses  
het <- mr_heterogeneity(dat) %>% select(!dplyr::starts_with("id"))
pleio <- mr_pleiotropy_test(dat) %>% 
  select(!dplyr::starts_with("id"),
         !dplyr::ends_with("y")) %>% 
  dplyr::mutate(method = "MR Egger")
sens <- full_join(het, pleio, by = c("method","exposure"))
all_sens <- rbind(all_sens, sens)


write.table(univar_results,file="univar_results.txt",sep="\t",col.names=T,row.names=F,quote=F)

write.table(all_sens,file="het_pleio_analyses.txt", sep="\t",col.names=T,row.names=F,quote=F)

nested_exps <- exp_dat %>% group_split(exposure)
f_stats <- data.frame()

#calculate and save F stats for reporting
for (i in 1:11){
  ex <- nested_exps[[i]]
  ex$F <- (ex$beta.exposure*ex$beta.exposure)/(ex$se.exposure*ex$se.exposure)
  f_stat <- mean(ex$F)
  f <- data.frame(exposure = ex$exposure[1],F = f_stat)
  f_stats <- rbind(f_stats,f)
}

write.table(f_stats,file = "univar_f_stats.txt",sep="\t)