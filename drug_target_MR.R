################################################################################

                            # Drug scores on AVS #

################################################################################

library(TwoSampleMR)
library(MendelianRandomization)
library(TwoSampleMR)
library(data.table)
library(stringr)
library(dplyr)

setwd("~/wd/")

load("AVS_GWAS.RData")

#This is for the BOLT LMM correction
u <- 4807/(4807+458120)
v <- u*(1-u)

################################################################################
#all scores in the same file 

all_scores <- fread('~/wd/drug_targets.txt') 

all_results <- data.frame()
results<-data.frame()

nested_scores <- all_scores %>% group_split(Phenotype)

for (i in 1:length(nested_scores)) {
  print(i)
  snp <- nested_scores[[i]]
  exp_dat <- format_data(snp, type='exposure')
  n=nrow(exp_dat)
  
  
  dat <- harmonise_data(
    exposure_dat = exp_dat,
    outcome_dat = out_dat,
    action=1
    )
    
  dat <- dat[dat$mr_keep,]
  
  #Conditional F stats  
  print(mean((dat$beta.exposure*dat$beta.exposure)/(dat$se.exposure*dat$se.exposure)))
    
  if(nrow(dat)<1) {next}
    
  p =nrow(dat)
  if(p<2){
    #Wald ratio if only 1 SNP
    res0 <- mr(dat, method_list=c('mr_wald_ratio')) %>% 
      select(-c(1,2)) %>% 
      mutate(across(c(b,se,pval), ~ as.double(.x)))
    print(res0)
    results <- rbind(results, res0)
  } else {
      
    dat1 <- dat[match(snp$SNP, dat$SNP),]
    dat2 <- dat1[which(!is.na(dat1$SNP)),] %>% 
      distinct()
      
    SNPs = dat2$SNP
      
    ## use default 1kg matrix for MR accoutning for LD correlation
    ld <- ld_matrix(dat2$SNP)
    if(is.null(dim(ld))) {next}
    ld1 <- unlist(strsplit(row.names(ld), "_"))
    ld2 <- ld1[seq(1, length(ld1), 3)]
    dat3 <- dat2[which(dat2$SNP %in% ld2),]
    if(nrow(dat3)<2) {
      next
    } else {
        
      dat4 <- dat_to_MRInput(dat3, get_correlations = T)
        
      res = MendelianRandomization::mr_ivw(dat4[[1]], correl=T)
      print(res)
      dfvector <- c(exposure = dat$exposure[1], outcome=dat$outcome[1], method= 'IVW', "b"=res@Estimate, 
                    "se"=res$StdError
                    "pval"=as.numeric(res@Pvalue), 
                    "nsnp"=res@SNPs)
      dfmatrix = as.matrix(dfvector)
      df = t(dfmatrix)
      all_results <- rbind(all_results, df)
    }
  }
}

#transform the betas into ORs, etc 
final_results<- all_results %>% 
  rename(Estimate = b) %>% 
  mutate(pval= as.numeric(pval),
         Estimate = as.numeric(Estimate), se = as.numeric(se)) %>%
  mutate(Estimate =  -Estimate/v, # Reversing direction to reflect drug action
    se = se/v) %>% 
  mutate(target = case_when(
    exposure %in% c("SLC12A3","ADRB1","NCC") ~ "Systolic blood pressure", 
    exposure %in% c("PCSK9","HMGCR","NPC1L1") ~ "LDL cholesterol",
    exposure %in% c("ANGPTL4","APOC3","LPL") ~ "Triglycerides",
    exposure %in% c("LPA") ~ "Lipoprotein(a)")) %>% 
  mutate(
    OR = exp(Estimate), UCI = exp(Estimate + 1.96 * se), 
    LCI = exp(Estimate - 1.96 * se),
    OR = case_when(
      str_detect(target, "Systolic") ~ exp(log(OR)*(1/19.32)), #converting estimates to mmHg
      !str_detect(target, "Systolic") ~ OR),
    UCI=case_when(
      str_detect(target, "Systolic") ~ exp(log(UCI)*(1/19.32)), #converting estimates to mmHg
      !str_detect(target, "Systolic") ~ UCI),
    LCI=case_when(
      str_detect(target, "Systolic") ~ exp(log(LCI)*(1/19.32)), #converting estimates to mmHg
      !str_detect(target, "Systolic") ~ LCI)) 


write.table(final_results, file="drug_mr.txt",sep = "\t",row.names = F, col.names = T)
