#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

##Libraries
library(tobitnet)
library(tidyverse)
library(doMC)

registerDoMC(cores = 32L)
print(args)
combinations <- expand_grid(c("B","all"),c("freq","conc","slow","all"),c(5,10,20)) %>% as.data.frame()

subtype <- combinations[as.numeric(args[1]),1]
phenox <- combinations[as.numeric(args[1]),2]
topx <- combinations[as.numeric(args[1]),3]

print(topx)
print(phenox)
print(subtype)

top_snp_dat_for_crossval <- function(select_snps,topx,phenox){
  
   select_snps <- select_snps %>% 
    filter(pheno == phenox) %>%
    group_by(position,region) %>%
     arrange(p.value, .by_group = TRUE) %>%
     slice(1) %>%
     ungroup() %>%
    arrange(p.value) %>%
    slice(1:topx)
  
  combi_dat <- foreach(i = 1:topx, .combine = "rbind") %dopar% {

    return(raw_gwas[[select_snps$region[i]]][[1]] %>%
             filter(position == paste0("pos",select_snps$position[i])) %>%
             mutate(region = select_snps$region[i])
    )
  }
  
  merged_top_snps <- combi_dat %>% mutate(region_pos = paste0(region,"_",position)) %>% 
    select(id,region_pos,aa) %>%
    pivot_wider(names_from = region_pos, values_from = aa) %>%
    mutate_at(c(2:(topx+1)), ~replace(., is.na(.), "X"))
  
  merged_top_snps <- merge(merged_top_snps,pheno_covar_data,by = "id", all = FALSE)
  
  merged_top_snps <- merge(merged_top_snps,pca,by = "id", by.y = "id", all = FALSE)
  
  merged_top_snps <- merged_top_snps %>% mutate(drugs = ifelse(coca == "z_coca" | hero == "z_hero" | cana_ni == "z_cana","z_drugs","a_drugs")) %>%
    mutate(dis = ifelse(dis == "z_dis" | inf_dis == "z_inf_dis" | neuro_assoc == "z_efv_assoc","z_dis","a_dis")) %>%
    mutate(ethnicity = ifelse(ethnicity == "eth9" | ethnicity == "ethNA" |  ethnicity == "eth0","eth4",ethnicity)) %>%
    mutate(PC1 = PC1 * 50,PC2 = PC2 * 50,PC3 = PC3 * 50,PC4 = PC4 * 50,PC5 = PC5 * 50,PC6 = PC6 * 50,PC7 = PC7 * 50,PC8 = PC8 * 50,
           PC9 = PC9 * 50,PC10 = PC10 * 50) %>% mutate(risk = case_when(risk == "risk0" ~ "risk3",
                                                                        risk == "risk1" ~ "risk1",
                                                                        risk == "risk2" ~ "risk2",
                                                                        risk == "risk3" ~ "risk3",
                                                                        risk == "risk4" ~ "risk3",
                                                                        risk == "risk5" ~ "risk3",
                                                                        risk == "risk6" ~ "risk3",
                                                                        risk == "risk7" ~ "risk3",
                                                                        risk == "risk9" ~ "risk3",
                                                                        TRUE ~ "risk3"))
  
  return(merged_top_snps)
}

crossval_tobit_top_snps <- function(gwas_res_data,topx,phenox){

  top_snp_dat <- top_snp_dat_for_crossval(gwas_res_data,topx,phenox)
  
  covars_cv <- c("pheno",colnames(top_snp_dat[c(2:(topx+1))]),"sex","age","education","risk","auc_rna","auc_cd4","ethnicity","AUC_efv",
                 "antidepress","depri","drugs","dis","hep_c","hep_b",
                 "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
  
  top_snp_dat <- top_snp_dat %>% 
    mutate(pheno = get(paste0("auc_",phenox,"_avg"))) %>%
    select(covars_cv)
  
  covars_cv_cat <- c("sex","education","risk","ethnicity","antidepress","depri",
                     "drugs","dis","hep_c","hep_b")
  
  top_snp_dat <- fastDummies::dummy_cols(top_snp_dat, select_columns = c(covars_cv_cat,colnames(top_snp_dat[c(2:(topx+1))])), remove_selected_columns = TRUE, remove_most_frequent_dummy = TRUE)
  top_snp_dat <- top_snp_dat %>% mutate_if(is.character,as.numeric)
  
  trains_idx <- sample(size = nrow(top_snp_dat)*0.8,seq(1:nrow(top_snp_dat)))
  
  train_tobit <- top_snp_dat[trains_idx,]
  test_tobit <- top_snp_dat[-trains_idx,]
  
  ##Covar only
  elastic_tobit_covars_cv <- cv.tobitscad(y = train_tobit$pheno, x = as.matrix(train_tobit[,-grep("_pos|pheno",colnames(train_tobit))]), nfolds = 10)
  elastic_tobit_covars <- tobitscad(y = train_tobit$pheno, x = as.matrix(train_tobit[,-grep("_pos|pheno",colnames(train_tobit))]),  lambda = c(elastic_tobit_covars_cv$lambda.min,elastic_tobit_covars_cv$lambda.1se))
  test_preds_covars <- t(as.matrix(elastic_tobit_covars$beta[])) %*% t(as.matrix(test_tobit[,-grep("_pos|pheno",colnames(test_tobit))])) + elastic_tobit_covars$b0
  test_preds_covars <- ifelse(test_preds_covars < 0,0,test_preds_covars)
  
  ##full
  elastic_tobit_full_cv <- cv.tobitscad(y = train_tobit$pheno, x = as.matrix(train_tobit[,-c(1)]), early.stop = FALSE, nfolds = 10)
  elastic_tobit_full <- tobitscad(y = train_tobit$pheno, x = as.matrix(train_tobit[,-c(1)]),  lambda = c(elastic_tobit_full_cv$lambda.min,elastic_tobit_full_cv$lambda.1se))
  test_preds_full <- t(as.matrix(elastic_tobit_full$beta[])) %*% t(as.matrix(test_tobit[,-c(1)])) + elastic_tobit_full$b0
  test_preds_full <- ifelse(test_preds_full < 0,0,test_preds_full)

  return(c(cor(t(test_preds_covars),test_tobit$pheno)**2,cor(t(test_preds_full),test_tobit$pheno)**2))
}

load("./meta_data.RData")

if(subtype == "B"){
  gwas_dat <- gwas_res_B_clean
  raw_gwas <- B_gwas
}else{
  gwas_dat <- gwas_res_all_clean
  raw_gwas <- all_gwas
}

print("start foreach")

foreach::foreach(i = 1:10000, .combine = "rbind", .errorhandling = "remove", .verbose = TRUE) %dopar% {
  set.seed(i+10000000)
  write.csv(crossval_tobit_top_snps(gwas_dat,topx,phenox),paste0("./crossval_res/",subtype,"_",phenox,"_",topx,"_crossval_",i,"_v3_",Sys.Date(),".csv"))	
}
