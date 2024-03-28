##Libraries
library(tidyverse)
library(stringr)
library(forcats)
library(AER)
library(broom)
library(lubridate)
library(svMisc)
library(ape)
library(phangorn)
library(ggrepel)
library(phytools)
library(qdapTools)
library(foreach)
library(doParallel)

##NGS meta data
shcs_ngs_docu <- read.csv("NGS_samples_list_all\ runs_NGS_R1R2.csv")

#Load matching table between SHCS- and ZPHI IDs
matching <- readxl::read_excel("info.xlsx", col_names = TRUE)

##Phenotype and covariables data
pheno_covar_data <- read.csv("eacs_pheno_covar.csv")
pheno_covar_data <- pheno_covar_data %>% group_by(id) %>% slice(1) 

##HIV Subtype data
shcs_sub <- read.csv("subtypes.txt")

#Match the ZPHI ID with an SHCS ID
for(i in 1:length(shcs_ngs_docu$SHCS.ID)){
  for(y in 1:length(matching$patnr)){
    if(!is.na(matching$cohortnr[y]) && !is.na(shcs_ngs_docu$ZPHI.ID[i]) && matching$patnr[y] == shcs_ngs_docu$ZPHI.ID[i]){
      shcs_ngs_docu$SHCS.ID[i]=matching$cohortnr[y]
    }
  }
}

##Function to get fasta file to csv data
reformat_seq_to_csv <- function(msa){
  jk <- as.data.frame(msa)
  jk$V1 <- msa@ranges@NAMES
  ##Rename column names
  names(jk) <- c("seq","id")
  jk <- jk %>% select(id,seq)
  ##Dataset for sequencing data
  maxfa <- jk %>% separate(seq, into = paste("V", 1:(max(nchar(jk$seq))+1), sep = ""), sep = "")
  maxfa <- maxfa %>%
    dplyr::select(-V1)
  colnames(maxfa)[2:ncol(maxfa)] <- paste0("V",1:(ncol(maxfa)-1))
  return(maxfa)
}

##Rename UUIDs unique for a sequence into respective Patient IDs
rename_uuids <- function(msa,region,shcs_ngs_docu){
  hxb2_region <- paste0("HXB2_",toupper(region))
  msa <- msa[msa@ranges@NAMES %in% shcs_ngs_docu$base_uuid | msa@ranges@NAMES == hxb2_region]
  
  seq_names <- as.data.frame(msa@ranges@NAMES)
  colnames(seq_names) <- "uuid"
  seq_names <- seq_names %>% merge(.,shcs_ngs_docu[,c("SHCS.ID","base_uuid","sample.date")], by.x = "uuid", by.y = "base_uuid", all.x = TRUE, all.y = FALSE)
  seq_names$SHCS.ID[seq_names$uuid == hxb2_region] <- hxb2_region
  seq_names <- seq_names[order(match(seq_names$uuid, msa@ranges@NAMES)),]
  msa@ranges@NAMES[msa@ranges@NAMES != hxb2_region] <- seq_names$SHCS.ID[seq_names$SHCS.ID != hxb2_region]
  msa <- msa[!duplicated(msa@ranges@NAMES)]
  
  return(msa)
}

##Read multiple sequence alignment file
read_aa_msa <- function(region,type,file_annote,subtype,kind,res){
  
  msa = Biostrings::readAAStringSet(paste0("/alignment/",region,kind,"/nci_",type,"_",subtype,"_",region,kind,"_aligned_",file_annote,".fa"))
  
  if((kind != "_resdb")) {msa = rename_uuids(msa,region,shcs_ngs_docu)}
  
  msa <- reformat_seq_to_csv(msa)

  return(msa)
  
}

##Extract HXB2 positions for referencing later
hxb2_ref_pos <- function(region,msa){
  hxb2_reg <- paste0("HXB2_",toupper(region))
  
  #remove non HXB2 positions for imv sequences or when HXB2 is still included
  seqs_to_keep_neuro <- as.matrix(msa %>%
                                   filter(id == hxb2_reg))
  
  seqs_to_keep_neuro <- reshape2::melt(seqs_to_keep_neuro)
  
  seqs_to_keep_neuro <- seqs_to_keep_neuro  %>%
    dplyr::select(-Var1) %>%
    filter(value != hxb2_reg) %>%
    filter(value != "-") %>%
    mutate(hxb2_index = as.numeric(rownames(.))) %>%
    mutate(hxb2_index_overall = as.numeric(hxb2_index) )
  
  seqs_to_keep_neuro$gwas_pos <- as.numeric(gsub("V", "", seqs_to_keep_neuro$Var2)) 
  seqs_to_keep_neuro$msa_pos <- seqs_to_keep_neuro$gwas_pos 
  
  return(seqs_to_keep_neuro)
  
}

##Remove HXB2 sequence 
remove_hxb2_seq <- function(region,msa){
  hxb2_reg <- paste0("HXB2_",toupper(region))
  #remove hxb2
  msa <- msa %>%
    filter(id != hxb2_reg) 
  
  colnames(msa)[2:(ncol(msa))] <- paste0("pos",seq(1,(ncol(msa)-1)))
  
  return(msa)
}

##Frequencies of amino acids for each position
vcf_allels_counts <- function(msa){
  
  vcf_allels_b_aa <- NULL
  vcf_allels_b_aa <- as.data.frame(matrix(ncol = 1, nrow = (ncol(msa)-1)))
  
  for(i in 2:(ncol(msa))){
    vcf_allels_b_aa$V1[i-1] = i - 1
    vcf_allels_b_aa$F[i-1] = sum(msa[,i] == "F")
    vcf_allels_b_aa$L[i-1] = sum(msa[,i] == "L")
    vcf_allels_b_aa$I[i-1] = sum(msa[,i] == "I")
    vcf_allels_b_aa$M[i-1] = sum(msa[,i] == "M")
    vcf_allels_b_aa$V[i-1] = sum(msa[,i] == "V")
    vcf_allels_b_aa$S[i-1] = sum(msa[,i] == "S")
    vcf_allels_b_aa$P[i-1] = sum(msa[,i] == "P")
    vcf_allels_b_aa$T[i-1] = sum(msa[,i] == "T")
    vcf_allels_b_aa$A[i-1] = sum(msa[,i] == "A")
    vcf_allels_b_aa$Y[i-1] = sum(msa[,i] == "Y")
    vcf_allels_b_aa$H[i-1] = sum(msa[,i] == "H")
    vcf_allels_b_aa$Q[i-1] = sum(msa[,i] == "Q")
    vcf_allels_b_aa$N[i-1] = sum(msa[,i] == "N")
    vcf_allels_b_aa$K[i-1] = sum(msa[,i] == "K")
    vcf_allels_b_aa$D[i-1] = sum(msa[,i] == "D")
    vcf_allels_b_aa$E[i-1] = sum(msa[,i] == "E")
    vcf_allels_b_aa$C[i-1] = sum(msa[,i] == "C")
    vcf_allels_b_aa$W[i-1] = sum(msa[,i] == "W")
    vcf_allels_b_aa$R[i-1] = sum(msa[,i] == "R")
    vcf_allels_b_aa$G[i-1] = sum(msa[,i] == "G")
    vcf_allels_b_aa$Z[i-1] = sum(msa[,i] == "-")
    vcf_allels_b_aa$X[i-1] = sum(msa[,i] == "X")
  }
  
  vcf_allels_b_aa$var <- apply(vcf_allels_b_aa[2:22],1, function(x) sum(x>1))
  
  ##gwas as bases as characters
  msa <- msa[,c(1,2:ncol(msa))]
  msa[,2:ncol(msa)][msa[,2:ncol(msa)] == "-"] <- "Z" 
  return(msa)
}

##Wrapper function for alignment formating
gwas_formatting_main <- function(region,type,file_annote,subtype,kind,impute,res){
  
  msa <- read_aa_msa(region,type,file_annote,subtype,kind,impute,res)
  
  seq_to_keep <- hxb2_ref_pos(region,msa)
  
  msa <- remove_hxb2_seq(region,msa)
  
  msa <- vcf_allels_counts(msa)
  
  return(list(msa,seq_to_keep))
}

##Build combined covariable/pheno data frame
pheno_covar_data_formatting = function(region,pheno_covar_data,gwas_reg,subtype,kind,res){
  
  if(kind == "_resdb") {
    source <- "resdb"}
  else{
    source <- "wga"
  }
  
  pca <- read.csv(paste0("/eigensoft/",region,res,"_",source,"/",region,"_",subtype,res,"_",source,".evec"), sep = "", header = FALSE)[-1,]
  
  colnames(pca) <- c("id","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","famid")
  
  data <- pheno_covar_data %>% group_by(id) %>% slice(1) %>%
    merge(.,pca[,-12], by = "id", all.x = FALSE, all.y = FALSE) %>% 
    filter(id %in% gwas_reg$id)
  
  ##covar variables reformatting
  data <- data %>% mutate(drugs = ifelse(coca == "z_coca" | hero == "z_hero" | cana_ni == "z_cana","z_drugs","a_drugs")) %>%
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
  return(data)
}

##Clean alignments
region_set <- function(gwas_b_char_aa_set,covar_set,minvar){
  
  gwas_single_test <- gwas_b_char_aa_set %>%  pivot_longer(!id, names_to = "position", values_to = "aa")
  gwas_single_test <- merge(gwas_single_test,covar_set, by.x = "id", by.y = "id", all = FALSE)
  
  gwas_single_test <- gwas_single_test %>% 
    filter(aa != "Z")
  mincount <- gwas_single_test %>% dplyr::select(position,aa) %>% tidyr::gather(position,aa) %>% dplyr::count(position,aa, .drop=FALSE) %>% filter(n >= minvar)
  refcount <- mincount %>% group_by(position) %>% arrange(desc(n), .by_group = TRUE) %>% slice(1) %>% select(position,ref = aa,n_ref = n)
  pos_to_take <- gwas_single_test %>% dplyr::select(position,aa) %>% tidyr::gather(position,aa) %>% dplyr::count(position,aa, .drop=FALSE) %>% filter(n >= minvar) %>% group_by(position) %>% mutate(n_aas = sum(n>1)) %>% filter(n_aas > 1)
  gwas_single_test <- gwas_single_test  %>% filter(paste0(position,aa) %in% paste0(pos_to_take$position,pos_to_take$aa))
  gwas_single_test <- merge(gwas_single_test,mincount,by = c("position","aa"), all = FALSE) %>% merge(.,refcount, by = "position") %>% filter(n_ref > 100)
  return(gwas_single_test)
}

##GWAS ready dataset wrapper function
gwas_ready_formatting_final <- function(region,type,pheno_covar_data,file_annote,maf = 20,subtype = "B",kind="",impute, res = ""){
  
  ##geno type formatting (genotype file and hxb2 psotion relation as two output tables)
  gwas_geno <- gwas_formatting_main(region,type,file_annote,subtype,kind,impute,res)
  
  ##combined dataset pheno covar
  covar_data <- pheno_covar_data_formatting(region,pheno_covar_data,gwas_geno[[1]],subtype,kind,res)
  
  ##combined dataset geno pheno covar
  gwas_final <- gwas_geno[[1]] %>% region_set(.,covar_data,maf)
  
  return(list(gwas_final,gwas_geno[[2]]))
}

##GWAS
gwas_multi_var <- function(gwas_single_test_set){
  
  ##filter out "dis" with single factor level
  gwas_single_test_set <- gwas_single_test_set %>% group_by(position) %>% na.omit() %>% mutate(n_dis = length(unique(dis))) %>%
    filter(n_dis > 1) 
  
  all_nci <- gwas_single_test_set %>%
    group_by(position) %>% 
    do(lm = broom::tidy(summary(AER::tobit(auc_all_avg ~ fct_infreq(aa)+sex+age+education+risk+auc_rna+auc_cd4+ethnicity+AUC_efv+antidepress+depri+drugs+dis+hep_c+hep_b+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = .))$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "all")
  
  slow_nci <- gwas_single_test_set %>%
    group_by(position) %>% 
    do(lm = broom::tidy(summary(AER::tobit(auc_slow_avg ~ fct_infreq(aa)+sex+age+education+risk+auc_rna+auc_cd4+ethnicity+AUC_efv+antidepress+depri+drugs+dis+hep_c+hep_b+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = .))$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "slow")
  
  conc_nci <- gwas_single_test_set %>%
    group_by(position) %>% 
    do(lm = broom::tidy(summary(AER::tobit(auc_conc_avg ~ fct_infreq(aa)+sex+age+education+risk+auc_rna+auc_cd4+ethnicity+AUC_efv+antidepress+depri+drugs+dis+hep_c+hep_b+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = .))$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "conc")
  
  freq_nci <- gwas_single_test_set %>%
    group_by(position) %>% 
    do(lm = broom::tidy(summary(AER::tobit(auc_freq_avg ~ fct_infreq(aa)+sex+age+education+risk+auc_rna+auc_cd4+ethnicity+AUC_efv+antidepress+depri+drugs+dis+hep_c+hep_b+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = .))$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "freq")
  
  all_nci <- rbind(all_nci,slow_nci,conc_nci,freq_nci)
  
  return(all_nci)
}

##GWAS adjusted only with PCs for population structure
gwas_multi_var_unadjusted <- function(gwas_single_test_set){
  
  ##filter out "dis" with single factor level
  gwas_single_test_set <- gwas_single_test_set %>% group_by(position) %>% na.omit() %>% mutate(n_dis = length(unique(dis))) %>%
    filter(n_dis > 1) 
  
  all_nci <- gwas_single_test_set %>%
    group_by(position) %>% 
    do(lm = broom::tidy(summary(AER::tobit(auc_all_avg ~ fct_infreq(aa)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = .))$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "all")
  
  slow_nci <- gwas_single_test_set %>%
    group_by(position) %>% 
    do(lm = broom::tidy(summary(AER::tobit(auc_slow_avg ~ fct_infreq(aa)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = .))$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "slow")
  
  conc_nci <- gwas_single_test_set %>%
    group_by(position) %>% 
    do(lm = broom::tidy(summary(AER::tobit(auc_conc_avg ~ fct_infreq(aa)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = .))$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "conc")
  
  freq_nci <- gwas_single_test_set %>%
    group_by(position) %>% 
    do(lm = broom::tidy(summary(AER::tobit(auc_freq_avg ~ fct_infreq(aa)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = .))$coefficients)) %>%
    unnest(lm) %>%
    mutate(pheno = "freq")
  
  all_nci <- rbind(all_nci,slow_nci,conc_nci,freq_nci)
  
  return(all_nci)
}

##Clean GWAS results
filter_aa_results <- function(results,raw_data,seqs_to_keep){
  
  results <- results %>% 
    filter(grepl("aa",term)) %>% 
    mutate(position = gsub("pos","",position)) %>%
    mutate(term = gsub("fct\\_infreq\\(aa\\)","",term)) 
  
  results <- results %>%
    mutate(position = as.numeric(position)) %>%
    merge(.,seqs_to_keep[,c("msa_pos","hxb2_index")], by.x = "position", by.y = "msa_pos", all = TRUE) %>%
    group_by(pheno) %>%  
    arrange(position, .by_group = TRUE) %>%
    tidyr::fill(hxb2_index, .direction = "down")
  
  results <- results[complete.cases(results), ]
  
  results <- raw_data %>% mutate(position = gsub("pos","",position)) %>% group_by(position,aa) %>%
    slice(1) %>% select(position,n,ref,n_ref) %>%
    merge(results,., by.x = c("position","term"), by.y = c("position","aa"), all.x = TRUE, all.y = FALSE) %>%
    
    return(results)
}

##Manhattan plot
manhattan_gg <- function(data,multiple_testing_threshold){

  data <- data %>% mutate(estimate = ifelse(estimate >= 0 ,"increase","decrease"))

  data$pheno[data$pheno == "all"] = "Combination\n "
  data$pheno[data$pheno == "conc"] = "Concentration\ndifficulties\n "
  data$pheno[data$pheno == "freq"] = "Frequent\nmemory loss\n "
  data$pheno[data$pheno == "slow"] = "Cognitive\nslowing\n "
  data$region[data$region == "tat"] = "Tat"
  data$region[data$region == "gag"] = "Gag"
  data$region[data$region == "env"] = "Env"
  data$region[data$region == "pol"] = "Pol"
  data$region[data$region == "vif"] = "Vif"
  data$region[data$region == "vpr"] = "Vpr"
  data$region[data$region == "vpu"] = "Vpu"
  data$region[data$region == "nef"] = "Nef"
  data$region[data$region == "rev"] = "Rev"
  data <- data %>% mutate(region = factor(region, levels = c("Gag","Pol","Vif","Vpr","Tat","Rev","Env","Vpu","Nef")))
  data <- data %>% mutate(pheno = factor(pheno, levels = c("Frequent\nmemory loss\n ","Concentration\ndifficulties\n ",
                                                          "Cognitive\nslowing\n ","Combination\n "))) 
  data$dna <- as.numeric(data$dna)

  data %>%
    mutate(p.value = -log10(p.value)) %>%
    mutate(size_p = (exp(p.value)/exp(max(p.value)))*3) %>%
    mutate(alpha_p = (p.value/max(p.value))*1) %>%
    ggplot(.) + geom_point(aes(x = dna, y = p.value, size = size_p, colour = region, shape = pheno)) +
    xlab(paste0("HIV genome [Base position]")) + ylab("Significance [-log10(pvalue)]") + 
    geom_hline(yintercept = multiple_testing_threshold) + 
    theme_minimal() +  labs(colour = "", shape = "") +  ylim(0,6.1) +
    guides(size = FALSE,colour = guide_legend(override.aes = list(size=5)),shape = guide_legend(override.aes = list(size=5))) +
    theme(text = element_text(size = 40)) +
    geom_text_repel(aes(x = dna, y = p.value,label = ifelse(p.value>3.5,paste0(ref,hxb2_index,term),'')),size = 6) 
}

##Effective test size for helper function for "outer"-function
get.V2 <- function(x, y){
  rcompanion::cramerV(x, y, bias.correct = TRUE)
}

##Effective test size
effective_test_size <- function(subtype,data,file_annote) {

  eff_tests <- foreach(region_indx = c("env","tat","gag","pol","rev","vif","vpr","vpu","nef"), .combine = "rbind") %dopar% {
    
    eff_tes_siz <- Biostrings::readAAStringSet(paste0("/alignment/",region_indx,"/nci_AA_",subtype,"_",region_indx,"_aligned_",file_annote,".fa"))
    
    for_cramer <- as.matrix(eff_tes_siz[-1,])
    
    colnames(for_cramer) <- paste0("pos",seq(1:ncol(for_cramer)))
    
    for_cramer <- for_cramer[,colnames(for_cramer) %in% unique(paste0("pos",data$position[data$region == region_indx]))]
    
    cramers_v <- v_outer(for_cramer, get.V2)
    
    eigen <- as.data.frame(eigen(cramers_v)$values)
    colnames(eigen) <- "eigen"
    
    effec_test <- eigen %>%
      mutate(eigen = abs(eigen)) %>%
      dplyr::arrange(desc(eigen)) %>%
      mutate(n = dplyr::row_number())%>%
      mutate(totaleigen = sum(eigen)) %>%
      mutate(eigenthreshold = totaleigen*0.995) %>%
      mutate(eigensum = cumsum(eigen)) %>%
      mutate(eigensum_perc = eigensum/totaleigen) %>%
      mutate(effect_tests = eigensum >= eigenthreshold) %>%
      filter(effect_tests == TRUE) %>% select(n) %>% min()
    
    eff_test = length(data$position[data$region == region_indx]) * (as.numeric(effec_test)/as.numeric(ncol(for_cramer)))
    
    return(eff_test)
  }
  
  return(eff_tests)
}

