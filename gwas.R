##Required functions
source("/gwas_functions.R")

##File name annotation
file_annote <- ""

##Genomic regions 
regions <- c("env","gag","pol","tat","nef","rev","vpu","vif","vpr")

##Prepare dataset for HIV subtype B
B_gwas <- lapply(regions,gwas_ready_formatting_final,pheno_covar_data = pheno_covar_data, type = "AA", subtype = "B", maf = 20, file_annote = "2023-04-03", impute = FALSE)
names(B_gwas) <- c(regions)

##Prepare dataset for HIV all HIV subtypes
all_gwas <- lapply(regions,gwas_ready_formatting_final,pheno_covar_data = pheno_covar_data, type = "AA", subtype = "all", maf = 20, file_annote = "2023-03-21", impute = FALSE)
names(all_gwas) <- c(regions)

##Run gwas
registerDoParallel(cores = 8L)
##Subtype B
gwas_res_B <- foreach(x = names(B_gwas)) %dopar% {
  
  gwas_multi_var(B_gwas[[x]][[1]])
  
}
names(gwas_res_B) <- names(B_gwas)

##All subtypes
gwas_res_all <- foreach(x = names(all_gwas)) %dopar% {
  
  gwas_multi_var(all_gwas[[x]][[1]])
}
names(gwas_res_all) <- names(all_gwas)

##Clean gwas results

##HIV genes and positions according to HXB2
position_hxb2_ref <- list(c(seq(6046,8794,3)),#env
                          c(seq(791,2291,3)),#gag
                          c(seq(2066,5095,3)),#pol
                          c(c(seq(5831,6045,3),seq(8379,8469,3))),#tat
                          c(seq(8797,9417,3)),#nef
                          c(c(seq(5970,6045,3),seq(8381,8653,3))),#rev
                          c(seq(6062,6310,3)),#vpu
                          c(seq(5041,5619,3)),#vif
                          c(c(seq(5559,5771,3),seq(5773,5850,3)))#vpr
)
names(position_hxb2_ref) <- regions

##Subtype B
gwas_res_B_clean <- foreach(x = names(gwas_res_B), .combine = "rbind") %dopar% {
  posis <- as.data.frame(cbind(position_hxb2_ref[[x]],1:length(position_hxb2_ref[[x]])))
  colnames(posis) <- c("dna","aa")
  filter_aa_results(gwas_res_B[[x]],B_gwas[[x]][[1]],B_gwas[[x]][[2]]) %>%
    merge(.,posis, by.x = "hxb2_index", by.y = "aa", all.x = FALSE, all.y = FALSE) %>%
    mutate(region = x)
}

##All subtypes
gwas_res_all_clean <- foreach(x = names(gwas_res_all), .combine = "rbind") %dopar% {
  posis <- as.data.frame(cbind(position_hxb2_ref[[x]],1:length(position_hxb2_ref[[x]])))
  colnames(posis) <- c("dna","aa")
  filter_aa_results(gwas_res_all[[x]],all_gwas[[x]][[1]],all_gwas[[x]][[2]]) %>%
    merge(.,posis, by.x = "hxb2_index", by.y = "aa", all.x = FALSE, all.y = FALSE) %>%
    mutate(region = x)
}

##Filter minor frequencies and minimum major frequencies
gwas_res_all_clean <- gwas_res_all_clean %>% filter(n_ref > 300 & n > 60)
gwas_res_B_clean <- gwas_res_B_clean %>% filter(n_ref > 300 & n > 60)

##Multiple testing threshold
registerDoParallel(cores = 8L)
effs <- effective_test_size("B",gwas_res_B_clean,"2023-04-03")
effs_all <- effective_test_size("all",gwas_res_all_clean,"2023-03-21")

multiple_testing_threshold_B <- -log10(0.05/(sum(effs)))
multiple_testing_threshold_all <- -log10(0.05/(sum(effs_all)))

##Manhattan plots
amino_manhattan_b <- manhattan_gg(gwas_res_B_clean,multiple_testing_threshold_B)
amino_manhattan_all <- manhattan_gg(gwas_res_all_clean,multiple_testing_threshold_all)

##Top 50 
gwas_res_all_clean %>% arrange(p.value) %>% slice(1:50) %>%
  select(Phenotype = pheno,Protein = region,Position = hxb2_index, wt = ref, mut = term, n_wt = n_ref, n_term = n, effect = estimate, se = std.error, p = p.value)
gwas_res_B_clean %>% arrange(p.value) %>% slice(1:50) %>%
  select(Phenotype = pheno,Protein = region,Position = hxb2_index, wt = ref, mut = term, n_wt = n_ref, n_term = n, effect = estimate, se = std.error, p = p.value)

##QQ plots
qqman::qq(gwas_res_all_clean$p.value)
qqman::qq(gwas_res_B_clean$p.value)

##Genomic lambda
median(qchisq(gwas_res_all_clean$p.value, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
median(qchisq(gwas_res_B_clean$p.value, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
