##Swiss HIV Cohort Data
##lower function for datatsets 
d_tolower <- function(data){
  colnames(data) <- tolower(colnames(data))
  return(data)
}

#Loading SHCS data tables
shcs_pat <- d_tolower(read.csv("pat.csv"))
shcs_fup <- d_tolower(read.csv("fup.csv"))
shcs_lab <- d_tolower(read.csv("lab.csv"))
shcs_drg <- d_tolower(read.csv("med_drug.csv"))
shcs_drstop <- d_tolower(read.csv("var_stopdrug.csv"))
shcs_clin <- d_tolower(read.csv("clinical.csv"))
shcs_dis <- d_tolower(read.csv("dis.csv"))
shcs_vdis <- d_tolower(read.csv("var_disease.csv"))
shcs_vclin <- d_tolower(read.csv("var_clinical.csv"))
shcs_adh <- d_tolower(read.csv("adherence.csv"))
shcs_art <- d_tolower(read.csv("modif_art.csv"))

##first reported eacs question
first_eacs <- shcs_fup %>% 
  filter(((!is.na(cog_conc)) & cog_conc != 9) | ((!is.na(cog_slow)) & cog_slow != 9) | ((!is.na(cog_freq)) & cog_freq != 9)) %>%
  group_by(id) %>% arrange(fupdate, .by_group = TRUE) %>% slice(1) %>%
  select(id,first_eacs_date = fupdate)

##cd4 auc before first reported eacs
cd4_auc <- shcs_lab %>%
  fill(cd4, .direction = "downup") %>%
  mutate(cd4 = sqrt(cd4)) %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(labdate < first_eacs_date) %>%
  group_by(id) %>%
  filter(n() > 1) %>%
  arrange(labdate, .by_group = TRUE) %>%
  mutate(fuptime = as.numeric(difftime(labdate,lag(labdate)))) %>%
  mutate(fuptime = ifelse(is.na(fuptime),0,fuptime)) %>%
  mutate(fuptime = cumsum(fuptime)) %>%
  mutate(AUC = (fuptime - lag(fuptime)) * ((cd4 + lag(cd4)) / 2)) %>%
  mutate(AUC = ifelse(is.na(AUC),0,AUC)) %>%
  mutate(auc_cd4 = sum(AUC)/max(fuptime)) %>%
  slice(1) %>%
  select(id,auc_cd4)

##rna auc before first reported eacs
rna_auc = shcs_lab %>%
  fill(rna, .direction = "downup") %>%
  mutate(rna = log10(rna+1)) %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(labdate < first_eacs_date) %>%
  group_by(id) %>%
  filter(n() > 1) %>%
  arrange(labdate, .by_group = TRUE) %>%
  mutate(fuptime = as.numeric(difftime(labdate,lag(labdate)))) %>%
  mutate(fuptime = ifelse(is.na(fuptime),0,fuptime)) %>%
  mutate(fuptime = cumsum(fuptime)) %>%
  mutate(AUC = (fuptime - lag(fuptime)) * ((rna + lag(rna)) / 2)) %>%
  mutate(AUC = ifelse(is.na(AUC),0,AUC)) %>%
  mutate(auc_rna = sum(AUC)/max(fuptime)) %>%
  slice(1) %>%
  select(id,auc_rna)

##Efavirenz - depression AUC
##Efavirenz neurocogniteve associated
efv_neuro_associated <- shcs_drg[which(shcs_drg$drug == "EFV"),] %>%
  merge(.,shcs_drstop,by.x = "stop_why", by.y = "stop_why", all.x = TRUE, all.y = FALSE) %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(stopdate < first_eacs_date) %>%
  filter(stop_why %in% c(6.2,6.1,6)) %>%
  mutate(neuro_assoc = 1) %>%
  group_by(id) %>% slice(1) %>%
  select(id,neuro_assoc)
##relevant drug stop reasons
# 6.2 "Toxicity, predominantly from nervous system" 6.1 "Toxicity - peripheral neuropathy" 6 "Toxicity - neuropsychiatric"

##auc of Efavirenz use (in years)
nci_efv <- shcs_drg[which(shcs_drg$drug == "EFV"),] %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(startdate < first_eacs_date) %>%
  merge(.,shcs_pat[,c("id","regdate")], by = "id", all = FALSE) %>%
  mutate(stopdate = ifelse(stopdate > first_eacs_date | stopdate == "" | is.na(stopdate),first_eacs_date,stopdate)) %>%
  mutate(AUC = as.numeric(difftime(stopdate,startdate,units = "days"))) %>%
  group_by(id) %>%  
  mutate(totalfollow = as.numeric(difftime(first_eacs_date,regdate,units = "days"))) %>%
  mutate(totalfollow = ifelse(totalfollow <= 0,1,totalfollow)) %>%
  mutate(AUC_efv = sum(AUC)/totalfollow/365) %>% slice(1) %>%
  select(id,AUC_efv)


#######################
##Depression reported before first reported eacs
nci_depr <- shcs_fup %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(fupdate < first_eacs_date) %>%
  select(id,diag_psy,depression,diag_other_phys,diag_shcs_phys) %>%
  mutate(depri = (diag_psy == 1 | depression == 1 | diag_other_phys == 1 | diag_shcs_phys == 1)) %>%
  filter(!is.na(depri)) %>%
  group_by(id) %>%
  slice(1) %>%
  select(id,depri) 

#antidepression last value carried forward
nci_antidepr <- shcs_fup %>%
  filter(!is.na(antidepress)) %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(fupdate < first_eacs_date) %>%
  filter(antidepress == 1) %>%
  group_by(id) %>% slice(1) %>%
  select(id, antidepress)

##neuro disease before first eacs report
nro_dis_ids <- shcs_clin %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(clin_date < first_eacs_date) %>%
  filter(clin_id %in% c("CEI","CEH")) %>%
  mutate(dis = 1) %>%
  group_by(id) %>% slice (1) %>%
  select(id,dis)

##neuro infectious disease before first eacs report
inf_dis_ids <- shcs_dis %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(newdate < first_eacs_date) %>%
  filter(disease %in% c("DEM","ILE","PML","PNP","TOX")) %>%
  mutate(inf_dis = 1) %>%
  group_by(id) %>% slice (1) %>%
  select(id,inf_dis)

##Hep B
shcs_hep_b <- shcs_lab %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(labdate < first_eacs_date) %>%
  mutate(hepb = case_when(aghbs == "N" & antihbc == "N" &  antihbs == "N" ~ 0,
                          aghbs == "N" & antihbc == "P" &  antihbs == "P" ~ 1,
                          aghbs == "N" & antihbc == "N" &  antihbs == "P" ~ 2,
                          aghbs == "P" & antihbc == "P" &  antihbs == "N" ~ 3,
                          aghbs == "N" & antihbc == "P" &  antihbs == "N" ~ 4))  %>%
  select(id,hepb,heb_date) %>%
  filter(hepb == 1 |  hepb == 3) %>%
  mutate(hep_b = 1) %>%
  group_by(id) %>% slice (1) %>%
  select(id,hep_b)


##Hep C
shcs_hep_c <- shcs_lab %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(labdate < first_eacs_date) %>%
  select(id,antihcv,hec_date) %>%
  filter(antihcv == "P") %>%
  mutate(hep_c = 1) %>%
  group_by(id) %>% slice (1) %>%
  select(id,hep_c)

##cocain before first eacs report
coca_ids <- shcs_fup %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(fupdate < first_eacs_date) %>%
  filter(coca_iv == 1 | coca_ni == 1) %>%
  mutate(coca = 1) %>%
  group_by(id) %>% slice (1) %>%
  select(id,coca)

##heroin before first eacs report
hero_ids <- shcs_fup %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(fupdate < first_eacs_date) %>%
  filter(hero_iv == 1 | hero_ni == 1) %>%
  mutate(hero = 1) %>%
  group_by(id) %>% slice (1) %>%
  select(id,hero)

##canabis before first eacs report
cana_ids = shcs_fup %>%
  merge(.,first_eacs, by = "id", all = FALSE) %>%
  filter(fupdate < first_eacs_date) %>%
  filter(cana_ni == 1 ) %>%
  group_by(id) %>% slice (1) %>%
  select(id,cana_ni)

##merge all covariables into into one file
eacs_covar <- shcs_pat %>%
  select(id,born,ethnicity,regdate,sex,risk,education) %>%
  mutate(ethnicity = paste0("eth",ethnicity)) %>%
  mutate(risk = paste0("risk",risk)) %>%
  mutate(sex = paste0("sex",sex)) %>%
  mutate(education = paste0("edu",education)) %>%
  mutate(riskgroup = case_when(risk == "risk0" ~ "risk0",
                               risk == "risk1" ~ "risk1",
                               risk == "risk2" ~ "risk2",
                               risk == "risk3" ~ "risk3",
                               risk == "risk4" ~ "risk3",
                               risk == "risk5" ~ "risk0",
                               risk == "risk6" ~ "risk0",
                               risk == "risk7" ~ "risk0",
                               risk == "risk9" ~ "risk0")) %>%
  mutate(riskgroup = factor(riskgroup, levels = c("risk1","risk2","risk3","risk0"))) %>%
  mutate(education = case_when(education == "edu1" ~ "edu1",
                               education == "edu2" ~ "edu2",
                               education == "edu3" ~ "edu3",
                               education == "edu4" ~ "edu3",
                               education == "edu5" ~ "edu3",
                               education == "edu6" ~ "edu3",
                               education == "edu7" ~ "edu3",
                               education == "edu0" ~ "edu4",
                               education == "edu9" ~ "edu4",
                               TRUE ~ "edu4")) %>%
  mutate(education = factor(education, levels = c("edu1","edu2","edu3","edu4"))) %>%
  merge(.,first_eacs,by = "id", all = FALSE) %>%
  mutate(age = (year(first_eacs_date)-born)/10) %>%
  merge(.,eacs_sop_non_nam,by.x = "id" ,by.y = "ID", all = FALSE) %>%
  merge(.,rna_auc,by = "id" , all.x = FALSE, all.y = FALSE) %>%
  merge(.,cd4_auc,by = "id" , all.x = FALSE, all.y = FALSE) %>%
  merge(.,efv_neuro_associated,by = "id" , all.x = TRUE, all.y = FALSE) %>%
  mutate(neuro_assoc = ifelse(is.na(neuro_assoc),paste0("a_efv_assoc"),paste0("z_efv_assoc"))) %>%
  merge(.,nci_efv,by = "id" , all.x = TRUE, all.y = FALSE) %>%
  mutate(AUC_efv = ifelse(is.na(AUC_efv),0,AUC_efv)) %>%
  merge(.,nci_depr,by = "id" , all.x = TRUE, all.y = FALSE) %>%
  mutate(depri = ifelse(is.na(depri),paste0("a_depr"),paste0("z_depr"))) %>%
  merge(.,nci_antidepr,by = "id" , all.x = TRUE, all.y = FALSE) %>%
  mutate(antidepress = ifelse(is.na(antidepress),paste0("a_anti_de"),paste0("z_anti_de"))) %>%
  merge(.,shcs_hep_c,by = "id" , all.x = TRUE, all.y = FALSE) %>%
  mutate(hep_c = ifelse(is.na(hep_c),paste0("a_hep_c"),paste0("z_hep_c"))) %>%
  merge(.,shcs_hep_b,by = "id" , all.x = TRUE, all.y = FALSE) %>%
  mutate(hep_b = ifelse(is.na(hep_b),paste0("a_hep_b"),paste0("z_hep_b"))) %>%
  merge(.,nro_dis_ids,by = "id" , all.x = TRUE, all.y = FALSE) %>%
  mutate(dis = ifelse(is.na(dis),paste0("a_dis"),paste0("z_dis"))) %>%
  merge(.,inf_dis_ids,by = "id" , all.x = TRUE, all.y = FALSE) %>%
  mutate(inf_dis = ifelse(is.na(inf_dis),paste0("a_inf_dis"),paste0("z_inf_dis"))) %>%
  merge(.,coca_ids,by = "id" , all.x = TRUE, all.y = FALSE) %>%
  mutate(coca = ifelse(is.na(coca),paste0("a_coca"),paste0("z_coca"))) %>%
  merge(.,hero_ids,by = "id" , all.x = TRUE, all.y = FALSE) %>%
  mutate(hero = ifelse(is.na(hero),paste0("a_hero"),paste0("z_hero"))) %>%
  merge(.,cana_ids,by = "id" , all.x = TRUE, all.y = FALSE) %>%
  mutate(cana_ni = ifelse(is.na(cana_ni),paste0("a_cana"),paste0("z_cana"))) 

write.csv(eacs_covar,"eacs_pheno_covar.csv")
