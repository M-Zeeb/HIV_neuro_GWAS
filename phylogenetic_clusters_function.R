##Generates data file including a cluster variable extracted from a phylogenetic tree

region = c("env","gag","pol","tat","nef","rev","vpu","vif","vpr")
subtype = c("all","B")
seqtype = c("_resdb","_wga")

xttobit_choices = as.data.frame(expand_grid(region,subtype,seqtype))

xttobit_choices$file_annote = "..."
registerDoParallel(cores = 3L)
foreach(x = 1:nrow(xttobit_choices)) %dopar% {
  
  xtt_data(xttobit_choices$region[x],
           xttobit_choices$subtype[x],
           xttobit_choices$seqtype[x],
           xttobit_choices$file_annote[x])
  
}

##Functions

cluster_for_xtt = function(subtres){
  
  ##calculated maxdistance and number of tips for each subtree
  cluster_max_dist = foreach(x = 1:length(subtres), .combine = "rbind") %do% {
    
    return(cbind(x,max(cophenetic.phylo(subtres[[x]])),subtres[[x]]$Ntip))
    
  }
  
  ##data formatting
  cluster_max_dist = as.data.frame(cluster_max_dist)
  colnames(cluster_max_dist) = c("subtree","max_dist","ntips")
  
  return(cluster_max_dist)
  
}

extract_dist_clusters = function(cluster_max_dist,dist_thres,subtres){
  ##filter to clusters with a within maximum distance of xx
  cluster_max_dist_for_filter = cluster_max_dist[cluster_max_dist$max_dist < dist_thres,]
  
  ##extract maximum cluster
  clusters = NULL
  for(i in 1:nrow(cluster_max_dist_for_filter)){
    
    #end loop when all clusters are checked
    if(nrow(cluster_max_dist_for_filter) == 0) {break}
    
    ##index of largest cluster (take first one if multiple with same size)
    largest_clus_index = which(cluster_max_dist_for_filter$ntips == max(cluster_max_dist_for_filter$ntips))[1]
    
    ##clustername of largest cluster (index of subtree)
    cluster_name_maxtips_temp = cluster_max_dist_for_filter$subtree[largest_clus_index]
    
    ##tip names in selected cluster
    tipsnames = subtres[[cluster_name_maxtips_temp]]$tip.label
    
    ##remove cluster from (no longer needed)
    cluster_max_dist_for_filter = cluster_max_dist_for_filter[-largest_clus_index,]
    
    ##tips already accounted for? 
    ##if yes, next start with next cluster
    ##possible can be more efficient with recursion to filter out all already accounted clusters directly
    ##but its fast enough I would say
    if(sum(tipsnames %in% clusters[,2]) > 0){next}
    
    ##bind cluster name and tipnames within cluster
    clusters = rbind(clusters,
                     cbind(cluster_name_maxtips_temp,tipsnames)
    )
    
  }
  
  clusters = as.data.frame(clusters) %>% group_by(cluster_name_maxtips_temp) %>% mutate(n = n())
  
  return(clusters)
}


xtt_data = function(region,subtype,seqtype,file_annote){

  lineage_tree = read.tree(paste0("....fa.treefile"))
  
  ##extract all subtrees
  subtres = subtrees(lineage_tree)
  
  tree_clusters = cluster_for_xtt(subtres)
  
  if(seqtype == "") {seqtype = "_wga"}
  
  for(dist_thres in as.numeric(c(seq(0.003,0.009,0.001),seq(0.01,0.09,0.01),seq(0.1,0.5,0.1)))){
    clusters_extracted = extract_dist_clusters(tree_clusters,dist_thres,subtres)
    
    if(seqtype == "_wga") {
      clusters_extracted = merge(clusters_extracted,shcs_ngs_docu[,c("SHCS.ID","base_uuid")], by.x = "tipsnames", by.y = "base_uuid", all = FALSE )
      clusters_extracted = clusters_extracted %>% select(cluster_name_maxtips_temp,tipsnames = SHCS.ID, n)
    }
    
    ##new phenos 
    pheno_covar_data = read.csv("eacs_pheno_covar.csv")
    pca = read.csv(paste0("path_to_pcs.evec"), sep = "", header = FALSE)[-1,]
    
    colnames(pca) = c("id","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","famid")
    pheno_covar_data = pheno_covar_data %>% group_by(id) %>% dplyr::slice(1) %>%
      merge(.,pca[,-12], by = "id", all.x = TRUE, all.y = FALSE) 
    
    pheno_covar_data = pheno_covar_data %>% mutate(drugs = ifelse(coca == "z_coca" | hero == "z_hero" | cana_ni == "z_cana","z_drugs","a_drugs")) %>%
      mutate(dis = ifelse(dis == "z_dis" | inf_dis == "z_inf_dis" | neuro_assoc == "z_efv_assoc","z_dis","a_dis")) %>%
      mutate(ethnicity = ifelse(ethnicity == "eth9" | ethnicity == "ethNA" |  ethnicity == "eth0","eth4",ethnicity)) %>%
      mutate(PC1 = PC1 * 50,PC2 = PC2 * 50,PC3 = PC3 * 50,PC4 = PC4 * 50,PC5 = PC5 * 50,PC6 = PC6 * 50,PC7 = PC7 * 50,PC8 = PC8 * 50,
             PC9 = PC9 * 50,PC10 = PC10 * 50) %>% 
      mutate(risk = case_when(risk == "risk0" ~ "risk3",
                              risk == "risk1" ~ "risk1",
                              risk == "risk2" ~ "risk2",
                              risk == "risk3" ~ "risk3",
                              risk == "risk4" ~ "risk3",
                              risk == "risk5" ~ "risk3",
                              risk == "risk6" ~ "risk3",
                              risk == "risk7" ~ "risk3",
                              risk == "risk9" ~ "risk3",
                              TRUE ~ "risk3"))
    
    clusters_dists = clusters_extracted %>% select(id = tipsnames , cluster = cluster_name_maxtips_temp)
    
    pheno_covar_data = pheno_covar_data %>% 
      merge(.,clusters_dists, by = "id", all = FALSE) %>%
      group_by(cluster) %>% mutate(n = n()) %>% filter(n > 1)
    
    pheno_covar_data = pheno_covar_data %>% 
      select(id,cluster,auc_slow_avg,auc_conc_avg,auc_all_avg,auc_freq_avg,sex,age,education,auc_rna,auc_cd4,ethnicity,risk,AUC_efv,antidepress,depri,drugs,dis,hep_c,hep_b,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)
    
    write.csv(pheno_covar_data,paste0("file_for_stata_xtt_analysis.csv"))
  }
}