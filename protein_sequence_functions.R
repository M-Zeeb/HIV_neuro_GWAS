##Libraries
library(devtools)
library(tidyverse)
library(foreach)
library(doParallel)

#devtools::install_github("mhahsler/rBLAST", force = TRUE)

##Blast (requires local blast, available from NCBI)
region_extract = function(uuids,region,out,ref_data_base) {
  
  if(!region %in% c("gag","GAG","env","ENV","POL","pol","rev_exon1","rev_exon2","tat_exon1","tat_exon2","VIF","vif","VPR","vpr","VPU","vpu","NEF","nef")) {
    stop("Region doesnt exist or is not implemented yet (only POL, GAG, tat exons, NEF, VIF, VPU, VPR, rev exons, and ENV)")
  }
  
  region <- tolower(region)
  ##BLAST database
  rBLAST::makeblastdb(paste0("./references_seqs/",ref_data_base,"_",region,".fa"),dbtype = "nucl")
  bl <- rBLAST::blast(paste0("./references_seqs/",ref_data_base,"_",region,".fa"), type = "blastn")
  blastmeta <- data.frame(QueryID = character(),
                          SubjectID = character(),  
                          Perc.Ident = numeric(),  
                          Alignment.Length = integer(),
                          Mismatches = integer(),   
                          Gap.Openings = integer(),  
                          Q.start  = integer(),      
                          Q.end   = integer(),        
                          S.start  = integer(),        
                          S.end   = integer(),         
                          E = numeric(),             
                          Bits = numeric(), 
                          stringsAsFactors = FALSE)
  y <- 1
  blastmeta <- foreach(i = uuids, .combine = "rbind") %dopar%{
    ##sequences
    
    seq <- Biostrings::readDNAStringSet(paste0(consensus_path,i,".fa"))
    seq_df <- as.data.frame(Biostrings::readDNAStringSet(paste0(consensus_path,i,".fa")))
    
    ##BLAST
    
    blast_seq <- tryCatch({predict(bl, seq, BLAST_args = "-max_target_seqs 2000 -evalue 4000")}, 
                          error = function(e) { 
                            print("blast error")
                          })
    
    
    if(nrow(blast_seq) == 0){
      
      blast_seq [nrow(blast_seq) + 1, ] <- NA
      blast_seq[1,1] <- i
      blastmeta <- rbind(blastmeta,blast_seq)
      #next 
      return(NA)
    }   
    
    ##seems to be the most practical to sort for length and then for E score
    blast_seq <- blast_seq %>% arrange(desc(Alignment.Length),E) %>% slice(1)
    
    ##MERGE BLAST with raw sequence
    seq_df %>%
      mutate(reg_seq <- substr(x, blast_seq$Q.start, blast_seq$Q.end))  %>%
      mutate(forcsv <- paste0(">",rownames(seq_df)[1],"\n",reg_seq)) %>% {
        write.table(.$forcsv,paste0("./blasted_seqs_",region,"_partial_check/",i,"_blasted.fa"), row.names = FALSE, quote = FALSE,col.names = FALSE)
      }
    return(blast_seq)
  }
  write.csv(blastmeta,paste0("./blastmeta/blastmeta_",region,"_",out,".csv"))
}

##Merge exons
exon_merge <- function(uuids,region){
  
  if(!region %in% c("rev","tat","TAT","REV")) {
    stop("only REV or TAT")
  }
  
  blast_ex1 <- read.csv(paste0("./blastmeta/blastmeta_",tolower(region),"_exon1_eacs.csv")) %>% filter(!is.na(SubjectID))
  blast_ex2 <- read.csv(paste0("./blastmeta/blastmeta_",tolower(region),"_exon2_eacs.csv")) %>% filter(!is.na(SubjectID))
  uuids <- uuids %>% mutate(ex1 = uuid %in% blast_ex1$QueryID) %>%
    mutate(ex2 = uuid %in% blast_ex2$QueryID) 
  
  for(i in uuids$uuid[uuids$ex1 == TRUE & uuids$ex2 == TRUE]){
    prot_comb <- xscat(readDNAStringSet(paste0("./blasted_seqs_",tolower(region),"_exon1/",i,"_blasted.fa")),
                       readDNAStringSet(paste0("./blasted_seqs_",tolower(region),"_exon2/",i,"_blasted.fa")))
    prot_comb@ranges@NAMES <- i
    Biostrings::writeXStringSet(prot_comb,paste0("./blasted_seqs_",tolower(region),"/",i,"_blasted.fa"))                     
  }
  
  for(i in uuids$uuid[uuids$ex1 == TRUE & uuids$ex2 == FALSE]){
    prot_comb <- xscat(readDNAStringSet(paste0("./blasted_seqs_",tolower(region),"_exon1/",i,"_blasted.fa")))
    prot_comb@ranges@NAMES <- i
    Biostrings::writeXStringSet(prot_comb,paste0("./blasted_seqs_",tolower(region),"/",i,"_blasted.fa"))                        
  }
  
  for(i in uuids$uuid[uuids$ex1 == FALSE & uuids$ex2 == TRUE]){
    prot_comb <- xscat(readDNAStringSet(paste0("./blasted_seqs_",tolower(region),"_exon2/",i,"_blasted.fa")))
    prot_comb@ranges@NAMES <- i
    Biostrings::writeXStringSet(prot_comb,paste0("./blasted_seqs_",tolower(region),"/",i,"_blasted.fa"))                       
  }
}

##Condon alignment with MACSE
condon_align <- function(uuids,region){
  
  for(i in uuids){
    system(paste0("java -jar ",macse_path," -prog alignSequences -seq ./references_seqs/hxb2_",region,".fa -seq_lr ",
                  "./blasted_seqs_",region,"/",i,"_blasted.fa -fs_lr 10 -stop_lr 15 -fs 100 -fs_term 100 ",
                  " -out_AA ./",region,"_codon_align/",i,"_",region,"_AA.fa",
                  " -out_NT ./",region,"_codon_align/",i,"_",region,"_NT.fa"))
  }
  
}
