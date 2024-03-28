##Libraries
library(tidyverse)
library(svMisc)
library(foreach)
library(doParallel)
library(Rcpp)
library(qdapTools)
library(glmnet)
library(fastDummies)

##NGS meta data  
shcs_ngs_docu = read.csv("NGS_samples_list_all\ runs_NGS_R1R2.csv")

#Load the matching table between SHCS and ZPHI
matching = readxl::read_excel("info.xlsx", col_names = TRUE)

##phenotype and covariables for nullmodel to generate residuals phenotype
pheno_covar_data = read.csv("eacs_pheno_covar.csv")
pheno_covar_data = pheno_covar_data %>% group_by(id) %>% dplyr::slice(1) 

##subtype
shcs_sub = read.csv("subtypes.txt")

#Match the ZPHI ID with an SHCS ID
for(i in 1:length(shcs_ngs_docu$SHCS.ID)){
  for(y in 1:length(matching$patnr)){
    if(!is.na(matching$cohortnr[y]) && !is.na(shcs_ngs_docu$ZPHI.ID[i]) && matching$patnr[y] == shcs_ngs_docu$ZPHI.ID[i]){
      shcs_ngs_docu$SHCS.ID[i]=matching$cohortnr[y]
    }
  }
}

##path to plink
path_to_plink = "/plink"

rename_uuids = function(msa,region,shcs_ngs_docu){
  hxb2_region = paste0("HXB2_",toupper(region))
  msa = msa[msa@ranges@NAMES %in% shcs_ngs_docu$base_uuid | msa@ranges@NAMES == hxb2_region]
  
  ##THIS REQUIRES HXB2 REF AS THE FIRST POSITION...CHANGE IT IDIOT!!!!!!
  
  seq_names = as.data.frame(msa@ranges@NAMES)
  colnames(seq_names) = "uuid"
  seq_names = seq_names %>% merge(.,shcs_ngs_docu[,c("SHCS.ID","base_uuid","sample.date")], by.x = "uuid", by.y = "base_uuid", all.x = TRUE, all.y = FALSE)
  seq_names$SHCS.ID[seq_names$uuid == hxb2_region] = hxb2_region
  seq_names <- seq_names[order(match(seq_names$uuid, msa@ranges@NAMES)),]
  msa@ranges@NAMES[msa@ranges@NAMES != hxb2_region] = seq_names$SHCS.ID[seq_names$SHCS.ID != hxb2_region]
  msa = msa[!duplicated(msa@ranges@NAMES)]
  
  return(msa)
}

##reformat sequence data

reformat_seq_to_csv <- function(msa){
  jk <- as.data.frame(msa)
  jk$V1 <- msa@ranges@NAMES
  ##Rename column names
  names(jk) <- c("seq","id")
  jk <- jk %>% select(id,seq)
  ##Dataset for sequencing data
  maxfa <- jk %>% separate(seq, into = paste("V", 1:(max(nchar(jk$seq))+1), sep = ""), sep = "") %>% select(-V1)
  return(maxfa)
}

read_alignment_and_format <- function(region,shcs_ngs_docu,file_annote,res,sequencingtype,subtype){
  ##load alignment
  if(sequencingtype == "_wga") {
    msa <- Biostrings::readDNAStringSet(paste0("alignment/",region,"/nci_NT",subtype,"_",region,"_aligned_",file_annote,".fa"))
  }
  if(sequencingtype == "_resdb") {
    msa <- Biostrings::readDNAStringSet(paste0("alignment/",region,sequencingtype,"/nci_NT",subtype,"_",region,sequencingtype,"_aligned_",file_annote,".fa"))
  }
  
  ##renaming from uuid to shcsid (possibly subject to change)
  if((sequencingtype != "_resdb")){
    msa <- rename_uuids(msa,region,shcs_ngs_docu)
  }
  ##reformat dnastring-format to dataframe
  msa <- reformat_seq_to_csv(msa)
  write.csv(msa,paste0("msa_dataframe_",region,subtype,res,sequencingtype,".csv"))
  
  return(msa)
}

##store baseposition in relation to reference for later seqs to 
ref_seq_pos <- function(msa,region,res,sequencingtype,subtype){
  hxb2_region <- paste0("HXB2_",toupper(region))
  
  seqs_to_keep_neuro <- as.matrix(msa %>%
                                   filter(id == hxb2_region))
  
  seqs_to_keep_neuro  <- reshape2::melt(seqs_to_keep_neuro)
  
  seqs_to_keep_neuro <- seqs_to_keep_neuro  %>%
    dplyr::select(-Var1) %>%
    filter(value != hxb2_region) %>%
    filter(value != "-") %>%
    mutate(hxb2_index_pol = as.numeric(rownames(.))) %>%
    mutate(hxb2_index_overall = as.numeric(hxb2_index_pol) )
  
  seqs_to_keep_neuro$gwas_pos <- as.numeric(gsub("V", "", seqs_to_keep_neuro$Var2)) 
  seqs_to_keep_neuro$msa_pos <- seqs_to_keep_neuro$gwas_pos 
  write.csv(seqs_to_keep_neuro,paste0("seqs_to_keep_",region,subtype,res,sequencingtype,".csv"))
  return(seqs_to_keep_neuro)
}

##resitances definitions for pol only
resistances_pos <- function(sequencingtype,region){
  ##resistance positiion
  resistances <- as.data.frame(c(89,90,91,95,96,97,98,99,100,137,138,139,140,141,142,143,144,145,149,
                                150,151,161,162,163,173,174,175,221,222,223,227,228,229,245,246,247,
                                248,249,250,251,252,253,263,264,265,269,270,271,419,420,421,482,483,
                                484,491,492,493,497,498,499,506,507,508,518,519,520,521,522,523,527,
                                528,529,596,597,598,599,600,601,605,606,607,614,615,616,620,621,622,
                                641,642,643,644,645,646,710,711,712,749,750,751,833,834,835,839,840,
                                841,848,849,850,860,861,862,866,867,868,926,927,928,941,942,943,953,
                                954,955,959,960,961,971,972,973,977,978,979,986,987,988))
  colnames(resistances) <- "respos"
  if((sequencingtype == "_wga") & (region != "whole")){
    resistances <- resistances %>% mutate(respos3 = respos * 3)  %>% mutate(respos2 = respos3 - 1) %>% mutate(respos1 = respos3 - 2)
  }
  if(sequencingtype == "_resdb"){
    resistances <- resistances %>%  mutate(respos3 = (respos * 3))  %>% mutate(respos2 = respos3 - 1) %>% mutate(respos1 = respos3 - 2)
  }
  if((sequencingtype == "_wga") & (region == "whole")){
    resistances <- resistances %>%  mutate(respos3 = (respos * 3) + 2083)  %>% mutate(respos2 = respos3 - 1) %>% mutate(respos1 = respos3 - 2)
  }
  
  return(resistances)
}

##plinkformatting
plink_formatting_vcf_alleles <- function(msa,n_snps,n_samp){
  
  msa <- msa %>% mutate(across(2:(n_snps), ~ replace(., which(is.na(.) | !. %in% c("A","C","G","T")),"N")))
  
  vcf_allels <- NULL
  vcf_allels <- as.data.frame(matrix(ncol = 1, nrow = n_snps-1))
  for(i in 2:n_snps){
    vcf_allels$V1[i-1] = i-1
    vcf_allels$A[i-1] = sum(msa[,i] == "A",na.rm = TRUE)
    vcf_allels$C[i-1] = sum(msa[,i] == "C",na.rm = TRUE)
    vcf_allels$G[i-1] = sum(msa[,i] == "G",na.rm = TRUE)
    vcf_allels$T[i-1] = sum(msa[,i] == "T",na.rm = TRUE)
    vcf_allels$N[i-1] = sum(is.na(msa[,i]) | (!msa[,i] %in% c("A","G","C","T")))
  }
  
  return(vcf_allels)
}

generate_vcf_empty <- function(vcf_allels,maf,n_snps,n_samp){
  ##vcf file
  gwas_vcf <- NULL
  gwas_vcf <- as.data.frame(matrix(ncol = (n_samp+9), nrow = (n_snps-1)))
  for(z in 1:(n_snps-1)){
    gwas_vcf[z,1] = 1
    gwas_vcf[z,2] = z 
    gwas_vcf[z,3] = "." 
    gwas_vcf[z,4] = colnames(vcf_allels)[which(vcf_allels[z,2:5] == max(vcf_allels[z,2:5]))[1]+1]
    gwas_vcf[z,5] = paste0(unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][1]),
                           ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][2]),
                           ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][3]),
                           ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][4]),
                           ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][5]))
    gwas_vcf[z,5] = gsub(pattern = paste0("NA,|NA|,NA"),"",gwas_vcf[z,5])
    gwas_vcf[z,5] = gsub(pattern = paste0("\\",gwas_vcf[z,4],","),"",gwas_vcf[z,5])
    gwas_vcf[z,5] = gsub(pattern = paste0(",\\",gwas_vcf[z,4]),"",gwas_vcf[z,5])
    gwas_vcf[z,5] = gsub(pattern = paste0("\\",gwas_vcf[z,4]),"",gwas_vcf[z,5])
    gwas_vcf[z,5][gwas_vcf[z,5] == ""] = "."
    gwas_vcf[z,6] = "."
    gwas_vcf[z,7] = "."
    gwas_vcf[z,8] = "."
    gwas_vcf[z,9] = "GT"
    
  }  
  colnames(gwas_vcf)[1:9] <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  
  return(gwas_vcf)
}

generate_vcf_fill <- function(msa,gwas_vcf,n_snps,n_samp){
  for(z in 1:nrow(msa)){
    colnames(gwas_vcf)[z+9] <- as.character(msa$id[z])
    chrs <- strsplit(as.character(paste0(gwas_vcf[1:(n_snps-1),4],",",gwas_vcf[1:(n_snps-1),5],",",msa[z,2:n_snps])), split = ",")
    gwas_vcf[,z+9] <- paste0(unlist(lapply(chrs,anyDuplicated, fromLast = TRUE))-1)
    
    progress(z,progress.bar = FALSE, init = (z == 1), max.value = n_samp)
  }
  
  gwas_vcf[,][gwas_vcf[,] == "-1"] <- "." 
  
  return(gwas_vcf) 
}

binarising <- function(gwas_vcf,altc,n_snps,n_samp){
  othercs <- c(1,2,3,4)[c(1,2,3,4) != altc] 
  gwas_vcf_alt2 <- gwas_vcf
  nseq <- n_samp + 9 
  gwas_vcf_alt2$ID <- NA
  gwas_vcf_alt2$ID[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc)] <- paste0(gwas_vcf_alt2$POS[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc)],gwas_vcf_alt2$REF[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc)],altc,sapply(strsplit(gwas_vcf_alt2[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc),5],","),"[[", altc))  
  gwas_vcf_alt2 <- gwas_vcf_alt2[which(!is.na(gwas_vcf_alt2$ID)),]
  
  ##to keep non-missing with other base as 0 and complete missing as .
  truemiss <- str_locate(pattern = "N", gsub(",","",gwas_vcf_alt2[,"ALT"]))[,1]
  truemiss[is.na(truemiss)] <- 99
  miss <- gwas_vcf_alt2[,10:nseq] == truemiss
  
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[1] | miss)] = "."
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[1] | (!miss))] = "0" 
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[2] | miss)] = "."
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[2] | (!miss))] ="0"
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[3] | miss)] = "."
  gwas_vcf_alt2[,10:nseq][(gwas_vcf_alt2[,10:nseq] != altc) & (gwas_vcf_alt2[,10:nseq] == othercs[3] | (!miss))] ="0"
  gwas_vcf_alt2[,10:nseq][gwas_vcf_alt2[,10:nseq] == altc] = 1
  gwas_vcf_alt2$ALT = sapply(strsplit(gwas_vcf_alt2[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc),5],","),"[[", altc)
  return(gwas_vcf_alt2)
}


final_vcf <- function(gwas_vcf_alt,seqs_to_keep_neuro_todel,region,res,sequencingtype,subtype){
  gwas_vcf_alt$`#CHROM` <- 1
  ##remove resistance position
  if((res != "") & ((region == "pol") | (region == "whole"))){
    gwas_vcf_alt <- gwas_vcf_alt %>% filter(POS %in% seqs_to_keep_neuro_todel$msa_pos)
  }
  gwas_vcf_alt <- gwas_vcf_alt %>% filter(ALT != "N" & ALT != ".")
  gwas_vcf_alt$FILTER <- "PASS"
  gwas_vcf_alt$POS <- c(1:nrow(gwas_vcf_alt))
  write.table(gwas_vcf_alt,paste0("plink2vcf_",region,subtype,res,sequencingtype,".vcf"), sep ="\t", quote = FALSE, row.names = FALSE)
  
  header <- paste0("##fileformat=VCFv4.2 
##fileDate=20220614
##source=R4.00
##contig=<ID=1,length=15462>
##INFO=<ID=PR,Number=0,Type=Flag,Description=\"",region,subtype,res,sequencingtype,"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
  
  write.table(header,paste0("vcfheader_",region,subtype,res,sequencingtype,".txt"),sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  ##merge vcf with fileheader
  system(paste0("cat vcfheader_",region,subtype,res,sequencingtype,".txt plink2vcf_",region,subtype,res,sequencingtype,".vcf > final_vcf_",region,subtype,res,sequencingtype,".vcf"))
  
}

final_vcf_x <- function(gwas_vcf_alt,seqs_to_keep_neuro_todel,region,res,sequencingtype,subtype){
  gwas_vcf_alt$`#CHROM` <- 1
  ##remove resistance position
  if((res != "") & ((region == "pol") | (region == "whole"))){
    gwas_vcf_alt = gwas_vcf_alt %>% filter(POS %in% seqs_to_keep_neuro_todel$msa_pos)
  }
  gwas_vcf_alt <- gwas_vcf_alt %>% filter(ALT != "N" & ALT != ".")
  gwas_vcf_alt$FILTER <- "PASS"
  gwas_vcf_alt$POS <- c(1:nrow(gwas_vcf_alt))
  gwas_vcf_alt[,10:ncol(gwas_vcf_alt)][gwas_vcf_alt[,10:ncol(gwas_vcf_alt)] == 1] <- 2
  write.table(gwas_vcf_alt,paste0("plink2vcf_",region,subtype,res,sequencingtype,"x.vcf"), sep ="\t", quote = FALSE, row.names = FALSE)
  
  header = paste0("##fileformat=VCFv4.2 
##fileDate=20220614
##source=R4.00
##contig=<ID=1,length=15462>
##INFO=<ID=PR,Number=0,Type=Flag,Description=\"",region,subtype,res,sequencingtype,"\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
  
  write.table(header,paste0("vcfheader_",region,subtype,res,sequencingtype,"x.txt"),sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  ##merge vcf with fileheader
  system(paste0("cat vcfheader_",region,subtype,res,sequencingtype,"x.txt plink2vcf_",region,subtype,res,sequencingtype,"x.vcf > final_vcf_",region,subtype,res,sequencingtype,"x.vcf"))
  
}


plink2_file_generation <- function(path_to_plink,region,res,sequencingtype,subtype){
  ##GWAS
  system(paste0(path_to_plink,"/plink2 --vcf final_vcf_",region,subtype,res,sequencingtype,".vcf --keep-allele-order --allow-extra-chr --chr-set -1 --make-pgen --sort-vars --out ",path_to_plink,"/temp/gwasp2"))
  system(paste0(path_to_plink,"/plink2 --pfile ",path_to_plink,"/temp/gwasp2 --pheno pheno_type_all.txt --nonfounders  --allow-extra-chr  --make-pgen --out ",path_to_plink,"/temp/hivp2_pheno1"))
  ##remove with missing phenotypes
  system(paste0(path_to_plink,"/plink2 --pfile ",path_to_plink,"/temp/hivp2_pheno1 --require-pheno auc_slow_avg --allow-extra-chr  --make-pgen --out ",path_to_plink,"/temp/hivp2_pheno2"))
  
  ##linkage correlation matrix hopefully
  system(paste0(path_to_plink,"/plink2 --pfile ",path_to_plink,"/temp/hivp2_pheno2 --geno 0.2 --nonfounders --make-bed --out ",region,subtype,res,sequencingtype,"_forcor"))
  
  ##linkage
  system(paste0(path_to_plink,"/plink2 --pfile ",path_to_plink,"/temp/hivp2_pheno2 --geno 0.2 --nonfounders --indep-pairwise 50 10 0.5 --out gcta_",region,subtype,res,sequencingtype,"_exclude"))
  ##save for gcta
  system(paste0(path_to_plink,"/plink2 --pfile ",path_to_plink,"/temp/hivp2_pheno2 --geno 0.2 --nonfounders  --exclude gcta_",region,subtype,res,sequencingtype,"_exclude.prune.out --make-bed --out gcta_",region,subtype,res,sequencingtype))
  ##further for eigensoft
  system(paste0(path_to_plink,"/plink2 --pfile ",path_to_plink,"/temp/hivp2_pheno2 --mind 0.9 --nonfounders --allow-extra-chr  --make-pgen --out ",path_to_plink,"/temp/hivp2_mind3"))
  system(paste0(path_to_plink,"/plink2 --pfile ",path_to_plink,"/temp/hivp2_mind3 --geno  0.99  --nonfounders --allow-extra-chr  --make-pgen --out ",path_to_plink,"/temp/hivp2_mind3"))
  
  ##make bed for eigensoft
  system(paste0(path_to_plink,"/plink2 --pfile ",path_to_plink,"/temp/hivp2_mind3  --nonfounders --allow-extra-chr  --export ped --out ",region,subtype,res,sequencingtype))
  
  ##plink cleanup
  system(paste0("rm -rfv ",path_to_plink,"/temp/*"))
  
}

##DUmmy coding for sex, required for PLINK
sex_update <- function(iddata,region,res,sequencingtype,subtype){
  sex_data <- as.data.frame(colnames(iddata)[10:ncol(iddata)])
  sex_data <- sex_data %>% mutate(V1 = 0, V3 = 1) %>% select(c(2,1,3))
  write.table(sex_data,paste0("gcta_",region,subtype,res,sequencingtype,"_sex_update"), sep= " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

##create parfile for convertf (ped to eigenstrat format)
create_parfile <- function(region,res,sequencingtype,subtype){
  
  ##pheno type in column six as 1 (errors otherwise with "0" or "-9" values)
  set_pheno_to_one <- read.csv(paste0(region,subtype,res,sequencingtype,".ped"), sep = "\t", header = FALSE)
  set_pheno_to_one$V6 <- 1
  write.table(set_pheno_to_one,paste0(region,subtype,res,sequencingtype,".ped"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  partext <- paste0("genotypename:    ",region,subtype,res,sequencingtype,".ped
snpname:         ",region,subtype,res,sequencingtype,".map 
indivname:       ",region,subtype,res,sequencingtype,".ped
outputformat:    EIGENSTRAT
genotypeoutname: ",region,subtype,res,sequencingtype,".eigenstratgeno
snpoutname:      ",region,subtype,res,sequencingtype,".snp
indivoutname:    ",region,subtype,res,sequencingtype,".ind
familynames:     NO
         ")
  
  write.table(partext,paste0("par.ped.eigenstrat"),sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}

##convert ped to eigenstrat format
convertf <- function(){
  
  system(paste0("/opt/anaconda3/bin/conda run -n eigensoft convertf -p par.ped.eigenstrat"))
  
}

eigensoft_uniform_population <- function(region,res,sequencingtype,subtype){
  
  ##Uniform population for eigensoft
  pop_remove_eigensoft <- read.csv(paste0(region,subtype,res,sequencingtype,".ind"), sep= "", header = FALSE)
  pop_remove_eigensoft$V3 <- 1
  write.table(pop_remove_eigensoft,paste0(region,subtype,res,sequencingtype,".ind"), sep= " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}

##create parfile for eigensoft pca
smartpca_parfile <- function(region,res,sequencingtype,subtype){
  
  partext <- paste0("genotypename:    ",region,subtype,res,sequencingtype,".eigenstratgeno
snpname:         ",region,subtype,res,sequencingtype,".snp
indivname:       ",region,subtype,res,sequencingtype,".ind
evecoutname:     ",region,subtype,res,sequencingtype,".evec
evaloutname:     ",region,subtype,res,sequencingtype,".eval
altnormstyle:    NO 
numoutevec:      10
familynames:     NO
grmoutname:      grmjunk")
  
  write.table(partext,paste0("parfile"),sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}

smartpca <- function(){
  
  system(paste0("/opt/anaconda3/bin/conda run -n eigensoft smartpca -p parfile > logfile_test_check"))

  }

##wrapper function
main_function <- function(genreg,res = "",sequencingtype,subtype,file_annote){
  
  print(paste0(genreg," start"))
  
  ##set working directory for gene region
  setwd(paste0("/Users/mariuszeeb/Documents/neuro/eigensoft/",genreg,res,sequencingtype,"/"))
  
  ##read and format
  msa <- read_alignment_and_format(genreg,shcs_ngs_docu,file_annote,res,sequencingtype,subtype)
  ##filter to subtype B via covar (quick and dirty here needs to be dynamics !!!TBD!!!)
  msa <- msa %>% filter((id %in% pheno_covar_data$id) | (id == paste0("HXB2_",toupper(genreg))))
  ##store ref seq pos
  seqs_to_keep <- ref_seq_pos(msa,genreg,res,sequencingtype,subtype)
  ##resistance definitions for pol
  resistances <- resistances_pos(sequencingtype,genreg)
  
  ##when region is pol
  seqs_to_keep_neuro_todel <- NULL
  if((genreg == "pol") | (genreg == "whole")){
    seqs_to_keep_neuro_todel <- seqs_to_keep %>% filter(!(hxb2_index_pol %in% resistances$respos1 | hxb2_index_pol %in% resistances$respos2 | hxb2_index_pol %in% resistances$respos3))
  }
  
  ##number of snps
  n_snps <- length(msa)
  #number of samples
  n_samp <- nrow(msa)
  
  ##ambguuious nucleotides as "N"
  msa <- msa %>% mutate(across(2:(n_snps), ~ replace(., which(is.na(.) | !. %in% c("A","C","G","T")),"N")))
  
  ##make table with nucleotides at each position across sequence samples 
  vcf_alleles <- plink_formatting_vcf_alleles(msa,n_snps = n_snps,n_samp = n_samp)
  
  ##create empty vcf file
  ##maf is minor allele frequency threshold
  gwas_vcf <- generate_vcf_empty(vcf_alleles,maf=5,n_snps = n_snps,n_samp = n_samp)
  
  ##file empty vcf file
  gwas_vcf <- generate_vcf_fill(msa,gwas_vcf,n_snps = n_snps,n_samp = n_samp)
  
  ##binarise vcf file
  ##altc for different possible minor alleles
  gwas_vcf_bin <- rbind(binarising(gwas_vcf = gwas_vcf, altc = 1,n_snps = n_snps,n_samp = n_samp),
                       binarising(gwas_vcf = gwas_vcf, altc = 2,n_snps = n_snps,n_samp = n_samp),
                       binarising(gwas_vcf = gwas_vcf, altc = 3,n_snps = n_snps,n_samp = n_samp),
                       binarising(gwas_vcf = gwas_vcf, altc = 4,n_snps = n_snps,n_samp = n_samp))   
  
  ##generate final vcf file for processing with plink
  final_vcf(gwas_vcf_bin,seqs_to_keep_neuro_todel,genreg,res,sequencingtype,subtype)
  final_vcf_x(gwas_vcf_bin,seqs_to_keep_neuro_todel,genreg,res,sequencingtype,subtype)
  
  ##sex update
  sex_update(gwas_vcf_bin,genreg,res,sequencingtype,subtype)
  
  ##formatt for eigensoft with plink2
  plink2_file_generation(path_to_plink,genreg,res,sequencingtype,subtype)
  
  print(paste0(genreg," plink done"))
  
  ##general note
  ##parfile (parameter file I assume?) is requireed for eigensoft-software.
  ##Necessary file names are stored there
  
  ##create parfile for convertf
  create_parfile(genreg,res,sequencingtype,subtype)
  
  ##convert to eigenstrat format
  convertf()
  
  ##uniform population for eigensoft (errors otherwise)
  eigensoft_uniform_population(genreg,res,sequencingtype,subtype)
  
  ##create parfile for pca 
  smartpca_parfile(genreg,res,sequencingtype,subtype)
  
  ##run smartpca
  smartpca()
  
  ##for keeping track within loop
  print(paste0(genreg," finished"))
  
}

