##Required functions
source("protein_sequence_functions.R")

##Sequence selection
##Specify path to fasta files
consensus_path <- "path_to_nucleotide_consensus_fasta-files"
##specify path to MACSE (codon alignment)
macse_path <- "path_to_macse_v2.05.jar"

##Define list of uuids
uuids <- list.files(path = consensus_path)
uuids <- gsub("\\.fa","",uuids)
uuids <- as.data.frame(uuids)
colnames(uuids) <- "uuid"

##Extract regions via blast (reference subtype panel) (see folder references_seqs)
registerDoParallel(cores = 8L)
region_extract(uuids$uuid,"POL","partial_pol_check","hxx")
region_extract(uuids$uuid,"ENV","hyper_05","hxx")
region_extract(uuids$uuid,"GAG","hyper_05","hxx")
region_extract(uuids$uuid,"VIF","eacs","hxx")
region_extract(uuids$uuid,"NEF","eacs","hxx")
region_extract(uuids$uuid,"VPU","eacs","hxx")
region_extract(uuids$uuid,"VPR","eacs","hxx")
region_extract(uuids$uuid,"tat_exon1","eacs","hxx")
region_extract(uuids$uuid,"tat_exon2","eacs","hxx")
region_extract(uuids$uuid,"rev_exon1","eacs","hxx")
region_extract(uuids$uuid,"rev_exon2","eacs","hxx")

##Extract regions via blast (reference hxb2) (see folder references_seqs)
# region_extract(uuids$uuid,"POL","eacs")
# region_extract(uuids$uuid,"env","eacs")
# region_extract(uuids$uuid,"gag","eacs")
# region_extract(uuids$uuid,"nef","eacs")
# region_extract(uuids$uuid,"tat_exon1","eacs")
# region_extract(uuids$uuid,"tat_exon2","eacs")
# region_extract(uuids$uuid,"vpu","eacs")
# region_extract(uuids$uuid,"vpr","eacs")
# region_extract(uuids$uuid,"vif","eacs")
# region_extract(uuids$uuid,"rev_exon1","eacs")
# region_extract(uuids$uuid,"rev_exon1","eacs")

##Merge exons from rev and tat
exon_merge(uuids,"rev")
exon_merge(uuids,"tat")

##Codon alignment 
##Parallelized
for(j in c("pol","env","gag","nef","tat","vpu","vpr","vif","rev")){
  uuids = read.csv(paste0("./blastmeta/blastmeta_",j,"_partial_pol_check.csv")) %>% filter(!is.na(Alignment.Length)) %>% select(QueryID)
  colnames(uuids) = "uuid"
  registerDoParallel(cores = 8L)
  foreach(x = unique(uuids$uuid)) %dopar% {
    
    condon_align(x,j)
    
  }
  
}
