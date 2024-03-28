This repository includes the core functions for the manuscript titled "Self-reported neurocognitive complaints in the Swiss HIV Cohort Study - a viral genome wide association study" published in Brain Communications.

The underlying data is not publicly available, as explained by the data sharing statement of the Swiss HIV Cohort Study: 

The individual level datasets generated or analyzed during the current study do not fulfill the requirements for open data access:

1) The SHCS informed consent states that sharing data outside the SHCS network is only permitted for specific studies on HIV infection and its complications, and to researchers who have signed an agreement detailing the use of the data and biological samples; 

and

2) the data is too dense and comprehensive to preserve patient privacy in persons living with HIV.

According to the Swiss law, data cannot be shared if data subjects have not agreed or data is too sensitive to share. 
Investigators with a request for selected data should send a proposal to the respective SHCS address (www.shcs.ch/contact). 
The provision of data will be considered by the Scientific Board of the SHCS and the study team and is subject to Swiss legal and ethical regulations, and is outlined in a material and data transfer agreement.




Workflow:

1. Gene sequence generation
	Script: protein_sequence_generation.R, protein_sequence_functions.R
	Requirements: rBLAST from NCBI (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)

2. Protein sequence generation
	Script: protein_sequence_generation.R, protein_sequence_functions.R
	Requirements: MACSE (https://www.agap-ge2pop.org/macse/)

3. Co-receptor tropism prediction
	Geno2Pheno using nucleotide sequences (https://coreceptor.geno2pheno.org/)

4. HIV-1 Subtype classification
	Requirements: rega (https://www.genomedetective.com/app/typingtool/hiv), comet (https://comet.lih.lu/), geno2pheno (https://coreceptor.geno2pheno.org/)
		We determined the HIV-1 subtype for each individual from near whole genome sequences using comet. 
		In case of multiple samples assigned to different subtypes, we used the subtype with highest certainty. 
		If near whole genome was unavailable or the subtype unassigned by comet, we determined it from partial 
		pol sequences using rega. If both near whole genome and partial pol sequences had unassigned subtypes, 
		we used the subtype from geno2pheno. 
		If no subtype from any method was assigned, we classified the individual as subtype “others”. 
	
5. Phenotypes
	Script: phenotypes.R

6. Multiple sequence alignment (with MAFFT)
	for each gene on nucleotide and protein level:
	mafft --localpair --maxiterate 1000  --thread 64 --add sequences.fa  hxb2_sequence.fa >  alignment.fa
	Requirements: MAFFT (https://mafft.cbrc.jp/alignment/software/)

7. Genomic population structure, Principal components
	Scripts: PCAs.R, PCAs_functions.R
	Requirements: PLINK2 (https://www.cog-genomics.org/plink/2.0/), eigensoft (https://reich.hms.harvard.edu/software)

8. Phylogenetic trees and cluster extraction
	Requirements: IQ-TREE (http://www.iqtree.org/)
	IQ-TREE command: for all combinations of Hiv-1 genes, HIV-1 subtype
		iqtree2 -s alignment.fa -alrt 1000 -bb 1000 -T AUTO

9. Phylogenetic cluster extraction and Heritability
	Script: phylogenetic_clusters_functions.R
	Requirements: stata	
	Stata command: for all combinations of Hiv-1 genes, HIV-1 subtype, and phylogenetic distance thresholds

		import delimited  gene_x_subtype_j_distance_threshold_y.csv
		xtset(cluster)
		encode sex, gen(nSex)
		encode depri , gen(nDepri)
		encode ethnicity, gen(nEthnicity)
		encode hep_c, gen(nhep_c)
		encode hep_b, gen(nhep_b)
		encode education, gen(neducation)
		encode risk, gen(nrisk)
		encode antidepress, gen(nantidepress)
		encode drugs, gen(ndrugs)
		encode dis, gen(ndis)
		xttobit auc_freq_avg i.nSex i.nrisk i.ndis i.ndrugs i.nantidepress i.nhep_c i.nhep_b i.neducation auc_efv age auc_rna auc_cd4 i.nEthnicity i.nDepri,ll(0) tobit
		xttobit auc_conc_avg i.nSex i.nrisk i.ndis i.ndrugs i.nantidepress i.nhep_c i.nhep_b i.neducation auc_efv age auc_rna auc_cd4 i.nEthnicity i.nDepri,ll(0) tobit
		xttobit auc_slow_avg i.nSex i.nrisk i.ndis i.ndrugs i.nantidepress i.nhep_c i.nhep_b i.neducation auc_efv age auc_rna auc_cd4 i.nEthnicity i.nDepri,ll(0) tobit
		xttobit auc_all_avg i.nSex i.nrisk i.ndis i.ndrugs i.nantidepress i.nhep_c i.nhep_b i.neducation auc_efv age auc_rna auc_cd4 i.nEthnicity i.nDepri,ll(0) tobit
		clear

10. GWAS
	Script: gwas.R, gwas_functions.R

11. GWAS crossvalidation
	Script: tobit_crossval_cluster.R, tobitnet_main_array.sh, tobitnet_array_wrapper.sh

    