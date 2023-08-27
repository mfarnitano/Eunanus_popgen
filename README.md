### Script and bioinformatics details for Eunanus popgen

This repository contains code used in the paper 'Strong postmating reproductive isolation in Mimulus section Eunanus', by Matthew Farnitano and Andrea Sweigart, accepted in the Journal of Evolutionary Biology 2023. Preprint available at https://doi.org/10.1101/2022.12.21.521469. 

I. Alignment of raw reads to reference genomes

Main script:
	fastq_align_wrapper.sh
Additional scripts:
	fastq_align_each.sh
Input files:
	genome reference fasta
	table listing samples and .fastq file locations (Brevipes_complete_fastq_table.txt)
	sample .fastq files
Outputs:
	filtered and processed .bam files for each sample
Description:
	Runs multistep pipeline to trim, align, filter, and de-duplicate reads.
	Wrapper scripts loops through all samples and calls an individual script for each sample.

II. SNP calling and filtering

Main scripts:
	genotype_wrapper.sh
	genotype_combine.sh
	filter_VCF.sh
Additional scripts:
	genotype_each.sh
Input files:
	genome reference fasta
	table listing samples and batch info (Brevipes_complete_fastq_table.txt)
	sample .bam files (from step I)
Outputs:
	filtered .vcf files with SNPs or SNPs+invariant sites
Description
	Uses GATK to filter .vcf files for a number of standard parameters such as total depth and number of called samples.

III. Subsetting VCFs to only four-fold degenerate sites

Main scripts:
	4D_sites.sh
Additional scripts:
	Identify_4D_Sites.pl
Input files:
	filtered .vcf file
	degenerate-sites codon table (codon_table.txt)
	genome annotation file (.gff)
Outputs:
	filtered .vcf file containing only four-fold degenerate sites
Description:
	Method and scripts adopted from Tim Sackton https://github.com/tsackton/linked-selection/blob/master/misc_scripts/Identify_4D_Sites.pl to subset VCF to only include sites that are four-fold degenerate, based on a genome annotation.

IV. Downsampling, aligment, and SNP calling

Main scripts:
	samplefastqs.sh 		##for pre-chosen group of samples
	samplefastqs_random.sh 	##for random group of samples
	genotype_wrapper.sh 	##change BATCH_NAME to match each downsampled batch
	genotype_combine.sh 	##change BATCH_NAME to match each downsampled batch
	filter_VCF_downsample.sh
Additional scripts:
	fastq_align_each.sh
	genotype_each.sh
Input files:
	genome reference fasta
	table listing samples and batch info (Brevipes_complete_fastq_table.txt, or subset if samples are pre-chosen)
	full .fastqs for each sample
Outputs:
	filtered .vcf files based on randomly downsampled reads from a pre-chosen or random subset of samples
Description:
	Randomly downsamples fastq files to the same number of reads, then aligns to the genome. Uses a subset of samples, either pre-chosen or randomly chosen from each species.

V. Phylogenetic trees using maximum likelihood and neighbor-joining

Main scripts:
	phylo.sh
Additional scripts:
	raxml_iter.sh
Input files:
	filtered .vcf file with variable sites only
Outputs:
	multisample .fasta containing variable sites only, with randomly selected allele at het sites (input for tree building)
	RAxML outputs: best tree from each iteration with bootstrap values, as well as bootstrap trees and likelihood info
Description:
	Prepares .fasta input from variable sites .vcf file, with Ns for missing data and randomly selected allele at het sites.
	Runs multiple iterations of RAxML using .fasta input, with ascertainment bias correction, and rapid bootstrapping.
	Must specify number of variable and invariable sites for ascertainment bias correction.
	Neighbor-joining tree can be produced using the same .fasta input, using R packages or another method.

VI. Gene tree analysis

Main scripts:
	trees_by_genes.sh
	twisst_bygenes.sh
	astral_by_genes.sh
Input files:
	filtered .vcf file for whole genome
	genome annotation file
	population grouping file for twisst (twisst_groups_main.txt)
Outputs:
	.fasta file of variable sites for each gene in the annotation
	maximum-likelihood tree for each gene in the annotation (if sufficient data)
	all trees with sufficient data combined into a single file
	TWISST outputs: quartet summaries for each tested quartet
	ASTRAL outputs: consensus tree with quartet or posterior-probability annotations
Description:
	Creates a .fasta file of variable sites for each gene separately from the genome annotation. Hets are replaced with a randomly selected allele.
	Uses RAxML to make a gene tree for each of these genes. Genes with insufficient data to run RAxML, or where RAxML collapses multiple samples due to a lack of informative sites, are removed from further analysis.
	Uses TWISST with specified species quartets to get summaries of the gene tree support for each possible quartet topology.
	Uses ASTRAL to get a consensus species tree based on gene trees, with quartet scores and posterior support values.

VII. Diversity and divergence statistics

Main scripts:
	pixy.sh
Input files:
	filtered .vcf with SNPs and invariant sites, complete or with only four-fold degenerate sites
	pixy populations file (pixy_popfile.txt)
	pixy individuals-as-populations file (pixy_popfile_inds.txt) for heterozygosity calcs
Outputs:
	estimates of pi, dxy, and fst between species, as well as pi for individuals (i.e. heterozygosity), in windows across the genome, with total counts so that genome-wide values can be easily calculated.
Description:
	Uses pixy to calculate pi, dxy, and fst between species for 10Kb and 1Mb window sizes across the genome. Also calculates heterozygosity by using the pi function with groups of 1.

VIII. Phylogenetic network models with TreeMix

Main scripts:
	TreeMix.sh
Additional scripts:
	vcf2treemix.sh
	plink2treemix.sh
	ldPruning.sh from https://github.com/joanam/scripts/raw/master/ldPruning.sh
Input files:
	filtered .vcf file
	species groupings file (treemix_popfile.clust)
Outputs:
	TreeMix networks for 10 iterations each of 1-8 migration edges, with likelihood values for each.
Description:
	Follows advice and scripts in https://speciationgenomics.github.io/Treemix/ to prepare TreeMix input files from VCF, including pruning to remove sites in high LD and removing sites with missing data.
	Runs TreeMix for 1 to 8 migration edges, 10 iterations each with randomly sampled SNPs for each iteration.

IX. Tests for introgression using D- and f-statistics

Main scripts:
	Dsuite.sh
	Dsuite_downsample.sh
Input files:
	filtered .vcf file
	Dsuite populations file (Dsuite_popfile_AURA.txt)
	Tree file showing population relationships (poplevel_treefile.tre)
Outputs:
	Dsuite outputs: D, f4, and related statistcs for each trio of species.
Description:
	Uses the Dsuite program to calculate D and f4 for each trio of species, using M. aurantiacus as the outgroup. Uses block-jack-knifing (implemented in Dsuite) with 100 blocks to estimate Z-score and p-value for each trio.
	Calculates the Fbranch statistics (implemented in Dsuite) based on f4 values for each trio.
	Run the same analyses using each of the downsampled datasets, using the Dsuite_downsample.sh script.

X. Window-based genome scan for introgression outliers

Main scripts:
	Dinvestigate.sh
Input files:
	filtered .vcf file
	Populations file for each trio of interest (Dinv_BCNA.txt, Dinv_CSFA.txt)
	Trio file for each trio of interest (Trio_BCN.txt, Trio_CSF.txt)
Outputs:
	D, d_f, and related statistics in 100-SNP windows (overlapping by 50 SNPs) across the genome, for each trio of interest.
Description:
	Uses the Dinvestigate function in Dsuite to calculate f-statistics in 100-SNP windows (overlapping by 50 SNPs) across the genome, for two tested trios with M. aurantiacus as the outgroup. 
