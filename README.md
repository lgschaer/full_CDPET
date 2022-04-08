# Full Analysis for Deconstructed Polyethylene Terephthalate and River Water Paper

# Analysis of Growth Data & Biodegradation
### Plotting growth curves and biodegradation bar chart

R Script: 04012022_DCPET_Growth_Stats.R

Input: DCPET_Experient_Data.csv

# Analysis of Chemistry Data
### Making table of estimated nutrient concentrations

R Script: 04042022_DCPET_Chemistry_Calculations.R

Input: DCPET_AQUA_Results_Formatted_R.xlsx

Output: seqtab.rds and taxa.rds

# Sequencing Data Analysis

## Dada2
### Making ASVs from demultiplexed sequencing files

R Script: dada2.R

Input: fastq files

Output: seqtab.rds and taxa.rds

## Phyloseq
### Alpha diversity (violin plot, Kruskal-Wallis, Dunn test)
### Beta diversity (ordination, PERMANOVA)
### Taxa plot

R Script: phyloseq.R

Input: seqtab.rds, taxa.rds, DCPET_Full_Metadata.csv

Output: dcpet_rarefied_nochloroplasts.rds

## DESeq2
### Differential abundance analysis

R Script: deseq2.R

Input: dcpet_rarefied_nochloroplasts.rds
