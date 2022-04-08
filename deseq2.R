#install DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("DESeq2", version = "3.10", dependencies = TRUE)

#packages used
library(tidyverse)
#install.packages("matrixStats")
library(matrixStats)
library(DESeq2)
library(phyloseq)
library(BiocManager)

#new approach

# Load phyloseq object
justbacteria <- readRDS("/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/phyloseq_output/dcpet_rarefied_nochloroplasts.rds")
head(sample_data(justbacteria))


#### Substrate Comparisons
# Start by pre-processing data
e2_only_ps <- subset_samples(justbacteria, Enrichment == "Emma2" & Media_Carbon == "BH_DCPET" | Media_Carbon == "BH_TPA")
l1_only_ps <- subset_samples(justbacteria, Enrichment == "Laura1" & Media_Carbon == "BH_DCPET" | Media_Carbon == "BH_TPA")
l2_only_ps <- subset_samples(justbacteria, Enrichment == "Laura2" & Media_Carbon == "BH_DCPET" | Media_Carbon == "BH_TPA")


# Convert phyloseq object to deseq2 format
ps.e2 <- phyloseq_to_deseq2(e2_only_ps, ~ Media_Carbon)
ps.l1 <- phyloseq_to_deseq2(l1_only_ps, ~ Media_Carbon)
ps.l2 <- phyloseq_to_deseq2(l2_only_ps, ~ Media_Carbon)

ps.e2$Media_Carbon<-relevel(ps.e2$Media_Carbon, "BH_DCPET")
ps.l1$Media_Carbon<-relevel(ps.l1$Media_Carbon, "BH_DCPET")
ps.l2$Media_Carbon<-relevel(ps.l2$Media_Carbon, "BH_DCPET")

# Run DESeq2 analysis (all taxa at once!)
dds_ps.e2 <- DESeq(ps.e2)
dds_ps.l1 <- DESeq(ps.l1)
dds_ps.l2 <- DESeq(ps.l2)

# Investigate results
resultsNames(dds_ps.e2)
resultsNames(dds_ps.l1)
resultsNames(dds_ps.l2)

E2 <- as.data.frame(results(dds_ps.e2, contrast=c("Media_Carbon","BH_TPA","BH_DCPET"))) %>% mutate(Comparison="Emma2\nTPA vs. DCPET") %>% rownames_to_column(var = "taxon")
L1 <- as.data.frame(results(dds_ps.l1, contrast=c("Media_Carbon","BH_TPA","BH_DCPET"))) %>% mutate(Comparison="Laura1\nTPA vs. DCPET") %>% rownames_to_column(var = "taxon")
L2 <- as.data.frame(results(dds_ps.l2, contrast=c("Media_Carbon","BH_TPA","BH_DCPET"))) %>% mutate(Comparison="Laura2\nTPA vs. DCPET") %>% rownames_to_column(var = "taxon")
head(L2)

enriched <- E2 %>% 
  full_join(L1) %>%
  full_join(L2) %>%
  mutate(
    threshold = ifelse(padj <= 0.05 & abs(log2FoldChange) >= 2, "Enriched", "Not_Enriched")
  ) 
head(enriched)

any(enriched$threshold == "Enriched")

# Save a csv of the results
write.csv(enriched,"/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/deseq2_output/enriched_within_enrichment_substrate_comparisons.csv")

# Add taxonomy to DESeq output
taxps.e2 = as(tax_table(justbacteria)[rownames(ps.e2), ], "matrix") %>% as.data.frame() %>% rownames_to_column(var = "taxon")
taxps.l1 = as(tax_table(justbacteria)[rownames(ps.l1), ], "matrix") %>% as.data.frame() %>% rownames_to_column(var = "taxon")
taxps.l2 = as(tax_table(justbacteria)[rownames(ps.l2), ], "matrix") %>% as.data.frame() %>% rownames_to_column(var = "taxon")
head(taxps.l2)

enriched_w_tax <- enriched %>%
  full_join(taxps.e2) %>%
  full_join(taxps.l1) %>%
  full_join(taxps.l2) %>%
  filter(!is.na(threshold)) %>%
  mutate(
    Enriched_Genus = ifelse(threshold == "Enriched", as.character(Genus), "Not_Enriched"),
    Enriched_Genus = ifelse(is.na(Enriched_Genus), "Unclassified_Genus", Enriched_Genus)
  ) %>%
  filter(!is.na(Enriched_Genus))
head(enriched_w_tax)

summary_table <- enriched_w_tax %>%
  filter(threshold == "Enriched") %>%
  dplyr::select(taxon, log2FoldChange, padj, Comparison, Phylum, Class, Order, Family, Genus) %>%
  mutate(Comparison = gsub("\n", " ", Comparison)) 
dim(summary_table)
View(summary_table)

write.csv(summary_table,"/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/deseq2_output/only_enriched_w_tax_within_enrichment_substrate_comparisons.csv")


enriched_counts <- enriched_w_tax %>%
  mutate(Category = ifelse(log2FoldChange > 2, "right", "not_sig"),
         Category = ifelse(log2FoldChange < -2 & Category == "not_sig", "left", Category)) %>%
  group_by(Comparison, Category) %>%
  summarise(Count = n()) %>%
  filter(Category != "not_sig")
enriched_counts

#any(is.na(enriched_w_tax))
#View(enriched_w_tax)

write.csv(enriched_w_tax,"/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/deseq2_output/enriched_w_tax_enriched_within_enrichment_substrate_comparisons.csv")

# Re-Order data to organize legend
sort(unique(enriched_w_tax$Enriched_Genus))
length(unique(enriched_w_tax$Enriched_Genus))

enriched_w_tax$Enriched_Genus <- factor(enriched_w_tax$Enriched_Genus, 
                                        levels = c("Achromobacter",      "Aminobacter",        "Ancylobacter",       "Aquamicrobium",      "Brevundimonas",      "Chelatococcus",     
                                                   "Chitinophaga",       "Flavobacterium",     "Hydrogenophaga",     "Hyphomonas",         "Legionella",         "Luteimonas",        
                                                   "Mesorhizobium",      "Millisia",           "Parapusillimonas",   "Parvibaculum",       "Pelagibacterium",   
                                                   "Persicitalea",       "Planktosalinus",     "Pseudaminobacter",   "Pseudolabrys",       "Pseudomonas",        "Pseudoxanthomonas", 
                                                   "Pusillimonas",       "Reyranella",         "Rhodobacter",        "Rhodococcus",        "Shinella",           "Sphingobacterium",  
                                                   "Tepidimonas",        "Verticiella",        "Youhaiella",   
                                                   "Unclassified_Genus",   "Not_Enriched"))
# Save color palette

cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
colors11 <- sample(cl, 54)

#colors11 <- c(
 # "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
  #"magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#  "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
 # "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
  #"magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","gray63","white"
#)

shapes <- c(
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24, 
  21, 22, 23, 24, 
  21, 22, 23, 24,
  21, 22
)
length(shapes)

colors11 <- c(
  "orangered",      "purple",        "green",           "cyan",          "orange",        "khaki4",             "mediumslateblue",
  "mediumpurple1",  "darkmagenta",   "darkgreen",       "wheat4",        "yellow",        "lawngreen",          "plum",  
  "royalblue",      "magenta",       "mediumseagreen",  "palegoldenrod", "grey47",        "chocolate4",         "darkorange3",        
  "lightblue",      "firebrick",     "yellowgreen",     "turquoise3",    "purple4",       "red",                "blue",            
  "lightcyan",       "coral1",       "cyan",            "goldenrod",     "black",         "white"   
) 
length(colors11)

#volcano plot:
ggplot(data=enriched_w_tax, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(fill=Enriched_Genus, shape = Enriched_Genus), color = "black", size=4,
  #shape = 21
  ) +
  facet_grid(cols = vars(Comparison))+
  xlim(c(-30, 30)) + 
  ylim(c(0, 17)) +
  scale_fill_manual(values=colors11) +
  scale_shape_manual(values=shapes) +
  labs(x = "log2 fold change", 
       y = "-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = -2, colour="#990000", linetype="dashed")+ 
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18))+
  guides(fill = guide_legend(override.aes = list(shape = shapes)))


#ggarrange(PlotA,PlotB, align = "hv", common.legend = FALSE)


#### Enrichment Comparisons
# Start by pre-processing data
enrichments_only_ps <- subset_samples(justbacteria, Enrichment == "Laura1"|Enrichment == "Laura2"|Enrichment == "Emma2")

# Convert phyloseq object to deseq2 format
ps1 <- phyloseq_to_deseq2(enrichments_only_ps, ~ Enrichment)
ps2 <- ps1
ps1$Enrichment<-relevel(ps1$Enrichment, "Emma2")
ps2$Enrichment<-relevel(ps2$Enrichment, "Laura1")

# Run DESeq2 analysis (all taxa at once!)
dds_ps1 <- DESeq(ps1)
dds_ps2 <- DESeq(ps2)

# Investigate results
resultsNames(dds_ps1)
resultsNames(dds_ps2)

E2L1 <- as.data.frame(results(dds_ps1, contrast=c("Enrichment","Emma2","Laura1"))) %>% mutate(Comparison="Laura1_vs_Emma2") %>% rownames_to_column(var = "taxon")
E2L2 <- as.data.frame(results(dds_ps1, contrast=c("Enrichment","Emma2","Laura2"))) %>% mutate(Comparison="Laura2_vs_Emma2") %>% rownames_to_column(var = "taxon")
L1L2 <- as.data.frame(results(dds_ps2, contrast=c("Enrichment","Laura1","Laura2"))) %>% mutate(Comparison="Laura2_vs_Laura1") %>% rownames_to_column(var = "taxon")
head(L1L2)

enriched <- E2L1 %>% 
  full_join(E2L2) %>%
  full_join(L1L2) %>%
  mutate(
    threshold = ifelse(padj <= 0.05 & abs(log2FoldChange) >= 2, "Enriched", "Not_Enriched")
  ) 
head(enriched)

any(enriched$threshold == "Enriched")

# Save a csv of the results
write.csv(enriched,"/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/deseq2_output/enriched_all_enrichment_comparisons.csv")

# Add taxonomy to DESeq output
taxps1 = as(tax_table(justbacteria)[rownames(ps1), ], "matrix") %>% as.data.frame() %>% rownames_to_column(var = "taxon")
taxps2 = as(tax_table(justbacteria)[rownames(ps2), ], "matrix") %>% as.data.frame() %>% rownames_to_column(var = "taxon")
head(taxps2)

enriched_w_tax <- enriched %>%
  full_join(taxps1) %>%
  full_join(taxps2) %>%
  filter(!is.na(threshold)) %>%
  mutate(
    Enriched_Genus = ifelse(threshold == "Enriched", as.character(Genus), "Not_Enriched"),
    Enriched_Genus = ifelse(is.na(Enriched_Genus), "Unclassified_Genus", Enriched_Genus)
         )
head(enriched_w_tax)
#View(enriched_w_tax)


enriched_counts <- enriched_w_tax %>%
  mutate(Category = ifelse(log2FoldChange > 2, "right", "not_sig"),
         Category = ifelse(log2FoldChange < -2 & Category == "not_sig", "left", Category)) %>%
  group_by(Comparison, Category) %>%
  summarise(Count = n()) %>%
  filter(Category != "not_sig")
enriched_counts

#any(is.na(enriched_w_tax))
#View(enriched_w_tax)

write.csv(enriched_w_tax,"/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/deseq2_output/enriched_w_tax_all_enrichment_comparisons.csv")

# Re-Order data to organize legend
sort(unique(enriched_w_tax$Enriched_Genus))
enriched_w_tax$Enriched_Genus <- factor(enriched_w_tax$Enriched_Genus, 
                                        levels = c("Achromobacter",        "Afipia",               "Aminobacter",          "Ancylobacter",         "Aquamicrobium",       
                                                   "Brevundimonas",        "Chelatococcus",        "Chitinophaga",         "Cohnella",             "Devosia",             
                                                   "Dokdonella",           "Eoetvoesia",           "Escherichia/Shigella", "Flavobacterium",       "Hydrogenophaga",      
                                                   "Luteimonas",           "Mesorhizobium",        "Millisia",             "Moheibacter",                 
                                                   "Ochrobactrum",         "Orrella",              "Paenibacillus",        "Parapusillimonas",     "Parvibaculum",        
                                                   "Pelagibacterium",      "Persicitalea",         "Planktosalinus",       "Pseudaminobacter",     "Pseudolabrys",        
                                                   "Pseudomonas",          "Pseudoxanthomonas",    "Pusillimonas",         "Reyranella",           "Rhodobacter",         
                                                   "Rhodococcus",          "Shinella",             "SN8",                  "Solitalea",            "Sphingobacterium",    
                                                   "Sphingopyxis",         "Tepidimonas",          "Variovorax",           "Verticiella",
                                                   "Unclassified_Genus",   "Not_Enriched"))
# Save color palette

cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
colors11 <- sample(cl, 9)

colors11 <- c(
  "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
  "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
  "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
  "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
  "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","gray63","white"
)

shapes <- c(
  21, 22, 23, 24, 25,
  21, 22, 23, 24, 25,
  21, 22, 23, 24, 25,
  21, 22, 23, 24, 25,
  21, 22, 23, 24, 25,
  21, 22, 23, 24, 25,
  21, 22, 23, 24, 25,
  21, 22, 23, 24, 25,
  21, 22, 23, 24, 25
)

colors11 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   "tan1", "purple1",
  "grey77",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "orange",  "darkblue",     "white",    "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue", "green4", "grey50",
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   "tan1", "purple1",
  "grey77",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "orange",  "darkblue",     "white",    "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue", "green4", "grey50"
) 

#volcano plot:
ggplot(data=enriched_w_tax, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(fill=Enriched_Genus, 
                 #shape = Enriched_Genus
                 ), color = "black", size=4,
             shape = 21
             ) +
  facet_grid(cols = vars(Comparison))+
  #xlim(c(-15, 15)) + 
  #ylim(c(0, 23)) +
  scale_fill_manual(values=colors11) +
  #scale_shape_manual(values=shapes) +
  labs(x = "log2 fold change", 
       y = "-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_hline(yintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = -2, colour="#990000", linetype="dashed")+ 
  theme(axis.text.y.left = element_text(size = 10),
      axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 0.5),
      axis.title.y = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.position = "bottom",
      title = element_text(size = 18))+
  guides(fill = guide_legend(override.aes = list(shape = 21))) #shapes)))


#ggarrange(PlotA,PlotB, align = "hv", common.legend = FALSE)
