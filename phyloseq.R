#-----MAKING A PHYLOSEQ OBJECT-----#

lib_loc <- "/home/lgschaer/R/x86_64-redhat-linux-gnu-library/3.6/"
lib_loc
to_install <- unname(installed.packages(lib.loc = lib_loc)[, "Package"])
to_install
install.packages(pkgs = to_install, dependencies = TRUE)

library(ggpubr)
library(tidyverse)
library(phyloseq)
library(csv)
library(pwr)
library(FSA)
library(tsnemicrobiota)
library(vegan)
library(ape)

# Load data into R for phyloseq analysis
# We will need (1) sample/meta data, (2) sequence table, and (3) taxa table
# (1) is a csv file with information about each sample (media, carbon, conditions, enrichment details, etc)
# (2) and (3) are output files from dada2.

#load sample data
sdata <- as.csv("/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/DCPET_Full_Metadata.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)

# reformat metadata (may or may not be necessary depending on how it is formatted)
sdata2 <- sdata %>% 
  rownames_to_column(var = "SampleID") %>%
  separate(Metadata_Info, into = c("Enrichment", "Carbon", "Media", "Replicate")) %>%
  tidyr::unite(Media_Carbon, Media, Carbon, sep = "_", remove = FALSE) %>%
  dplyr::select(c("SampleID", "Sample_Number", "Enrichment", "Carbon", "Media", "Replicate", "Media_Carbon"))
head(sdata2)                                                       #view data to make sure everything is OK

sdata3 <- sdata2 %>%                                               #make data to use in phyloseq object with Sample_ID as rownames
  column_to_rownames(var = "SampleID")
head(sdata3)                                                       #view data to make sure everything is OK

#load sequence table
sequence_table <- readRDS("/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/dada2_output/seqtab.rds")
colnames(sequence_table) <- NULL                                   #remove column names

#make nonzero subset to remove all columns with zero taxa counts
sequence_table <- as.matrix(sequence_table)                                #change to matrix format
m <- (colSums(sequence_table, na.rm=TRUE) != 0)                        #T if colSum is not 0, F otherwise
nonzero <- sequence_table[, m]                                         #all the non-zero columns

#load taxa table
taxa_table <- readRDS("/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/dada2_output/taxa.rds")
taxa_table <- as.matrix(taxa_table)                                #change to matrix format
taxa_table[1:5,1:5]                                                #view to make sure everything looks good

#make phyloseq object
samdata = sample_data(sdata3)                                      #define sample data
colnames(nonzero) <- NULL                                          #remove column names from "nonzero"
seqtab = otu_table(nonzero, taxa_are_rows = FALSE)                 #define sequence table
taxtab = tax_table(taxa_table)                                     #define taxa table
rownames(taxtab) <- NULL                                           #remove rownames from taxa table


length(sample_names(samdata))
length(sample_names(seqtab))

phyloseq_object_all = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
phyloseq_object_all

sample_sums(phyloseq_object_all)

#sample counts before rarefying
sample_counts <- sample_data(phyloseq_object_all) %>%
  group_by(Enrichment, Carbon, Media) %>%
  mutate(Count = 1) %>%
  summarise(SumCount = sum(Count))
sample_counts

#normalize data
#Delete samples with a mean of less than 1000
samplesover1000_all <- subset_samples(phyloseq_object_all, sample_sums(phyloseq_object_all) > 1000)
min(sample_sums(samplesover1000_all))

#Check if there are OTUs with no counts, if so how many?
any(taxa_sums(samplesover1000_all) == 0)
sum(taxa_sums(samplesover1000_all) == 0)

#Prune OTUs with no counts 
prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)
any(taxa_sums(prune_samplesover1000_all) == 0)

#make sure seed is set the same each time, set to 81 here
rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))

#number rarefied at:
min(sample_sums(prune_samplesover1000_all))

#filter out eukaryotes and mitochondria
head(sample_data(rarefy_samplesover1000_all))

justbacteria <- rarefy_samplesover1000_all %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  ) %>%
  subset_samples(Sample_Number != "26")
justbacteria

#saveRDS file
saveRDS(justbacteria, file = "/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/phyloseq_output/dcpet_rarefied_nochloroplasts.rds")

#-----ALPHA DIVERISTY-----#

#Violin plot of alpha diversity Observed OTUs and Shannon Diversity (with color)
colors <- c("PositiveControl" = "lightblue", "Emma2" = "purple", "Laura1" = "orange", "Laura2" = "firebrick")

PlotA <- justbacteria %>%                                                     #phyloseq object
  plot_richness(
    x = "Enrichment",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = Enrichment), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = colors)+                         #set fill colors
  #ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position
PlotA

colors2 <- c("Emma2" = "purple", "Laura1" = "orange", "Laura2" = "firebrick")

justbacteria2 <- subset_samples(justbacteria, Enrichment!="PositiveControl")
justbacteria2

justbacteria2 %>%                                                     #phyloseq object
  plot_richness(
    x = "Media_Carbon",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  #geom_violin(aes(fill = Media_Carbon), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=4) +                                          #add boxplot, set width
  geom_jitter(aes(fill = Enrichment),color = "black", size = 5, shape = 21)+
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = colors2)+                         #set fill colors
  #ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

#-----ALPHA DIVERISTY STATISTICS-----#

#power analysis to determime whether we have sufficent power use an ANOVA with non-normally distributed data.
#count samples
sample_counts <- sample_data(justbacteria2) %>%
  group_by(Media_Carbon) %>%
  mutate(Count = 1)%>%
  summarise(SumCount = sum(Count)) %>%
  mutate(countpercategory = mean(SumCount))
sample_counts
head(sample_data(justbacteria))


#where k = number of groups, f = effect size, sig.level = type 1 error probability, power (1-probability of type 2 error)
pwr.anova.test(k=4, f=0.23, sig.level=.05, power=.8)
# in the output the number given for n must be smaller than or equal to the number of samples in the smallest group
#this data set is too small to meet the assumptions of an ANOVA

#add alpha diversity data to a data frame
richness <- justbacteria2 %>%
  estimate_richness(measures = c("Observed", "Shannon")) %>%           #specify which measures
  rownames_to_column(var = "SampleID") %>%                             #add column name to SampleID column
  as_tibble() 
head(richness)

alphadiv <- richness %>%
  left_join(sdata2, by = "SampleID")
head(alphadiv)

#Kruskal-Wallis Test
set.seed(81)

##BY ENRICHMENT
#Observed
kruskal.test(Observed ~ Enrichment, data = alphadiv) 

#Shannon
kruskal.test(Shannon ~ Enrichment, data = alphadiv)

#Dunn test (post hoc)

##Shannon
dunnS <- dunnTest(Shannon ~ Enrichment,
                  data=alphadiv,
                  method="bh")
dunnS

##BY MEDIA/CARBON & ENRICHMENT
head(alphadiv)

L1_alphadiv <- filter(alphadiv, Enrichment == "Laura1")
L2_alphadiv <- filter(alphadiv, Enrichment == "Laura2")
E2_alphadiv <- filter(alphadiv, Enrichment == "Emma2")

#Observed
kruskal.test(Observed ~ Media_Carbon, data = L1_alphadiv) 
kruskal.test(Observed ~ Media_Carbon, data = L2_alphadiv)
kruskal.test(Observed ~ Media_Carbon, data = E2_alphadiv) 

#Shannon
kruskal.test(Shannon ~ Media_Carbon, data = L1_alphadiv) 
kruskal.test(Shannon ~ Media_Carbon, data = L2_alphadiv) 
kruskal.test(Shannon ~ Media_Carbon, data = E2_alphadiv) 


#-----BETA DIVERSITY-----#

#adding a phylogenetic tree to phyloseq object using ape library

random_tree = rtree(ntaxa(justbacteria2), rooted=TRUE, tip.label=taxa_names(justbacteria2))
plot(random_tree)

justbacteria3 = merge_phyloseq(justbacteria2, samdata, random_tree)
justbacteria3

#ordination
distance <- ordinate(
  physeq = justbacteria3, 
  method = "PCoA", 
  distance = "unifrac"
)
#summary(distance)
#distance

#plot
PlotB <- plot_ordination(
  physeq = justbacteria3,                                                          #phyloseq object
  ordination = distance)+                                                #ordination
  geom_point(aes(fill = Enrichment, shape = Media_Carbon), size = 6) +                         #sets fill color to sampletype
  scale_fill_manual(values = colors2) +
  scale_shape_manual(values = c(21, 22, 23, 24))+
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"),  #adds black boarder around legend
    legend.position = c(0.15,0.2),
    legend.text = element_text(size = 12, face = "bold"),                                 
    axis.text.y.left = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))          #fills legend points based on the fill command
PlotB

ggarrange(PlotA,PlotB, common.legend = FALSE, labels = c("A", "B"))

# I made the next two plots, but they were not included in the paper.
#t-SNE plot
head(sdata3)

tsne <- tsne_phyloseq(justbacteria2, distance = "bray", perplexity = 15, dimensions = 2,
                      precomputed_distance = NULL, pseudocounts = 1, verbose = 1,
                      rng_seed = 81, philr_options = list(), control = list())
summary(tsne)

justbacteria2

#tSNE Plot
plot_tsne_phyloseq(justbacteria2, tsne, color = "Enrichment", shape = "Carbon") +
  geom_point(aes(color = Enrichment, fill = Enrichment), color = "black", size = 5, show.legend = TRUE) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = colors2) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(color = FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

#PCOA Plot

#ordination
all_pcoa <- ordinate(
  physeq = justbacteria2, 
  method = "PCoA", 
  distance = "bray"
)

colors3 <- c("green", "blue", "lightblue")

#plot
PlotB <- plot_ordination(
  physeq = justbacteria2,                                                          #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(fill = Carbon, shape = Enrichment), size = 6) +                         #sets fill color to sampletype
  scale_fill_manual(values = colors3) +
  scale_shape_manual(values = c(21, 22, 23))+
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank())+                                      #removes legend title
    #legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.position = "bottom",
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
PlotB

#-----BETA DIVERSTIY STATISTICS-----#

#PERMANOVA
set.seed(81)
head(sample_data(justbacteria2))

#looking at differences between enrichments
#subset phyloseq object
E2_L1 <- subset_samples(justbacteria2, Enrichment %in% c("Emma2", "Laura1"))
E2_L2 <- subset_samples(justbacteria2, Enrichment %in% c("Emma2", "Laura2"))
L1_L2 <- subset_samples(justbacteria2, Enrichment %in% c("Laura1", "Laura2"))

# Calculate bray curtis distance matrix, all samples
E2_L1bray1 <- phyloseq::distance(E2_L1, method = "bray")
E2_L2bray2 <- phyloseq::distance(E2_L2, method = "bray")
L1_L2bray3 <- phyloseq::distance(L1_L2, method = "bray")

# make a data frame from the sample_data, all samples
E2_L1sam1 <- data.frame(sample_data(E2_L1))
E2_L2sam2 <- data.frame(sample_data(E2_L2))
L1_L2sam3 <- data.frame(sample_data(L1_L2))

# Adonis test, all samples
adonis(E2_L1bray1 ~ Enrichment, data = E2_L1sam1)
adonis(E2_L2bray2 ~ Enrichment, data = E2_L2sam2)
adonis(L1_L2bray3 ~ Enrichment, data = L1_L2sam3)

#looking at differences between media/carbon type within each enrichment
#subset phyloseq object
E2<-subset_samples(justbacteria2, Enrichment == "Emma2")
L1<-subset_samples(justbacteria2, Enrichment == "Laura1")
L2<-subset_samples(justbacteria2, Enrichment == "Laura2")

e2one <- subset_samples(E2, Media_Carbon %in% c("ARW_DCPET", "AMW_DCPET"))
e2two <- subset_samples(E2, Media_Carbon %in% c("ARW_DCPET", "BH_DCPET"))
e2three <- subset_samples(E2, Media_Carbon %in% c("ARW_DCPET", "BH_TPA"))
e2four <- subset_samples(E2, Media_Carbon %in% c("AMW_DCPET", "BH_DCPET"))
e2five <- subset_samples(E2, Media_Carbon %in% c("AMW_DCPET", "BH_TPA"))
e2six <- subset_samples(E2, Media_Carbon %in% c("BH_DCPET", "BH_TPA"))

l1one <- subset_samples(L1, Media_Carbon %in% c("ARW_DCPET", "AMW_DCPET"))
l1two <- subset_samples(L1, Media_Carbon %in% c("ARW_DCPET", "BH_DCPET"))
l1three <- subset_samples(L1, Media_Carbon %in% c("ARW_DCPET", "BH_TPA"))
l1four <- subset_samples(L1, Media_Carbon %in% c("AMW_DCPET", "BH_DCPET"))
l1five <- subset_samples(L1, Media_Carbon %in% c("AMW_DCPET", "BH_TPA"))
l1six <- subset_samples(L1, Media_Carbon %in% c("BH_DCPET", "BH_TPA"))

l2one <- subset_samples(L2, Media_Carbon %in% c("ARW_DCPET", "AMW_DCPET"))
l2two <- subset_samples(L2, Media_Carbon %in% c("ARW_DCPET", "BH_DCPET"))
l2three <- subset_samples(L2, Media_Carbon %in% c("ARW_DCPET", "BH_TPA"))
l2four <- subset_samples(L2, Media_Carbon %in% c("AMW_DCPET", "BH_DCPET"))
l2five <- subset_samples(L2, Media_Carbon %in% c("AMW_DCPET", "BH_TPA"))
l2six <- subset_samples(L2, Media_Carbon %in% c("BH_DCPET", "BH_TPA"))


# Calculate bray curtis distance matrix, all samples
E2bray1 <- phyloseq::distance(e2one, method = "bray")
E2bray2 <- phyloseq::distance(e2two, method = "bray")
E2bray3 <- phyloseq::distance(e2three, method = "bray")
E2bray4 <- phyloseq::distance(e2four, method = "bray")
E2bray5 <- phyloseq::distance(e2five, method = "bray")
E2bray6 <- phyloseq::distance(e2six, method = "bray")

L1bray1 <- phyloseq::distance(l1one, method = "bray")
L1bray2 <- phyloseq::distance(l1two, method = "bray")
L1bray3 <- phyloseq::distance(l1three, method = "bray")
L1bray4 <- phyloseq::distance(l1four, method = "bray")
L1bray5 <- phyloseq::distance(l1five, method = "bray")
L1bray6 <- phyloseq::distance(l1six, method = "bray")

L2bray1 <- phyloseq::distance(l2one, method = "bray")
L2bray2 <- phyloseq::distance(l2two, method = "bray")
L2bray3 <- phyloseq::distance(l2three, method = "bray")
L2bray4 <- phyloseq::distance(l2four, method = "bray")
L2bray5 <- phyloseq::distance(l2five, method = "bray")
L2bray6 <- phyloseq::distance(l2six, method = "bray")

# make a data frame from the sample_data, all samples
E2sam1 <- data.frame(sample_data(e2one))
E2sam2 <- data.frame(sample_data(e2two))
E2sam3 <- data.frame(sample_data(e2three))
E2sam4 <- data.frame(sample_data(e2four))
E2sam5 <- data.frame(sample_data(e2five))
E2sam6 <- data.frame(sample_data(e2six))

L1sam1 <- data.frame(sample_data(l1one))
L1sam2 <- data.frame(sample_data(l1two))
L1sam3 <- data.frame(sample_data(l1three))
L1sam4 <- data.frame(sample_data(l1four))
L1sam5 <- data.frame(sample_data(l1five))
L1sam6 <- data.frame(sample_data(l1six))

L2sam1 <- data.frame(sample_data(l2one))
L2sam2 <- data.frame(sample_data(l2two))
L2sam3 <- data.frame(sample_data(l2three))
L2sam4 <- data.frame(sample_data(l2four))
L2sam5 <- data.frame(sample_data(l2five))
L2sam6 <- data.frame(sample_data(l2six))

# Adonis test, all samples
adonis(E2bray1 ~ Media_Carbon, data = E2sam1)
adonis(E2bray2 ~ Media_Carbon, data = E2sam2)
adonis(E2bray3 ~ Media_Carbon, data = E2sam3)
adonis(E2bray4 ~ Media_Carbon, data = E2sam4)
adonis(E2bray5 ~ Media_Carbon, data = E2sam5)
adonis(E2bray6 ~ Media_Carbon, data = E2sam6)

adonis(L1bray1 ~ Media_Carbon, data = L1sam1)
adonis(L1bray2 ~ Media_Carbon, data = L1sam2)
adonis(L1bray3 ~ Media_Carbon, data = L1sam3)
adonis(L1bray4 ~ Media_Carbon, data = L1sam4)
adonis(L1bray5 ~ Media_Carbon, data = L1sam5)
adonis(L1bray6 ~ Media_Carbon, data = L1sam6)

adonis(L2bray1 ~ Media_Carbon, data = L2sam1)
adonis(L2bray2 ~ Media_Carbon, data = L2sam2)
adonis(L2bray3 ~ Media_Carbon, data = L2sam3)
adonis(L2bray4 ~ Media_Carbon, data = L2sam4)
adonis(L2bray5 ~ Media_Carbon, data = L2sam5)
adonis(L2bray6 ~ Media_Carbon, data = L2sam6)



#EXPLORING TAXA

phyloseq_object_all
justbacteria2

#Summarize abundance of each class
genusabundance <- justbacteria2 %>%
  tax_glom(taxrank = "Genus") %>%                      # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 
head(genusabundance)

#just the positive control
##pos <- genusabundance %>%
  #dplyr::select(Phylum, Class, Family, Genus, Sample, Abundance, Enrichment, Media_Carbon, Replicate) %>%
#  filter(Abundance != 0) %>%
 # filter(Enrichment == "PositiveControl") %>%
  #mutate(
   # Phylum = as.character(Phylum),
    #Class = as.character(Class),
#    Family = as.character(Family),
 #   Genus = as.character(Genus)) %>%
  #group_by(Enrichment, Media_Carbon, Genus)%>%
#  summarise(Abundance = sum(Abundance)) %>%
 # mutate(Genus.2p = ifelse(Abundance < 0.02, "<2%", Genus))#%>%
  #arrange(desc(Abundance, Sample))                           #Arrange with descending abundances
#head(pos)

#View(pos)

#save color palatte
colors10 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   
  "grey47",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue",    "white", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue"
) 

#ggplot(pos)+
 # geom_col(mapping = aes(x = Media_Carbon, y = Abundance, fill = Genus), color = "black", position = "fill", show.legend = TRUE)+
  #facet_grid(rows = vars(Enrichment))+
#  ylab("Proportion of Community") +
 # scale_fill_manual(values = colors10) +
  #xlab(NULL)+
#  theme_minimal()+
 # theme(axis.text.y.left = element_text(size = 20),
  #      axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0.5),
   #     axis.title.y = element_text(size = 20),
    #    legend.text = element_text(size = 20),
     #   title = element_text(size = 25))


#Select and summarize necessary variables
all <- genusabundance %>%
  select(Phylum, Class, Family, Genus, Sample, Abundance, Enrichment, Media_Carbon, Replicate) %>%
  filter(Abundance != 0) %>%
  filter(Enrichment != "PositiveControl") %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus)
    )
head(all)

phylum <- all %>%
  dplyr::group_by(Enrichment, Media_Carbon, Phylum)%>%
  summarise(Abundance = sum(Abundance)/n()) %>%
  mutate(Phylum.1p = ifelse(Abundance < 0.01, "<1%", Phylum))
head(phylum)

range(phylum$Abundance)

class <- all %>%
  group_by(Enrichment, Media_Carbon, Class)%>%
  summarise(Abundance = sum(Abundance)/n()) #%>%
  #mutate(Class = ifelse(Abundance < 0.01, "<1%", Class))
head(class)

family <- all %>%
  group_by(Enrichment, Media_Carbon, Family)%>%
  summarise(Abundance = sum(Abundance)/n()) %>%
  mutate(Family = ifelse(Abundance < 0.01, "<1%", Family))
head(family)

genus <- all %>%
  group_by(Enrichment, Media_Carbon, Genus)%>%
  summarise(Abundance = sum(Abundance)/n()) %>%
  mutate(Genus = ifelse(Abundance < 0.02, "<2%", Genus))
head(genus)


#MAKING A TAXA PLOT 

#save color palatte
colors9 <- c(
  "orchid1",     "darkcyan",     "green",       "blue",   
  "grey47",     "cyan",         "coral1",      "darkgreen",   "palegoldenrod",    
  "grey77",     "darkblue",     "orange",      "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue", "firebrick",    "yellowgreen", "magenta",     "green", "red", "orchid", "lightblue"
)  


colors10 <- c(
  "black",      "darkcyan",     "orchid1",     "green",       "blue",   
  "grey47",     "cyan",         "coral1",      "darkgreen",   "palegoldenrod",    
  "grey77",     "darkblue",     "orange",      "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue", "firebrick",    "yellowgreen", "magenta",     "green", "red", "orchid", "lightblue"
)  

length(colors10)

#plot
unique(phylum$Media_Carbon)
Media_Carbon_Labels <- c("AMW_DCPET" = "AMW\nDCPET", "ARW_DCPET" = "ARW\nDCPET", "BH_DCPET" = "BH\nDCPET" , "BH_TPA" = "BH\nTPA")

phy <- ggplot(phylum)+
  geom_col(mapping = aes(x = Media_Carbon, y = Abundance, fill = Phylum), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(cols = vars(Enrichment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors9) +
  scale_x_discrete(labels = Media_Carbon_Labels)+
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 12, face = "bold", angle = 0),
    #    legend.title = element_blank(),
        title = element_text(size = 18))

cla <- ggplot(class)+
  geom_col(mapping = aes(x = Media_Carbon, y = Abundance, fill = Class), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(cols = vars(Enrichment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors9) +
  scale_x_discrete(labels = Media_Carbon_Labels)+
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 12, face = "bold", angle = 0),
     #   legend.title = element_blank(),
        title = element_text(size = 18))

fam <- ggplot(family)+
  geom_col(mapping = aes(x = Media_Carbon, y = Abundance, fill = Family), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(cols = vars(Enrichment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  scale_x_discrete(labels = Media_Carbon_Labels)+
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 12, face = "bold", angle = 0),
       # legend.title = element_blank(),
        title = element_text(size = 18))

gen <- ggplot(genus)+
  geom_col(mapping = aes(x = Media_Carbon, y = Abundance, fill = Genus), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(cols = vars(Enrichment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  scale_x_discrete(labels = Media_Carbon_Labels)+
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 12, face = "bold", angle = 0),
      #  legend.title = element_blank(),
        title = element_text(size = 18))

ggarrange(phy, cla, fam, gen, 
          align = "hv", 
          legend = "bottom", 
          ncol = 2, nrow = 2,
         # hjust = -6, vjust = 0.7,
          labels = c("A", "B", "C", "D"))


#summaries

sample_counts <- sample_data(justbacteria2) %>%
  group_by(Enrichment) %>%
  mutate(Count = 1)%>%
  summarise(SumCount = sum(Count))# %>%
#mutate(countpercategory = mean(SumCount))
sample_counts


phy.summary <- phylum  %>%
  mutate(totalSum = sum(Abundance),
         RelAb = Abundance/totalSum) %>%
  filter(Phylum =="Actinobacteriota" | Phylum == "Proteobacteria") %>%
  ungroup() %>%
  dplyr::group_by(Enrichment, Phylum) %>%
  summarise(
    maxRelAb = max(RelAb),
    minRelAb = min(RelAb)
    ) %>%
  unique()
View(phy.summary)

class.summary <- class  %>%
  mutate(totalSum = sum(Abundance),
         RelAb = Abundance/totalSum) %>%
  ungroup() %>%
  dplyr::group_by(Enrichment, Class) %>%
  summarise(
    maxRelAb = format(max(RelAb), scientific = F, digits = 3),
    minRelAb = format(min(RelAb), scientific = F, digits = 3)
  ) 
View(class.summary)


family.genus.summary <- all  %>%
  group_by(Enrichment, Family, Genus)%>%
  summarise(maxRelAb = format(max(Abundance), scientific = FALSE, digits = 3),
            meanRelAb = format(mean(Abundance), scientific = FALSE, digits = 3),
            minRelAb = format(min(Abundance), scientific = FALSE, digits = 3)) 
View(family.genus.summary)
write_csv(family.genus.summary, "/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/phyloseq_output/max_mean_min_relative_abundance_family_genus_level.csv")


#genus by replicate taxa plot
genus2 <- all %>%
  group_by(Enrichment, Media_Carbon, Genus, Replicate)%>%
  summarise(Abundance = sum(Abundance)) %>%
  mutate(Genus.2p = ifelse(Abundance < 0.02, "<2%", Genus))#%>%
  #arrange(desc(Abundance, Sample))                           #Arrange with descending abundances
head(genus2)

ggplot(genus2)+
  geom_col(mapping = aes(x = Replicate, y = Abundance, fill = Genus.2p), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(rows = vars(Enrichment), cols = vars(Media_Carbon))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 25))

# A closer look at Laura2 at the genus level

L2_genus <- all %>%
    mutate(Genus.2p = ifelse(Abundance < 0.02, "<2%", Genus)) %>%
    filter(Enrichment == "Laura2")
head(L2_genus)

ggplot(L2_genus)+
  geom_col(mapping = aes(x = Replicate, y = Abundance, fill = Genus.2p), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(cols = vars(Media_Carbon))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  scale_x_discrete(labels = Media_Carbon_Labels)+
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 12, face = "bold", angle = 0),
        #  legend.title = element_blank(),
        title = element_text(size = 18))
