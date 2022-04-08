#Load required packages
library(tidyverse)
library(csv)

#Using an ANCOVA to compare growth rates in DCPET treatments
#http://r-eco-evo.blogspot.com/2011/08/comparing-two-regression-slopes-by.html

growth <- as.csv("/Users/lgschaer/Desktop/MTU_Projects/DCPET/06072021_DCPET_Experiment/DCPET_Experiment_Data.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(growth)

#Start with Laura1
l1_ancova_data <- growth %>%
  select(Enrichment, Carbon, Media, Hours, Absorbance, Time_Point) %>%
  filter(Time_Point=="T0"|Time_Point=="T1"|Time_Point=="T2"|Time_Point=="T3"|Time_Point=="T4") %>%
  filter(Enrichment == "Laura1" & Carbon == "DCPET") %>%
  unite(Group, Enrichment, Carbon, Media, sep = "_")
head(l1_ancova_data)

#To test for significant differences in slope between groups
l1_ancova_mod1 <- aov(Absorbance~Hours*Group, data=l1_ancova_data) #test for significant interaction between two variables
summary(l1_ancova_mod1)
l1_ancova_mod2 <- aov(Absorbance~Hours+Group, data=l1_ancova_data) #test for significant differences between slope (growth rates)
summary(l1_ancova_mod2)
#ANOVA to compare the two models and see if the interaction term adds important information to the model. 
##If p-value is non-significant, interaction term does not add information.
anova(l1_ancova_mod1,l1_ancova_mod2)   

#Repeat for Laura2
l2_ancova_data <- growth %>%
  select(Enrichment, Carbon, Media, Hours, Absorbance, Time_Point) %>%
  filter(Time_Point=="T0"|Time_Point=="T1"|Time_Point=="T2"|Time_Point=="T3"|Time_Point=="T4") %>%
  filter(Enrichment == "Laura2" & Carbon == "DCPET") %>%
  unite(Group, Enrichment, Carbon, Media, sep = "_")
head(l2_ancova_data)

l2_ancova_mod1 <- aov(Absorbance~Hours*Group, data=l2_ancova_data) 
summary(l2_ancova_mod1)
l2_ancova_mod2 <- aov(Absorbance~Hours+Group, data=l2_ancova_data) 
summary(l2_ancova_mod2)
anova(l2_ancova_mod1,l2_ancova_mod2)   

#Repeat for Emma2
e2_ancova_data <- growth %>%
  select(Enrichment, Carbon, Media, Hours, Absorbance, Time_Point) %>%
  filter(Time_Point=="T0"|Time_Point=="T1"|Time_Point=="T2"|Time_Point=="T3"|Time_Point=="T4") %>%
  filter(Enrichment == "Emma2" & Carbon == "DCPET") %>%
  unite(Group, Enrichment, Carbon, Media, sep = "_")
head(e2_ancova_data)

e2_ancova_mod1 <- aov(Absorbance~Hours*Group, data=e2_ancova_data) #test for significant interaction between two variables
summary(e2_ancova_mod1)
e2_ancova_mod2 <- aov(Absorbance~Hours+Group, data=e2_ancova_data) #test for significant differences between slope (growth rates)
summary(e2_ancova_mod2)
anova(e2_ancova_mod1,e2_ancova_mod2)   

#TukeyHSD post-hoc test
#for Emma2 DCPET treatments shows significant difference between ARW and AMW/BH treatments
TukeyHSD(e2_ancova_mod2)
#TukeyHSD test is not recommended as a post hoc test for ANCOVA since it does not appropriately correct p-values

#Pairwise t-test post hoc test with bonferroni p-value adjustment
#performed for Emma2 to determine between which groups there is a statistically significant difference.
pairwise.t.test(e2_ancova_data$Absorbance, e2_ancova_data$Group, p.adj = "none")
#Adjusted p-values all non-significant

