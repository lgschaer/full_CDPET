#---------------Chemistry Analysis---------------#

#Packages used
library(csv)
library(tidyverse)
#devtools::install_github("ianmoran11/mmtable2")
library(mmtable2)
library(gt)
#install.packages("reactable")
library(reactable)


#Load data
##Your data should be formatted in a CSV file
chem <- as.csv("/Users/lgschaer/Desktop/DCPET_AQUA_Results_forR.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(chem)

chem2 <- chem %>%
  mutate(
    SRP = SRP_mg_P_L*Dilution_Factor,
    NO3NO2 = NO3_NO2_mg_N_L*Dilution_Factor,
    NH4 = NH4_mg_N_L*Dilution_Factor,
    DOC = DOC_mg_L*Dilution_Factor,
    TDN = TDN_mg_L*Dilution_Factor
  ) %>%
  select(SRP, NO3NO2, NH4, DOC, TDN) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Measurement_Type") %>%
  pivot_longer(cols = c(MilliQ, Pilgrim_River, BH_1_1000), names_to = "Media", values_to = "Value") %>%
  mutate(
    DCPET_Conc_low = (DCPET_1_10to5th*.0025)*50,
    DCPET_Conc_high = (DCPET_1_10to6th*.0025)*50, 
    with_DCPET_low = Value + DCPET_Conc_low,
    with_DCPET_low = format(round(with_DCPET_low, 0), nsmall = 0, format="d", big.mark=","),
    with_DCPET_high = Value + DCPET_Conc_high,
    with_DCPET_high = format(round(with_DCPET_high, 0), nsmall = 0, format="d", big.mark=","),
    without_DCPET_low = format(round(Value, 0), nsmall = 0, format="d", big.mark=","),
    without_DCPET_high = format(round(Value, 0), nsmall = 0, format="d", big.mark=","),
    with_DCPET = paste(with_DCPET_low,"to",with_DCPET_high),
    without_DCPET = paste("approx.", without_DCPET_low)
  ) %>%
  select(Measurement_Type, Media, with_DCPET, without_DCPET) %>%
  pivot_longer(cols = c(with_DCPET, without_DCPET), names_to = "Type", values_to = "Value", names_repair = "minimal") %>%
  mutate(
    Type = ifelse(Type=="with_DCPET", "With DCPET", "No DCPET"),
    Media = case_when(Media=="Pilgrim_River"~"Pilgrim River",
                      Media=="MilliQ"~"MilliQ Water",
                      Media=="BH_1_1000"~"Bushnell Haas"),
    Measurement_Type_Units = case_when(Measurement_Type=="NH4"~"(N mg/L)",
                                       Measurement_Type=="SRP"~"(P mg/L)",
                                       Measurement_Type=="TDN"~"(mg/L)",
                                       Measurement_Type=="DOC"~"(mg/L)", 
                                       Measurement_Type=="NO3NO2"~"(N mg/L)"),
    Measurement_Type = case_when(Measurement_Type=="NH4"~"Ammonium",
                                 Measurement_Type=="SRP"~"Soluble Reactive Phosphorus",
                                 Measurement_Type=="TDN"~"Total Dissolved Nitrogen",
                                 Measurement_Type=="DOC"~"Dissolved Organic Carbon", 
                                 Measurement_Type=="NO3NO2"~"Nitrate & Nitrite")
  )
head(chem2)
#View(chem2)

chem_table <- 
  chem2 %>% 
  mmtable(cells = Value) +
  header_left(Type) +
  header_top(Measurement_Type_Units) +
  header_top(Measurement_Type) +
  header_left_top(Media)  +
  header_format(header = Measurement_Type_Units, style = list(cell_text(align = "center", weight = "normal"),
                                                              cell_borders(sides = "bottom",color = "grey")))+
  header_format(header = Measurement_Type, style = list(cell_text(align = "center"),
                                                        cell_borders(sides = "top",color = "grey")))+
  cells_format(style = list(cell_text(align = "center")))
  #header_top_left(continent) + 
  #header_format(var, scope = "table", style = style_list)
chem_table
