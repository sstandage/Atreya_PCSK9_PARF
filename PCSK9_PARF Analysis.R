# Analysis Script

library(tidyverse)
library(readxl)

#########################
## Reading in the Data ##
#########################

PARF <- read_excel("BALI Clinical data_Atreya_Dahmer.xlsx") %>% 
  rename(Subject)

Geno <- read_excel("BALI Genotyping results_Atreya_Dahmer.xlsx") %>% 
  slice(1:462) %>% 
  rename(baliid = SampleID)

Combined <- PARF %>% left_join(Geno, by = "")