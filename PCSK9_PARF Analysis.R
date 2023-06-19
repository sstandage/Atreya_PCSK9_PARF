# Analysis Script

library(tidyverse)
library(readxl)

#########################
## Reading in the Data ##
#########################

PARF <- read_excel("BALI Clinical data_Atreya_Dahmer.xlsx") %>% 
  rename(prism = prismscore,
         age = ageyrs,
         d1ARDS = day_firstARDS2_1_OIOSI_criteria,
         ARDS = everARDS2_1_OIOSI_criteria,
         cvfail = ever_cvfail,
         dialysis = ever_dialysis_c,
         ecmo = ever_ecmo) %>% 
  mutate(across(gender, ~factor(., levels = c("Male", "Female")))) %>% 
  mutate(across(race, ~case_match(., "Black/African American" ~ "Black",
                                  "American Indian or Alaskan Native" ~ "NatAmerican",
                                  "Native Hawaiian or Other Pacific Islander" ~ "Islander",
                                  "Multiracial/More than one race" ~ "Multiracial",
                                  c("Declined", "Unknown/Unavailable") ~ "Unknown",
                                  .default = race) %>% 
                  fct() %>% fct_infreq())) %>% 
  mutate(across(ethnic, ~case_match(., "Not Hispanic or Latino" ~ "NotHispanic",
                                    "Hispanic or Latino" ~ "Hispanic",
                                    "Unknown/Unavailable" ~ "Unknown") %>% 
                  fct(levels = c("NotHispanic", "Hispanic", "Unknown")))) %>%
  mutate(across(primary, ~case_match(., "Pneumonia (any organism)" ~ "Pneumonia",
                                     "Laryngotracheobronchitis (croup/tracheitis)" ~ "CroupLTB",
                                     "Asthma or reactive airway disease" ~ "Asthma",
                                     "Acute respiratory failure post BMT" ~ "BMT", 
                                     "Acute respiratory failure related to sepsis" ~ "Sepsis",
                                     "Acute chest syndrome/sickle cell disease" ~ "ACS",
                                     "Pulmonary edema" ~ "PulmEdema",
                                     "Thoracic trauma: pulmonary contusion or inhalation burns" ~ "TraumaBurn",
                                     "Pulmonary hemorrhage" ~ "PulmHemorrhage",
                                     "Aspiration pneumonia" ~ "AspPNA",
                                     "Chronic lung disease: cystic fibrosis or BPD" ~ "CLD",
                                     "Other, specify:" ~ "Other",
                                     "Acute respiratory failure related to multiple blood transfusions" ~ "BloodTrnsfsn",
                                     .default = primary) %>% 
                  fct() %>% fct_infreq())) %>% 
  mutate(across(c("primdeath", "secdeath"), ~case_match(., "Respiratory failure" ~ "RespFailure",
                                       "Sepsis/septic shock" ~ "Sepsis",
                                       "Other, specify:" ~ "Other",
                                       "Multisystem organ failure" ~ "MSOF",
                                       "Cancer" ~ "Cancer") %>% 
                  fct(levels = c("RespFailure", "Sepsis", "MSOF", "Cancer", "Other")))) %>% 
  mutate(premature = if_else(pastmhx_premature_c == 1, TRUE, FALSE),
         asthma = if_else(pastmhx_asthma_c == 1, TRUE, FALSE),
         bpd = if_else(pastmhx_dysplasia_c == 1, TRUE, FALSE),
         cf = if_else(pastmhx_cf_c == 1, TRUE, FALSE),
         t2dm = if_else(pastmhx_insulindiab_c == 1, TRUE, FALSE),
         immunodef = if_else(pastmhx_immunodef_c == 1, TRUE, FALSE),
         asppna = if_else(pastmhx_aspiration_c == 1, TRUE, FALSE),
         seizure = if_else(pastmhx_seizures_c == 1, TRUE, FALSE),
         cancerpmh = if_else(pastmhx_cancerrx_c == 1, TRUE, FALSE),
         scd = if_else(pastmhx_sicklecell_c == 1, TRUE, FALSE),
         cancer = if_else(cancer_c == 1, TRUE, FALSE),
         chromabn = if_else(chromabn_c == 1, TRUE, FALSE),
         dead90 = if_else(death_c == 1, TRUE, FALSE)) %>% 
  mutate(across(c("d1ARDS", "ARDS", "cvfail", "dialysis", "ecmo"), ~if_else(. == 1, TRUE, FALSE))) %>% 
  select(!contains("_c"), -matches("thrbm|IL_8|IL_1RA")) %>% 
  rename_with(.cols = everything(), ~str_replace(., "_03", ""))

Labs <- read_excel("BALI Clinical data_Atreya_Dahmer.xlsx") %>% 
  select(baliid, matches("thrbm|IL_8|IL_1RA")) %>% 
  rename_with(.cols = matches("thrbm|IL_8|IL_1RA"), ~str_replace(., "__pg_mL_|ngml", "")) %>% 
  pivot_longer(cols = 2:ncol(.), 
               names_to = c("Analyte", "Day"),
               names_sep = "_day",
               values_to = "Value",
               values_drop_na = TRUE) %>% 
  mutate(across(Analyte, ~fct(.))) %>% 
  mutate(across(Day, ~fct(., levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8")))) %>% 
  arrange(baliid, Analyte, Day)

LabsSum <- Labs %>% 
  group_by(baliid, Analyte) %>% 
  summarize(Count = n(),
            Max = max(Value),
            Min = min(Value),
            Mean = mean(Value),
            Median = median(Value),
            First = first(Value)) %>% 
  ungroup()

LabSum_wide <- LabsSum %>% 
  select(-Count) %>% 
  pivot_wider(names_from = Analyte, values_from = c("Max", "Min", "Mean", "Median", "First"), names_glue = "{Analyte}_{.value}", names_sort = TRUE) %>% 
  select(baliid, sort(colnames(.)))

Geno <- read_excel("BALI Genotyping results_Atreya_Dahmer.xlsx") %>%
  rename(baliid = SampleID) %>% 
  slice(1:462) %>% 
  select(where(~ !(all(is.na(.)) | all(. == "")))) %>% 
  mutate(LDLR_homo = if_else(LDLR_rs688 == "TT", TRUE, FALSE),
         PCSK9_LOF = if_else(PCSK9_rs562556 == "AG" | PCSK9_rs562556 == "GG" | PCSK9_rs11591147 == "GT" | PCSK9_rs11591147 == "TT" | PCSK9_rs11583680 == "CT" | PCSK9_rs11583680 == "TT", TRUE, FALSE),
         HMGCR_LOF = if_else(HMGCR_rs12916 == "CT" | HMGCR_rs12916 == "CC" | HMGCR_rs17238484 == "GT" | HMGCR_rs17238484 == "TT", TRUE, FALSE)) %>% 
  mutate(across(PCSK9_rs11583680, ~case_match(., "NA" ~ NA_character_, #Removing the one "NA" written in there.
                                              .default = PCSK9_rs11583680)))

Comorbidities <- colnames(select(PARF, premature:chromabn))
Severity <- colnames(select(PARF, d1ARDS:ecmo, highestOI:mods))
Outcomes <- colnames(select(PARF, dead90, primdeath, secdeath, durmv28, piculos, hosplos))

Combined <- PARF %>% 
  left_join(Geno, by = "baliid") %>% 
  left_join(LabSum_wide, by = "baliid") %>% 
  relocate(age:primary, .after = restoreid) %>% 
  relocate(all_of(Comorbidities), .after = ethnic) %>% 
  relocate(all_of(Severity), .after = prism) %>% 
  relocate(all_of(Outcomes), .after = mods) %>% 
  mutate(across(contains("_rs"), ~fct(.) %>% fct_infreq()))

SNPs <- colnames(select(Combined, contains("_rs")))
Cytokines <- colnames(select(Combined, matches("IL_|thrbm")))
