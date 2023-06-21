# Analysis Script

library(tidyverse)
library(readxl)
library(broom)

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
  mutate(across(mods, ~if_else(. > 0, TRUE, FALSE))) %>% 
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
Severity <- colnames(select(PARF, d1ARDS:ecmo, highestOI, mods))
Outcomes <- colnames(select(PARF, dead90, durmv28, piculos, hosplos))

Combined <- PARF %>% 
  left_join(Geno, by = "baliid") %>% 
  left_join(LabSum_wide, by = "baliid") %>% 
  relocate(age:primary, .after = restoreid) %>% 
  relocate(all_of(Comorbidities), .after = ethnic) %>% 
  relocate(all_of(Severity), .after = prism) %>% 
  relocate(all_of(Outcomes), .after = mods) %>% 
  mutate(across(contains("_rs"), ~fct(.) %>% fct_infreq())) %>% 
  mutate(PG = if_else(PCSK9_LOF == FALSE | is.na(PCSK9_LOF), "Other", "LOF"),
         female = if_else(gender == "Female", TRUE, FALSE)) %>% 
  mutate(across(PG, ~fct(., levels = c("Other", "LOF")))) %>% 
  relocate(PG, .after = HMGCR_LOF) %>% 
  relocate(female, .after = gender)
  

SNPs <- colnames(select(Combined, contains("_rs")))
Cytokines <- colnames(select(Combined, matches("IL_|thrbm")))


#####################################
## Building the demographic tables ##
#####################################

DemoVars1 <- c("female", Comorbidities)
DemoTab1 <- Combined %>% 
  group_by(PG) %>% 
  summarize(across(all_of(DemoVars1), list(
    N = ~ n(),
    Count = ~sum(. == TRUE, na.rm = TRUE)),
    .names = "{col}.{fn}")) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Vars", values_to = "Vals") %>% 
  pivot_wider(names_from = PG, values_from = Vals) %>% 
  mutate(across(where(is.numeric), ~round(., 1))) %>% 
  mutate(across("Vars", ~replace(., Vars == "female.N", "N") %>% str_replace("\\.Count", ""))) %>% 
  filter(!str_detect(Vars, "\\.N"))

DemoProp <- NULL
for (i in 2:nrow(DemoTab1)) {
  tmp <- prop.test(as.numeric(DemoTab1[i,2:3]), as.numeric(DemoTab1[1,2:3]))
  
  Props <- tmp$estimate %>% 
    as_tibble() %>% 
    t() %>% 
    `colnames<-`(colnames(DemoTab1[,2:3]))
  DemoProp <- tibble(Vars = as.character(DemoTab1[i,1]),
                     p.value = tmp$p.value) %>% 
    bind_cols(Props) %>% 
    bind_rows(DemoProp, .) %>% 
    relocate(p.value, .after = last_col())
}

DemoProp_f <- DemoProp %>% 
  mutate(across(c("Other", "LOF"), ~round(.*100, 0))) %>% 
  mutate(across(p.value, ~round(., 2)))

Partial <- left_join(DemoTab1, DemoProp_f, by = "Vars") %>% 
  mutate(Other = if_else(is.na(p.value), paste0(Other.x, " (", round(Other.x/nrow(Combined)*100, 0), ")"), paste0(Other.x, " (", Other.y, ")")),
         LOF = if_else(is.na(p.value), paste0(LOF.x, " (", round(LOF.x/nrow(Combined)*100, 0), ")"), paste0(LOF.x, " (", LOF.y, ")"))) %>% 
  dplyr::select(Vars, Other, LOF, p.value) %>% 
  mutate(across(Vars, ~if_else(. == "N", paste0("n = ", nrow(Combined), ", n (%)"), paste0(., ", n (%)"))))

DemoVars2 <- c("age", "prism")
DemoTab2 <- Combined %>% 
  group_by(PG) %>% 
  summarise(across(all_of(DemoVars2), list(
    Median = ~median(., na.rm = TRUE),
    q25 = ~quantile(., 0.25, na.rm = TRUE),
    q75 = ~quantile(., 0.75, na.rm = TRUE)),
    .names = "{col}.{fn}")) %>% 
  pivot_longer(cols = 2:ncol(.), names_to = "Vars", values_to = "Vals") %>% 
  pivot_wider(names_from = PG, values_from = Vals) %>% 
  mutate(across(where(is.numeric), ~round(., 1))) %>% 
  separate(col = Vars, sep = "\\.", into = c("Vars", "Term"))

SumStats <- NULL
for (i in DemoVars2) {
  tmp <- Combined %>% 
    group_by(PG) %>% 
    summarise(Mean = mean(get(i), na.rm = TRUE),
              SD = sd(get(i), na.rm = TRUE),
              Median = median(get(i), na.rm = TRUE),
              q25 = quantile(get(i), probs = 0.25, na.rm = TRUE),
              q75 = round(quantile(get(i), probs = 0.75, na.rm = TRUE), 0),
              Count = sum(!is.na(get(i))),
              SWilk = shapiro.test(get(i))$p.value)
  
  mod <- wilcox.test(get(i) ~ PG, data = Combined) %>% 
    tidy()
  
  SumStats <- tibble(Vars = rep(i, 2),
                     tmp,
                     p.value = rep(round(mod[1,2], 2), 2)) %>% 
    bind_rows(SumStats, .)
}

SumStats_f <- SumStats %>% 
  mutate(across(c("Median", "q25", "q75"), ~if_else(Vars %in% c("age", "highestOI", "durmv28", "piculos"), round(., 1), .))) %>% 
  mutate(Final = paste0(Median, " (", q25, ", ", q75, ")")) %>% 
  pivot_wider(id_cols = c("Vars", "p.value"), names_from = "PG", values_from = "Final") %>% 
  relocate(p.value, .after = LOF) %>% 
  unnest(cols = p.value) %>% 
  mutate(across(Vars, ~paste0(., ", median (IQR)")))


Full <- bind_rows(Partial, SumStats_f)

write_csv(Full, "DemographicTable.csv")

rm(tmp, mod, Props, Partial, SumStats, SumStats_f, list = ls(pattern = "Demo"))


##################################
## Exploratory Outcome Analysis ##
##################################

Genotypes <- colnames(select(Combined, LDLR_rs688:PG))
OutcomesLog <- c(Severity, Outcomes) %>% 
  .[! . %in% c("highestOI", "durmv28", "piculos", "hosplos")]

OutcomesLin <- c(Severity, Outcomes) %>% 
  .[. %in% c("highestOI", "durmv28", "piculos", "hosplos")]

UniReg <- NULL
for (i in Genotypes) {
  for (j in OutcomesLog) {
    tmpform <- as.formula(paste0(j, " ~ ", i))
    tmp <- glm(tmpform, data = Combined, family = "binomial")
    res <- tidy(tmp) %>% 
      mutate(OR = exp(estimate)) %>% 
      bind_cols(exp(confint(tmp))) #%>% 
    #mutate(across(5:8, ~round(., 3)))
    UniReg <- tibble(Var1 = j,
                     Var2 = i,
                     nObs = nobs(tmp),
                     res[2,]) %>% 
      select(-term) %>% 
      bind_rows(UniReg, .)
  }
}

UniReg_NoLDLR <- NULL
for (i in Genotypes) {
  for (j in OutcomesLog) {
    tmpform <- as.formula(paste0(j, " ~ ", i))
    tmp <- glm(tmpform, data = filter(Combined, LDLR_homo == FALSE), family = "binomial")
    res <- tidy(tmp) %>% 
      mutate(OR = exp(estimate)) %>% 
      bind_cols(exp(confint(tmp))) #%>% 
    #mutate(across(5:8, ~round(., 3)))
    UniReg_NoLDLR <- tibble(Var1 = j,
                     Var2 = i,
                     nObs = nobs(tmp),
                     res[2,]) %>% 
      select(-term) %>% 
      bind_rows(UniReg_NoLDLR, .)
  }
}

# Looking at other outcome variables (the continuous).

UniRegLin <- NULL
for (i in Genotypes) {
  for (j in OutcomesLin) {
    tmpform <- as.formula(paste0(j, " ~ ", i))
    tmp <- lm(tmpform, data = Combined)
    res <- tidy(tmp) %>% 
    mutate(across(2:5, ~round(., 3)))
    UniRegLin <- tibble(Var1 = j,
                     Var2 = i,
                     nObs = nobs(tmp),
                     res) %>% 
      bind_rows(UniRegLin, .)
  }
}

UniRegLin_NoLDLR <- NULL
for (i in Genotypes) {
  for (j in OutcomesLin) {
    tmpform <- as.formula(paste0(j, " ~ ", i))
    tmp <- lm(tmpform, data = filter(Combined, LDLR_homo == FALSE))
    res <- tidy(tmp) %>% 
      mutate(across(2:5, ~round(., 3)))
    UniRegLin_NoLDLR <- tibble(Var1 = j,
                        Var2 = i,
                        nObs = nobs(tmp),
                        res) %>% 
      bind_rows(UniRegLin_NoLDLR, .)
  }
}

# Controlling for prism score.

MultiReg_prism <- NULL
for (i in Genotypes) {
  for (j in OutcomesLog) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + prism"))
    tmp <- glm(tmpform, data = Combined, family = "binomial")
    res <- tidy(tmp) %>% 
      mutate(OR = exp(estimate)) %>% 
      bind_cols(exp(confint(tmp))) #%>% 
    #mutate(across(5:8, ~round(., 3)))
    MultiReg_prism <- tibble(Var1 = j,
                     Var2 = i,
                     nObs = nobs(tmp),
                     res[2,]) %>% 
      select(-term) %>% 
      bind_rows(MultiReg_prism, .)
  }
}

MultiReg_NoLDLR_prism <- NULL
for (i in Genotypes) {
  for (j in OutcomesLog) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + prism"))
    tmp <- glm(tmpform, data = filter(Combined, LDLR_homo == FALSE), family = "binomial")
    res <- tidy(tmp) %>% 
      mutate(OR = exp(estimate)) %>% 
      bind_cols(exp(confint(tmp))) #%>% 
    #mutate(across(5:8, ~round(., 3)))
    MultiReg_NoLDLR_prism <- tibble(Var1 = j,
                            Var2 = i,
                            nObs = nobs(tmp),
                            res[2,]) %>% 
      select(-term) %>% 
      bind_rows(MultiReg_NoLDLR_prism, .)
  }
}

MultiRegLin_prism <- NULL
for (i in Genotypes) {
  for (j in OutcomesLin) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + prism"))
    tmp <- lm(tmpform, data = Combined)
    res <- tidy(tmp) %>% 
      mutate(across(2:5, ~round(., 3)))
    MultiRegLin_prism <- tibble(Var1 = j,
                        Var2 = i,
                        nObs = nobs(tmp),
                        res) %>% 
      bind_rows(MultiRegLin_prism, .)
  }
}

MultiRegLin_NoLDLR_prism <- NULL
for (i in Genotypes) {
  for (j in OutcomesLin) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + prism"))
    tmp <- lm(tmpform, data = filter(Combined, LDLR_homo == FALSE))
    res <- tidy(tmp) %>% 
      mutate(across(2:5, ~round(., 3)))
    MultiRegLin_NoLDLR_prism <- tibble(Var1 = j,
                               Var2 = i,
                               nObs = nobs(tmp),
                               res) %>% 
      bind_rows(MultiRegLin_NoLDLR_prism, .)
  }
}

# Controlling for age.

MultiReg_age <- NULL
for (i in Genotypes) {
  for (j in OutcomesLog) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + age"))
    tmp <- glm(tmpform, data = Combined, family = "binomial")
    res <- tidy(tmp) %>% 
      mutate(OR = exp(estimate)) %>% 
      bind_cols(exp(confint(tmp))) #%>% 
    #mutate(across(5:8, ~round(., 3)))
    MultiReg_age <- tibble(Var1 = j,
                             Var2 = i,
                             nObs = nobs(tmp),
                             res[2,]) %>% 
      select(-term) %>% 
      bind_rows(MultiReg_age, .)
  }
}

MultiReg_NoLDLR_age <- NULL
for (i in Genotypes) {
  for (j in OutcomesLog) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + age"))
    tmp <- glm(tmpform, data = filter(Combined, LDLR_homo == FALSE), family = "binomial")
    res <- tidy(tmp) %>% 
      mutate(OR = exp(estimate)) %>% 
      bind_cols(exp(confint(tmp))) #%>% 
    #mutate(across(5:8, ~round(., 3)))
    MultiReg_NoLDLR_age <- tibble(Var1 = j,
                                    Var2 = i,
                                    nObs = nobs(tmp),
                                    res[2,]) %>% 
      select(-term) %>% 
      bind_rows(MultiReg_NoLDLR_age, .)
  }
}

MultiRegLin_age <- NULL
for (i in Genotypes) {
  for (j in OutcomesLin) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + age"))
    tmp <- lm(tmpform, data = Combined)
    res <- tidy(tmp) %>% 
      mutate(across(2:5, ~round(., 3)))
    MultiRegLin_age <- tibble(Var1 = j,
                                Var2 = i,
                                nObs = nobs(tmp),
                                res) %>% 
      bind_rows(MultiRegLin_age, .)
  }
}

MultiRegLin_NoLDLR_age <- NULL
for (i in Genotypes) {
  for (j in OutcomesLin) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + age"))
    tmp <- lm(tmpform, data = filter(Combined, LDLR_homo == FALSE))
    res <- tidy(tmp) %>% 
      mutate(across(2:5, ~round(., 3)))
    MultiRegLin_NoLDLR_age <- tibble(Var1 = j,
                                       Var2 = i,
                                       nObs = nobs(tmp),
                                       res) %>% 
      bind_rows(MultiRegLin_NoLDLR_age, .)
  }
}

# Controlling for gender.

MultiReg_gender <- NULL
for (i in Genotypes) {
  for (j in OutcomesLog) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + gender"))
    tmp <- glm(tmpform, data = Combined, family = "binomial")
    res <- tidy(tmp) %>% 
      mutate(OR = exp(estimate)) %>% 
      bind_cols(exp(confint(tmp))) #%>% 
    #mutate(across(5:8, ~round(., 3)))
    MultiReg_gender <- tibble(Var1 = j,
                           Var2 = i,
                           nObs = nobs(tmp),
                           res[2,]) %>% 
      select(-term) %>% 
      bind_rows(MultiReg_gender, .)
  }
}

MultiReg_NoLDLR_gender <- NULL
for (i in Genotypes) {
  for (j in OutcomesLog) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + gender"))
    tmp <- glm(tmpform, data = filter(Combined, LDLR_homo == FALSE), family = "binomial")
    res <- tidy(tmp) %>% 
      mutate(OR = exp(estimate)) %>% 
      bind_cols(exp(confint(tmp))) #%>% 
    #mutate(across(5:8, ~round(., 3)))
    MultiReg_NoLDLR_gender <- tibble(Var1 = j,
                                  Var2 = i,
                                  nObs = nobs(tmp),
                                  res[2,]) %>% 
      select(-term) %>% 
      bind_rows(MultiReg_NoLDLR_gender, .)
  }
}

MultiRegLin_gender <- NULL
for (i in Genotypes) {
  for (j in OutcomesLin) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + gender"))
    tmp <- lm(tmpform, data = Combined)
    res <- tidy(tmp) %>% 
      mutate(across(2:5, ~round(., 3)))
    MultiRegLin_gender <- tibble(Var1 = j,
                              Var2 = i,
                              nObs = nobs(tmp),
                              res) %>% 
      bind_rows(MultiRegLin_gender, .)
  }
}

MultiRegLin_NoLDLR_gender <- NULL
for (i in Genotypes) {
  for (j in OutcomesLin) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + gender"))
    tmp <- lm(tmpform, data = filter(Combined, LDLR_homo == FALSE))
    res <- tidy(tmp) %>% 
      mutate(across(2:5, ~round(., 3)))
    MultiRegLin_NoLDLR_gender <- tibble(Var1 = j,
                                     Var2 = i,
                                     nObs = nobs(tmp),
                                     res) %>% 
      bind_rows(MultiRegLin_NoLDLR_gender, .)
  }
}

# Controlling for all 3 together!

MultiReg_all <- NULL
for (i in Genotypes) {
  for (j in OutcomesLog) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + prism + age + gender"))
    tmp <- glm(tmpform, data = Combined, family = "binomial")
    res <- tidy(tmp) %>% 
      mutate(OR = exp(estimate)) %>% 
      bind_cols(exp(confint(tmp))) #%>% 
    #mutate(across(5:8, ~round(., 3)))
    MultiReg_all <- tibble(Var1 = j,
                              Var2 = i,
                              nObs = nobs(tmp),
                              res[2,]) %>% 
      select(-term) %>% 
      bind_rows(MultiReg_all, .)
  }
}

MultiReg_NoLDLR_all <- NULL
for (i in Genotypes) {
  for (j in OutcomesLog) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + prism + age + gender"))
    tmp <- glm(tmpform, data = filter(Combined, LDLR_homo == FALSE), family = "binomial")
    res <- tidy(tmp) %>% 
      mutate(OR = exp(estimate)) %>% 
      bind_cols(exp(confint(tmp))) #%>% 
    #mutate(across(5:8, ~round(., 3)))
    MultiReg_NoLDLR_all <- tibble(Var1 = j,
                                     Var2 = i,
                                     nObs = nobs(tmp),
                                     res[2,]) %>% 
      select(-term) %>% 
      bind_rows(MultiReg_NoLDLR_all, .)
  }
}


MultiRegLin_all <- NULL
for (i in Genotypes) {
  for (j in OutcomesLin) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + prism + age + gender"))
    tmp <- lm(tmpform, data = Combined)
    res <- tidy(tmp) %>% 
      mutate(across(2:5, ~round(., 3)))
    MultiRegLin_all <- tibble(Var1 = j,
                                 Var2 = i,
                                 nObs = nobs(tmp),
                                 res) %>% 
      bind_rows(MultiRegLin_all, .)
  }
}

MultiRegLin_NoLDLR_all <- NULL
for (i in Genotypes) {
  for (j in OutcomesLin) {
    tmpform <- as.formula(paste0(j, " ~ ", i, " + prism + age + gender"))
    tmp <- lm(tmpform, data = filter(Combined, LDLR_homo == FALSE))
    res <- tidy(tmp) %>% 
      mutate(across(2:5, ~round(., 3)))
    MultiRegLin_NoLDLR_all <- tibble(Var1 = j,
                                        Var2 = i,
                                        nObs = nobs(tmp),
                                        res) %>% 
      bind_rows(MultiRegLin_NoLDLR_all, .)
  }
}

for (i in ls(pattern = "Uni|Multi")) {
  get(i) %>% 
    mutate(fdr.p = p.adjust(p = p.value, method = "fdr")) %>% 
    relocate(fdr.p, .after = p.value) %>% 
    arrange(fdr.p) %>% 
    assign(i, value = ., envir = .GlobalEnv) %>% 
    write_csv(paste0("./RegressionOutput/", i, ".csv"))
}


Counts <- Combined %>% 
  group_by(PCSK9_LOF) %>% 
  summarise(Count = n())

IL8MEd <- Combined %>% 
  group_by(PCSK9_LOF) %>% 
  summarise(IL8max = median(IL_8_Max, na.rm = TRUE),
            ILmed = median(IL_8_Median, na.rm = TRUE))

TMMEd <- Combined %>% 
  group_by(PCSK9_LOF) %>% 
  summarise(TMmax = median(thrbm_Max, na.rm = TRUE),
          TMmed = median(thrbm_Median, na.rm = TRUE))

Combined %>% 
  filter(!LDLR_homo == TRUE) %>% 
  ggplot(aes(x = PCSK9_LOF, y = log(thrbm_Max), fill = PCSK9_LOF)) +
  geom_violin() +
  geom_point(position = position_jitter(width = 0.2))
