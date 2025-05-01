# Load required libraries
install.packages("mice")

library(tidyverse)
library(tibble)
library(mice)

setwd('/Users/user')

# Read in metabolite data which also includes sex, age, bmi 
df <- read_csv("cohort_data6.csv")

# Define biomarkers and weights from Deelen et al., 2019
biomarkers <- tribble(
  ~column,                    ~label,
  "participant.eid",          "eid",
  "participant.p21003_i0",    "Age",
  "participant.p31",          "Sex",
  "participant.p21001_i0",    "BMI",
  "participant.p22009_a1",    "PC1",
  "participant.p22009_a2",    "PC2",
  "participant.p22009_a3",    "PC3",
  "participant.p22009_a4",    "PC4",
  "participant.p22009_a5",    "PC5",
  "participant.p22009_a6",    "PC6",
  "participant.p22009_a7",    "PC7",
  "participant.p22009_a8",    "PC8",
  "participant.p22009_a9",    "PC9",
  "participant.p22009_a10",   "PC10",
  "participant.p23470_i0",    "Glc",
  "participant.p23471_i0",    "Lac",
  "participant.p23463_i0",    "His",
  "participant.p23465_i0",    "Ile",
  "participant.p23466_i0",    "Leu",
  "participant.p23467_i0",    "Val",
  "participant.p23468_i0",    "Phe",
  "participant.p23476_i0",    "AcAce",
  "participant.p30600_i0",    "Alb",
  "participant.p23480_i0",    "GlycA",
  "participant.p23473_i0",    "Cit",
  "participant.p23453_i0",    "PUFA_FA",
  "participant.p23482_i0",    "XXL_VLDL_L",
  "participant.p23573_i0",    "S_HDL_L",
  "participant.p23431_i0",    "VLDL_D"
)

# Rename the columns based on the biomarkers table
df <- df %>%
  rename_with(~ biomarkers$label[match(.x, biomarkers$column)], .cols = everything())

# Remove the BMI column - justification for this https://www.nature.com/articles/s41467-022-29143-5
df <- df %>%
  select(-BMI)

# Remove the Cit column - appears in the paper but has no weight 
df <- df %>%
  select(-Cit)

# View the first few rows to confirm
head(df)

# Remove all those that have NA values in the 10 PC columns 
# Remove rows that have NA values in the 10 PC columns (PC1 to PC10)
df <- df %>%
  filter(
    !is.na(PC1) & !is.na(PC2) & !is.na(PC3) & !is.na(PC4) & 
      !is.na(PC5) & !is.na(PC6) & !is.na(PC7) & !is.na(PC8) & 
      !is.na(PC9) & !is.na(PC10)
  )

# View the first few rows to confirm the rows are removed
head(df)


# Define the metabolite columns (from the biomarkers table)
metabolites <- c("Glc", "Lac", "His", "Ile", "Leu", "Val", "Phe", "AcAce", "Alb", 
                 "GlycA", "PUFA_FA", "XXL_VLDL_L", "S_HDL_L", "VLDL_D")

# Add 1 to any 0 values (this is the preprocessing step)
df_adjusted <- df %>%
  mutate(
    across(
      all_of(metabolites),
      ~ ifelse(. == 0, . + 1, .),
      .names = "adjusted_{.col}"
    )
  )

# Check missingness and zeros for each metabolite
metabolite_summary <- df %>%
  summarise(across(
    all_of(metabolites),
    list(
      missing_pct = ~ mean(is.na(.)) * 100,
      zero_pct = ~ mean(. == 0, na.rm = TRUE) * 100
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("Metabolite", "Metric"),
    names_pattern = "^(.*)_(missing_pct|zero_pct)$"
  ) %>%
  pivot_wider(names_from = Metric, values_from = value) %>%
  arrange(desc(missing_pct))

print(metabolite_summary)


# Log-transform the adjusted values (this is the actual log transformation)
df_log <- df_adjusted %>%
  mutate(
    across(
      starts_with("adjusted_"),
      ~ log(.),
      .names = "log_{.col}"
    )
  )

# View the new dataframe with log-transformed metabolite values
head(df_log)

# Standardize the log-transformed columns
df_scaled <- df_log %>%
  mutate(across(starts_with("log_"), scale, .names = "z_{.col}"))

# Create a new dataframe that starts with df_log
df_extended <- df_scaled %>%
  # For each metabolite in the weights list, create a new column with the weight value
  mutate(
    weight_XXL_VLDL_L = log(0.80),
    weight_S_HDL_L    = log(0.87),
    weight_VLDL_D     = log(0.85),
    weight_PUFA_FA    = log(0.78),
    weight_Glc        = log(1.16),
    weight_Lac        = log(1.06),
    weight_His        = log(0.93),
    weight_Ile        = log(1.23),
    weight_Leu        = log(0.82),
    weight_Val        = log(0.87),
    weight_Phe        = log(1.13),
    weight_AcAce      = log(1.08),
    weight_Alb        = log(0.89),
    weight_GlycA      = log(1.32)
  )

# View the new dataframe with log-transformed metabolite values
head(df_extended)
head(df_log)

# Create Z-score columns: multiply log_metabolite by weight and standardize 
df_extended <- df_extended %>%
  mutate(
    MH_Glc        = (scale(log_adjusted_Glc))        * (weight_Glc),
    MH_Lac        = (scale(log_adjusted_Lac))        * (weight_Lac),
    MH_His        = (scale(log_adjusted_His))        * (weight_His),
    MH_Ile        = (scale(log_adjusted_Ile))        * (weight_Ile),
    MH_Leu        = (scale(log_adjusted_Leu))        * (weight_Leu),
    MH_Val        = (scale(log_adjusted_Val))        * (weight_Val),
    MH_Phe        = (scale(log_adjusted_Phe))        * (weight_Phe),
    MH_AcAce      = (scale(log_adjusted_AcAce))      * (weight_AcAce),
    MH_Alb        = (scale(log_adjusted_Alb))        * (weight_Alb),
    MH_GlycA      = (scale(log_adjusted_GlycA))      * (weight_GlycA),
    MH_PUFA_FA    = (scale(log_adjusted_PUFA_FA))    * (weight_PUFA_FA),
    MH_XXL_VLDL_L = (scale(log_adjusted_XXL_VLDL_L)) * (weight_XXL_VLDL_L),
    MH_S_HDL_L    = (scale(log_adjusted_S_HDL_L))    * (weight_S_HDL_L),
    MH_VLDL_D     = (scale(log_adjusted_VLDL_D))     * (weight_VLDL_D)
  )

# Add MetaboHealth_Score by summing all MH_ (metabohealth score) columns
df_extended <- df_extended %>%
  mutate(
    MetaboHealth_Score = rowSums(select(., starts_with("MH_")), na.rm = TRUE)
  )
