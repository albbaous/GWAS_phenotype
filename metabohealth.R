# Install packages if they are not already installed
install.packages("mice")

# Load required libraries
library(tidyverse)
library(tibble)
library(mice)
library(jsonlite)
library(tidyr)

setwd('/Users/user')

# Read in metabolite data which also includes sex, age, bmi 
df <- read_csv("cohort_data2.csv")

# Convert to date properly (only if it's numeric!)
df$`participant.p42018` <- as.Date(df$`participant.p42018`, origin = "1970-01-01")

# Define biomarkers and weights from Deelen et al., 2019
biomarkers <- tribble(
  ~column,                    ~label,
  "participant.eid",          "eid",
  "participant.p21003_i0",    "Age",
  "participant.p31",          "Sex",
  "participant.p42018",       "Date of dementia diagnosis",
  "participant.p42019",       "Source of all cause dementia",
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
  "participant.p23479_i0",    "Alb",
  "participant.p23480_i0",    "GlycA",
  "participant.p23453_i0",    "PUFA_FA",
  "participant.p23482_i0",    "XXL_VLDL_L",
  "participant.p23573_i0",    "S_HDL_L",
  "participant.p23431_i0",    "VLDL_D"
)

df <- df %>%
  rename_with(~ biomarkers$label[match(.x, biomarkers$column)], .cols = everything())

# -------------------------------
# Initial sample size
# -------------------------------
cat("Number of individuals before all filtering:", nrow(df), "\n")

# -------------------------------
# Filter: remove rows with NA in any of the 10 PC columns
# -------------------------------
df <- df %>%
  filter(
    !is.na(PC1) & !is.na(PC2) & !is.na(PC3) & !is.na(PC4) & 
      !is.na(PC5) & !is.na(PC6) & !is.na(PC7) & !is.na(PC8) & 
      !is.na(PC9) & !is.na(PC10)
  )

# -------------------------------
# Define metabolite columns and filter for missingness
# -------------------------------
metabolites <- c("Glc", "Lac", "His", "Ile", "Leu", "Val", "Phe", "AcAce", "Alb", 
                 "GlycA", "PUFA_FA", "XXL_VLDL_L", "S_HDL_L", "VLDL_D")

df <- df %>%
  filter(
    rowSums(is.na(df[metabolites])) < length(metabolites) & 
      rowSums(is.na(df[setdiff(metabolites, "Alb")])) < (length(metabolites) - 1)
  )


# -------------------------------
# Final sample size
# -------------------------------
cat("Number of individuals after all filtering:", nrow(df), "\n")

# -------------------------------
# Missing summary for metabolites
# -------------------------------
missing_summary_post_filter <- df %>%
  summarise(across(
    all_of(metabolites),
    list(
      missing_pct = ~ mean(is.na(.)) * 100,
      missing_count = ~ sum(is.na(.))
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("Metabolite", "Metric"),
    names_pattern = "^(.*)_(missing_pct|missing_count)$"
  ) %>%
  pivot_wider(names_from = Metric, values_from = value) %>%
  arrange(desc(missing_pct))

# Print updated missing summary
print(missing_summary_post_filter)

# -------------------------------
# Imputation
# -------------------------------
# Prepare the data for imputation (Age, Sex, and metabolites)
df_imputation <- df %>%
  select(c("Age", "Sex", all_of(metabolites)))

# Create the predictor matrix for 'mice' (1 means predictor, 0 means not predictor)
predictor_matrix <- matrix(1, nrow = ncol(df_imputation), ncol = ncol(df_imputation))
colnames(predictor_matrix) <- colnames(df_imputation)
rownames(predictor_matrix) <- colnames(df_imputation)

# Exclude metabolites as predictors for each other (set to 0)
predictor_matrix[which(colnames(predictor_matrix) %in% metabolites), which(colnames(predictor_matrix) %in% metabolites)] <- 0

# Perform imputation using the 'mice' function
imputed_data <- mice(df_imputation, method = 'rf', m = 5, predictorMatrix = predictor_matrix, verbose = FALSE)

# Extract the imputed data (choose the first imputed dataset)
df_imputed <- complete(imputed_data, action = 1)

# Clean column names to avoid renaming issues after imputation (like Age...2, Sex...3)
colnames(df_imputed) <- gsub("\\.\\.\\d+$", "", colnames(df_imputed))

# Combine the imputed data back with the rest of the dataframe (excluding metabolites in the original df)
df <- df %>%
  select(-all_of(metabolites)) %>%
  bind_cols(df_imputed)

# Analysis as defined in Deelen et al., 2019 after clearing and imputing data 
# Add 1 to any 0 values (this is the preprocessing step)
df_adjusted <- df %>%
  mutate(
    across(
      all_of(metabolites),
      ~ ifelse(. == 0, . + 1, .),
      .names = "adjusted_{.col}"
    )
  )

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
    weight_XXL_VLDL_L = 0.80,
    weight_S_HDL_L    = 0.87,
    weight_VLDL_D     = 0.85,
    weight_PUFA_FA    = 0.78,
    weight_Glc        = 1.16,
    weight_Lac        = 1.06,
    weight_His        = 0.93,
    weight_Ile        = 1.23,
    weight_Leu        = 0.82,
    weight_Val        = 0.87,
    weight_Phe        = 1.13,
    weight_AcAce      = 1.08,
    weight_Alb        = 0.89,
    weight_GlycA      = 1.32
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

# Clean up column names in df_extended
df_extended2 <- df_extended %>%
  # Rename columns Age...4, Sex...5, to Age, Sex, and BMI
  rename(
    Age = `Age...4`,
    Sex = `Sex...5`
  ) %>%
  # Remove the unwanted columns like BMI...15, Age...16, Sex...17
  select(-`Age...16`, -`Sex...17`)

# Create final dataframe with selected columns in the specified order
df_final <- df_extended2 %>%
  select(eid,                    # Participant ID
         Age,                     # Age
         Sex,                     # Sex
         starts_with("PC"),       # Principal components PC1–PC10
         all_of(metabolites),     # Original metabolite values
         MetaboHealth_Score)      # Final score

# Create PHENOTYPE file '.pheno'
df_phen <- df_final %>%
  # Rename 'eid' to 'IID'
  rename(IID = eid) %>%
  # Add 'FID' column with values same as 'IID' (i.e., both FID and IID will be the same)
  mutate(FID = IID) %>%
  # Select the desired columns in the specified order
  select(FID, IID, MetaboHealth_Score)

# View the final phenotype dataframe
head(df_phen)

# Save the result to a file (e.g., tab-delimited .pheno)
write.table(df_phen, "ukb_phenotype_data.pheno", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Create COVARIATE file '.cov'
df_cov <- df_final %>%
  select(eid, Age, Sex, starts_with("PC")) %>%
  # Rename 'eid' to 'IID'
  rename(IID = eid) %>%
  # Add 'FID' column with values same as 'IID' (i.e., both FID and IID will be the same)
  mutate(FID = IID) %>%
  # Select the desired columns in the specified order
  select(FID, IID, everything())

# View the final covariate dataframe
head(df_cov)

# Save the result to a file (e.g., tab-delimited .cov)
write.table(df_cov, "ukb_covariates.cov", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

