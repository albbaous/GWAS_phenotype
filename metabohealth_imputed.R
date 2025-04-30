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

# Define the metabolite columns (from the biomarkers table)
metabolites <- c("Glc", "Lac", "His", "Ile", "Leu", "Val", "Phe", "AcAce", "Alb", 
                 "GlycA", "Cit", "PUFA_FA", "XXL_VLDL_L", "S_HDL_L", "VLDL_D")

# Prepare the data for imputation (BMI, Age, Sex, and metabolites)
df_imputation <- df %>%
  select(c("BMI", "Age", "Sex", all_of(metabolites)))

# Create the predictor matrix for 'mice' (1 means predictor, 0 means not predictor)
predictor_matrix <- matrix(1, nrow = ncol(df_imputation), ncol = ncol(df_imputation))
colnames(predictor_matrix) <- colnames(df_imputation)
rownames(predictor_matrix) <- colnames(df_imputation)

# Exclude metabolites as predictors for each other (set to 0)
predictor_matrix[which(colnames(predictor_matrix) %in% metabolites), which(colnames(predictor_matrix) %in% metabolites)] <- 0

# Perform imputation using the 'mice' function
imputed_data <- mice(df_imputation, method = 'pmm', m = 5, predictorMatrix = predictor_matrix, verbose = FALSE)

# Extract the imputed data (choose the first imputed dataset)
df_imputed <- complete(imputed_data, action = 1)

# Clean column names to avoid renaming issues after imputation (like Age...2, Sex...3)
colnames(df_imputed) <- gsub("\\.\\.\\d+$", "", colnames(df_imputed))

# Combine the imputed data back with the rest of the dataframe (excluding metabolites in the original df)
df <- df %>%
  select(-all_of(metabolites)) %>%
  bind_cols(df_imputed)

# Create a new dataframe with log-transformed values in new columns
df_log <- df %>%
  mutate(
    across(all_of(metabolites), log, .names = "log_{.col}")  # Create new columns with log-transformed metabolite values
  )

# View the new dataframe with log-transformed metabolite values
head(df_log)

# Create a new dataframe that starts with df_log
df_extended <- df_log %>%
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
    Z_Glc        = scale(log_Glc        * weight_Glc),
    Z_Lac        = scale(log_Lac        * weight_Lac),
    Z_His        = scale(log_His        * weight_His),
    Z_Ile        = scale(log_Ile        * weight_Ile),
    Z_Leu        = scale(log_Leu        * weight_Leu),
    Z_Val        = scale(log_Val        * weight_Val),
    Z_Phe        = scale(log_Phe        * weight_Phe),
    Z_AcAce      = scale(log_AcAce      * weight_AcAce),
    Z_Alb        = scale(log_Alb        * weight_Alb),
    Z_GlycA      = scale(log_GlycA      * weight_GlycA),
    Z_Cit        = scale(log_Cit),  # No weight provided, so just standardize the log value
    Z_PUFA_FA    = scale(log_PUFA_FA    * weight_PUFA_FA),
    Z_XXL_VLDL_L = scale(log_XXL_VLDL_L * weight_XXL_VLDL_L),
    Z_S_HDL_L    = scale(log_S_HDL_L    * weight_S_HDL_L),
    Z_VLDL_D     = scale(log_VLDL_D     * weight_VLDL_D)
  )

# Add MetaboHealth_Score by summing all Z_ columns
df_extended <- df_extended %>%
  mutate(
    MetaboHealth_Score = rowSums(select(., starts_with("Z_")), na.rm = TRUE)
  )

# Clean up column names in df_extended
df_extended2 <- df_extended %>%
  # Rename columns Age...2, Sex...3, and BMI...4 to Age, Sex, and BMI
  rename(
    Age = `Age...2`,
    Sex = `Sex...3`,
    BMI = `BMI...4`
  ) %>%
  # Remove the unwanted columns like BMI...15, Age...16, Sex...17
  select(-`BMI...15`, -`Age...16`, -`Sex...17`)

# Create final dataframe with selected columns in the specified order
df_final <- df_extended2 %>%
  select(eid,                    # Participant ID
         Age,                     # Age
         Sex,                     # Sex
         BMI,                     # BMI
         starts_with("PC"),       # Principal components PC1–PC10
         all_of(metabolites),     # Original metabolite values
         MetaboHealth_Score)      # Final score

# Filter out rows where MetaboHealth_Score is 0 and create a new dataframe
df_no_zero_score <- df_final %>%
  filter(MetaboHealth_Score != 0)

# Print the number of rows in df_final and df_no_zero_score
cat("Number of rows in df_final: ", nrow(df_final), "\n")
cat("Number of rows in df_no_zero_score: ", nrow(df_no_zero_score), "\n")
