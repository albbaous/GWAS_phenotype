# GWAS_prep
A simple pipeline for extracting my phenotype — all the way through to running my GWAS. This is MetaboHealth as defined in Deelen et al., 2019

# Phenotype
This is a breakdown of the logic behind the script titled metabohealth.R

## Step 1 — Extracting Metabolites for Phenotype
To get the correct data from **auth.dnanexus.com** (UK Biobank Research Analysis Platform), you need to use the `dx extract_dataset` command. 

### What’s Extracted:
- **14 metabolites** listed in the paper by *Deelen et al., 2019* (see [Supplement, page 30](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11311-9/MediaObjects/41467_2019_11311_MOESM1_ESM.pdf)):
- Covariates:
  - Age (`21003`)
  - Sex (`31`)
  - 10 Genetic PCs (`22009_a1` to `22009_a10`)
  - ICD-10 codes (`41202`, `41204`)
  - Covariates in `.fam` files are Family ID (FIID), Individual ID (IID),  paternal ID, maternal ID and phenotype.

```bash
dx extract_dataset project-Gzyb0j8JQYbBjQxYqfb4xJYX:record-GzyfX70Jfj0bvy8YfvYQ302v --fields participant.eid,participant.p41202,participant.p41204,participant.p21003_i0,participant.p31,participant.p22009_a1,participant.p22009_a2,participant.p22009_a3,participant.p22009_a4,participant.p22009_a5,participant.p22009_a6,participant.p22009_a7,participant.p22009_a8,participant.p22009_a9,participant.p22009_a10,participant.p23470_i0,participant.p23471_i0,participant.p23463_i0,participant.p23465_i0,participant.p23466_i0,participant.p23467_i0,participant.p23468_i0,participant.p23476_i0,participant.p30600_i0,participant.p23480_i0,participant.p23453_i0,participant.p23482_i0,participant.p23573_i0,participant.p23431_i0 -o cohort_data.csv
```

> ⚠️ **Note**:  
>   
> Also retrieves ICD-10 diagnoses:
> - `41202`: main diagnosis
> - `41204`: secondary diagnosis  
> This matches methods in [BMJ Mental Health 2023](https://pmc.ncbi.nlm.nih.gov/articles/instance/10577770/bin/bmjment-2023-300719supp001.pdf)
>
> In that paper, they also included any **medications associated with dementia** to help define cases. However, I am *not* doing this, as many dementia medications are also used for **other indications** (e.g., antipsychotics or stimulants like **modafinil**), which may **skew results**.  

### dx Command Explanation

- `dx extract_dataset`: DNAnexus CLI command for extracting a dataset from RAP  
- `project-...`: The ID of your RAP project, extracted using:
```
project = os.getenv('DX_PROJECT_CONTEXT_ID')
project
```
- `record-...`: The file record (metabolite dataset), extracted using:
```
record = os.popen("dx find data --type Dataset --delimiter ',' | awk -F ',' '{print $5}'").read().rstrip()
record
```
- `--fields`: List of fields to extract (in `participant.p[FIELD_ID]_i0` format)  
- `-o`: Specifies the output CSV file  

> ⚠️ **Note**: You can check values by mirroring/building the same cohort on UKB RAP, selecting a participant and seeing if the values are the same.
  - i.e, `grep '100010' cohort2.csv` and check they are all the same as column names in UKB rap have metabolite name
---

### Step 2 (just notes)— Exclude Dementia Diagnoses and pick filters

We exclude individuals with these ICD-10 codes (described in `icd10_meaning`) related to dementia:

```
A810, F00, F000, F001, F002, F009, F01, F010, F011, F012, F013, F018, F019,
F02, F020, F021, F022, F023, F024, F028, F03, F051, F106, G30, G300, G301,
G308, G309, G310, G311, G318, I673
```

Focus is only on **hospital-diagnosed cases**, not medications.

- Suggestions have been made to filter on metabolic syndrome.
  
**Metabolic syndrome is defined as:**
Metabolic syndrome is a group of conditions that increase the risk of heart disease, stroke and type 2 diabetes. These conditions include high blood pressure, high blood sugar, too much fat around the waist, and high cholesterol or triglyceride levels.
https://www.mayoclinic.org/diseases-conditions/metabolic-syndrome/symptoms-causes/syc-20351916

> ⚠️ **Note** These are some things you can do but with imputation and a diff focus this is no longer NEEDED:
>
> Also extract features for metabolic syndrome i.e., heart disease, stroke, diabetes and cancer - removing the unhealthy lot as theres a lot of evidence that they may skew the results i.e., https://pmc.ncbi.nlm.nih.gov/articles/PMC8504077/
>
> Then clean data to remove N/A values and to remove unhealthy (so those with the above diseases)
>
> ❌ Do NOT remove:
> - Missing data → instead, **impute**
> - Unhealthy individuals → **include all**
> - Relatives → remove **later** if needed
---
### Step 3 (on local R script) — Add weights and make score 
For each biomarker value, the following steps are applied:

---
#### 1. Handle Zero Values

If the raw observed value is 0, add 1 to make it compatible with log-transformation

#### 2. Log Transformation

#### 3. Standardization (Z-score Scaling)

Scale to standard deviation units (mean = 0, sd = 1)

#### 4.Add weights from Deelen et al., paper to get score

---

#### **Run the R script saved here as `metabohealth.R`** 
- This does each step by step so I can see the resulting columns and then multiply them by each other
- this creates a file phenotype dataframe at the end whereby we have all the metabolites, characteristics and then a score - this is not in the phenotype format i need for GWAS yet as certain values need to be imputed in original data to avoid scores of 0
- Towards the end it might look like `Alb` is missing at an alarming rate but it is oly 12%
  
---
### Step 4 (on local R script)— Imputation using Age, and Sex as predictors 
- Followed this lovely tutorial: https://libguides.princeton.edu/R-Missingdata
- This is all in 'metabohealth.R'

**Explanation of the Parameters:**

- `method`: The imputation method used (in this case, pmm stands for predictive mean matching and can alternatively use rf which is random forest).

- `m`: The number of imputed datasets (typically 5 or more).

---
### Step 5 (on local R script)— GWAS format 
- The `.fam` file in UKB has the following columns which correspond to Family ID (FIID), Individual ID (IID),  paternal ID, maternal ID and phenotype. 

```
==> ukb22418_c9_b0_v2.fam <==
4459327 4459327 0 0 1 Batch_b001
```

- `.pheno` file needs to be made in the same format (code for this in lines 273-286 in `metabohealth_imputed.R`), whereby phenotype is not just a batch number, but our phenotype score. It looks like this: 

```
FID	IID	Age	Sex	MetaboHealth_Score
1000073	1000073	55	1	0.120481033803025
```


- `.cov` file also needs to be made in the same format (code for this in lines 289-297 in `metabohealth_imputed.R`). It looks like this: 

```
FID	IID	Age	Sex	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10
1000073	1000073	55	1	-12.5122	2.52381	-1.68065	1.15155-3.75587	-1.52478	-0.271479	-0.655341	-0.848979	0.933165
```
### Step 6 (on UKB RAP) — Running GWAS on subset of data 
- On UKB Rap in Jupyter notebooks, I opened up the Terminal and ran the following:
```
conda install -c bioconda plink2 -y
```
```
mkdir -p /mnt/project/plink_output
```
```
plink2 \
  --bfile "/mnt/project/Bulk/Genotype Results/Genotype calls/ukb22418_c22_b0_v2" \
  --pheno /mnt/project/ukb_phenotype_data.pheno \
  --pheno-name MetaboHealth_Score \
  --covar /mnt/project/ukb_covariates.cov \
  --covar-name Age,Sex,PC1-PC10 \
  --glm \
  --thin 0.05 \
  --out /mnt/project/plink_output/metabo_analysis

```
- `--thin` is to pick out 0.05% of the chrom 22 data so i am testing on only a subset
- this produces `metabo_analysis.log` and `metabo_analysis.MetaboHealth_Score.glm.linear`
---

## EXTRA - IGNORE 
Just saving these here as it is the UKB IDs mapped to the names of variables in the Deelen paper: 
```
biomarkers <- tribble(
  ~column,                    ~label,
  "participant.eid",          "eid",
  "participant.p21003_i0",    "Age",
  "participant.p31",          "Sex",
  "participant.p41202",       "ICD10_Main",
  "participant.p41204",       "ICD10_Secondary",
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
  "participant.p23453_i0",    "PUFA_FA",
  "participant.p23482_i0",    "XXL_VLDL_L",
  "participant.p23573_i0",    "S_HDL_L",
  "participant.p23431_i0",    "VLDL_D"
)
```
I got the mappings of what they come up as in UKB from the table in `ukb_names` - which I got from a paper i cant remember the exact name of
