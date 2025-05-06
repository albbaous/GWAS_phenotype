# GWAS_prep
A simple pipeline for extracting my phenotype — all the way through to running my GWAS. This is MetaboHealth as defined in Deelen et al., 2019

# Phenotype
## Step 1 — Extracting Metabolites for Phenotype

To get the correct data from **auth.dnanexus.com** (UK Biobank Research Analysis Platform), you need to use the following command. This extracts all **14 metabolites** listed in the paper by *Deelen et al., 2019* (see [Supplement, page 30](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11311-9/MediaObjects/41467_2019_11311_MOESM1_ESM.pdf)):

```bash
dx extract_dataset \
  project-Gzyb0j8JQYbBjQxYqfb4xJYX:record-GzyfX70Jfj0bvy8YfvYQ302v \
  --fields participant.eid,participant.p23470_i0,participant.p23471_i0,participant.p23463_i0,participant.p23465_i0,participant.p23466_i0,participant.p23467_i0,participant.p23468_i0,participant.p23476_i0,participant.p30600_i0,participant.p23480_i0,participant.p23473_i0,participant.p23453_i0,participant.p23482_i0,participant.p23573_i0,participant.p23431_i0 \
  -o cohort_data2.csv
```
> ⚠️ **Note**: This command also retrieves citrate (`participant.p23473_i0`) which is listed in the paper but has no weight so i remove IT later on 

### Explanation

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

> ⚠️ **Note**: This command only retrieves metabolite values.  
> You will also need to extract baseline characteristics like age, sex, and BMI to complete phenotype mapping
> You can check values by mirroring/building the same cohort on UKB RAP, selecting a participant and seeing if the values are the same.
  - i.e, `grep '100010' cohort2.csv` and check they are all the same as column names in UKB rap have metabolite name
---

### Step 2 — Adding some baseline characteristics/covariates to the command and proceeding to filter on those

Once you've extracted the metabolite data, the next step is to map participants to those on the UoM db
#### Actions

1. **Add and extract additional covariates** from UK Biobank
   - Age (`21003`)
   - Sex (`31`)
   - 10 Principal Components (`22009_a1`) - all the way to 10

```bash
dx extract_dataset project-Gzyb0j8JQYbBjQxYqfb4xJYX:record-GzyfX70Jfj0bvy8YfvYQ302v --fields participant.eid,participant.p21003_i0,participant.p31,participant.p21001_i0,participant.p22009_a1,participant.p22009_a2,participant.p22009_a3,participant.p22009_a4,participant.p22009_a5,participant.p22009_a6,participant.p22009_a7,participant.p22009_a8,participant.p22009_a9,participant.p22009_a10,participant.p23470_i0,participant.p23471_i0,participant.p23463_i0,participant.p23465_i0,participant.p23466_i0,participant.p23467_i0,participant.p23468_i0,participant.p23476_i0,participant.p30600_i0,participant.p23480_i0,participant.p23473_i0,participant.p23453_i0,participant.p23482_i0,participant.p23573_i0,participant.p23431_i0 -o cohort_data6.csv
```
> ⚠️ **Note**: This command also retrieves BMI (`participant.p21001_i0`) which I remove later on as it does not appear much in GWAS 
- Covariates in `.fam` files are Family ID (FIID), Individual ID (IID),  paternal ID, maternal ID and phenotype.

**Metabolic syndrome is defined as:**
Metabolic syndrome is a group of conditions that increase the risk of heart disease, stroke and type 2 diabetes. These conditions include high blood pressure, high blood sugar, too much fat around the waist, and high cholesterol or triglyceride levels.
https://www.mayoclinic.org/diseases-conditions/metabolic-syndrome/symptoms-causes/syc-20351916

> ⚠️ **Note** These are some things you can do but with imputation and a diff focus this is no longer NEEDED: 
> Also extract features for metabolic syndrome i.e., heart disease, stroke, diabetes and cancer - removing the unhealthy lot as theres a lot of evidence that they may skew the results i.e., https://pmc.ncbi.nlm.nih.gov/articles/PMC8504077/
> Then clean data to remove N/A values and to remove unhealthy (so those with the above diseases)

> **Instead**
> DO NOT remove NAs- impute them
> DO NOT remove unhealthy - focus on everyone
> DO NOT remove relatives - do later 
---

### Step 3 — Add weights and make score
For each biomarker value, the following steps are applied:
#### 1. Handle Zero Values

If the raw observed value is 0, add 1 to make it compatible with log-transformation

#### 2. Log Transformation

#### 3. Standardization (Z-score Scaling)

Scale to standard deviation units (mean = 0, sd = 1)

#### 4.Add weights from Deelen et al., paper to get score
It has alreadt been calculated in van Holstein et al., 2024 using this 

```
MetaboHealth = (((Z(ln[XXL_VLDL_L]))*ln(0.80)) + ((Z(ln[S_HDL_L]))*ln(0.87)) + ((Z(ln[VLDL-D]))*ln(0.85)) + ((Z(ln[PUFA/FA]))*ln(0.78)) + ((Z(ln[Glucose]))*ln(1.16)) + 
((Z(ln[Lactate]))*ln(1.06)) + ((Z(ln[Histidine]))*ln(0.93)) + ((Z(ln[Isoleucine]))*ln(1.23)) + ((Z(ln[Leucine]))*ln(0.82)) + ((Z(ln[Valine]))*ln(0.87)) + ((Z(ln[Phenylalanine]))*ln(1.13)) + ((Z(ln[Acetoacetate]))*ln(1.08)) + ((Z(ln[Albumin]))*ln(0.89)) + ((Z(ln[Glycoprotein_acetyls]))*ln(1.32))).
Z states for z-scaling and ln states for natural logarithm.
```
**Run the R script saved here as `metabohealth.R`** 
- This does each step by step so I can see the resulting columns and then multiply them by each other
- this creates a file phenotype dataframe at the end whereby we have all the metabolites, characteristics and then a score - this is not in the phenotype format i need for GWAS yet as certain values need to be imputed in original data to avoid scores of 0
- Towards the end it might look like `Alb` is missing at an alarming rate but it is oly 12%
---
### Step 4 — Imputation using Age, and Sex as predictors
- This was a bit more sophisticated so I wanted to get scores before diving into it but **run the script saved as `metabohealth_imputed.R` to get imputed results**
- followed this lovely tutorial: https://libguides.princeton.edu/R-Missingdata

**Explanation of the Parameters:**

- `method`: The imputation method used (in this case, pmm stands for predictive mean matching and can alternatively use rf which is random forest).

- `m`: The number of imputed datasets (typically 5 or more).

---
### Step 5 — GWAS format 
- The `.fam` file in UKB has the following columns which correspond to Family ID (FIID), Individual ID (IID),  paternal ID, maternal ID and phenotype. 

```
==> ukb22418_c9_b0_v2.fam <==
4459327 4459327 0 0 1 Batch_b001
```

- `.pheno` file needs to be made in the same format (code for this in lines 273-286 in `metabohealth_imputed.R`), whereby phenotype is not just a batch number, but our phenotype score. It looks like this: 

```
FID	IID	Age	Sex	MetaboHealth_Score
0	1000073	55	1	0.120481033803025
0	1000312	62	0	-0.163868412744039
```


- `.cov` file also needs to be made in the same format (code for this in lines 289-297 in `metabohealth_imputed.R`). It looks like this: 

```
FID	IID	Age	Sex	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10
0	1000073	55	1	-12.5122	2.52381	-1.68065	1.15155-3.75587	-1.52478	-0.271479	-0.655341	-0.848979	0.933165
```

---

## EXTRA - IGNORE 
Just saving these here as it is the UKB IDs mapped to the names of variables in the Deelen paper: 
```
library(tibble)

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
```
I got the mappings of what they come up as in UKB from the table in `ukb_names` - which I got from a paper i cant remember the exact name of
