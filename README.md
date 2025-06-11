# GWAS Phenotype 
A simple pipeline for extracting my phenotype — all the way through to running my GWAS. This is MetaboHealth as defined in Deelen et al., 2019

## Step 1 — Extracting Metabolites for Phenotype
To get the correct data from **auth.dnanexus.com** (UK Biobank Research Analysis Platform), you need to use the `dx extract_dataset` command. 

### What’s Extracted:
- **14 metabolites** listed in the paper by *Deelen et al., 2019* (see [Supplement, page 30](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11311-9/MediaObjects/41467_2019_11311_MOESM1_ESM.pdf)):
- Covariates:
  - Age (`21003`)
  - Sex (`31`)
  - 10 Genetic PCs (`22009_a1` to `22009_a10`)
  - Covariates placed in `.fam` files are Family ID (FIID), Individual ID (IID),  paternal ID, maternal ID and phenotype.

 - Others
    - Algorithmically-defined dementia codes (source of dementia:`41202`, date of dementia diagnosis: `41204`) - this is just to check stats and unnecessary step - this is entirely omited from final files anyway. 

```bash
dx extract_dataset project-Gzyb0j8JQYbBjQxYqfb4xJYX:record-GzyfX70Jfj0bvy8YfvYQ302v --fields participant.eid,participant.p42019,participant.p42018,participant.p21003_i0,participant.p31,participant.p22009_a1,participant.p22009_a2,participant.p22009_a3,participant.p22009_a4,participant.p22009_a5,participant.p22009_a6,participant.p22009_a7,participant.p22009_a8,participant.p22009_a9,participant.p22009_a10,participant.p23470_i0,participant.p23471_i0,participant.p23463_i0,participant.p23465_i0,participant.p23466_i0,participant.p23467_i0,participant.p23468_i0,participant.p23476_i0,participant.p30600_i0,participant.p23480_i0,participant.p23453_i0,participant.p23482_i0,participant.p23573_i0,participant.p23431_i0 -o cohort_data2.csv

```


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
### Step 2 (on local R script) — Add weights and make score 
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

- `method`: The imputation method used (in this case, rf stands for random forest).

- `m`: The number of imputed datasets (typically 5 or more).

---
### Step 5 (on local R script)— GWAS format 
- The `.fam` file in UKB has the following columns which correspond to Family ID (FIID), Individual ID (IID),  paternal ID, maternal ID and phenotype. 

```
==> ukb22418_c9_b0_v2.fam <==
4459327 4459327 0 0 1 Batch_b001
```

- `.pheno` file needs to be made in the same format (code for this in lines 254-269 in `metabohealth_imputed.R`), whereby phenotype is not just a batch number, but our phenotype score. It looks like this: 

```
FID	IID	Age	Sex	MetaboHealth_Score
1000073	1000073	55	1	0.120481033803025
```


- `.cov` file also needs to be made in the same format (code for this in lines 260 onwards in `metabohealth_imputed.R`). It looks like this: 

```
FID	IID	Age	Sex	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10
1000073	1000073	55	1	-12.5122	2.52381	-1.68065	1.15155-3.75587	-1.52478	-0.271479	-0.655341	-0.848979	0.933165
```

### Step 6 (on UKB RAP) — Running GWAS on subset of data - this was done using just array data as a test
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
