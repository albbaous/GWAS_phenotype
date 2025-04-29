# GWAS_prep
A simple pipeline for extracting my phenotype — all the way through to running my GWAS.

---

## Step 1 — Extracting Metabolites for Phenotype

To get the correct data from **auth.dnanexus.com** (UK Biobank Research Analysis Platform), you need to use the following command. This extracts all **14 metabolites** listed in the paper by *Deelen et al., 2019* (see [Supplement, page 30](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11311-9/MediaObjects/41467_2019_11311_MOESM1_ESM.pdf)):

```bash
dx extract_dataset \
  project-Gzyb0j8JQYbBjQxYqfb4xJYX:record-GzyfX70Jfj0bvy8YfvYQ302v \
  --fields participant.eid,participant.p23470_i0,participant.p23471_i0,participant.p23463_i0,participant.p23465_i0,participant.p23466_i0,participant.p23467_i0,participant.p23468_i0,participant.p23476_i0,participant.p30600_i0,participant.p23480_i0,participant.p23473_i0,participant.p23453_i0,participant.p23482_i0,participant.p23573_i0,participant.p23431_i0 \
  -o cohort_data2.csv
```

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

### Step 2 — Mapping for CSF db

Once you've extracted the metabolite data, the next step is to map participants to those on the UoM db
#### Actions

1. **Extract additional baseline characteristics** from UK Biobank and add those to the cohort file and command from Step 1:
   - Age (`21003`)
   - Sex (`31`)
   - BMI (`21001`)
   - Smoking (`20116`)
   - Other health status variables as needed

2. **Integrate phenotype definitions**:
   - Map participants to CSF data
   - Use identifying characteristics for merging datasets

---

### Step 3 — Add weights
Add weights from Deelen et al., paper to get score 

Just saving these here as it is the UKB IDs mapped to the names of variables in UKB: 
# Define biomarkers and weights from Deelen et al., 2019
# The "column" field matches exactly with the headers in your CSV
```
biomarkers <- tribble(
  ~column,                   ~label,       ~lnHR,
  "participant.p23470_i0",   "Glc",         0.246,
  "participant.p23471_i0",   "Lac",         0.130,
  "participant.p23463_i0",   "His",        -0.067,
  "participant.p23465_i0",   "Ile",         0.056,
  "participant.p23466_i0",   "Leu",        -0.318,
  "participant.p23467_i0",   "Val",        -0.035,
  "participant.p23468_i0",   "Phe",         0.115,
  "participant.p23476_i0",   "AcAce",       0.100,
  "participant.p30600_i0",   "Alb",        -0.149,
  "participant.p23480_i0",   "GlycA",       0.272,
  "participant.p23473_i0",   "Cit",         0.289,
  "participant.p23453_i0",   "PUFA_FA",    -0.253,
  "participant.p23482_i0",   "XXL_VLDL_L",  0.022,
  "participant.p23573_i0",   "S_HDL_L",    -0.123,
  "participant.p23431_i0",   "VLDL_D",     -0.245
)
```
