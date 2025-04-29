# GWAS_prep
A simple pipeline for extracting my phenotype — all the way through to running my GWAS.

---

## 📌 Step 1 — Extracting Metabolites for Phenotype

To get the correct data from **auth.dnanexus.com** (UK Biobank Research Analysis Platform), you need to use the following command. This extracts all **14 metabolites** listed in the paper by *Deelen et al., 2019* (see [Supplement, page 30](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11311-9/MediaObjects/41467_2019_11311_MOESM1_ESM.pdf)):

```bash
dx extract_dataset \
  project-Gzyb0j8JQYbBjQxYqfb4xJYX:record-GzyfX70Jfj0bvy8YfvYQ302v \
  --fields participant.eid,participant.p23470_i0,participant.p23471_i0,participant.p23463_i0,participant.p23465_i0,participant.p23466_i0,participant.p23467_i0,participant.p23468_i0,participant.p23476_i0,participant.p30600_i0,participant.p23480_i0,participant.p23473_i0,participant.p23453_i0,participant.p23482_i0,participant.p23573_i0,participant.p23431_i0 \
  -o cohort_data2.csv

### Explanation

- `dx extract_dataset`: DNAnexus CLI command for extracting a dataset from RAP  
- `project-...`: The ID of your RAP project, extracted using:
project = os.getenv('DX_PROJECT_CONTEXT_ID')
project

project-Gzyb0j8JQYbBjQxYqfb4xJYX

- `record-...`: The file record (metabolite dataset), extracted using:

record = os.popen("dx find data --type Dataset --delimiter ',' | awk -F ',' '{print $5}'").read().rstrip()
record

'record-GzyfX70Jfj0bvy8YfvYQ302v'
- `--fields`: List of fields to extract (in `participant.p[FIELD_ID]_i0` format)  
- `-o`: Specifies the output CSV file  

> ⚠️ **Note**: This command only retrieves metabolite values.  
> You will also need to extract baseline characteristics like age, sex, and BMI to complete phenotype mapping.

---

### Step 2 — Mapping for CSF db

Once you've extracted the metabolite data, the next step is to map participants to those on the UoM db
#### ✅ Actions

1. **Extract baseline characteristics** from UK Biobank:
   - Age (`21003`)
   - Sex (`31`)
   - BMI (`21001`)
   - Smoking (`20116`)
   - Other health status variables as needed

2. **Integrate phenotype definitions**:
   - Map participants to CSF data, case/control status, or other outcomes
   - Use `eid` (participant ID) for merging datasets

3. **Perform quality control**:
   - Remove individuals with missing values in critical variables
   - Filter by ancestry if necessary (e.g. using `21000` for ethnicity)

4. **Generate final phenotype file**:
   - Merge all features into a single file:
     - `phenotype_data.csv`
   - Include: `eid`, phenotype (case/control), covariates, and selected biomarkers
