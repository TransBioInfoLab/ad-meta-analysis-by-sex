# Sex-specific analysis of DNA methylation changes implicates new genes in Alzheimer’s Disease pathology
Lanyu Zhang, Juan I. Young, Lissette Gomez, Tiago C. Silva, Michael A. Schmidt, Jesse Cai, Xi Chen, Eden R. Martin, Lily Wang

### Description

Sex is an important factor that contributes to the clinical and biological heterogeneities in Alzheimer’s disease. DNA methylation is a major epigenetic modification that is known to be involved in AD. 

We studied DNA methylation changes across different AD Braak stages in over 1000 postmortem prefrontal cortex brain samples by performing meta-analyses on four brain studies. To identify sex-specific methylation changes, we employed two complementary approaches, a sex-stratified analysis that examined methylation-Braak stage associations in male and female samples separately, and a sex-by-Braak stage interaction analysis that compared the magnitude of these associations between different sexes.  

### Single cohort analysis

| File                 | Dataset | Link |
|----------------------|-------------|-------------|
| single_cohort_analysis_gasparoni.Rmd        |   Gasparoni (Gasparoni, 2018) | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/single_cohort_analysis_gasparoni.Rmd) |
| single_cohort_analysis_london.Rmd           |   London (Lunnon, 2014)    | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/single_cohort_analysis_london.Rmd) |
| single_cohort_analysis_mtsinai.Rmd          |   Mt. Sinai (Smith, 2018)  | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/single_cohort_analysis_mtsinai.Rmd) |
| single_cohort_analysis_rosmap.Rmd           |   ROSMAP (PMID: 29865057)    | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/single_cohort_analysis_rosmap.Rmd) |


### Meta-analysis 

| File                 | Link |
|----------------------|-------------|
| meta_analysis_single_cpg.Rmd | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/meta_analysis_single_cpg.rmd) |
| meta_analysis_region_based.Rmd | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/meta_analysis_region_based.rmd) |


### Meta-analysis results

| File                 | Link |
|----------------------|-------------|
| single_cpg_female_meta_df.csv         | [Link to the result](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/single_cpg_female_meta_df.csv) |
| single_cpg_male_meta_df.csv           | [Link to the result](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/single_cpg_male_meta_df.csv) |
| region_based_female_meta_df.csv       | [Link to the result](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/region_based_female_meta_df.csv) |
| region_based_male_meta_df.csv         | [Link to the result](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/region_based_male_meta_df.csv) |


### Enrichment analysis of significant DNA methylation changes 

| File                 | Link |
|----------------------|-------------|
| enrichment_analysis.rmd | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/enrichment_analysis.rmd) |

### Correlation of significant DNA methylation changes with expressions of nearby genes

| File                 | Link |
|----------------------|-------------|
| eqtm_analysis_single_cpg.R | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/eqtm_analysis_single_cpg.R) |
| eqtm_analysis_region_based.R | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/eqtm_analysis_region_based.R) |

### Correlation of significant DNA methylation changes with genetic variants

| File                 | Link |
|----------------------|-------------|
| mqtl_analysis_PCbatch_CpGs_female.R | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/mqtl_analysis_PCbatch_CpGs_female.R) |
| mqtl_analysis_PCbatch_CpGs_male.R | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/mqtl_analysis_PCbatch_CpGs_male.R) |
| mqtl_analysis_PCbatch_intCpGs_female.R | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/mqtl_analysis_PCbatch_intCpGs_female.R) |
| mqtl_analysis_PCbatch_intCpGs_male.R | [Link to the script](https://github.com/TransBioInfoLab/ADMetaBySex/blob/master/mqtl_analysis_PCbatch_intCpGs_male.R) |

### Acknowledgement
All datasets used in this study are publicly available. The Mt. Sinai, London, Gasparoni and ROSMAP datasets were obtained from GEO (accessions GSE80970, GSE59685, GSE66351) and Synapse (accession syn3157275). The ROSMAP study data were provided by the Rush Alzheimer’s Disease Center, Rush University Medical Center, Chicago. Data collection was supported through funding by NIA grants P30AG10161, R01AG15819, R01AG17917, R01AG30146, R01AG36836, U01AG32984, U01AG46152, the Illinois Department of Public Health, and the Translational Genomics Research Institute. 


