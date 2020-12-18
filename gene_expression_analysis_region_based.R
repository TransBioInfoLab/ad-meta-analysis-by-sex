library(rGREAT)
library(dplyr)

# Regions
readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_region_based_results/summary_table_for_region_based_analysis.csv")  %>% 
  dplyr::filter(female_fdr_sig == 1) %>% dim()

# checking numbers
readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_region_based_results/summary_table_for_region_based_analysis.csv")  %>% 
  dplyr::filter(male_fdr_sig == 1) %>% dim()

meta.by.sex.region <- readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_region_based_results/summary_table_for_region_based_analysis.csv")  %>% 
  dplyr::filter((male_fdr_sig | female_fdr_sig) == 1)

meta.by.sex.region.female <- readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_region_based_results/summary_table_for_region_based_analysis.csv")  %>% 
  dplyr::filter((female_fdr_sig) == 1)

meta.by.sex.region.male <- readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_region_based_results/summary_table_for_region_based_analysis.csv")  %>% 
  dplyr::filter((male_fdr_sig) == 1)

job_regions <- submitGreatJob(
  gr = meta.by.sex.region %>% makeGRangesFromDataFrame(),
  species = "hg19" # or "hg38".
)
regionsToGenes <- data.frame(plotRegionGeneAssociationGraphs(job_regions))
regionsToGenes$type <- "Promoter"
regionsToGenes$type[abs(regionsToGenes$distTSS) > 2000] <- "Distal"
regionsToGenes$GREAT_annotation <- ifelse(
  regionsToGenes$distTSS > 0,
  paste0(regionsToGenes$gene, "(+", regionsToGenes$distTSS, ")"),
  paste0(regionsToGenes$gene, "(", regionsToGenes$distTSS, ")"))

#-------------------------------------------------------------------------------
# (3) the rest is the same as we did before, 
# as in this file https://www.dropbox.com/s/cgp2qoi650fhaht/5_gene-expression-analysis-for-all-not-by-case-control.Rmd?dl=0
# only difference is that here we need one analysis that fit model to both cases and controls in female samples (adjust for Braak stage)
# another analysis that fit model to both cases and controls in male samples (adjust for braak stage)
#-------------------------------------------------------------------------------
load("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/Tiago/data_male_female_resid_dmrs.rda")

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# Analysis                                                                     |
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# sex specific                                                                     |
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+
# Male                  |
#=+=+=+=+=+=+=+=+=+=+=+=+

gene.info <- TCGAbiolinks::get.GRCh.bioMart(genome = "hg19")
regionsToGenes$target <- gene.info$ensembl_gene_id[match(regionsToGenes$gene,gene.info$external_gene_name)]

regionsToGenes$regionID <-  paste0(regionsToGenes$seqnames,":",regionsToGenes$start,"-",regionsToGenes$end)
regionsToGenes.male <- regionsToGenes[regionsToGenes$regionID %in%  meta.by.sex.region.male$inputRegion,]
regionsToGenes.male$regionID %>% unique %>% length # 76
regionsToGenes.female <- regionsToGenes[regionsToGenes$regionID %in% meta.by.sex.region.female$inputRegion,]
regionsToGenes.female$regionID %>% unique %>% length # 381

regionsToGenes.female <- regionsToGenes.female[regionsToGenes.female$target %in% rownames(exp.female.resid),]
regionsToGenes.male <- regionsToGenes.male[regionsToGenes.male$target %in% rownames(exp.male.resid),]

resid_exp <- exp.male.resid
resid_met <- dnam.male.resid
metada.samples <- metada.samples.male

doParallel::registerDoParallel(cores = 2)
results.male <- plyr::adply(
  regionsToGenes.male,
  .margins = 1,
  .fun = function(row) {
    tryCatch({
      rna.target <- resid_exp[rownames(resid_exp) == row$target, , drop = FALSE]
      met.residual <- resid_met[rownames(resid_met) == as.character(row$regionID), ]
      
      df <- data.frame(
        rna.residual = rna.target %>% as.numeric,
        met.residual = met.residual %>% as.numeric,
        Braak_stage = metada.samples$braaksc %>% as.numeric
      )
      
      # fit linear model:
      results.all <-  lm(
        rna.residual ~ met.residual + Braak_stage, data = df
      )
      results.all.pval <- summary(results.all)$coefficients[
        2, "Pr(>|t|)", drop = F] %>% 
        t %>% as.data.frame()
      results.all.estimate <- summary(results.all)$coefficients[
        2, "Estimate", drop = F] %>% 
        t %>% as.data.frame()
      colnames(results.all.pval) <- paste0(
        "all_pval_", colnames(results.all.pval))
      colnames(results.all.estimate) <- paste0(
        "all_estimate_", colnames(results.all.estimate))
      
      return(
        data.frame(
          cbind(results.all.pval, results.all.estimate),
          row.names = NULL,
          stringsAsFactors = FALSE
        )
      )
    }, error = function(e) {
      print(row)
      return()
    })
  },
  .id = NULL,
  .progress = "time",
  .parallel = TRUE,
  .inform = FALSE
)
# Apply BH to distal, then to promoter
results.male <- results.male[!is.na(results.male$all_pval_met.residual),]
results.male <- results.male %>%
  group_by(type) %>% 
  mutate(pval.adj = p.adjust (all_pval_met.residual, method = 'BH'))
meta.by.sex.region$regionID <- meta.by.sex.region$inputRegion
results.male <- merge(results.male,meta.by.sex.region[,c( "regionID" , "UCSC_RefGene_Group",  "Relation_to_Island", "state")])

writexl::write_xlsx(results.male, path = "results/results.dmrs.male_GREAT.xlsx")


#=+=+=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Female
#=+=+=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
resid_exp <- exp.female.resid
resid_met <- dnam.female.resid
metada.samples <- metada.samples.female
doParallel::registerDoParallel(cores = 2)
results.female <- plyr::adply(
  regionsToGenes.female,
  .margins = 1,
  .fun = function(row) {
    tryCatch({
      rna.target <- resid_exp[rownames(resid_exp) == row$target, , drop = FALSE]
      met.residual <- resid_met[rownames(resid_met) == as.character(row$regionID), ]
      
      df <- data.frame(
        rna.residual = rna.target %>% as.numeric,
        met.residual = met.residual %>% as.numeric,
        Braak_stage = metada.samples$braaksc %>% as.numeric
      )
      
      # fit linear model:
      results.all <-  lm(
        rna.residual ~ met.residual + Braak_stage, data = df
      )
      results.all.pval <- summary(results.all)$coefficients[
        2, "Pr(>|t|)", drop = F] %>% 
        t %>% as.data.frame()
      results.all.estimate <- summary(results.all)$coefficients[
        2, "Estimate", drop = F] %>% 
        t %>% as.data.frame()
      colnames(results.all.pval) <- paste0(
        "all_pval_", colnames(results.all.pval))
      colnames(results.all.estimate) <- paste0(
        "all_estimate_", colnames(results.all.estimate))
      
      return(
        data.frame(
          cbind(results.all.pval, results.all.estimate),
          row.names = NULL,
          stringsAsFactors = FALSE
        )
      )
    }, error = function(e) {
      print(row)
      return()
    })
  },
  .id = NULL,
  .progress = "time",
  .parallel = TRUE,
  .inform = FALSE
)
results.female <- results.female[!is.na(results.female$all_pval_met.residual),]
results.female <- results.female %>%
  group_by(type) %>% 
  mutate(pval.adj = p.adjust (all_pval_met.residual, method = 'BH'))
results.female <- merge(results.female,meta.by.sex.region[,c( "regionID" , "UCSC_RefGene_Group",  "Relation_to_Island", "state")])

writexl::write_xlsx(results.female, path = "results/results.dmrs.female_GREAT.xlsx")