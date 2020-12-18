library(rGREAT)

readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/summary_table_for_single_cpg_analysis.csv") %>% 
  dplyr::filter(male_fdr_sig == 1) %>% dim()

# checking numbers
readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/summary_table_for_single_cpg_analysis.csv") %>% 
  dplyr::filter(female_fdr_sig == 1) %>% dim()

meta.by.sex.cpg.female <- readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/summary_table_for_single_cpg_analysis.csv") %>% 
  dplyr::filter((female_fdr_sig) == 1)

meta.by.sex.cpg.male <- readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/summary_table_for_single_cpg_analysis.csv") %>% 
  dplyr::filter((male_fdr_sig) == 1)

meta.by.sex.cpg <- readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/summary_table_for_single_cpg_analysis.csv") %>% 
  dplyr::filter((male_fdr_sig | female_fdr_sig) == 1)

meta.by.sex.cpg$regionID <- paste0(meta.by.sex.cpg$chr,":", meta.by.sex.cpg$start,"-", meta.by.sex.cpg$end)

# checking numbers
readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/summary_table_for_single_cpg_analysis.csv") %>% 
  dplyr::filter(intxn_stageR_sig == 1) %>% dim()

meta.by.sex.cpg.intxn_stageR_sig <- readr::read_csv("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/summary_table_for_single_cpg_analysis.csv") %>% 
  dplyr::filter((intxn_stageR_sig) == 1)

hm450.hg19 <- sesameData::sesameDataGet("HM450.hg19.manifest")
hm450.hg19.sex <- hm450.hg19[c(meta.by.sex.cpg.intxn_stageR_sig$cpg,meta.by.sex.cpg$cpg) %>% unique,]

job_regions <- submitGreatJob(
  gr = hm450.hg19.sex,
  species = "hg19" # or "hg38".
)
prebes.info <- as.data.frame(hm450.hg19.sex)[,c("seqnames","start","end")]
prebes.info$cpg <- rownames(prebes.info)
regionsToGenes <- data.frame(plotRegionGeneAssociationGraphs(job_regions))
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
load("~/Dropbox (BBSR)/AD_metaAnalysis_bySex/Tiago/data_male_female_resid.rda")

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# Analysis                                                                     |
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# sex specific                                                                     |
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+
# Male                  |
#=+=+=+=+=+=+=+=+=+=+=+=+
regionsToGenes$type <- "Promoter"
regionsToGenes$type[abs(regionsToGenes$distTSS) > 2000] <- "Distal"
gene.info <- TCGAbiolinks::get.GRCh.bioMart(genome = "hg19")
regionsToGenes$target <- gene.info$ensembl_gene_id[match(regionsToGenes$gene,gene.info$external_gene_name)]
cpg.target <- merge(prebes.info,regionsToGenes)

cpg.target$regionID <-  paste0(cpg.target$seqnames,":",cpg.target$start,"-",cpg.target$end)
cpg.target.male <- cpg.target[cpg.target$regionID %in% paste0(meta.by.sex.cpg.male$chr,":",cpg.target$start,"-",cpg.target$end),]
cpg.target.male$regionID %>% unique %>% length # 76
cpg.target.female <- cpg.target[cpg.target$regionID %in% paste0(meta.by.sex.cpg.female$chr,":",meta.by.sex.cpg.female$start,"-",meta.by.sex.cpg.female$end),]
cpg.target.female$regionID %>% unique %>% length # 381

cpg.target.female <- cpg.target.female[cpg.target.female$target %in% rownames(resid_exp),]
cpg.target.male <- cpg.target.male[cpg.target.male$target %in% rownames(resid_exp),]

resid_exp <- exp.male.resid
resid_met <- dnam.male.resid
metada.samples <- metada.samples.male

doParallel::registerDoParallel(cores = 2)
results.male <- plyr::adply(
  cpg.target.male,
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
results.male <- merge(results.male,meta.by.sex.cpg[,c( "regionID" , "UCSC_RefGene_Group",  "Relation_to_Island", "state","cpg")])
results.male <- merge(results.male,meta.by.sex.cpg[,c( "regionID" ,"cpg")])

writexl::write_xlsx(results.male, path = "results/results.cpgs.male_GREAT.xlsx")


#=+=+=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Female
#=+=+=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
resid_exp <- exp.female.resid
resid_met <- dnam.female.resid
metada.samples <- metada.samples.female
doParallel::registerDoParallel(cores = 2)
results.female <- plyr::adply(
  cpg.target.female,
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
results.female <- merge(results.female,meta.by.sex.cpg[,c( "regionID" , "UCSC_RefGene_Group",  "Relation_to_Island", "state","cpg")])

writexl::write_xlsx(results.female, path = "results/results.cpgs.female_GREAT.xlsx")


#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# interaction specific                                                                     |
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# Male                  
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
cpg.intxn.target <- cpg.target[cpg.target$regionID %in% paste0(meta.by.sex.cpg.intxn_stageR_sig$chr,":",meta.by.sex.cpg.intxn_stageR_sig$start,"-",meta.by.sex.cpg.intxn_stageR_sig$end),]
cpg.intxn.target <- cpg.intxn.target[cpg.intxn.target$target %in% rownames(resid_exp),]

resid_exp <- exp.male.resid
resid_met <- dnam.male.resid
metada.samples <- metada.samples.male

doParallel::registerDoParallel(cores = 2)
results.intxn.male <- plyr::adply(
  cpg.intxn.target,
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
results.intxn.male <- results.intxn.male[!is.na(results.intxn.male$all_pval_met.residual),]
results.intxn.male <- results.intxn.male %>%
  group_by(type) %>% 
  mutate(pval.adj = p.adjust (all_pval_met.residual, method = 'BH'))
results.intxn.male <- merge(results.intxn.male,meta.by.sex.cpg[,c( "regionID" , "UCSC_RefGene_Group",  "Relation_to_Island", "state","cpg")])
#results.intxn.male <- merge(results.intxn.male,meta.by.sex.cpg[,c( "regionID" ,"cpg")])

writexl::write_xlsx(results.intxn.male, path = "results/results.intxn.cpgs.male_GREAT.xlsx")


#=+=+=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+
# Female
#=+=+=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+=+=+=+=++=+=+=+=+=+=+
resid_exp <- exp.female.resid
resid_met <- dnam.female.resid
metada.samples <- metada.samples.female
doParallel::registerDoParallel(cores = 2)
results.intxn.female <- plyr::adply(
  cpg.intxn.target,
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
results.intxn.female <- results.intxn.female[!is.na(results.intxn.female$all_pval_met.residual),]
results.intxn.female <- results.intxn.female %>%
  group_by(type) %>% 
  mutate(pval.adj = p.adjust (all_pval_met.residual, method = 'BH'))
results.intxn.female <- merge(results.intxn.female,meta.by.sex.cpg[,c( "regionID" , "UCSC_RefGene_Group",  "Relation_to_Island", "state","cpg")])

writexl::write_xlsx(results.intxn.female, path = "results/results.intxn.cpgs.female_GREAT.xlsx")

