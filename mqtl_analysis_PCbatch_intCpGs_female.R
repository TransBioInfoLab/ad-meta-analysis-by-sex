setwd("M:/AD/analysis/projects/Sex_strat/Aim3/Meta-analysis/Lissette/RosmapImputedData/mQTLanalysis_bySex/intCpGs/Female")

## Read pheno file
RosmapCov <- read.table("../covAgeBatchPCs.txt", header = TRUE)


## Read nominal significant SNPs dosage file
#Read sample file
SNPfileSamples <- read.table("rosmap12.fam")

#AD SNPs that overlap with AD CpGs in female
SNPfileF <- read.table("rosmap12Merged_CpGoverlap_MAF1_dosage_female_allChr_signSNPs.gen")
colnames(SNPfileF) <- c("SNP", SNPfileSamples$V1)

## Read Rosmap residuals file
resRosmap <- load("C:/Users/lxg255/Dropbox (BBSR)/AD_metaAnalysis_bySex/Lissette/ADmethyBySex/Rosmap/ROSMAP_mval_noCovs_and_pheno.rda")

## Read CpG/SNP overlap file

#female
overlapF <- read.csv("../CpGintSNPsOverlap.csv")
overlapF$SNP <- paste0(overlapF$chr_SNP, ":", overlapF$pos_SNP)
overlapF_ADsignSNPs <- overlapF[which(overlapF$SNP %in% SNPfileF$SNP), ]
saveRDS(overlapF_ADsignSNPs, "overlap_ADsignSNPs_female.rds")
overlapF_ADsignSNPs <- readRDS("overlap_ADsignSNPs_female.rds")
CpGsF_ADsign <- unique(overlapF_ADsignSNPs$cpg)

## Read sample key file
sampleKey <- read.table("../../sampleKey.txt", header = TRUE)
sampleKeyF <- sampleKey[which(sampleKey$SNPsampleID %in% SNPfileSamples$V1),]

## Filter residual matrix for the 434 females with genotype data and 
#the sign CpGs

#female
# resCpGsF_ADsignSNPs <- mval_noCovs_mat[which(rownames(mval_noCovs_mat) %in% CpGsF_ADsign),
#                                  which(colnames(mval_noCovs_mat) %in% sampleKeyF$methySampleID)]
# saveRDS(resCpGsF_ADsignSNPs, "residuals_CpGsF_434females_ADsignSNPs.rds")
resCpGsF_ADsign <- readRDS("residuals_CpGsF_434females_ADsignSNPs.rds")


## Match Ids in methy and SNP datasets

methyIDs <- data.frame(IDs = colnames(resCpGsF_ADsign))
methyIDsKey <- merge(methyIDs, 
                     sampleKeyF,
                     by.x = "IDs",
                     by.y = "methySampleID",
                     sort = FALSE)
identical(methyIDs$IDs, methyIDsKey$IDs)


#female
colnames(resCpGsF_ADsign) <- methyIDsKey$SNPsampleID
RosmapCovF <- RosmapCov[which(RosmapCov$FID %in% SNPfileSamples$V1),]
resMethyCpGsF_ADsign <- resCpGsF_ADsign[, as.character(RosmapCovF$FID)]

## Final list of CpGs
#female
CpGsF <- rownames(resMethyCpGsF_ADsign)

identical(as.vector(RosmapCovF$FID), colnames(SNPfileF)[-1])
identical(as.vector(RosmapCovF$FID), colnames(resMethyCpGsF_ADsign))

## mQTL analysis
## Methylation_resid ~ SNP dosage + batch + PC1 + PC2 + PC3
lmodel1 <- function(SNP) {
  
  f <- lm(CpGmethy ~ SNP + RosmapCovF$batch + RosmapCovF$PC1 + RosmapCovF$PC2 + RosmapCovF$PC3)
  
  res <- t(coef(summary(f))[2, c(1, 2, 4)])
  
  res
  
}


#female
allF_m1 <- data.frame(matrix(ncol=6,nrow=0))
for (i in CpGsF[1:length(CpGsF)]){
  
  CpGmethy <- as.numeric(as.vector(resMethyCpGsF_ADsign[i,]))
  CpG_SNPs <- subset(overlapF_ADsignSNPs, cpg == i, select = SNP)
  CpG_SNPs <- CpG_SNPs$SNP
  SNPdosage <- SNPfileF[which(SNPfileF$SNP %in% CpG_SNPs), ]
  CpG_SNP <- cbind(paste0(i, "-", SNPdosage$SNP), i, SNPdosage$SNP)
  
  m1 <- cbind(CpG_SNP, as.data.frame(t(apply(SNPdosage[,-1], 1, lmodel1))))
  allF_m1 <- rbind(allF_m1, m1)
  
}


colnames(allF_m1) <- c("CpG-SNP", "CpG", "SNP", "Estimate_SNP", "StdError_SNP", "P_SNP")
allF_m1$FDR <- p.adjust(allF_m1$P_SNP, "fdr")
write.csv(allF_m1, "ADmqtlAnalysisResults_intCpGs_female.csv")
saveRDS(allF_m1, "ADmqtlAnalysisResults_intCpGs_female.rds")

SNPglmResults <- read.table(
  "glmRosmap12_CpGoverlap_MAF1_female/glmCogdx_ResultsFreqAll_CpGs_female.txt",
  header = TRUE
  )

results <- merge(allF_m1, SNPglmResults, by.x = "SNP", by.y = "ID")

results$A1_FREQS_cases <- results$ALT_FREQS_cases
results$A1_FREQS_cases[results$ALT != results$A1] <- (1-results$ALT_FREQS_cases[results$ALT != results$A1])

results$A1_FREQS_controls <- results$ALT_FREQS_controls
results$A1_FREQS_controls[results$ALT != results$A1] <- (1-results$ALT_FREQS_controls[results$ALT != results$A1])


write.csv(results, "ADmqtlAnalysisResults_intCpGs_female_SNPglm.csv")
