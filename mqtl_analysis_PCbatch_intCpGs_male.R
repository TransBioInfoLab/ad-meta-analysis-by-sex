#author: Lissette Gomez, Lily Wang

## Read pheno file
RosmapCov <- read.table("../covAgeBatchPCs.txt", header = TRUE)


## Read nominal significant SNPs dosage file
#Read sample file
SNPfileSamples <- read.table("rosmap12.fam")

#AD SNPs that overlap with AD CpGs in male
SNPfileM <- read.table("rosmap12Merged_CpGoverlap_MAF1_dosage_male_allChr_signSNPs.gen")
colnames(SNPfileM) <- c("SNP", SNPfileSamples$V1)

## Read Rosmap residuals file
resRosmap <- load("C:/Users/lxg255/Dropbox (BBSR)/AD_metaAnalysis_bySex/Lissette/ADmethyBySex/Rosmap/ROSMAP_mval_noCovs_and_pheno.rda")

## Read CpG/SNP overlap file


#male
overlapM <- read.csv("../CpGintSNPsOverlap.csv")
overlapM$SNP <- paste0(overlapM$chr_SNP, ":", overlapM$pos_SNP)
overlapM_ADsignSNPs <- overlapM[which(overlapM$SNP %in% SNPfileM$SNP), ]
saveRDS(overlapM_ADsignSNPs, "overlap_ADsignSNPs_male.rds")
overlapM_ADsignSNPs <- readRDS("overlap_ADsignSNPs_male.rds")
CpGsM_ADsign <- unique(overlapM_ADsignSNPs$cpg)

## Read sample key file
sampleKey <- read.table("../../sampleKey.txt", header = TRUE)
sampleKeyM <- sampleKey[which(sampleKey$SNPsampleID %in% SNPfileSamples$V1),]

## Filter residual matrix for the 434 females with genotype data and 
#the sign CpGs

#male
resCpGsM_ADsignSNPs <- mval_noCovs_mat[which(rownames(mval_noCovs_mat) %in% CpGsM_ADsign),
                                 which(colnames(mval_noCovs_mat) %in% sampleKeyM$methySampleID)]
saveRDS(resCpGsM_ADsignSNPs, "residuals_CpGsF_254males_ADsignSNPs.rds")
resCpGsM_ADsign <- readRDS("residuals_CpGsF_254males_ADsignSNPs.rds")


## Match Ids in methy and SNP datasets

methyIDs <- data.frame(IDs = colnames(resCpGsM_ADsign))
methyIDsKey <- merge(methyIDs, 
                     sampleKeyM,
                     by.x = "IDs",
                     by.y = "methySampleID",
                     sort = FALSE)
identical(methyIDs$IDs, methyIDsKey$IDs)

#male
colnames(resCpGsM_ADsign) <- methyIDsKey$SNPsampleID
RosmapCovM <- RosmapCov[which(RosmapCov$FID %in% SNPfileSamples$V1),]
resMethyCpGsM_ADsign <- resCpGsM_ADsign[, as.character(RosmapCovM$FID)]

## Final list of CpGs
#male
CpGsM <- rownames(resMethyCpGsM_ADsign)

identical(as.vector(RosmapCovM$FID), colnames(SNPfileM)[-1])
identical(as.vector(RosmapCovM$FID), colnames(resMethyCpGsM_ADsign))

## mQTL analysis
## Methylation_resid ~ SNP dosage + batch + PC1 + PC2 + PC3
lmodel1 <- function(SNP) {
  
  f <- lm(CpGmethy ~ SNP + RosmapCovM$batch + RosmapCovM$PC1 + RosmapCovM$PC2 + RosmapCovM$PC3)
  
  res <- t(coef(summary(f))[2, c(1, 2, 4)])
  
  res
  
}



#male
allM_m1 <- data.frame(matrix(ncol=6,nrow=0))
for (i in CpGsM[1:length(CpGsM)]){
  
  CpGmethy <- as.numeric(as.vector(resMethyCpGsM_ADsign[i,]))
  CpG_SNPs <- subset(overlapM_ADsignSNPs, cpg == i, select = SNP)
  CpG_SNPs <- CpG_SNPs$SNP
  SNPdosage <- SNPfileM[which(SNPfileM$SNP %in% CpG_SNPs), ]
  CpG_SNP <- cbind(paste0(i, "-", SNPdosage$SNP), i, SNPdosage$SNP)
  
  m1 <- cbind(CpG_SNP, as.data.frame(t(apply(SNPdosage[,-1], 1, lmodel1))))
  allM_m1 <- rbind(allM_m1, m1)
  
}


colnames(allM_m1) <- c("CpG-SNP", "CpG", "SNP", "Estimate_SNP", "StdError_SNP", "P_SNP")
allM_m1$FDR <- p.adjust(allM_m1$P_SNP, "fdr")
write.csv(allM_m1, "ADmqtlAnalysisResults_CpGs_male.csv")
saveRDS(allM_m1, "ADmqtlAnalysisResults_CpGs_male.rds")

SNPglmResults <- read.table(
  "glmRosmap12_CpGoverlap_MAF1_male/glmCogdx_ResultsFreqAll_CpGs_male.txt",
  header = TRUE
  )

results <- merge(allM_m1, SNPglmResults, by.x = "SNP", by.y = "ID")

results$A1_FREQS_cases <- results$ALT_FREQS_cases
results$A1_FREQS_cases[results$ALT != results$A1] <- (1-results$ALT_FREQS_cases[results$ALT != results$A1])

results$A1_FREQS_controls <- results$ALT_FREQS_controls
results$A1_FREQS_controls[results$ALT != results$A1] <- (1-results$ALT_FREQS_controls[results$ALT != results$A1])


write.csv(results, "ADmqtlAnalysisResults_CpGs_male_SNPglm.csv")
