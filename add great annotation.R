### add GREAT annotation
job_regions <- submitGreatJob(dat[, c("chr", "start", "end")], version = "3.0")
regionsToGenes <- data.frame(plotRegionGeneAssociationGraphs(job_regions))
regionsToGenes$GREAT_annotation <- ifelse(
  regionsToGenes$distTSS > 0,
  paste0(regionsToGenes$gene, " (+", regionsToGenes$distTSS, ")"),
  paste0(regionsToGenes$gene, " (", regionsToGenes$distTSS, ")"))
regionsToGenes <- regionsToGenes[
  ,c("seqnames", "start", "end", "GREAT_annotation")]

great <- regionsToGenes %>%
  group_by(seqnames, start, end) %>%
  mutate(time = paste0("time", row_number())) %>%
  tidyr::spread(time, GREAT_annotation) %>%
  ungroup()

great_final <- plyr::adply(great, 1, function(row){
  great_chr <- as.character(row[, paste0("time", 1:(ncol(row) - 3))])
  row$GREAT_annotation <- paste(great_chr, collapse = "; ")
  row})

great_final <- great_final[, c("seqnames", "start", "end", "GREAT_annotation")