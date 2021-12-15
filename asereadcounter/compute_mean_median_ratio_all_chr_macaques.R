args = commandArgs(trailingOnly=TRUE)

id=args[1]
tissue=args[2]
chr_n = args[3]

setwd(file.path("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/", id, sep=""))

# chrX
chrX = read.csv(paste(id, "_", tissue, "_Xchr.tsv", sep=""), sep = '\t')
colnames(chrX) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')

chrX$refRatio = chrX$refCount/chrX$totalCount
chrX$altRatio = chrX$altCount/chrX$totalCount

chrX_ratio = c()
for (i in 1:nrow(chrX)) {
  row = chrX[i,]
  if (row$refRatio > row$altRatio) {
    chrX_ratio = c(chrX_ratio, row$refRatio)
  }
  else (chrX_ratio = c(chrX_ratio, row$altRatio))
}

chrX_data = data.frame(chr="chrX", mean=mean(chrX_ratio), median=median(chrX_ratio))
write.table(chrX_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chrX.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr1
chr1 = read.csv(paste(id, "_", tissue, "_chr1.tsv", sep=""), sep = '\t')
colnames(chr1) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr1$refRatio = chr1$refCount/chr1$totalCount
chr1$altRatio = chr1$altCount/chr1$totalCount

chr1_ratio = c()
for (i in 1:nrow(chr1)) {
  row = chr1[i,]
  if (row$refRatio > row$altRatio) {
    chr1_ratio = c(chr1_ratio, row$refRatio)
  }
  else (chr1_ratio = c(chr1_ratio, row$altRatio))
}
chr1_data = data.frame(chr="chr1", mean=mean(chr1_ratio), median=median(chr1_ratio))
write.table(chr1_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr1.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr2
chr2 = read.csv(paste(id, "_", tissue, "_chr2.tsv", sep=""), sep = '\t')
colnames(chr2) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr2$refRatio = chr2$refCount/chr2$totalCount
chr2$altRatio = chr2$altCount/chr2$totalCount

chr2_ratio = c()
for (i in 1:nrow(chr2)) {
  row = chr2[i,]
  if (row$refRatio > row$altRatio) {
    chr2_ratio = c(chr2_ratio, row$refRatio)
  }
  else (chr2_ratio = c(chr2_ratio, row$altRatio))
}
chr2_data = data.frame(chr="chr2", mean=mean(chr2_ratio), median=median(chr2_ratio))
write.table(chr2_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr2.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr3
chr3 = read.csv(paste(id, "_", tissue, "_chr3.tsv", sep=""), sep = '\t')
colnames(chr3) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr3$refRatio = chr3$refCount/chr3$totalCount
chr3$altRatio = chr3$altCount/chr3$totalCount

chr3_ratio = c()
for (i in 1:nrow(chr3)) {
  row = chr3[i,]
  if (row$refRatio > row$altRatio) {
    chr3_ratio = c(chr3_ratio, row$refRatio)
  }
  else (chr3_ratio = c(chr3_ratio, row$altRatio))
}
chr3_data = data.frame(chr="chr3", mean=mean(chr3_ratio), median=median(chr3_ratio))
write.table(chr3_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr3.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr4
chr4 = read.csv(paste(id, "_", tissue, "_chr4.tsv", sep=""), sep = '\t')
colnames(chr4) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr4$refRatio = chr4$refCount/chr4$totalCount
chr4$altRatio = chr4$altCount/chr4$totalCount

chr4_ratio = c()
for (i in 1:nrow(chr4)) {
  row = chr4[i,]
  if (row$refRatio > row$altRatio) {
    chr4_ratio = c(chr4_ratio, row$refRatio)
  }
  else (chr4_ratio = c(chr4_ratio, row$altRatio))
}
chr4_data = data.frame(chr="chr4", mean=mean(chr4_ratio), median=median(chr4_ratio))
write.table(chr4_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr4.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr5
chr5 = read.csv(paste(id, "_", tissue, "_chr5.tsv", sep=""), sep = '\t')
colnames(chr5) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr5$refRatio = chr5$refCount/chr5$totalCount
chr5$altRatio = chr5$altCount/chr5$totalCount

chr5_ratio = c()
for (i in 1:nrow(chr5)) {
  row = chr5[i,]
  if (row$refRatio > row$altRatio) {
    chr5_ratio = c(chr5_ratio, row$refRatio)
  }
  else (chr5_ratio = c(chr5_ratio, row$altRatio))
}
chr5_data = data.frame(chr="chr5", mean=mean(chr5_ratio), median=median(chr5_ratio))
write.table(chr5_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr5.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr6
chr6 = read.csv(paste(id, "_", tissue, "_chr6.tsv", sep=""), sep = '\t')
colnames(chr6) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr6$refRatio = chr6$refCount/chr6$totalCount
chr6$altRatio = chr6$altCount/chr6$totalCount

chr6_ratio = c()
for (i in 1:nrow(chr6)) {
  row = chr6[i,]
  if (row$refRatio > row$altRatio) {
    chr6_ratio = c(chr6_ratio, row$refRatio)
  }
  else (chr6_ratio = c(chr6_ratio, row$altRatio))
}
chr6_data = data.frame(chr="chr6", mean=mean(chr6_ratio), median=median(chr6_ratio))
write.table(chr6_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr6.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr7
chr7 = read.csv(paste(id, "_", tissue, "_chr7.tsv", sep=""), sep = '\t')
colnames(chr7) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr7$refRatio = chr7$refCount/chr7$totalCount
chr7$altRatio = chr7$altCount/chr7$totalCount

chr7_ratio = c()
for (i in 1:nrow(chr7)) {
  row = chr7[i,]
  if (row$refRatio > row$altRatio) {
    chr7_ratio = c(chr7_ratio, row$refRatio)
  }
  else (chr7_ratio = c(chr7_ratio, row$altRatio))
}
chr7_data = data.frame(chr="chr7", mean=mean(chr7_ratio), median=median(chr7_ratio))
write.table(chr7_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr7.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr8
chr8 = read.csv(paste(id, "_", tissue, "_chr8.tsv", sep=""), sep = '\t')
colnames(chr8) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr8$refRatio = chr8$refCount/chr8$totalCount
chr8$altRatio = chr8$altCount/chr8$totalCount

chr8_ratio = c()
for (i in 1:nrow(chr8)) {
  row = chr8[i,]
  if (row$refRatio > row$altRatio) {
    chr8_ratio = c(chr8_ratio, row$refRatio)
  }
  else (chr8_ratio = c(chr8_ratio, row$altRatio))
}
chr8_data = data.frame(chr="chr8", mean=mean(chr8_ratio), median=median(chr8_ratio))
write.table(chr8_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr8.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr9
chr9 = read.csv(paste(id, "_", tissue, "_chr9.tsv", sep=""), sep = '\t')
colnames(chr9) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr9$refRatio = chr9$refCount/chr9$totalCount
chr9$altRatio = chr9$altCount/chr9$totalCount

chr9_ratio = c()
for (i in 1:nrow(chr9)) {
  row = chr9[i,]
  if (row$refRatio > row$altRatio) {
    chr9_ratio = c(chr9_ratio, row$refRatio)
  }
  else (chr9_ratio = c(chr9_ratio, row$altRatio))
}
chr9_data = data.frame(chr="chr9", mean=mean(chr9_ratio), median=median(chr9_ratio))
write.table(chr9_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr9.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr10
chr10 = read.csv(paste(id, "_", tissue, "_chr10.tsv", sep=""), sep = '\t')
colnames(chr10) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr10$refRatio = chr10$refCount/chr10$totalCount
chr10$altRatio = chr10$altCount/chr10$totalCount

chr10_ratio = c()
for (i in 1:nrow(chr10)) {
  row = chr10[i,]
  if (row$refRatio > row$altRatio) {
    chr10_ratio = c(chr10_ratio, row$refRatio)
  }
  else (chr10_ratio = c(chr10_ratio, row$altRatio))
}
chr10_data = data.frame(chr="chr10", mean=mean(chr10_ratio), median=median(chr10_ratio))
write.table(chr10_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr10.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr11
chr11 = read.csv(paste(id, "_", tissue, "_chr11.tsv", sep=""), sep = '\t')
colnames(chr11) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr11$refRatio = chr11$refCount/chr11$totalCount
chr11$altRatio = chr11$altCount/chr11$totalCount

chr11_ratio = c()
for (i in 1:nrow(chr11)) {
  row = chr11[i,]
  if (row$refRatio > row$altRatio) {
    chr11_ratio = c(chr11_ratio, row$refRatio)
  }
  else (chr11_ratio = c(chr11_ratio, row$altRatio))
}
chr11_data = data.frame(chr="chr11", mean=mean(chr11_ratio), median=median(chr11_ratio))
write.table(chr11_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr11.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr12
chr12 = read.csv(paste(id, "_", tissue, "_chr12.tsv", sep=""), sep = '\t')
colnames(chr12) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr12$refRatio = chr12$refCount/chr12$totalCount
chr12$altRatio = chr12$altCount/chr12$totalCount

chr12_ratio = c()
for (i in 1:nrow(chr12)) {
  row = chr12[i,]
  if (row$refRatio > row$altRatio) {
    chr12_ratio = c(chr12_ratio, row$refRatio)
  }
  else (chr12_ratio = c(chr12_ratio, row$altRatio))
}
chr12_data = data.frame(chr="chr12", mean=mean(chr12_ratio), median=median(chr12_ratio))
write.table(chr12_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr12.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr13
chr13 = read.csv(paste(id, "_", tissue, "_chr13.tsv", sep=""), sep = '\t')
colnames(chr13) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr13$refRatio = chr13$refCount/chr13$totalCount
chr13$altRatio = chr13$altCount/chr13$totalCount

chr13_ratio = c()
for (i in 1:nrow(chr13)) {
  row = chr13[i,]
  if (row$refRatio > row$altRatio) {
    chr13_ratio = c(chr13_ratio, row$refRatio)
  }
  else (chr13_ratio = c(chr13_ratio, row$altRatio))
}
chr13_data = data.frame(chr="chr13", mean=mean(chr13_ratio), median=median(chr13_ratio))
write.table(chr13_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr13.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr14
chr14 = read.csv(paste(id, "_", tissue, "_chr14.tsv", sep=""), sep = '\t')
colnames(chr14) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr14$refRatio = chr14$refCount/chr14$totalCount
chr14$altRatio = chr14$altCount/chr14$totalCount

chr14_ratio = c()
for (i in 1:nrow(chr14)) {
  row = chr14[i,]
  if (row$refRatio > row$altRatio) {
    chr14_ratio = c(chr14_ratio, row$refRatio)
  }
  else (chr14_ratio = c(chr14_ratio, row$altRatio))
}
chr14_data = data.frame(chr="chr14", mean=mean(chr14_ratio), median=median(chr14_ratio))
write.table(chr14_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr14.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr15
chr15 = read.csv(paste(id, "_", tissue, "_chr15.tsv", sep=""), sep = '\t')
colnames(chr15) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr15$refRatio = chr15$refCount/chr15$totalCount
chr15$altRatio = chr15$altCount/chr15$totalCount

chr15_ratio = c()
for (i in 1:nrow(chr15)) {
  row = chr15[i,]
  if (row$refRatio > row$altRatio) {
    chr15_ratio = c(chr15_ratio, row$refRatio)
  }
  else (chr15_ratio = c(chr15_ratio, row$altRatio))
}
chr15_data = data.frame(chr="chr15", mean=mean(chr15_ratio), median=median(chr15_ratio))
write.table(chr15_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr15.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr16
chr16 = read.csv(paste(id, "_", tissue, "_chr16.tsv", sep=""), sep = '\t')
colnames(chr16) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr16$refRatio = chr16$refCount/chr16$totalCount
chr16$altRatio = chr16$altCount/chr16$totalCount

chr16_ratio = c()
for (i in 1:nrow(chr16)) {
  row = chr16[i,]
  if (row$refRatio > row$altRatio) {
    chr16_ratio = c(chr16_ratio, row$refRatio)
  }
  else (chr16_ratio = c(chr16_ratio, row$altRatio))
}
chr16_data = data.frame(chr="chr16", mean=mean(chr16_ratio), median=median(chr16_ratio))
write.table(chr16_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr16.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr17
chr17 = read.csv(paste(id, "_", tissue, "_chr17.tsv", sep=""), sep = '\t')
colnames(chr17) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr17$refRatio = chr17$refCount/chr17$totalCount
chr17$altRatio = chr17$altCount/chr17$totalCount

chr17_ratio = c()
for (i in 1:nrow(chr17)) {
  row = chr17[i,]
  if (row$refRatio > row$altRatio) {
    chr17_ratio = c(chr17_ratio, row$refRatio)
  }
  else (chr17_ratio = c(chr17_ratio, row$altRatio))
}
chr17_data = data.frame(chr="chr17", mean=mean(chr17_ratio), median=median(chr17_ratio))
write.table(chr17_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr17.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr18
chr18 = read.csv(paste(id, "_", tissue, "_chr18.tsv", sep=""), sep = '\t')
colnames(chr18) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr18$refRatio = chr18$refCount/chr18$totalCount
chr18$altRatio = chr18$altCount/chr18$totalCount

chr18_ratio = c()
for (i in 1:nrow(chr18)) {
  row = chr18[i,]
  if (row$refRatio > row$altRatio) {
    chr18_ratio = c(chr18_ratio, row$refRatio)
  }
  else (chr18_ratio = c(chr18_ratio, row$altRatio))
}
chr18_data = data.frame(chr="chr18", mean=mean(chr18_ratio), median=median(chr18_ratio))
write.table(chr18_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr18.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr19
chr19 = read.csv(paste(id, "_", tissue, "_chr19.tsv", sep=""), sep = '\t')
colnames(chr19) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr19$refRatio = chr19$refCount/chr19$totalCount
chr19$altRatio = chr19$altCount/chr19$totalCount

chr19_ratio = c()
for (i in 1:nrow(chr19)) {
  row = chr19[i,]
  if (row$refRatio > row$altRatio) {
    chr19_ratio = c(chr19_ratio, row$refRatio)
  }
  else (chr19_ratio = c(chr19_ratio, row$altRatio))
}
chr19_data = data.frame(chr="chr19", mean=mean(chr19_ratio), median=median(chr19_ratio))
write.table(chr19_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr19.tsv", sep=""), quote = F, sep = '\t', row.names = F)

# chr20
chr20 = read.csv(paste(id, "_", tissue, "_chr20.tsv", sep=""), sep = '\t')
colnames(chr20) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
chr20$refRatio = chr20$refCount/chr20$totalCount
chr20$altRatio = chr20$altCount/chr20$totalCount

chr20_ratio = c()
for (i in 1:nrow(chr20)) {
  row = chr20[i,]
  if (row$refRatio > row$altRatio) {
    chr20_ratio = c(chr20_ratio, row$refRatio)
  }
  else (chr20_ratio = c(chr20_ratio, row$altRatio))
}
chr20_data = data.frame(chr="chr20", mean=mean(chr20_ratio), median=median(chr20_ratio))
write.table(chr20_data, paste("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/", id, "_", tissue, "_mean_median_ratio_chr20.tsv", sep=""), quote = F, sep = '\t', row.names = F)