
#setwd("C:/Users/brend/Desktop/snakemake_workshop/gene_analysis/activation_status")
#install.packages("plyr", "dplyr")
library(dplyr)
library(plyr)

'%!in%' <- function(x,y)!('%in%'(x,y))

GTEX <- c("GTEX_adrenal", "GTEX_gonads", "GTEX_heart", "GTEX_liver", "GTEX_lung", "GTEX_pituitary")
PT <- c("PT_adrenal", "PT_gonads", "PT_heart", "PT_liver", "PT_lung", "PT_pituitary")

GTEX_adrenal    <- read.delim("GTEX_adrenal_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)
GTEX_gonads     <- read.delim("GTEX_gonads_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)
GTEX_heart      <- read.delim("GTEX_heart_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)
GTEX_liver      <- read.delim("GTEX_liver_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)
GTEX_lung       <- read.delim("GTEX_lung_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)
GTEX_pituitary  <- read.delim("GTEX_pituitary_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)
PT_adrenal      <- read.delim("PT_adrenal_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)
PT_gonads       <- read.delim("PT_gonads_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)
PT_heart        <- read.delim("PT_heart_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)
PT_liver        <- read.delim("PT_liver_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)
PT_lung         <- read.delim("PT_lung_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)
PT_pituitary    <- read.delim("PT_pituitary_allele_balance_counts.tsv", header=F, stringsAsFactors = TRUE)

GTEX_adrenal    <- transform(GTEX_adrenal, inactive = V8 >= 5)
GTEX_gonads     <- transform(GTEX_gonads, inactive = V8 >= 5)
GTEX_heart      <- transform(GTEX_heart, inactive = V8 >= 5)
GTEX_liver      <- transform(GTEX_liver, inactive = V8 >= 5)
GTEX_lung       <- transform(GTEX_lung, inactive = V8 >= 5)
GTEX_pituitary  <- transform(GTEX_pituitary, inactive = V8 >= 5)
PT_adrenal      <- transform(PT_adrenal, inactive = V16 >= 10)
PT_gonads       <- transform(PT_gonads, inactive = V14 >= 10)
PT_heart        <- transform(PT_heart, inactive = V16 >= 10)
PT_liver        <- transform(PT_liver, inactive = V16 >= 10)
PT_lung         <- transform(PT_lung, inactive = V14 >= 10)
PT_pituitary    <- transform(PT_pituitary,  inactive = V16 >= 10)

GTEX_adrenal    <- transform(GTEX_adrenal, escape = V8 <= 2)
GTEX_gonads     <- transform(GTEX_gonads, escape = V8 <= 2)
GTEX_heart      <- transform(GTEX_heart, escape = V8 <= 2)
GTEX_liver      <- transform(GTEX_liver, escape = V8 <= 2)
GTEX_lung       <- transform(GTEX_lung, escape = V8 <= 2)
GTEX_pituitary  <- transform(GTEX_pituitary, escape = V8 <= 2)
PT_adrenal      <- transform(PT_adrenal, escape = V16 <= 4)
PT_gonads       <- transform(PT_gonads, escape = V14 <= 4)          #missing 2 individuals
PT_heart        <- transform(PT_heart, escape = V16 <= 4)
PT_liver        <- transform(PT_liver, escape = V16 <= 4)
PT_lung         <- transform(PT_lung, escape = V14 <= 4)            #missing 2 individuals
PT_pituitary    <- transform(PT_pituitary,  escape = V16 <= 4)

write.table(GTEX_adrenal, file = "R_output/GTEX_adrenal.tsv", col.names=T, row.names = F, sep="\t")
write.table(GTEX_gonads, file = "R_output/GTEX_gonads.tsv", col.names=T, row.names = F, sep="\t")
write.table(GTEX_heart, file = "R_output/GTEX_heart.tsv", col.names=T, row.names = F, sep="\t")
write.table(GTEX_liver, file = "R_output/GTEX_liver.tsv", col.names=T, row.names = F, sep="\t")
write.table(GTEX_lung, file = "R_output/GTEX_lung.tsv", col.names=T, row.names = F, sep="\t")
write.table(GTEX_pituitary, file = "R_output/GTEX_pituitary.tsv", col.names=T, row.names = F, sep="\t")

write.table(PT_adrenal, file = "R_output/Cayo_adrenal.tsv", col.names=T, row.names = F, sep="\t")
write.table(PT_gonads, file = "R_output/Cayo_gonads.tsv", col.names=T, row.names = F, sep="\t")
write.table(PT_heart, file = "R_output/Cayo_heart.tsv", col.names=T, row.names = F, sep="\t")
write.table(PT_liver, file = "R_output/Cayo_liver.tsv", col.names=T, row.names = F, sep="\t")
write.table(PT_lung, file = "R_output/Cayo_lung.tsv", col.names=T, row.names = F, sep="\t")
write.table(PT_pituitary, file = "R_output/Cayo_pituitary.tsv", col.names=T, row.names = F, sep="\t")
