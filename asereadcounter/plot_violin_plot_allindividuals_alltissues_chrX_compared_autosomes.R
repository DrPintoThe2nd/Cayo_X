library(ggplot2)

setwd("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/")

# ---------
# autosomes
# ---------
# heart
heart_PT104219_chr13 = read.table("compare_X_A/results/PT104219_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT104219_chr13) = c("chr", "mean", "median")
heart_PT104219_chr13$tissue='heart'
heart_PT104219_chr13$id='PT104219'
heart_PT104219_chr13$chr='chr13'

heart_PT105290_chr13 = read.table("compare_X_A/results/PT105290_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT105290_chr13) = c("chr", "mean", "median")
heart_PT105290_chr13$tissue='heart'
heart_PT105290_chr13$id='PT105290'
heart_PT105290_chr13$chr='chr13'

heart_PT102485_chr13 = read.table("compare_X_A/results/PT102485_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT102485_chr13) = c("chr", "mean", "median")
heart_PT102485_chr13$tissue='heart'
heart_PT102485_chr13$id='PT102485'
heart_PT102485_chr13$chr='chr13'

heart_PT101159_chr13 = read.table("compare_X_A/results/PT101159_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT101159_chr13) = c("chr", "mean", "median")
heart_PT101159_chr13$tissue='heart'
heart_PT101159_chr13$id='PT101159'
heart_PT101159_chr13$chr='chr13'

heart_PT105239_chr13 = read.table("compare_X_A/results/PT105239_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT105239_chr13) = c("chr", "mean", "median")
heart_PT105239_chr13$tissue='heart'
heart_PT105239_chr13$id='PT105239'
heart_PT105239_chr13$chr='chr13'

heart_PT104883_chr13 = read.table("compare_X_A/results/PT104883_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT104883_chr13) = c("chr", "mean", "median")
heart_PT104883_chr13$tissue='heart'
heart_PT104883_chr13$id='PT104883'
heart_PT104883_chr13$chr='chr13'

heart_PT103811_chr13 = read.table("compare_X_A/results/PT103811_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT103811_chr13) = c("chr", "mean", "median")
heart_PT103811_chr13$tissue='heart'
heart_PT103811_chr13$id='PT103811'
heart_PT103811_chr13$chr='chr13'

heart_PT103352_chr13 = read.table("compare_X_A/results/PT103352_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT103352_chr13) = c("chr", "mean", "median")
heart_PT103352_chr13$tissue='heart'
heart_PT103352_chr13$id='PT103352'
heart_PT103352_chr13$chr='chr13'

heart_PT102842_chr13 = read.table("compare_X_A/results/PT102842_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT102842_chr13) = c("chr", "mean", "median")
heart_PT102842_chr13$tissue='heart'
heart_PT102842_chr13$id='PT102842'
heart_PT102842_chr13$chr='chr13'

heart_PT101210_chr13 = read.table("compare_X_A/results/PT101210_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT101210_chr13) = c("chr", "mean", "median")
heart_PT101210_chr13$tissue='heart'
heart_PT101210_chr13$id='PT101210'
heart_PT101210_chr13$chr='chr13'

heart_PT103046_chr13 = read.table("compare_X_A/results/PT103046_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT103046_chr13) = c("chr", "mean", "median")
heart_PT103046_chr13$tissue='heart'
heart_PT103046_chr13$id='PT103046'
heart_PT103046_chr13$chr='chr13'

heart_PT104733_chr13 = read.table("compare_X_A/results/PT104733_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT104733_chr13) = c("chr", "mean", "median")
heart_PT104733_chr13$tissue='heart'
heart_PT104733_chr13$id='PT104733'
heart_PT104733_chr13$chr='chr13'

heart_PT103760_chr13 = read.table("compare_X_A/results/PT103760_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT103760_chr13) = c("chr", "mean", "median")
heart_PT103760_chr13$tissue='heart'
heart_PT103760_chr13$id='PT103760'
heart_PT103760_chr13$chr='chr13'

heart_PT104270_chr13 = read.table("compare_X_A/results/PT104270_heart_mean_median_ratio_chr13.tsv", header=T)
colnames(heart_PT104270_chr13) = c("chr", "mean", "median")
heart_PT104270_chr13$tissue='heart'
heart_PT104270_chr13$id='PT104270'
heart_PT104270_chr13$chr='chr13'

# lung
lung_PT104219_chr13 = read.table("compare_X_A/results/PT104219_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT104219_chr13) = c("chr", "mean", "median")
lung_PT104219_chr13$tissue='lung'
lung_PT104219_chr13$id='PT104219'
lung_PT104219_chr13$chr='chr13'

lung_PT105290_chr13 = read.table("compare_X_A/results/PT105290_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT105290_chr13) = c("chr", "mean", "median")
lung_PT105290_chr13$tissue='lung'
lung_PT105290_chr13$id='PT105290'
lung_PT105290_chr13$chr='chr13'

lung_PT102485_chr13 = read.table("compare_X_A/results/PT102485_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT102485_chr13) = c("chr", "mean", "median")
lung_PT102485_chr13$tissue='lung'
lung_PT102485_chr13$id='PT102485'
lung_PT102485_chr13$chr='chr13'

lung_PT101159_chr13 = read.table("compare_X_A/results/PT101159_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT101159_chr13) = c("chr", "mean", "median")
lung_PT101159_chr13$tissue='lung'
lung_PT101159_chr13$id='PT101159'
lung_PT101159_chr13$chr='chr13'

lung_PT105239_chr13 = read.table("compare_X_A/results/PT105239_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT105239_chr13) = c("chr", "mean", "median")
lung_PT105239_chr13$tissue='lung'
lung_PT105239_chr13$id='PT105239'
lung_PT105239_chr13$chr='chr13'

lung_PT104883_chr13 = read.table("compare_X_A/results/PT104883_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT104883_chr13) = c("chr", "mean", "median")
lung_PT104883_chr13$tissue='lung'
lung_PT104883_chr13$id='PT104883'
lung_PT104883_chr13$chr='chr13'

lung_PT103811_chr13 = read.table("compare_X_A/results/PT103811_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT103811_chr13) = c("chr", "mean", "median")
lung_PT103811_chr13$tissue='lung'
lung_PT103811_chr13$id='PT103811'
lung_PT103811_chr13$chr='chr13'

lung_PT103352_chr13 = read.table("compare_X_A/results/PT103352_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT103352_chr13) = c("chr", "mean", "median")
lung_PT103352_chr13$tissue='lung'
lung_PT103352_chr13$id='PT103352'
lung_PT103352_chr13$chr='chr13'

lung_PT102842_chr13 = read.table("compare_X_A/results/PT102842_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT102842_chr13) = c("chr", "mean", "median")
lung_PT102842_chr13$tissue='lung'
lung_PT102842_chr13$id='PT102842'
lung_PT102842_chr13$chr='chr13'

lung_PT101210_chr13 = read.table("compare_X_A/results/PT101210_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT101210_chr13) = c("chr", "mean", "median")
lung_PT101210_chr13$tissue='lung'
lung_PT101210_chr13$id='PT101210'
lung_PT101210_chr13$chr='chr13'

lung_PT103760_chr13 = read.table("compare_X_A/results/PT103760_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT103760_chr13) = c("chr", "mean", "median")
lung_PT103760_chr13$tissue='lung'
lung_PT103760_chr13$id='PT103760'
lung_PT103760_chr13$chr='chr13'

lung_PT104270_chr13 = read.table("compare_X_A/results/PT104270_lung_mean_median_ratio_chr13.tsv", header=T)
colnames(lung_PT104270_chr13) = c("chr", "mean", "median")
lung_PT104270_chr13$tissue='lung'
lung_PT104270_chr13$id='PT104270'
lung_PT104270_chr13$chr='chr13'

# liver
liver_PT104219_chr13 = read.table("compare_X_A/results/PT104219_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT104219_chr13) = c("chr", "mean", "median")
liver_PT104219_chr13$tissue='liver'
liver_PT104219_chr13$id='PT104219'
liver_PT104219_chr13$chr='chr13'

liver_PT105290_chr13 = read.table("compare_X_A/results/PT105290_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT105290_chr13) = c("chr", "mean", "median")
liver_PT105290_chr13$tissue='liver'
liver_PT105290_chr13$id='PT105290'
liver_PT105290_chr13$chr='chr13'

liver_PT102485_chr13 = read.table("compare_X_A/results/PT102485_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT102485_chr13) = c("chr", "mean", "median")
liver_PT102485_chr13$tissue='liver'
liver_PT102485_chr13$id='PT102485'
liver_PT102485_chr13$chr='chr13'

liver_PT101159_chr13 = read.table("compare_X_A/results/PT101159_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT101159_chr13) = c("chr", "mean", "median")
liver_PT101159_chr13$tissue='liver'
liver_PT101159_chr13$id='PT101159'
liver_PT101159_chr13$chr='chr13'

liver_PT105239_chr13 = read.table("compare_X_A/results/PT105239_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT105239_chr13) = c("chr", "mean", "median")
liver_PT105239_chr13$tissue='liver'
liver_PT105239_chr13$id='PT105239'
liver_PT105239_chr13$chr='chr13'

liver_PT104883_chr13 = read.table("compare_X_A/results/PT104883_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT104883_chr13) = c("chr", "mean", "median")
liver_PT104883_chr13$tissue='liver'
liver_PT104883_chr13$id='PT104883'
liver_PT104883_chr13$chr='chr13'

liver_PT103811_chr13 = read.table("compare_X_A/results/PT103811_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT103811_chr13) = c("chr", "mean", "median")
liver_PT103811_chr13$tissue='liver'
liver_PT103811_chr13$id='PT103811'
liver_PT103811_chr13$chr='chr13'

liver_PT103352_chr13 = read.table("compare_X_A/results/PT103352_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT103352_chr13) = c("chr", "mean", "median")
liver_PT103352_chr13$tissue='liver'
liver_PT103352_chr13$id='PT103352'
liver_PT103352_chr13$chr='chr13'

liver_PT102842_chr13 = read.table("compare_X_A/results/PT102842_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT102842_chr13) = c("chr", "mean", "median")
liver_PT102842_chr13$tissue='liver'
liver_PT102842_chr13$id='PT102842'
liver_PT102842_chr13$chr='chr13'

liver_PT101210_chr13 = read.table("compare_X_A/results/PT101210_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT101210_chr13) = c("chr", "mean", "median")
liver_PT101210_chr13$tissue='liver'
liver_PT101210_chr13$id='PT101210'
liver_PT101210_chr13$chr='chr13'

liver_PT103046_chr13 = read.table("compare_X_A/results/PT103046_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT103046_chr13) = c("chr", "mean", "median")
liver_PT103046_chr13$tissue='liver'
liver_PT103046_chr13$id='PT103046'
liver_PT103046_chr13$chr='chr13'

liver_PT104733_chr13 = read.table("compare_X_A/results/PT104733_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT104733_chr13) = c("chr", "mean", "median")
liver_PT104733_chr13$tissue='liver'
liver_PT104733_chr13$id='PT104733'
liver_PT104733_chr13$chr='chr13'

liver_PT103760_chr13 = read.table("compare_X_A/results/PT103760_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT103760_chr13) = c("chr", "mean", "median")
liver_PT103760_chr13$tissue='liver'
liver_PT103760_chr13$id='PT103760'
liver_PT103760_chr13$chr='chr13'

liver_PT104270_chr13 = read.table("compare_X_A/results/PT104270_liver_mean_median_ratio_chr13.tsv", header=T)
colnames(liver_PT104270_chr13) = c("chr", "mean", "median")
liver_PT104270_chr13$tissue='liver'
liver_PT104270_chr13$id='PT104270'
liver_PT104270_chr13$chr='chr13'

# gonads
gonads_PT104219_chr13 = read.table("compare_X_A/results/PT104219_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT104219_chr13) = c("chr", "mean", "median")
gonads_PT104219_chr13$tissue='ovary'
gonads_PT104219_chr13$id='PT104219'
gonads_PT104219_chr13$chr='chr13'

gonads_PT105290_chr13 = read.table("compare_X_A/results/PT105290_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT105290_chr13) = c("chr", "mean", "median")
gonads_PT105290_chr13$tissue='ovary'
gonads_PT105290_chr13$id='PT105290'
gonads_PT105290_chr13$chr='chr13'

gonads_PT102485_chr13 = read.table("compare_X_A/results/PT102485_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT102485_chr13) = c("chr", "mean", "median")
gonads_PT102485_chr13$tissue='ovary'
gonads_PT102485_chr13$id='PT102485'
gonads_PT102485_chr13$chr='chr13'

gonads_PT101159_chr13 = read.table("compare_X_A/results/PT101159_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT101159_chr13) = c("chr", "mean", "median")
gonads_PT101159_chr13$tissue='ovary'
gonads_PT101159_chr13$id='PT101159'
gonads_PT101159_chr13$chr='chr13'

gonads_PT105239_chr13 = read.table("compare_X_A/results/PT105239_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT105239_chr13) = c("chr", "mean", "median")
gonads_PT105239_chr13$tissue='ovary'
gonads_PT105239_chr13$id='PT105239'
gonads_PT105239_chr13$chr='chr13'

gonads_PT104883_chr13 = read.table("compare_X_A/results/PT104883_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT104883_chr13) = c("chr", "mean", "median")
gonads_PT104883_chr13$tissue='ovary'
gonads_PT104883_chr13$id='PT104883'
gonads_PT104883_chr13$chr='chr13'

gonads_PT103811_chr13 = read.table("compare_X_A/results/PT103811_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT103811_chr13) = c("chr", "mean", "median")
gonads_PT103811_chr13$tissue='ovary'
gonads_PT103811_chr13$id='PT103811'
gonads_PT103811_chr13$chr='chr13'

gonads_PT103352_chr13 = read.table("compare_X_A/results/PT103352_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT103352_chr13) = c("chr", "mean", "median")
gonads_PT103352_chr13$tissue='ovary'
gonads_PT103352_chr13$id='PT103352'
gonads_PT103352_chr13$chr='chr13'

gonads_PT102842_chr13 = read.table("compare_X_A/results/PT102842_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT102842_chr13) = c("chr", "mean", "median")
gonads_PT102842_chr13$tissue='ovary'
gonads_PT102842_chr13$id='PT102842'
gonads_PT102842_chr13$chr='chr13'

gonads_PT101210_chr13 = read.table("compare_X_A/results/PT101210_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT101210_chr13) = c("chr", "mean", "median")
gonads_PT101210_chr13$tissue='ovary'
gonads_PT101210_chr13$id='PT101210'
gonads_PT101210_chr13$chr='chr13'

gonads_PT103046_chr13 = read.table("compare_X_A/results/PT103046_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT103046_chr13) = c("chr", "mean", "median")
gonads_PT103046_chr13$tissue='ovary'
gonads_PT103046_chr13$id='PT103046'
gonads_PT103046_chr13$chr='chr13'

gonads_PT104733_chr13 = read.table("compare_X_A/results/PT104733_gonads_mean_median_ratio_chr13.tsv", header=T)
colnames(gonads_PT104733_chr13) = c("chr", "mean", "median")
gonads_PT104733_chr13$tissue='ovary'
gonads_PT104733_chr13$id='PT104733'
gonads_PT104733_chr13$chr='chr13'

# adrenal
adrenal_PT104219_chr13 = read.table("compare_X_A/results/PT104219_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT104219_chr13) = c("chr", "mean", "median")
adrenal_PT104219_chr13$tissue='adrenal'
adrenal_PT104219_chr13$id='PT104219'
adrenal_PT104219_chr13$chr='chr13'

adrenal_PT105290_chr13 = read.table("compare_X_A/results/PT105290_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT105290_chr13) = c("chr", "mean", "median")
adrenal_PT105290_chr13$tissue='adrenal'
adrenal_PT105290_chr13$id='PT105290'
adrenal_PT105290_chr13$chr='chr13'

adrenal_PT102485_chr13 = read.table("compare_X_A/results/PT102485_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT102485_chr13) = c("chr", "mean", "median")
adrenal_PT102485_chr13$tissue='adrenal'
adrenal_PT102485_chr13$id='PT102485'
adrenal_PT102485_chr13$chr='chr13'

adrenal_PT101159_chr13 = read.table("compare_X_A/results/PT101159_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT101159_chr13) = c("chr", "mean", "median")
adrenal_PT101159_chr13$tissue='adrenal'
adrenal_PT101159_chr13$id='PT101159'
adrenal_PT101159_chr13$chr='chr13'

adrenal_PT105239_chr13 = read.table("compare_X_A/results/PT105239_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT105239_chr13) = c("chr", "mean", "median")
adrenal_PT105239_chr13$tissue='adrenal'
adrenal_PT105239_chr13$id='PT105239'
adrenal_PT105239_chr13$chr='chr13'

adrenal_PT104883_chr13 = read.table("compare_X_A/results/PT104883_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT104883_chr13) = c("chr", "mean", "median")
adrenal_PT104883_chr13$tissue='adrenal'
adrenal_PT104883_chr13$id='PT104883'
adrenal_PT104883_chr13$chr='chr13'

adrenal_PT103811_chr13 = read.table("compare_X_A/results/PT103811_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT103811_chr13) = c("chr", "mean", "median")
adrenal_PT103811_chr13$tissue='adrenal'
adrenal_PT103811_chr13$id='PT103811'
adrenal_PT103811_chr13$chr='chr13'

adrenal_PT103352_chr13 = read.table("compare_X_A/results/PT103352_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT103352_chr13) = c("chr", "mean", "median")
adrenal_PT103352_chr13$tissue='adrenal'
adrenal_PT103352_chr13$id='PT103352'
adrenal_PT103352_chr13$chr='chr13'

adrenal_PT102842_chr13 = read.table("compare_X_A/results/PT102842_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT102842_chr13) = c("chr", "mean", "median")
adrenal_PT102842_chr13$tissue='adrenal'
adrenal_PT102842_chr13$id='PT102842'
adrenal_PT102842_chr13$chr='chr13'

adrenal_PT101210_chr13 = read.table("compare_X_A/results/PT101210_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT101210_chr13) = c("chr", "mean", "median")
adrenal_PT101210_chr13$tissue='adrenal'
adrenal_PT101210_chr13$id='PT101210'
adrenal_PT101210_chr13$chr='chr13'

adrenal_PT103046_chr13 = read.table("compare_X_A/results/PT103046_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT103046_chr13) = c("chr", "mean", "median")
adrenal_PT103046_chr13$tissue='adrenal'
adrenal_PT103046_chr13$id='PT103046'
adrenal_PT103046_chr13$chr='chr13'

adrenal_PT104733_chr13 = read.table("compare_X_A/results/PT104733_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT104733_chr13) = c("chr", "mean", "median")
adrenal_PT104733_chr13$tissue='adrenal'
adrenal_PT104733_chr13$id='PT104733'
adrenal_PT104733_chr13$chr='chr13'

adrenal_PT103760_chr13 = read.table("compare_X_A/results/PT103760_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT103760_chr13) = c("chr", "mean", "median")
adrenal_PT103760_chr13$tissue='adrenal'
adrenal_PT103760_chr13$id='PT103760'
adrenal_PT103760_chr13$chr='chr13'

adrenal_PT104270_chr13 = read.table("compare_X_A/results/PT104270_adrenal_mean_median_ratio_chr13.tsv", header=T)
colnames(adrenal_PT104270_chr13) = c("chr", "mean", "median")
adrenal_PT104270_chr13$tissue='adrenal'
adrenal_PT104270_chr13$id='PT104270'
adrenal_PT104270_chr13$chr='chr13'

# pituitary
pituitary_PT104219_chr13 = read.table("compare_X_A/results/PT104219_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT104219_chr13) = c("chr", "mean", "median")
pituitary_PT104219_chr13$tissue='pituitary'
pituitary_PT104219_chr13$id='PT104219'
pituitary_PT104219_chr13$chr='chr13'

pituitary_PT105290_chr13 = read.table("compare_X_A/results/PT105290_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT105290_chr13) = c("chr", "mean", "median")
pituitary_PT105290_chr13$tissue='pituitary'
pituitary_PT105290_chr13$id='PT105290'
pituitary_PT105290_chr13$chr='chr13'

pituitary_PT102485_chr13 = read.table("compare_X_A/results/PT102485_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT102485_chr13) = c("chr", "mean", "median")
pituitary_PT102485_chr13$tissue='pituitary'
pituitary_PT102485_chr13$id='PT102485'
pituitary_PT102485_chr13$chr='chr13'

pituitary_PT101159_chr13 = read.table("compare_X_A/results/PT101159_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT101159_chr13) = c("chr", "mean", "median")
pituitary_PT101159_chr13$tissue='pituitary'
pituitary_PT101159_chr13$id='PT101159'
pituitary_PT101159_chr13$chr='chr13'

pituitary_PT105239_chr13 = read.table("compare_X_A/results/PT105239_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT105239_chr13) = c("chr", "mean", "median")
pituitary_PT105239_chr13$tissue='pituitary'
pituitary_PT105239_chr13$id='PT105239'
pituitary_PT105239_chr13$chr='chr13'

pituitary_PT104883_chr13 = read.table("compare_X_A/results/PT104883_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT104883_chr13) = c("chr", "mean", "median")
pituitary_PT104883_chr13$tissue='pituitary'
pituitary_PT104883_chr13$id='PT104883'
pituitary_PT104883_chr13$chr='chr13'

pituitary_PT103811_chr13 = read.table("compare_X_A/results/PT103811_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT103811_chr13) = c("chr", "mean", "median")
pituitary_PT103811_chr13$tissue='pituitary'
pituitary_PT103811_chr13$id='PT103811'
pituitary_PT103811_chr13$chr='chr13'

pituitary_PT103352_chr13 = read.table("compare_X_A/results/PT103352_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT103352_chr13) = c("chr", "mean", "median")
pituitary_PT103352_chr13$tissue='pituitary'
pituitary_PT103352_chr13$id='PT103352'
pituitary_PT103352_chr13$chr='chr13'

pituitary_PT102842_chr13 = read.table("compare_X_A/results/PT102842_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT102842_chr13) = c("chr", "mean", "median")
pituitary_PT102842_chr13$tissue='pituitary'
pituitary_PT102842_chr13$id='PT102842'
pituitary_PT102842_chr13$chr='chr13'

pituitary_PT101210_chr13 = read.table("compare_X_A/results/PT101210_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT101210_chr13) = c("chr", "mean", "median")
pituitary_PT101210_chr13$tissue='pituitary'
pituitary_PT101210_chr13$id='PT101210'
pituitary_PT101210_chr13$chr='chr13'

pituitary_PT103046_chr13 = read.table("compare_X_A/results/PT103046_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT103046_chr13) = c("chr", "mean", "median")
pituitary_PT103046_chr13$tissue='pituitary'
pituitary_PT103046_chr13$id='PT103046'
pituitary_PT103046_chr13$chr='chr13'

pituitary_PT104733_chr13 = read.table("compare_X_A/results/PT104733_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT104733_chr13) = c("chr", "mean", "median")
pituitary_PT104733_chr13$tissue='pituitary'
pituitary_PT104733_chr13$id='PT104733'
pituitary_PT104733_chr13$chr='chr13'

pituitary_PT103760_chr13 = read.table("compare_X_A/results/PT103760_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT103760_chr13) = c("chr", "mean", "median")
pituitary_PT103760_chr13$tissue='pituitary'
pituitary_PT103760_chr13$id='PT103760'
pituitary_PT103760_chr13$chr='chr13'

pituitary_PT104270_chr13 = read.table("compare_X_A/results/PT104270_pituitary_mean_median_ratio_chr13.tsv", header=T)
colnames(pituitary_PT104270_chr13) = c("chr", "mean", "median")
pituitary_PT104270_chr13$tissue='pituitary'
pituitary_PT104270_chr13$id='PT104270'
pituitary_PT104270_chr13$chr='chr13'


# combine_data = rbind(heart_PT104219, heart_PT105290, heart_PT102485, heart_PT101159, heart_PT105239, heart_PT104883, heart_PT103811, heart_PT103352, heart_PT102842, heart_PT101210, heart_PT103046, heart_PT104733, heart_PT103760, heart_PT104270, lung_PT104219, lung_PT105290, lung_PT102485, lung_PT101159, lung_PT105239, lung_PT104883, lung_PT103811, lung_PT103352, lung_PT102842, lung_PT101210, lung_PT103760, lung_PT104270, liver_PT104219, liver_PT105290, liver_PT102485, liver_PT101159, liver_PT105239, liver_PT104883, liver_PT103811, liver_PT103352, liver_PT102842, liver_PT101210, liver_PT103046, liver_PT104733, liver_PT103760, liver_PT104270, gonads_PT104219, gonads_PT105290, gonads_PT102485, gonads_PT101159, gonads_PT105239, gonads_PT104883, gonads_PT103811, gonads_PT103352, gonads_PT102842, gonads_PT101210, gonads_PT103046, gonads_PT104733, adrenal_PT104219, adrenal_PT105290, adrenal_PT102485, adrenal_PT101159, adrenal_PT105239, adrenal_PT104883, adrenal_PT103811, adrenal_PT103352, adrenal_PT102842, adrenal_PT101210, adrenal_PT103046, adrenal_PT104733, adrenal_PT103760, adrenal_PT104270, pituitary_PT104219, pituitary_PT105290, pituitary_PT102485, pituitary_PT101159, pituitary_PT105239, pituitary_PT104883, pituitary_PT103811, pituitary_PT103352, pituitary_PT102842, pituitary_PT101210, pituitary_PT103046, pituitary_PT104733, pituitary_PT103760, pituitary_PT104270)
#
# p = ggplot(combine_data, aes(x=tissue, y=median)) +
#   geom_violin(color="darkgray") +
#   geom_boxplot(width=0.1, color="darkgray") +
#   theme_bw() +
#   coord_cartesian(ylim=c(0.5,1)) +
#   labs(x='Tissue', y="Allele balance") +
#   theme(plot.title = element_text(hjust = 0.5, size=18), legend.position = "none", axis.text = element_text(size=16), axis.title = element_text(size=18))

# -----
# chrX
# -----
# heart
heart_PT104219_chrX = read.table("compare_X_A/results/PT104219_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT104219_chrX$tissue='heart'
heart_PT104219_chrX$id='PT104219'
heart_PT104219_chrX$chr='chrX'

heart_PT105290_chrX = read.table("compare_X_A/results/PT105290_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT105290_chrX$tissue='heart'
heart_PT105290_chrX$id='PT105290'
heart_PT105290_chrX$chr='chrX'

heart_PT102485_chrX = read.table("compare_X_A/results/PT102485_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT102485_chrX$tissue='heart'
heart_PT102485_chrX$id='PT102485'
heart_PT102485_chrX$chr='chrX'

heart_PT101159_chrX = read.table("compare_X_A/results/PT101159_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT101159_chrX$tissue='heart'
heart_PT101159_chrX$id='PT101159'
heart_PT101159_chrX$chr='chrX'

heart_PT105239_chrX = read.table("compare_X_A/results/PT105239_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT105239_chrX$tissue='heart'
heart_PT105239_chrX$id='PT105239'
heart_PT105239_chrX$chr='chrX'

heart_PT104883_chrX = read.table("compare_X_A/results/PT104883_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT104883_chrX$tissue='heart'
heart_PT104883_chrX$id='PT104883'
heart_PT104883_chrX$chr='chrX'

heart_PT103811_chrX = read.table("compare_X_A/results/PT103811_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT103811_chrX$tissue='heart'
heart_PT103811_chrX$id='PT103811'
heart_PT103811_chrX$chr='chrX'

heart_PT103352_chrX = read.table("compare_X_A/results/PT103352_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT103352_chrX$tissue='heart'
heart_PT103352_chrX$id='PT103352'
heart_PT103352_chrX$chr='chrX'

heart_PT102842_chrX = read.table("compare_X_A/results/PT102842_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT102842_chrX$tissue='heart'
heart_PT102842_chrX$id='PT102842'
heart_PT102842_chrX$chr='chrX'

heart_PT101210_chrX = read.table("compare_X_A/results/PT101210_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT101210_chrX$tissue='heart'
heart_PT101210_chrX$id='PT101210'
heart_PT101210_chrX$chr='chrX'

heart_PT103046_chrX = read.table("compare_X_A/results/PT103046_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT103046_chrX$tissue='heart'
heart_PT103046_chrX$id='PT103046'
heart_PT103046_chrX$chr='chrX'

heart_PT104733_chrX = read.table("compare_X_A/results/PT104733_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT104733_chrX$tissue='heart'
heart_PT104733_chrX$id='PT104733'
heart_PT104733_chrX$chr='chrX'

heart_PT103760_chrX = read.table("compare_X_A/results/PT103760_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT103760_chrX$tissue='heart'
heart_PT103760_chrX$id='PT103760'
heart_PT103760_chrX$chr='chrX'

heart_PT104270_chrX = read.table("compare_X_A/results/PT104270_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT104270_chrX$tissue='heart'
heart_PT104270_chrX$id='PT104270'
heart_PT104270_chrX$chr='chrX'

# lung
lung_PT104219_chrX = read.table("compare_X_A/results/PT104219_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT104219_chrX$tissue='lung'
lung_PT104219_chrX$id='PT104219'
lung_PT104219_chrX$chr='chrX'

lung_PT105290_chrX = read.table("compare_X_A/results/PT105290_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT105290_chrX$tissue='lung'
lung_PT105290_chrX$id='PT105290'
lung_PT105290_chrX$chr='chrX'

lung_PT102485_chrX = read.table("compare_X_A/results/PT102485_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT102485_chrX$tissue='lung'
lung_PT102485_chrX$id='PT102485'
lung_PT102485_chrX$chr='chrX'

lung_PT101159_chrX = read.table("compare_X_A/results/PT101159_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT101159_chrX$tissue='lung'
lung_PT101159_chrX$id='PT101159'
lung_PT101159_chrX$chr='chrX'

lung_PT105239_chrX = read.table("compare_X_A/results/PT105239_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT105239_chrX$tissue='lung'
lung_PT105239_chrX$id='PT105239'
lung_PT105239_chrX$chr='chrX'

lung_PT104883_chrX = read.table("compare_X_A/results/PT104883_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT104883_chrX$tissue='lung'
lung_PT104883_chrX$id='PT104883'
lung_PT104883_chrX$chr='chrX'

lung_PT103811_chrX = read.table("compare_X_A/results/PT103811_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT103811_chrX$tissue='lung'
lung_PT103811_chrX$id='PT103811'
lung_PT103811_chrX$chr='chrX'

lung_PT103352_chrX = read.table("compare_X_A/results/PT103352_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT103352_chrX$tissue='lung'
lung_PT103352_chrX$id='PT103352'
lung_PT103352_chrX$chr='chrX'

lung_PT102842_chrX = read.table("compare_X_A/results/PT102842_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT102842_chrX$tissue='lung'
lung_PT102842_chrX$id='PT102842'
lung_PT102842_chrX$chr='chrX'

lung_PT101210_chrX = read.table("compare_X_A/results/PT101210_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT101210_chrX$tissue='lung'
lung_PT101210_chrX$id='PT101210'
lung_PT101210_chrX$chr='chrX'

lung_PT103760_chrX = read.table("compare_X_A/results/PT103760_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT103760_chrX$tissue='lung'
lung_PT103760_chrX$id='PT103760'
lung_PT103760_chrX$chr='chrX'

lung_PT104270_chrX = read.table("compare_X_A/results/PT104270_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT104270_chrX$tissue='lung'
lung_PT104270_chrX$id='PT104270'
lung_PT104270_chrX$chr='chrX'

# liver
liver_PT104219_chrX = read.table("compare_X_A/results/PT104219_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT104219_chrX$tissue='liver'
liver_PT104219_chrX$id='PT104219'
liver_PT104219_chrX$chr='chrX'

liver_PT105290_chrX = read.table("compare_X_A/results/PT105290_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT105290_chrX$tissue='liver'
liver_PT105290_chrX$id='PT105290'
liver_PT105290_chrX$chr='chrX'

liver_PT102485_chrX = read.table("compare_X_A/results/PT102485_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT102485_chrX$tissue='liver'
liver_PT102485_chrX$id='PT102485'
liver_PT102485_chrX$chr='chrX'

liver_PT101159_chrX = read.table("compare_X_A/results/PT101159_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT101159_chrX$tissue='liver'
liver_PT101159_chrX$id='PT101159'
liver_PT101159_chrX$chr='chrX'

liver_PT105239_chrX = read.table("compare_X_A/results/PT105239_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT105239_chrX$tissue='liver'
liver_PT105239_chrX$id='PT105239'
liver_PT105239_chrX$chr='chrX'

liver_PT104883_chrX = read.table("compare_X_A/results/PT104883_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT104883_chrX$tissue='liver'
liver_PT104883_chrX$id='PT104883'
liver_PT104883_chrX$chr='chrX'

liver_PT103811_chrX = read.table("compare_X_A/results/PT103811_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT103811_chrX$tissue='liver'
liver_PT103811_chrX$id='PT103811'
liver_PT103811_chrX$chr='chrX'

liver_PT103352_chrX = read.table("compare_X_A/results/PT103352_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT103352_chrX$tissue='liver'
liver_PT103352_chrX$id='PT103352'
liver_PT103352_chrX$chr='chrX'

liver_PT102842_chrX = read.table("compare_X_A/results/PT102842_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT102842_chrX$tissue='liver'
liver_PT102842_chrX$id='PT102842'
liver_PT102842_chrX$chr='chrX'

liver_PT101210_chrX = read.table("compare_X_A/results/PT101210_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT101210_chrX$tissue='liver'
liver_PT101210_chrX$id='PT101210'
liver_PT101210_chrX$chr='chrX'

liver_PT103046_chrX = read.table("compare_X_A/results/PT103046_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT103046_chrX$tissue='liver'
liver_PT103046_chrX$id='PT103046'
liver_PT103046_chrX$chr='chrX'

liver_PT104733_chrX = read.table("compare_X_A/results/PT104733_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT104733_chrX$tissue='liver'
liver_PT104733_chrX$id='PT104733'
liver_PT104733_chrX$chr='chrX'

liver_PT103760_chrX = read.table("compare_X_A/results/PT103760_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT103760_chrX$tissue='liver'
liver_PT103760_chrX$id='PT103760'
liver_PT103760_chrX$chr='chrX'

liver_PT104270_chrX = read.table("compare_X_A/results/PT104270_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT104270_chrX$tissue='liver'
liver_PT104270_chrX$id='PT104270'
liver_PT104270_chrX$chr='chrX'

# gonads
gonads_PT104219_chrX = read.table("compare_X_A/results/PT104219_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT104219_chrX$tissue='ovary'
gonads_PT104219_chrX$id='PT104219'
gonads_PT104219_chrX$chr='chrX'

gonads_PT105290_chrX = read.table("compare_X_A/results/PT105290_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT105290_chrX$tissue='ovary'
gonads_PT105290_chrX$id='PT105290'
gonads_PT105290_chrX$chr='chrX'

gonads_PT102485_chrX = read.table("compare_X_A/results/PT102485_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT102485_chrX$tissue='ovary'
gonads_PT102485_chrX$id='PT102485'
gonads_PT102485_chrX$chr='chrX'

gonads_PT101159_chrX = read.table("compare_X_A/results/PT101159_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT101159_chrX$tissue='ovary'
gonads_PT101159_chrX$id='PT101159'
gonads_PT101159_chrX$chr='chrX'

gonads_PT105239_chrX = read.table("compare_X_A/results/PT105239_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT105239_chrX$tissue='ovary'
gonads_PT105239_chrX$id='PT105239'
gonads_PT105239_chrX$chr='chrX'

gonads_PT104883_chrX = read.table("compare_X_A/results/PT104883_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT104883_chrX$tissue='ovary'
gonads_PT104883_chrX$id='PT104883'
gonads_PT104883_chrX$chr='chrX'

gonads_PT103811_chrX= read.table("compare_X_A/results/PT103811_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT103811_chrX$tissue='ovary'
gonads_PT103811_chrX$id='PT103811'
gonads_PT103811_chrX$chr='chrX'

gonads_PT103352_chrX = read.table("compare_X_A/results/PT103352_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT103352_chrX$tissue='ovary'
gonads_PT103352_chrX$id='PT103352'
gonads_PT103352_chrX$chr='chrX'

gonads_PT102842_chrX = read.table("compare_X_A/results/PT102842_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT102842_chrX$tissue='ovary'
gonads_PT102842_chrX$id='PT102842'
gonads_PT102842_chrX$chr='chrX'

gonads_PT101210_chrX = read.table("compare_X_A/results/PT101210_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT101210_chrX$tissue='ovary'
gonads_PT101210_chrX$id='PT101210'
gonads_PT101210_chrX$chr='chrX'

gonads_PT103046_chrX = read.table("compare_X_A/results/PT103046_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT103046_chrX$tissue='ovary'
gonads_PT103046_chrX$id='PT103046'
gonads_PT103046_chrX$chr='chrX'

gonads_PT104733_chrX = read.table("compare_X_A/results/PT104733_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT104733_chrX$tissue='ovary'
gonads_PT104733_chrX$id='PT104733'
gonads_PT104733_chrX$chr='chrX'

# adrenal
adrenal_PT104219_chrX = read.table("compare_X_A/results/PT104219_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT104219_chrX$tissue='adrenal'
adrenal_PT104219_chrX$id='PT104219'
adrenal_PT104219_chrX$chr='chrX'

adrenal_PT105290_chrX = read.table("compare_X_A/results/PT105290_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT105290_chrX$tissue='adrenal'
adrenal_PT105290_chrX$id='PT105290'
adrenal_PT105290_chrX$chr='chrX'

adrenal_PT102485_chrX = read.table("compare_X_A/results/PT102485_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT102485_chrX$tissue='adrenal'
adrenal_PT102485_chrX$id='PT102485'
adrenal_PT102485_chrX$chr='chrX'

adrenal_PT101159_chrX = read.table("compare_X_A/results/PT101159_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT101159_chrX$tissue='adrenal'
adrenal_PT101159_chrX$id='PT101159'
adrenal_PT101159_chrX$chr='chrX'

adrenal_PT105239_chrX = read.table("compare_X_A/results/PT105239_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT105239_chrX$tissue='adrenal'
adrenal_PT105239_chrX$id='PT105239'
adrenal_PT105239_chrX$chr='chrX'

adrenal_PT104883_chrX = read.table("compare_X_A/results/PT104883_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT104883_chrX$tissue='adrenal'
adrenal_PT104883_chrX$id='PT104883'
adrenal_PT104883_chrX$chr='chrX'

adrenal_PT103811_chrX = read.table("compare_X_A/results/PT103811_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT103811_chrX$tissue='adrenal'
adrenal_PT103811_chrX$id='PT103811'
adrenal_PT103811_chrX$chr='chrX'

adrenal_PT103352_chrX = read.table("compare_X_A/results/PT103352_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT103352_chrX$tissue='adrenal'
adrenal_PT103352_chrX$id='PT103352'
adrenal_PT103352_chrX$chr='chrX'

adrenal_PT102842_chrX = read.table("compare_X_A/results/PT102842_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT102842_chrX$tissue='adrenal'
adrenal_PT102842_chrX$id='PT102842'
adrenal_PT102842_chrX$chr='chrX'

adrenal_PT101210_chrX = read.table("compare_X_A/results/PT101210_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT101210_chrX$tissue='adrenal'
adrenal_PT101210_chrX$id='PT101210'
adrenal_PT101210_chrX$chr='chrX'

adrenal_PT103046_chrX = read.table("compare_X_A/results/PT103046_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT103046_chrX$tissue='adrenal'
adrenal_PT103046_chrX$id='PT103046'
adrenal_PT103046_chrX$chr='chrX'

adrenal_PT104733_chrX = read.table("compare_X_A/results/PT104733_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT104733_chrX$tissue='adrenal'
adrenal_PT104733_chrX$id='PT104733'
adrenal_PT104733_chrX$chr='chrX'

adrenal_PT103760_chrX = read.table("compare_X_A/results/PT103760_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT103760_chrX$tissue='adrenal'
adrenal_PT103760_chrX$id='PT103760'
adrenal_PT103760_chrX$chr='chrX'

adrenal_PT104270_chrX = read.table("compare_X_A/results/PT104270_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT104270_chrX$tissue='adrenal'
adrenal_PT104270_chrX$id='PT104270'
adrenal_PT104270_chrX$chr='chrX'

# pituitary
pituitary_PT104219_chrX = read.table("compare_X_A/results/PT104219_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT104219_chrX$tissue='pituitary'
pituitary_PT104219_chrX$id='PT104219'
pituitary_PT104219_chrX$chr='chrX'

pituitary_PT105290_chrX = read.table("compare_X_A/results/PT105290_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT105290_chrX$tissue='pituitary'
pituitary_PT105290_chrX$id='PT105290'
pituitary_PT105290_chrX$chr='chrX'

pituitary_PT102485_chrX = read.table("compare_X_A/results/PT102485_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT102485_chrX$tissue='pituitary'
pituitary_PT102485_chrX$id='PT102485'
pituitary_PT102485_chrX$chr='chrX'

pituitary_PT101159_chrX = read.table("compare_X_A/results/PT101159_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT101159_chrX$tissue='pituitary'
pituitary_PT101159_chrX$id='PT101159'
pituitary_PT101159_chrX$chr='chrX'

pituitary_PT105239_chrX = read.table("compare_X_A/results/PT105239_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT105239_chrX$tissue='pituitary'
pituitary_PT105239_chrX$id='PT105239'
pituitary_PT105239_chrX$chr='chrX'

pituitary_PT104883_chrX = read.table("compare_X_A/results/PT104883_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT104883_chrX$tissue='pituitary'
pituitary_PT104883_chrX$id='PT104883'
pituitary_PT104883_chrX$chr='chrX'

pituitary_PT103811_chrX = read.table("compare_X_A/results/PT103811_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT103811_chrX$tissue='pituitary'
pituitary_PT103811_chrX$id='PT103811'
pituitary_PT103811_chrX$chr='chrX'

pituitary_PT103352_chrX = read.table("compare_X_A/results/PT103352_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT103352_chrX$tissue='pituitary'
pituitary_PT103352_chrX$id='PT103352'
pituitary_PT103352_chrX$chr='chrX'

pituitary_PT102842_chrX = read.table("compare_X_A/results/PT102842_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT102842_chrX$tissue='pituitary'
pituitary_PT102842_chrX$id='PT102842'
pituitary_PT102842_chrX$chr='chrX'

pituitary_PT101210_chrX = read.table("compare_X_A/results/PT101210_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT101210_chrX$tissue='pituitary'
pituitary_PT101210_chrX$id='PT101210'
pituitary_PT101210_chrX$chr='chrX'

pituitary_PT103046_chrX = read.table("compare_X_A/results/PT103046_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT103046_chrX$tissue='pituitary'
pituitary_PT103046_chrX$id='PT103046'
pituitary_PT103046_chrX$chr='chrX'

pituitary_PT104733_chrX = read.table("compare_X_A/results/PT104733_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT104733_chrX$tissue='pituitary'
pituitary_PT104733_chrX$id='PT104733'
pituitary_PT104733_chrX$chr='chrX'

pituitary_PT103760_chrX = read.table("compare_X_A/results/PT103760_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT103760_chrX$tissue='pituitary'
pituitary_PT103760_chrX$id='PT103760'
pituitary_PT103760_chrX$chr='chrX'

pituitary_PT104270_chrX = read.table("compare_X_A/results/PT104270_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT104270_chrX$tissue='pituitary'
pituitary_PT104270_chrX$id='PT104270'
pituitary_PT104270_chrX$chr='chrX'

combine_data <- rbind(heart_PT104219_chr13, heart_PT105290_chr13, heart_PT102485_chr13, heart_PT101159_chr13, heart_PT105239_chr13, heart_PT104883_chr13, heart_PT103811_chr13, heart_PT103352_chr13, heart_PT102842_chr13, heart_PT101210_chr13, heart_PT103046_chr13, heart_PT104733_chr13, heart_PT103760_chr13, heart_PT104270_chr13, lung_PT104219_chr13, lung_PT105290_chr13, lung_PT102485_chr13, lung_PT101159_chr13, lung_PT105239_chr13, lung_PT104883_chr13, lung_PT103811_chr13, lung_PT103352_chr13, lung_PT102842_chr13, lung_PT101210_chr13, lung_PT103760_chr13, lung_PT104270_chr13, liver_PT104219_chr13, liver_PT105290_chr13, liver_PT102485_chr13, liver_PT101159_chr13, liver_PT105239_chr13, liver_PT104883_chr13, liver_PT103811_chr13, liver_PT103352_chr13, liver_PT102842_chr13, liver_PT101210_chr13, liver_PT103046_chr13, liver_PT104733_chr13, liver_PT103760_chr13, liver_PT104270_chr13, gonads_PT104219_chr13, gonads_PT105290_chr13, gonads_PT102485_chr13, gonads_PT101159_chr13, gonads_PT105239_chr13, gonads_PT104883_chr13, gonads_PT103811_chr13, gonads_PT103352_chr13, gonads_PT102842_chr13, gonads_PT101210_chr13, gonads_PT103046_chr13, gonads_PT104733_chr13, adrenal_PT104219_chr13, adrenal_PT105290_chr13, adrenal_PT102485_chr13, adrenal_PT101159_chr13, adrenal_PT105239_chr13, adrenal_PT104883_chr13, adrenal_PT103811_chr13, adrenal_PT103352_chr13, adrenal_PT102842_chr13, adrenal_PT101210_chr13, adrenal_PT103046_chr13, adrenal_PT104733_chr13, adrenal_PT103760_chr13, adrenal_PT104270_chr13, pituitary_PT104219_chr13, pituitary_PT105290_chr13, pituitary_PT102485_chr13, pituitary_PT101159_chr13, pituitary_PT105239_chr13, pituitary_PT104883_chr13, pituitary_PT103811_chr13, pituitary_PT103352_chr13, pituitary_PT102842_chr13, pituitary_PT101210_chr13, pituitary_PT103046_chr13, pituitary_PT104733_chr13, pituitary_PT103760_chr13, pituitary_PT104270_chr13, heart_PT104219_chrX, heart_PT105290_chrX, heart_PT102485_chrX, heart_PT101159_chrX, heart_PT105239_chrX, heart_PT104883_chrX, heart_PT103811_chrX, heart_PT103352_chrX, heart_PT102842_chrX, heart_PT101210_chrX, heart_PT103046_chrX, heart_PT104733_chrX, heart_PT103760_chrX, heart_PT104270_chrX, lung_PT104219_chrX, lung_PT105290_chrX, lung_PT102485_chrX, lung_PT101159_chrX, lung_PT105239_chrX, lung_PT104883_chrX, lung_PT103811_chrX, lung_PT103352_chrX, lung_PT102842_chrX, lung_PT101210_chrX, lung_PT103760_chrX, lung_PT104270_chrX, liver_PT104219_chrX, liver_PT105290_chrX, liver_PT102485_chrX, liver_PT101159_chrX, liver_PT105239_chrX, liver_PT104883_chrX, liver_PT103811_chrX, liver_PT103352_chrX, liver_PT102842_chrX, liver_PT101210_chrX, liver_PT103046_chrX, liver_PT104733_chrX, liver_PT103760_chrX, liver_PT104270_chrX, gonads_PT104219_chrX, gonads_PT105290_chrX, gonads_PT102485_chrX, gonads_PT101159_chrX, gonads_PT105239_chrX, gonads_PT104883_chrX, gonads_PT103811_chrX, gonads_PT103352_chrX, gonads_PT102842_chrX, gonads_PT101210_chrX, gonads_PT103046_chrX, gonads_PT104733_chrX, adrenal_PT104219_chrX, adrenal_PT105290_chrX, adrenal_PT102485_chrX, adrenal_PT101159_chrX, adrenal_PT105239_chrX, adrenal_PT104883_chrX, adrenal_PT103811_chrX, adrenal_PT103352_chrX, adrenal_PT102842_chrX, adrenal_PT101210_chrX, adrenal_PT103046_chrX, adrenal_PT104733_chrX, adrenal_PT103760_chrX, adrenal_PT104270_chrX, pituitary_PT104219_chrX, pituitary_PT105290_chrX, pituitary_PT102485_chrX, pituitary_PT101159_chrX, pituitary_PT105239_chrX, pituitary_PT104883_chrX, pituitary_PT103811_chrX, pituitary_PT103352_chrX, pituitary_PT102842_chrX, pituitary_PT101210_chrX, pituitary_PT103046_chrX, pituitary_PT104733_chrX, pituitary_PT103760_chrX, pituitary_PT104270_chrX)

png("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/plots/plot_violin_plot_allindidivuals_alltissues_chrX_comparedautosomes.png", width = 9, height=4, units = "in", res = 600)
ggplot(combine_data, aes(x=tissue, y=median, fill=chr)) +
  geom_violin() +
  coord_cartesian(ylim=c(0.5, 1)) +
  theme_bw() + 
  scale_fill_manual(values = c("black", "#40b0a6")) +
  labs(x="Tissue", y = "Allele balance", title = "Rhesus Macaque") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18), plot.title = element_text(hjust = 0.5, size = 22), legend.title = element_blank(), legend.position="top", legend.text = element_text(margin = margin(r = 30, unit = "pt"), size=16), legend.spacing.x = unit(0.2, 'cm'))

dev.off()
