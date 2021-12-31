library(ggplot2)

setwd("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/")

# heart
heart_PT101159 = read.table("compare_X_A/results/PT101159_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT101159$tissue='heart'
heart_PT101159$id='PT101159'

heart_PT105290 = read.table("compare_X_A/results/PT105290_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT105290$tissue='heart'
heart_PT105290$id='PT105290'

heart_PT102485 = read.table("compare_X_A/results/PT102485_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT102485$tissue='heart'
heart_PT102485$id='PT102485'

heart_PT101159 = read.table("compare_X_A/results/PT101159_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT101159$tissue='heart'
heart_PT101159$id='PT101159'

heart_PT105239 = read.table("compare_X_A/results/PT105239_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT105239$tissue='heart'
heart_PT105239$id='PT105239'

heart_PT104883 = read.table("compare_X_A/results/PT104883_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT104883$tissue='heart'
heart_PT104883$id='PT104883'

heart_PT103811 = read.table("compare_X_A/results/PT103811_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT103811$tissue='heart'
heart_PT103811$id='PT103811'

heart_PT103352 = read.table("compare_X_A/results/PT103352_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT103352$tissue='heart'
heart_PT103352$id='PT103352'

heart_PT102842 = read.table("compare_X_A/results/PT102842_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT102842$tissue='heart'
heart_PT102842$id='PT102842'

heart_PT101210 = read.table("compare_X_A/results/PT101210_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT101210$tissue='heart'
heart_PT101210$id='PT101210'

heart_PT103046 = read.table("compare_X_A/results/PT103046_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT103046$tissue='heart'
heart_PT103046$id='PT103046'

heart_PT104733 = read.table("compare_X_A/results/PT104733_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT104733$tissue='heart'
heart_PT104733$id='PT104733'

heart_PT103760 = read.table("compare_X_A/results/PT103760_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT103760$tissue='heart'
heart_PT103760$id='PT103760'

heart_PT104270 = read.table("compare_X_A/results/PT104270_heart_mean_median_ratio_chrX.tsv", header=T)
heart_PT104270$tissue='heart'
heart_PT104270$id='PT104270'

# lung
lung_PT101159 = read.table("compare_X_A/results/PT101159_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT101159$tissue='lung'
lung_PT101159$id='PT101159'

lung_PT105290 = read.table("compare_X_A/results/PT105290_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT105290$tissue='lung'
lung_PT105290$id='PT105290'

lung_PT102485 = read.table("compare_X_A/results/PT102485_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT102485$tissue='lung'
lung_PT102485$id='PT102485'

lung_PT101159 = read.table("compare_X_A/results/PT101159_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT101159$tissue='lung'
lung_PT101159$id='PT101159'

lung_PT105239 = read.table("compare_X_A/results/PT105239_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT105239$tissue='lung'
lung_PT105239$id='PT105239'

lung_PT104883 = read.table("compare_X_A/results/PT104883_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT104883$tissue='lung'
lung_PT104883$id='PT104883'

lung_PT103811 = read.table("compare_X_A/results/PT103811_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT103811$tissue='lung'
lung_PT103811$id='PT103811'

lung_PT103352 = read.table("compare_X_A/results/PT103352_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT103352$tissue='lung'
lung_PT103352$id='PT103352'

lung_PT102842 = read.table("compare_X_A/results/PT102842_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT102842$tissue='lung'
lung_PT102842$id='PT102842'

lung_PT101210 = read.table("compare_X_A/results/PT101210_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT101210$tissue='lung'
lung_PT101210$id='PT101210'

lung_PT103760 = read.table("compare_X_A/results/PT103760_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT103760$tissue='lung'
lung_PT103760$id='PT103760'

lung_PT104270 = read.table("compare_X_A/results/PT104270_lung_mean_median_ratio_chrX.tsv", header=T)
lung_PT104270$tissue='lung'
lung_PT104270$id='PT104270'

# liver
liver_PT101159 = read.table("compare_X_A/results/PT101159_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT101159$tissue='liver'
liver_PT101159$id='PT101159'

liver_PT105290 = read.table("compare_X_A/results/PT105290_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT105290$tissue='liver'
liver_PT105290$id='PT105290'

liver_PT102485 = read.table("compare_X_A/results/PT102485_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT102485$tissue='liver'
liver_PT102485$id='PT102485'

liver_PT101159 = read.table("compare_X_A/results/PT101159_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT101159$tissue='liver'
liver_PT101159$id='PT101159'

liver_PT105239 = read.table("compare_X_A/results/PT105239_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT105239$tissue='liver'
liver_PT105239$id='PT105239'

liver_PT104883 = read.table("compare_X_A/results/PT104883_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT104883$tissue='liver'
liver_PT104883$id='PT104883'

liver_PT103811 = read.table("compare_X_A/results/PT103811_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT103811$tissue='liver'
liver_PT103811$id='PT103811'

liver_PT103352 = read.table("compare_X_A/results/PT103352_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT103352$tissue='liver'
liver_PT103352$id='PT103352'

liver_PT102842 = read.table("compare_X_A/results/PT102842_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT102842$tissue='liver'
liver_PT102842$id='PT102842'

liver_PT101210 = read.table("compare_X_A/results/PT101210_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT101210$tissue='liver'
liver_PT101210$id='PT101210'

liver_PT103046 = read.table("compare_X_A/results/PT103046_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT103046$tissue='liver'
liver_PT103046$id='PT103046'

liver_PT104733 = read.table("compare_X_A/results/PT104733_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT104733$tissue='liver'
liver_PT104733$id='PT104733'

liver_PT103760 = read.table("compare_X_A/results/PT103760_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT103760$tissue='liver'
liver_PT103760$id='PT103760'

liver_PT104270 = read.table("compare_X_A/results/PT104270_liver_mean_median_ratio_chrX.tsv", header=T)
liver_PT104270$tissue='liver'
liver_PT104270$id='PT104270'

# gonads
gonads_PT101159 = read.table("compare_X_A/results/PT101159_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT101159$tissue='gonads'
gonads_PT101159$id='PT101159'

gonads_PT105290 = read.table("compare_X_A/results/PT105290_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT105290$tissue='gonads'
gonads_PT105290$id='PT105290'

gonads_PT102485 = read.table("compare_X_A/results/PT102485_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT102485$tissue='gonads'
gonads_PT102485$id='PT102485'

gonads_PT101159 = read.table("compare_X_A/results/PT101159_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT101159$tissue='gonads'
gonads_PT101159$id='PT101159'

gonads_PT105239 = read.table("compare_X_A/results/PT105239_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT105239$tissue='gonads'
gonads_PT105239$id='PT105239'

gonads_PT104883 = read.table("compare_X_A/results/PT104883_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT104883$tissue='gonads'
gonads_PT104883$id='PT104883'

gonads_PT103811 = read.table("compare_X_A/results/PT103811_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT103811$tissue='gonads'
gonads_PT103811$id='PT103811'

gonads_PT103352 = read.table("compare_X_A/results/PT103352_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT103352$tissue='gonads'
gonads_PT103352$id='PT103352'

gonads_PT102842 = read.table("compare_X_A/results/PT102842_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT102842$tissue='gonads'
gonads_PT102842$id='PT102842'

gonads_PT101210 = read.table("compare_X_A/results/PT101210_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT101210$tissue='gonads'
gonads_PT101210$id='PT101210'

gonads_PT103046 = read.table("compare_X_A/results/PT103046_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT103046$tissue='gonads'
gonads_PT103046$id='PT103046'

gonads_PT104733 = read.table("compare_X_A/results/PT104733_gonads_mean_median_ratio_chrX.tsv", header=T)
gonads_PT104733$tissue='gonads'
gonads_PT104733$id='PT104733'

# adrenal
adrenal_PT101159 = read.table("compare_X_A/results/PT101159_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT101159$tissue='adrenal'
adrenal_PT101159$id='PT101159'

adrenal_PT105290 = read.table("compare_X_A/results/PT105290_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT105290$tissue='adrenal'
adrenal_PT105290$id='PT105290'

adrenal_PT102485 = read.table("compare_X_A/results/PT102485_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT102485$tissue='adrenal'
adrenal_PT102485$id='PT102485'

adrenal_PT101159 = read.table("compare_X_A/results/PT101159_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT101159$tissue='adrenal'
adrenal_PT101159$id='PT101159'

adrenal_PT105239 = read.table("compare_X_A/results/PT105239_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT105239$tissue='adrenal'
adrenal_PT105239$id='PT105239'

adrenal_PT104883 = read.table("compare_X_A/results/PT104883_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT104883$tissue='adrenal'
adrenal_PT104883$id='PT104883'

adrenal_PT103811 = read.table("compare_X_A/results/PT103811_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT103811$tissue='adrenal'
adrenal_PT103811$id='PT103811'

adrenal_PT103352 = read.table("compare_X_A/results/PT103352_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT103352$tissue='adrenal'
adrenal_PT103352$id='PT103352'

adrenal_PT102842 = read.table("compare_X_A/results/PT102842_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT102842$tissue='adrenal'
adrenal_PT102842$id='PT102842'

adrenal_PT101210 = read.table("compare_X_A/results/PT101210_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT101210$tissue='adrenal'
adrenal_PT101210$id='PT101210'

adrenal_PT103046 = read.table("compare_X_A/results/PT103046_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT103046$tissue='adrenal'
adrenal_PT103046$id='PT103046'

adrenal_PT104733 = read.table("compare_X_A/results/PT104733_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT104733$tissue='adrenal'
adrenal_PT104733$id='PT104733'

adrenal_PT103760 = read.table("compare_X_A/results/PT103760_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT103760$tissue='adrenal'
adrenal_PT103760$id='PT103760'

adrenal_PT104270 = read.table("compare_X_A/results/PT104270_adrenal_mean_median_ratio_chrX.tsv", header=T)
adrenal_PT104270$tissue='adrenal'
adrenal_PT104270$id='PT104270'

# pituitary
pituitary_PT101159 = read.table("compare_X_A/results/PT101159_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT101159$tissue='pituitary'
pituitary_PT101159$id='PT101159'

pituitary_PT105290 = read.table("compare_X_A/results/PT105290_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT105290$tissue='pituitary'
pituitary_PT105290$id='PT105290'

pituitary_PT102485 = read.table("compare_X_A/results/PT102485_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT102485$tissue='pituitary'
pituitary_PT102485$id='PT102485'

pituitary_PT101159 = read.table("compare_X_A/results/PT101159_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT101159$tissue='pituitary'
pituitary_PT101159$id='PT101159'

pituitary_PT105239 = read.table("compare_X_A/results/PT105239_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT105239$tissue='pituitary'
pituitary_PT105239$id='PT105239'

pituitary_PT104883 = read.table("compare_X_A/results/PT104883_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT104883$tissue='pituitary'
pituitary_PT104883$id='PT104883'

pituitary_PT103811 = read.table("compare_X_A/results/PT103811_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT103811$tissue='pituitary'
pituitary_PT103811$id='PT103811'

pituitary_PT103352 = read.table("compare_X_A/results/PT103352_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT103352$tissue='pituitary'
pituitary_PT103352$id='PT103352'

pituitary_PT102842 = read.table("compare_X_A/results/PT102842_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT102842$tissue='pituitary'
pituitary_PT102842$id='PT102842'

pituitary_PT101210 = read.table("compare_X_A/results/PT101210_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT101210$tissue='pituitary'
pituitary_PT101210$id='PT101210'

pituitary_PT103046 = read.table("compare_X_A/results/PT103046_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT103046$tissue='pituitary'
pituitary_PT103046$id='PT103046'

pituitary_PT104733 = read.table("compare_X_A/results/PT104733_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT104733$tissue='pituitary'
pituitary_PT104733$id='PT104733'

pituitary_PT103760 = read.table("compare_X_A/results/PT103760_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT103760$tissue='pituitary'
pituitary_PT103760$id='PT103760'

pituitary_PT104270 = read.table("compare_X_A/results/PT104270_pituitary_mean_median_ratio_chrX.tsv", header=T)
pituitary_PT104270$tissue='pituitary'
pituitary_PT104270$id='PT104270'


combine_data = rbind(heart_PT101159, heart_PT105290, heart_PT102485, heart_PT101159, heart_PT105239, heart_PT104883, heart_PT103811, heart_PT103352, heart_PT102842, heart_PT101210, heart_PT103046, heart_PT104733, heart_PT103760, heart_PT104270, lung_PT101159, lung_PT105290, lung_PT102485, lung_PT101159, lung_PT105239, lung_PT104883, lung_PT103811, lung_PT103352, lung_PT102842, lung_PT101210, lung_PT103760, lung_PT104270, liver_PT101159, liver_PT105290, liver_PT102485, liver_PT101159, liver_PT105239, liver_PT104883, liver_PT103811, liver_PT103352, liver_PT102842, liver_PT101210, liver_PT103046, liver_PT104733, liver_PT103760, liver_PT104270, gonads_PT101159, gonads_PT105290, gonads_PT102485, gonads_PT101159, gonads_PT105239, gonads_PT104883, gonads_PT103811, gonads_PT103352, gonads_PT102842, gonads_PT101210, gonads_PT103046, gonads_PT104733, adrenal_PT101159, adrenal_PT105290, adrenal_PT102485, adrenal_PT101159, adrenal_PT105239, adrenal_PT104883, adrenal_PT103811, adrenal_PT103352, adrenal_PT102842, adrenal_PT101210, adrenal_PT103046, adrenal_PT104733, adrenal_PT103760, adrenal_PT104270, pituitary_PT101159, pituitary_PT105290, pituitary_PT102485, pituitary_PT101159, pituitary_PT105239, pituitary_PT104883, pituitary_PT103811, pituitary_PT103352, pituitary_PT102842, pituitary_PT101210, pituitary_PT103046, pituitary_PT104733, pituitary_PT103760, pituitary_PT104270)

png("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/plots/plot_median_chrX_allindividuals_alltissues.png", width = 8, height = 5, units = "in", res = 600)
ggplot(combine_data, aes(x=tissue, y=median)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_bw() +
  coord_cartesian(ylim=c(0.5,1)) +
  labs(x='Tissue', y="Allele balance", title="X chromosome") +
  theme(plot.title = element_text(hjust = 0.5, size=18), legend.position = "none", axis.text = element_text(size=16), axis.title = element_text(size=18))
dev.off()
