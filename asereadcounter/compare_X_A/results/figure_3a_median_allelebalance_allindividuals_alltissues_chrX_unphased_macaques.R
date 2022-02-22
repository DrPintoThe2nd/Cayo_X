library(ggplot2)
setwd("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/results/")

PT101159_heart = read.table("PT101159_heart_mean_median_ratio_chrX.tsv", header = T)
PT101159_lung = read.table("PT101159_lung_mean_median_ratio_chrX.tsv", header = T)
PT101159_liver = read.table("PT101159_liver_mean_median_ratio_chrX.tsv", header = T)
PT101159_adrenal = read.table("PT101159_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT101159_pituitary = read.table("PT101159_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT101159_ovary = read.table("PT101159_gonads_mean_median_ratio_chrX.tsv", header = T)

PT101210_heart = read.table("PT101210_heart_mean_median_ratio_chrX.tsv", header = T)
PT101210_lung = read.table("PT101210_lung_mean_median_ratio_chrX.tsv", header = T)
PT101210_liver = read.table("PT101210_liver_mean_median_ratio_chrX.tsv", header = T)
PT101210_adrenal = read.table("PT101210_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT101210_pituitary = read.table("PT101210_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT101210_ovary = read.table("PT101210_gonads_mean_median_ratio_chrX.tsv", header = T)

PT102485_heart = read.table("PT102485_heart_mean_median_ratio_chrX.tsv", header = T)
PT102485_lung = read.table("PT102485_lung_mean_median_ratio_chrX.tsv", header = T)
PT102485_liver = read.table("PT102485_liver_mean_median_ratio_chrX.tsv", header = T)
PT102485_adrenal = read.table("PT102485_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT102485_pituitary = read.table("PT102485_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT102485_ovary = read.table("PT102485_gonads_mean_median_ratio_chrX.tsv", header = T)

PT102842_heart = read.table("PT102842_heart_mean_median_ratio_chrX.tsv", header = T)
PT102842_lung = read.table("PT102842_lung_mean_median_ratio_chrX.tsv", header = T)
PT102842_liver = read.table("PT102842_liver_mean_median_ratio_chrX.tsv", header = T)
PT102842_adrenal = read.table("PT102842_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT102842_pituitary = read.table("PT102842_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT102842_ovary = read.table("PT102842_gonads_mean_median_ratio_chrX.tsv", header = T)

PT103352_heart = read.table("PT103352_heart_mean_median_ratio_chrX.tsv", header = T)
PT103352_lung = read.table("PT103352_lung_mean_median_ratio_chrX.tsv", header = T)
PT103352_liver = read.table("PT103352_liver_mean_median_ratio_chrX.tsv", header = T)
PT103352_adrenal = read.table("PT103352_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT103352_pituitary = read.table("PT103352_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT103352_ovary = read.table("PT103352_gonads_mean_median_ratio_chrX.tsv", header = T)

PT103811_heart = read.table("PT103811_heart_mean_median_ratio_chrX.tsv", header = T)
PT103811_lung = read.table("PT103811_lung_mean_median_ratio_chrX.tsv", header = T)
PT103811_liver = read.table("PT103811_liver_mean_median_ratio_chrX.tsv", header = T)
PT103811_adrenal = read.table("PT103811_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT103811_pituitary = read.table("PT103811_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT103811_ovary = read.table("PT103811_gonads_mean_median_ratio_chrX.tsv", header = T)

PT104219_heart = read.table("PT104219_heart_mean_median_ratio_chrX.tsv", header = T)
PT104219_lung = read.table("PT104219_lung_mean_median_ratio_chrX.tsv", header = T)
PT104219_liver = read.table("PT104219_liver_mean_median_ratio_chrX.tsv", header = T)
PT104219_adrenal = read.table("PT104219_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT104219_pituitary = read.table("PT104219_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT104219_ovary = read.table("PT104219_gonads_mean_median_ratio_chrX.tsv", header = T)

PT104883_heart = read.table("PT104883_heart_mean_median_ratio_chrX.tsv", header = T)
PT104883_lung = read.table("PT104883_lung_mean_median_ratio_chrX.tsv", header = T)
PT104883_liver = read.table("PT104883_liver_mean_median_ratio_chrX.tsv", header = T)
PT104883_adrenal = read.table("PT104883_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT104883_pituitary = read.table("PT104883_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT104883_ovary = read.table("PT104883_gonads_mean_median_ratio_chrX.tsv", header = T)

PT105239_heart = read.table("PT105239_heart_mean_median_ratio_chrX.tsv", header = T)
PT105239_lung = read.table("PT105239_lung_mean_median_ratio_chrX.tsv", header = T)
PT105239_liver = read.table("PT105239_liver_mean_median_ratio_chrX.tsv", header = T)
PT105239_adrenal = read.table("PT105239_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT105239_pituitary = read.table("PT105239_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT105239_ovary = read.table("PT105239_gonads_mean_median_ratio_chrX.tsv", header = T)

PT105290_heart = read.table("PT105290_heart_mean_median_ratio_chrX.tsv", header = T)
PT105290_lung = read.table("PT105290_lung_mean_median_ratio_chrX.tsv", header = T)
PT105290_liver = read.table("PT105290_liver_mean_median_ratio_chrX.tsv", header = T)
PT105290_adrenal = read.table("PT105290_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT105290_pituitary = read.table("PT105290_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT105290_ovary = read.table("PT105290_gonads_mean_median_ratio_chrX.tsv", header = T)

PT103046_heart = read.table("PT103046_heart_mean_median_ratio_chrX.tsv", header = T)
PT103046_liver = read.table("PT103046_liver_mean_median_ratio_chrX.tsv", header = T)
PT103046_adrenal = read.table("PT103046_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT103046_pituitary = read.table("PT103046_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT103046_ovary = read.table("PT103046_gonads_mean_median_ratio_chrX.tsv", header = T)

PT104733_heart = read.table("PT104733_heart_mean_median_ratio_chrX.tsv", header = T)
PT104733_liver = read.table("PT104733_liver_mean_median_ratio_chrX.tsv", header = T)
PT104733_adrenal = read.table("PT104733_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT104733_pituitary = read.table("PT104733_pituitary_mean_median_ratio_chrX.tsv", header = T)
PT104733_ovary = read.table("PT104733_gonads_mean_median_ratio_chrX.tsv", header = T)

PT103760_heart = read.table("PT103760_heart_mean_median_ratio_chrX.tsv", header = T)
PT103760_lung = read.table("PT103760_lung_mean_median_ratio_chrX.tsv", header = T)
PT103760_liver = read.table("PT103760_liver_mean_median_ratio_chrX.tsv", header = T)
PT103760_adrenal = read.table("PT103760_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT103760_pituitary = read.table("PT103760_pituitary_mean_median_ratio_chrX.tsv", header = T)

PT104270_heart = read.table("PT104270_heart_mean_median_ratio_chrX.tsv", header = T)
PT104270_lung = read.table("PT104270_lung_mean_median_ratio_chrX.tsv", header = T)
PT104270_liver = read.table("PT104270_liver_mean_median_ratio_chrX.tsv", header = T)
PT104270_adrenal = read.table("PT104270_adrenal_mean_median_ratio_chrX.tsv", header = T)
PT104270_pituitary = read.table("PT104270_pituitary_mean_median_ratio_chrX.tsv", header = T)

df_plot = data.frame(vals = c(PT101159_heart$median, PT101159_lung$median, PT101159_liver$median, PT101159_adrenal$median, PT101159_pituitary$median, PT101159_ovary$median, PT101210_heart$median, PT101210_lung$median, PT101210_liver$median, PT101210_adrenal$median, PT101210_pituitary$median, PT101210_ovary$median, PT102485_heart$median, PT102485_lung$median, PT102485_liver$median, PT102485_adrenal$median, PT102485_pituitary$median, PT102485_ovary$median, PT102842_heart$median, PT102842_lung$median, PT102842_liver$median, PT102842_adrenal$median, PT102842_pituitary$median, PT102842_ovary$median, PT103352_heart$median, PT103352_lung$median, PT103352_liver$median, PT103352_adrenal$median, PT103352_pituitary$median, PT103352_ovary$median, PT103811_heart$median, PT103811_lung$median, PT103811_liver$median, PT103811_adrenal$median, PT103811_pituitary$median, PT103811_ovary$median, PT104219_heart$median, PT104219_lung$median, PT104219_liver$median, PT104219_adrenal$median, PT104219_pituitary$median, PT104219_ovary$median, PT104883_heart$median, PT104883_lung$median, PT104883_liver$median, PT104883_adrenal$median, PT104883_pituitary$median, PT104883_ovary$median, PT105239_heart$median, PT105239_lung$median, PT105239_liver$median, PT105239_adrenal$median, PT105239_pituitary$median, PT105239_ovary$median, PT105290_heart$median, PT105290_lung$median, PT105290_liver$median, PT105290_adrenal$median, PT105290_pituitary$median, PT105290_ovary$median, PT103046_heart$median, NA, PT103046_liver$median, PT103046_adrenal$median, PT103046_pituitary$median, PT103046_ovary$median, PT104733_heart$median, NA, PT104733_liver$median, PT104733_adrenal$median, PT104733_pituitary$median, PT104733_ovary$median, PT103760_heart$median, PT103760_lung$median, PT103760_liver$median, PT103760_adrenal$median, PT103760_pituitary$median, NA, PT104270_heart$median, PT104270_lung$median, PT104270_liver$median, PT104270_adrenal$median, PT104270_pituitary$median, NA), Sites = c("Heart", "Lung", "Liver", "Adrenal", "Pituitary", "Ovary"), samples = c("PT101159", "PT101159", "PT101159", "PT101159", "PT101159", "PT101159", "PT101210", "PT101210", "PT101210", "PT101210", "PT101210", "PT101210", "PT102485", "PT102485", "PT102485", "PT102485", "PT102485", "PT102485", "PT102842", "PT102842", "PT102842", "PT102842", "PT102842", "PT102842", "PT103046", "PT103046", "PT103046", "PT103046", "PT103046", "PT103046", "PT103352", "PT103352", "PT103352", "PT103352", "PT103352", "PT103352", "PT103760", "PT103760", "PT103760", "PT103760", "PT103760", "PT103760", "PT103811", "PT103811", "PT103811", "PT103811", "PT103811", "PT103811", "PT104219", "PT104219", "PT104219", "PT104219", "PT104219", "PT104219", "PT104270", "PT104270", "PT104270", "PT104270", "PT104270", "PT104270", "PT104733", "PT104733", "PT104733", "PT104733", "PT104733", "PT104733", "PT104883", "PT104883", "PT104883", "PT104883", "PT104883", "PT104883", "PT105239", "PT105239", "PT105239", "PT105239", "PT105239", "PT105239", "PT105290", "PT105290", "PT105290", "PT105290", "PT105290", "PT105290"))

df_plot$Sites <- as.character(df_plot$Sites)
# df_plot$Sites <- factor(df_plot$Sites, levels=unique(df_plot$Sites))
df_plot$Sites <- factor(df_plot$Sites, levels=c("Heart", "Ovary", "Lung", "Liver", "Adrenal", "Pituitary"))

df_plot$samples <- as.character(df_plot$samples)
df_plot$samples <- factor(df_plot$samples, levels=c("PT103811", "PT104733", "PT101159", "PT101210", "PT103352", "PT104270", "PT105290", "PT103046", "PT104219", "PT105239", "PT102485", "PT103760", "PT104883", "PT102842"))

png("/data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/compare_X_A/plots/figure_3a_median_allelebalance_allindividuals_alltissues_chrX_unphased_macaques.png", width = 13, height = 6, units = "in", res = 300)
ggplot(df_plot, aes(x=samples, y=vals, color=Sites, shape=Sites, size=Sites, fill=Sites)) +
  geom_point() +
  theme_bw() +
  coord_cartesian(ylim=c(0.5, 1)) +
  geom_hline(yintercept = 0.8, color="darkgray", size=0.75, linetype=2) +
  scale_colour_manual(values=c("black", "black", "darkgray", "darkgray", "black", "black")) +
  scale_shape_manual(values=c(21, 24, 21, 24, 24, 21)) +
  scale_size_manual(values=c(5, 5, 5, 5, 5, 5)) +
  scale_fill_manual(values=c(NA, NA, "darkgray", "darkgray", "black", "black")) +
  theme(axis.text.x = element_text(size=14, colour="black"), axis.text.y = element_text(size=14, colour = "black"), axis.title= element_text(size=18), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=14), plot.title = element_text(hjust = 0.5, size = 22)) +
  labs(x="Samples", y="Allele balance", title = "Macaque")

dev.off()