library(ggplot2)
setwd("/scratch/bpinto2/Cayo_GTEx/asereadcounter/compare_X_A/results/")

GTEX_131YS_0003_heart = read.table("GTEX_131YS_heart_mean_median_ratio_chrX.tsv", header = T)
GTEX_131YS_0003_lung = read.table("GTEX_131YS_lung_mean_median_ratio_chrX.tsv", header = T)
GTEX_131YS_0003_liver = read.table("GTEX_131YS_liver_mean_median_ratio_chrX.tsv", header = T)
GTEX_131YS_0003_adrenal = read.table("GTEX_131YS_adrenal_mean_median_ratio_chrX.tsv", header = T)
GTEX_131YS_0003_pituitary = read.table("GTEX_131YS_pituitary_mean_median_ratio_chrX.tsv", header = T)
GTEX_131YS_0003_ovary = read.table("GTEX_131YS_heart_mean_median_ratio_chrX.tsv", header = T)

GTEX_1A32A_0001_heart = read.table("GTEX_1A32A_heart_mean_median_ratio_chrX.tsv", header = T)
GTEX_1A32A_0001_lung = read.table("GTEX_1A32A_lung_mean_median_ratio_chrX.tsv", header = T)
GTEX_1A32A_0001_liver = read.table("GTEX_1A32A_liver_mean_median_ratio_chrX.tsv", header = T)
GTEX_1A32A_0001_adrenal = read.table("GTEX_1A32A_adrenal_mean_median_ratio_chrX.tsv", header = T)
GTEX_1A32A_0001_pituitary = read.table("GTEX_1A32A_pituitary_mean_median_ratio_chrX.tsv", header = T)
GTEX_1A32A_0001_ovary = read.table("GTEX_1A32A_heart_mean_median_ratio_chrX.tsv", header = T)

GTEX_14DAQ_0004_heart = read.table("GTEX_14DAQ_heart_mean_median_ratio_chrX.tsv", header = T)
GTEX_14DAQ_0004_lung = read.table("GTEX_14DAQ_lung_mean_median_ratio_chrX.tsv", header = T)
GTEX_14DAQ_0004_liver = read.table("GTEX_14DAQ_liver_mean_median_ratio_chrX.tsv", header = T)
GTEX_14DAQ_0004_adrenal = read.table("GTEX_14DAQ_adrenal_mean_median_ratio_chrX.tsv", header = T)
GTEX_14DAQ_0004_pituitary = read.table("GTEX_14DAQ_pituitary_mean_median_ratio_chrX.tsv", header = T)
GTEX_14DAQ_0004_ovary = read.table("GTEX_14DAQ_heart_mean_median_ratio_chrX.tsv", header = T)

GTEX_12WSD_0002_heart = read.table("GTEX_12WSD_heart_mean_median_ratio_chrX.tsv", header = T)
GTEX_12WSD_0002_lung = read.table("GTEX_12WSD_lung_mean_median_ratio_chrX.tsv", header = T)
GTEX_12WSD_0002_liver = read.table("GTEX_12WSD_liver_mean_median_ratio_chrX.tsv", header = T)
GTEX_12WSD_0002_adrenal = read.table("GTEX_12WSD_adrenal_mean_median_ratio_chrX.tsv", header = T)
GTEX_12WSD_0002_pituitary = read.table("GTEX_12WSD_pituitary_mean_median_ratio_chrX.tsv", header = T)
GTEX_12WSD_0002_ovary = read.table("GTEX_12WSD_heart_mean_median_ratio_chrX.tsv", header = T)

#GTEX_1A8FM_0004_heart = read.table("GTEX_1A8FM_heart_mean_median_ratio_chrX.tsv", header = T)
#GTEX_1A8FM_0004_lung = read.table("GTEX_1A8FM_lung_mean_median_ratio_chrX.tsv", header = T)
#GTEX_1A8FM_0004_liver = read.table("GTEX_1A8FM_liver_mean_median_ratio_chrX.tsv", header = T)
#GTEX_1A8FM_0004_adrenal = read.table("GTEX_1A8FM_adrenal_mean_median_ratio_chrX.tsv", header = T)
#GTEX_1A8FM_0004_pituitary = read.table("GTEX_1A8FM_pituitary_mean_median_ratio_chrX.tsv", header = T)
#GTEX_1A8FM_0004_ovary = read.table("GTEX_1A8FM_heart_mean_median_ratio_chrX.tsv", header = T)

GTEX_15ER7_0002_heart = read.table("GTEX_15ER7_heart_mean_median_ratio_chrX.tsv", header = T)
GTEX_15ER7_0002_lung = read.table("GTEX_15ER7_lung_mean_median_ratio_chrX.tsv", header = T)
GTEX_15ER7_0002_liver = read.table("GTEX_15ER7_liver_mean_median_ratio_chrX.tsv", header = T)
GTEX_15ER7_0002_adrenal = read.table("GTEX_15ER7_adrenal_mean_median_ratio_chrX.tsv", header = T)
GTEX_15ER7_0002_pituitary = read.table("GTEX_15ER7_pituitary_mean_median_ratio_chrX.tsv", header = T)
GTEX_15ER7_0002_ovary = read.table("GTEX_15ER7_heart_mean_median_ratio_chrX.tsv", header = T)

GTEX_11GSP_0004_heart = read.table("GTEX_11GSP_heart_mean_median_ratio_chrX.tsv", header = T)
GTEX_11GSP_0004_lung = read.table("GTEX_11GSP_lung_mean_median_ratio_chrX.tsv", header = T)
GTEX_11GSP_0004_liver = read.table("GTEX_11GSP_liver_mean_median_ratio_chrX.tsv", header = T)
GTEX_11GSP_0004_adrenal = read.table("GTEX_11GSP_adrenal_mean_median_ratio_chrX.tsv", header = T)
GTEX_11GSP_0004_pituitary = read.table("GTEX_11GSP_pituitary_mean_median_ratio_chrX.tsv", header = T)
GTEX_11GSP_0004_ovary = read.table("GTEX_11GSP_heart_mean_median_ratio_chrX.tsv", header = T)

df_plot = data.frame(vals = c(GTEX_131YS_0003_heart$median, GTEX_131YS_0003_lung$median, GTEX_131YS_0003_liver$median, GTEX_131YS_0003_adrenal$median, GTEX_131YS_0003_pituitary$median, GTEX_131YS_0003_ovary$median, GTEX_1A32A_0001_heart$median, GTEX_1A32A_0001_lung$median, GTEX_1A32A_0001_liver$median, GTEX_1A32A_0001_adrenal$median, GTEX_1A32A_0001_pituitary$median, GTEX_1A32A_0001_ovary$median, GTEX_14DAQ_0004_heart$median, GTEX_14DAQ_0004_lung$median, GTEX_14DAQ_0004_liver$median, GTEX_14DAQ_0004_adrenal$median, GTEX_14DAQ_0004_pituitary$median, GTEX_14DAQ_0004_ovary$median, GTEX_12WSD_0002_heart$median, GTEX_12WSD_0002_lung$median, GTEX_12WSD_0002_liver$median, GTEX_12WSD_0002_adrenal$median, GTEX_12WSD_0002_pituitary$median, GTEX_12WSD_0002_ovary$median, GTEX_15ER7_0002_heart$median, GTEX_15ER7_0002_lung$median, GTEX_15ER7_0002_liver$median, GTEX_15ER7_0002_adrenal$median, GTEX_15ER7_0002_pituitary$median, GTEX_15ER7_0002_ovary$median, GTEX_11GSP_0004_heart$median, GTEX_11GSP_0004_lung$median, GTEX_11GSP_0004_liver$median, GTEX_11GSP_0004_adrenal$median, GTEX_11GSP_0004_pituitary$median, GTEX_11GSP_0004_ovary$median), Sites = c("Heart", "Lung", "Liver", "Adrenal", "Pituitary", "Ovary"), samples = c("GTEX_131YS", "GTEX_131YS", "GTEX_131YS", "GTEX_131YS", "GTEX_131YS", "GTEX_131YS", "GTEX_1A32A", "GTEX_1A32A", "GTEX_1A32A", "GTEX_1A32A", "GTEX_1A32A", "GTEX_1A32A", "GTEX_14DAQ", "GTEX_14DAQ", "GTEX_14DAQ", "GTEX_14DAQ", "GTEX_14DAQ", "GTEX_14DAQ", "GTEX_12WSD", "GTEX_12WSD", "GTEX_12WSD", "GTEX_12WSD", "GTEX_12WSD", "GTEX_12WSD", "GTEX_15ER7", "GTEX_15ER7", "GTEX_15ER7", "GTEX_15ER7", "GTEX_15ER7", "GTEX_15ER7", "GTEX_11GSP", "GTEX_11GSP", "GTEX_11GSP", "GTEX_11GSP", "GTEX_11GSP", "GTEX_11GSP"))

df_plot$Sites <- as.character(df_plot$Sites)
df_plot$Sites <- factor(df_plot$Sites, levels=c("Heart", "Ovary", "Lung", "Liver", "Adrenal", "Pituitary"))

df_plot$samples <- as.character(df_plot$samples)

png("/scratch/bpinto2/Cayo_GTEx/asereadcounter/compare_X_A/plots/test_01272020_GTEx.png", width = 13, height = 6, units = "in", res = 300)
ggplot(df_plot, aes(x=samples, y=vals, color=Sites, shape=Sites, size=Sites, fill=Sites)) +
  geom_point() +
  theme_bw() +
  coord_cartesian(ylim=c(0.5, 1)) +
  geom_hline(yintercept = 0.8, color="darkgray", size=0.75, linetype=2) +
  scale_colour_manual(values=c("black", "black", "darkgray", "darkgray", "black", "black")) +
  scale_shape_manual(values=c(21, 24, 21, 24, 24, 21)) +
  scale_size_manual(values=c(5, 5, 5, 5, 5, 5)) +
  scale_fill_manual(values=c(NA, NA, "darkgray", "darkgray", "black", "black")) +
  theme(axis.text.x = element_text(size=14, colour="black", angle = 30, hjust = 1), axis.text.y = element_text(size=14, colour = "black"), axis.title= element_text(size=18), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=14), plot.title = element_text(hjust = 0.5, size = 22)) +
  labs(x="Samples", y="Allele balance", title = "Human")

dev.off()