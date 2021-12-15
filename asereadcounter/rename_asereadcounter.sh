#!/bin/bash
# For convenience, rename rna id to tissue name

cd /data/CEM/wilsonlab/lab_generated/cayo/asereadcounter/
# ----
# chrX
# ----
cp PT101159_PT101114.tsv PT101159_heart.tsv
cp PT101159_PT101115.tsv PT101159_lung.tsv
cp PT101159_PT101118.tsv PT101159_liver.tsv
cp PT101159_PT101120.tsv PT101159_adrenal.tsv
cp PT101159_PT101121.tsv PT101159_pituitary.tsv
cp PT101159_PT101124.tsv PT101159_gonads.tsv

cp PT101210_PT101165.tsv PT101210_heart.tsv
cp PT101210_PT101166.tsv PT101210_lung.tsv
cp PT101210_PT101169.tsv PT101210_liver.tsv
cp PT101210_PT101171.tsv PT101210_adrenal.tsv
cp PT101210_PT101172.tsv PT101210_pituitary.tsv
cp PT101210_PT101175.tsv PT101210_gonads.tsv

cp PT102485_PT102440.tsv PT102485_heart.tsv
cp PT102485_PT102441.tsv PT102485_lung.tsv
cp PT102485_PT102444.tsv PT102485_liver.tsv
cp PT102485_PT102446.tsv PT102485_adrenal.tsv
cp PT102485_PT102447.tsv PT102485_pituitary.tsv
cp PT102485_PT102450.tsv PT102485_gonads.tsv

cp PT102842_PT102797.tsv PT102842_heart.tsv
cp PT102842_PT102798.tsv PT102842_lung.tsv
cp PT102842_PT102801.tsv PT102842_liver.tsv
cp PT102842_PT102803.tsv PT102842_adrenal.tsv
cp PT102842_PT102804.tsv PT102842_pituitary.tsv
cp PT102842_PT102807.tsv PT102842_gonads.tsv

cp PT103352_PT103307.tsv PT103352_heart.tsv
cp PT103352_PT103308.tsv PT103352_lung.tsv
cp PT103352_PT103311.tsv PT103352_liver.tsv
cp PT103352_PT103313.tsv PT103352_adrenal.tsv
cp PT103352_PT103314.tsv PT103352_pituitary.tsv
cp PT103352_PT103317.tsv PT103352_gonads.tsv

cp PT103811_PT103766.tsv PT103811_heart.tsv
cp PT103811_PT103767.tsv PT103811_lung.tsv
cp PT103811_PT103770.tsv PT103811_liver.tsv
cp PT103811_PT103772.tsv PT103811_adrenal.tsv
cp PT103811_PT103773.tsv PT103811_pituitary.tsv
cp PT103811_PT103776.tsv PT103811_gonads.tsv

cp PT104219_PT104174.tsv PT104219_heart.tsv
cp PT104219_PT104175.tsv PT104219_lung.tsv
cp PT104219_PT104178.tsv PT104219_liver.tsv
cp PT104219_PT104180.tsv PT104219_adrenal.tsv
cp PT104219_PT104181.tsv PT104219_pituitary.tsv
cp PT104219_PT104184.tsv PT104219_gonads.tsv

cp PT104883_PT104837.tsv PT104883_heart.tsv
cp PT104883_PT104838.tsv PT104883_lung.tsv
cp PT104883_PT104841.tsv PT104883_liver.tsv
cp PT104883_PT104843.tsv PT104883_adrenal.tsv
cp PT104883_PT104844.tsv PT104883_pituitary.tsv
cp PT104883_PT104847.tsv PT104883_gonads.tsv

cp PT105239_PT105194.tsv PT105239_heart.tsv
cp PT105239_PT105195.tsv PT105239_lung.tsv
cp PT105239_PT105198.tsv PT105239_liver.tsv
cp PT105239_PT105200.tsv PT105239_adrenal.tsv
cp PT105239_PT105201.tsv PT105239_pituitary.tsv
cp PT105239_PT105204.tsv PT105239_gonads.tsv

cp PT105290_PT105245.tsv PT105290_heart.tsv
cp PT105290_PT105246.tsv PT105290_lung.tsv
cp PT105290_PT105249.tsv PT105290_liver.tsv
cp PT105290_PT105251.tsv PT105290_adrenal.tsv
cp PT105290_PT105252.tsv PT105290_pituitary.tsv
cp PT105290_PT105255.tsv PT105290_gonads.tsv

cp PT103046_PT103001.tsv PT103046_heart.tsv
cp PT103046_PT103005.tsv PT103046_liver.tsv
cp PT103046_PT103007.tsv PT103046_adrenal.tsv
cp PT103046_PT103008.tsv PT103046_pituitary.tsv
cp PT103046_PT103011.tsv PT103046_gonads.tsv

cp PT103760_PT103715.tsv PT103760_heart.tsv
cp PT103760_PT103716.tsv PT103760_lung.tsv
cp PT103760_PT103719.tsv PT103760_liver.tsv
cp PT103760_PT103721.tsv PT103760_adrenal.tsv
cp PT103760_PT103722.tsv PT103760_pituitary.tsv

cp PT104270_PT104225.tsv PT104270_heart.tsv
cp PT104270_PT104226.tsv PT104270_lung.tsv
cp PT104270_PT104229.tsv PT104270_liver.tsv
cp PT104270_PT104231.tsv PT104270_adrenal.tsv
cp PT104270_PT104232.tsv PT104270_pituitary.tsv

cp PT104733_PT104888.tsv PT104733_heart.tsv
cp PT104733_PT104892.tsv PT104733_liver.tsv
cp PT104733_PT104894.tsv PT104733_adrenal.tsv
cp PT104733_PT104895.tsv PT104733_pituitary.tsv
cp PT104733_PT104898.tsv PT104733_gonads.tsv

# ---------
# autosomes
# ---------

for i in {1..20}
do
# PT101159
cp PT101159_PT101114_chr${i}.tsv PT101159_heart_chr${i}.tsv
cp PT101159_PT101115_chr${i}.tsv PT101159_lung_chr${i}.tsv
cp PT101159_PT101118_chr${i}.tsv PT101159_liver_chr${i}.tsv
cp PT101159_PT101120_chr${i}.tsv PT101159_adrenal_chr${i}.tsv
cp PT101159_PT101121_chr${i}.tsv PT101159_pituitary_chr${i}.tsv
cp PT101159_PT101124_chr${i}.tsv PT101159_gonads_chr${i}.tsv

#PT101210
cp PT101210_PT101165_chr${i}.tsv PT101210_heart_chr${i}.tsv
cp PT101210_PT101166_chr${i}.tsv PT101210_lung_chr${i}.tsv
cp PT101210_PT101169_chr${i}.tsv PT101210_liver_chr${i}.tsv
cp PT101210_PT101171_chr${i}.tsv PT101210_adrenal_chr${i}.tsv
cp PT101210_PT101172_chr${i}.tsv PT101210_pituitary_chr${i}.tsv
cp PT101210_PT101175_chr${i}.tsv PT101210_gonads_chr${i}.tsv

#PT102485
cp PT102485_PT102440_chr${i}.tsv PT102485_heart_chr${i}.tsv
cp PT102485_PT102441_chr${i}.tsv PT102485_lung_chr${i}.tsv
cp PT102485_PT102444_chr${i}.tsv PT102485_liver_chr${i}.tsv
cp PT102485_PT102446_chr${i}.tsv PT102485_adrenal_chr${i}.tsv
cp PT102485_PT102447_chr${i}.tsv PT102485_pituitary_chr${i}.tsv
cp PT102485_PT102450_chr${i}.tsv PT102485_gonads_chr${i}.tsv

#PT102842
cp PT102842_PT102797_chr${i}.tsv PT102842_heart_chr${i}.tsv
cp PT102842_PT102798_chr${i}.tsv PT102842_lung_chr${i}.tsv
cp PT102842_PT102801_chr${i}.tsv PT102842_liver_chr${i}.tsv
cp PT102842_PT102803_chr${i}.tsv PT102842_adrenal_chr${i}.tsv
cp PT102842_PT102804_chr${i}.tsv PT102842_pituitary_chr${i}.tsv
cp PT102842_PT102807_chr${i}.tsv PT102842_gonads_chr${i}.tsv

#PT103352
cp PT103352_PT103307_chr${i}.tsv PT103352_heart_chr${i}.tsv
cp PT103352_PT103308_chr${i}.tsv PT103352_lung_chr${i}.tsv
cp PT103352_PT103311_chr${i}.tsv PT103352_liver_chr${i}.tsv
cp PT103352_PT103313_chr${i}.tsv PT103352_adrenal_chr${i}.tsv
cp PT103352_PT103314_chr${i}.tsv PT103352_pituitary_chr${i}.tsv
cp PT103352_PT103317_chr${i}.tsv PT103352_gonads_chr${i}.tsv

#PT103811
cp PT103811_PT103766_chr${i}.tsv PT103811_heart_chr${i}.tsv
cp PT103811_PT103767_chr${i}.tsv PT103811_lung_chr${i}.tsv
cp PT103811_PT103770_chr${i}.tsv PT103811_liver_chr${i}.tsv
cp PT103811_PT103772_chr${i}.tsv PT103811_adrenal_chr${i}.tsv
cp PT103811_PT103773_chr${i}.tsv PT103811_pituitary_chr${i}.tsv
cp PT103811_PT103776_chr${i}.tsv PT103811_gonads_chr${i}.tsv

#PT104219
cp PT104219_PT104174_chr${i}.tsv PT104219_heart_chr${i}.tsv
cp PT104219_PT104175_chr${i}.tsv PT104219_lung_chr${i}.tsv
cp PT104219_PT104178_chr${i}.tsv PT104219_liver_chr${i}.tsv
cp PT104219_PT104180_chr${i}.tsv PT104219_adrenal_chr${i}.tsv
cp PT104219_PT104181_chr${i}.tsv PT104219_pituitary_chr${i}.tsv
cp PT104219_PT104184_chr${i}.tsv PT104219_gonads_chr${i}.tsv

#PT104883
cp PT104883_PT104837_chr${i}.tsv PT104883_heart_chr${i}.tsv
cp PT104883_PT104838_chr${i}.tsv PT104883_lung_chr${i}.tsv
cp PT104883_PT104841_chr${i}.tsv PT104883_liver_chr${i}.tsv
cp PT104883_PT104843_chr${i}.tsv PT104883_adrenal_chr${i}.tsv
cp PT104883_PT104844_chr${i}.tsv PT104883_pituitary_chr${i}.tsv
cp PT104883_PT104847_chr${i}.tsv PT104883_gonads_chr${i}.tsv

#PT105239
cp PT105239_PT105194_chr${i}.tsv PT105239_heart_chr${i}.tsv
cp PT105239_PT105195_chr${i}.tsv PT105239_lung_chr${i}.tsv
cp PT105239_PT105198_chr${i}.tsv PT105239_liver_chr${i}.tsv
cp PT105239_PT105200_chr${i}.tsv PT105239_adrenal_chr${i}.tsv
cp PT105239_PT105201_chr${i}.tsv PT105239_pituitary_chr${i}.tsv
cp PT105239_PT105204_chr${i}.tsv PT105239_gonads_chr${i}.tsv

#PT105290
cp PT105290_PT105245_chr${i}.tsv PT105290_heart_chr${i}.tsv
cp PT105290_PT105246_chr${i}.tsv PT105290_lung_chr${i}.tsv
cp PT105290_PT105249_chr${i}.tsv PT105290_liver_chr${i}.tsv
cp PT105290_PT105251_chr${i}.tsv PT105290_adrenal_chr${i}.tsv
cp PT105290_PT105252_chr${i}.tsv PT105290_pituitary_chr${i}.tsv
cp PT105290_PT105255_chr${i}.tsv PT105290_gonads_chr${i}.tsv

#PT103046
cp PT103046_PT103001_chr${i}.tsv PT103046_heart_chr${i}.tsv
cp PT103046_PT103005_chr${i}.tsv PT103046_liver_chr${i}.tsv
cp PT103046_PT103007_chr${i}.tsv PT103046_adrenal_chr${i}.tsv
cp PT103046_PT103008_chr${i}.tsv PT103046_pituitary_chr${i}.tsv
cp PT103046_PT103011_chr${i}.tsv PT103046_gonads_chr${i}.tsv

#PT103760
cp PT103760_PT103715_chr${i}.tsv PT103760_heart_chr${i}.tsv
cp PT103760_PT103716_chr${i}.tsv PT103760_lung_chr${i}.tsv
cp PT103760_PT103719_chr${i}.tsv PT103760_liver_chr${i}.tsv
cp PT103760_PT103721_chr${i}.tsv PT103760_adrenal_chr${i}.tsv
cp PT103760_PT103722_chr${i}.tsv PT103760_pituitary_chr${i}.tsv

#PT104270
cp PT104270_PT104225_chr${i}.tsv PT104270_heart_chr${i}.tsv
cp PT104270_PT104226_chr${i}.tsv PT104270_lung_chr${i}.tsv
cp PT104270_PT104229_chr${i}.tsv PT104270_liver_chr${i}.tsv
cp PT104270_PT104231_chr${i}.tsv PT104270_adrenal_chr${i}.tsv
cp PT104270_PT104232_chr${i}.tsv PT104270_pituitary_chr${i}.tsv

#PT104733
cp PT104733_PT104888_chr${i}.tsv PT104733_heart_chr${i}.tsv
cp PT104733_PT104892_chr${i}.tsv PT104733_liver_chr${i}.tsv
cp PT104733_PT104894_chr${i}.tsv PT104733_adrenal_chr${i}.tsv
cp PT104733_PT104895_chr${i}.tsv PT104733_pituitary_chr${i}.tsv
cp PT104733_PT104898_chr${i}.tsv PT104733_gonads_chr${i}.tsv

done
