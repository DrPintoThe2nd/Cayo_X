#!/bin/bash
# macaque

cd asereadcounter/

for id in PT104219 PT105290 PT102485 PT101159 PT105239 PT104883 PT103811 PT103352 PT102842 PT101210
do

for t in heart lung liver gonads adrenal pituitary
do

Rscript compute_mean_median_ratio_all_chr_macaques.R ${id} ${t}

done
done

#
for id in PT103046 PT104733
do

for t in heart liver gonads adrenal pituitary
do

Rscript compute_mean_median_ratio_all_chr_macaques.R ${id} ${t}

done
done

#
for id in PT103760 PT104270
do

for t in heart lung liver adrenal pituitary
do

Rscript compute_mean_median_ratio_all_chr_macaques.R ${id} ${t}

done
done
