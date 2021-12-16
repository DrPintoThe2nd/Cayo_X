# In this script, we want to find for each tissue, how many individuals have mean allele balance on X that is outside of the range of autosomes

import os
# # For tissues where there are all 14 individuals in each tissues (adrenal, heart, liver, pituitary)
individuals = ['PT104219', 'PT105290', 'PT102485', 'PT101159', 'PT105239', 'PT104883', 'PT103811', 'PT103352', 'PT102842', 'PT101210', 'PT103046', 'PT104733', 'PT103760', 'PT104270']
tissues = ['adrenal', 'heart', 'liver', 'pituitary']

# # For lung
lungs = ['PT104219', 'PT105290', 'PT102485', 'PT101159', 'PT105239', 'PT104883', 'PT103811', 'PT103352', 'PT102842', 'PT101210', 'PT103760', 'PT104270']

# For gonads
gonads = ['PT104219', 'PT105290', 'PT102485', 'PT101159', 'PT105239', 'PT104883', 'PT103811', 'PT103352', 'PT102842', 'PT101210', 'PT103046', 'PT104733']

count = 0
for i in individuals:
	for z in tissues:
		autosomes_mean_all = []
		autosomes_fn = os.path.join('asereadcounter/', i, i + '_' + z + '_mean_median_ratio_autosomes.tsv')
		with open(autosomes_fn, 'r') as f:
			for line in f:
				items = line.rstrip('\n').split('\t')
				if items[0] != 'chr':
					autosomes_mean_all.append(float(items[1]))
		chrX_fn = os.path.join('asereadcounter/compare_X_A/results/', i + '_' + z + '_mean_median_ratio_chrX.tsv')
		with open(chrX_fn, 'r') as f:
			for line in f:
				items = line.rstrip('\n').split('\t')
				if items[0] != 'chr':
					chrX_mean = float(items[1])
					if chrX_mean <= max(autosomes_mean_all) and chrX_mean >= min(autosomes_mean_all):
						count += 0
					else:
						count += 1
print (count)

count = 0
for i in gonads:
    autosomes_mean_all = []
    autosomes_fn = os.path.join('asereadcounter/', i, i + '_gonads_mean_median_ratio_autosomes.tsv')
    with open(autosomes_fn, 'r') as f:
        for line in f:
            items = line.rstrip('\n').split('\t')
            if items[0] != 'chr':
                autosomes_mean_all.append(float(items[1]))
    chrX_fn = os.path.join('asereadcounter/compare_X_A/results/', i + '_gonads_mean_median_ratio_chrX.tsv')
    with open(chrX_fn, 'r') as f:
        for line in f:
            items = line.rstrip('\n').split('\t')
            if items[0] != 'chr':
                chrX_mean = float(items[1])
                if chrX_mean <= max(autosomes_mean_all) and chrX_mean >= min(autosomes_mean_all):
                    count += 0
                else:
                    count += 1
print (count)

count = 0
for i in lungs:
    autosomes_mean_all = []
    autosomes_fn = os.path.join('asereadcounter/', i, i + '_lung_mean_median_ratio_autosomes.tsv')
    with open(autosomes_fn, 'r') as f:
        for line in f:
            items = line.rstrip('\n').split('\t')
            if items[0] != 'chr':
                autosomes_mean_all.append(float(items[1]))
    chrX_fn = os.path.join('asereadcounter/compare_X_A/results/', i + '_lung_mean_median_ratio_chrX.tsv')
    with open(chrX_fn, 'r') as f:
        for line in f:
            items = line.rstrip('\n').split('\t')
            if items[0] != 'chr':
                chrX_mean = float(items[1])
                if chrX_mean <= max(autosomes_mean_all) and chrX_mean >= min(autosomes_mean_all):
                    count += 0
                else:
                    count += 1
print (count)
