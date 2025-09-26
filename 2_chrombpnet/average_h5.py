
# Averages predicted contribution scores in counts_scores.h5 files across folds for a given cell type

import deepdish
import numpy as np

import os
import sys


celltype = sys.argv[1]
base_dir = sys.argv[2] #f"/oak/stanford/groups/wjg/skim/projects/LDA/2_chrombpnet/markerpeak_motifs_pred/contribution_scores"

avg_projected_shap = None
raw = None
avg_shap = None

for fold in range(5):
	
	print(f"\tfold {fold}")
	fold_dir = f"{base_dir}/{celltype}.fold_{fold}"
	deepshap = deepdish.io.load(f"{fold_dir}/{celltype}.fold_{fold}.counts_scores.h5")

	if raw is None:
		avg_projected_shap = deepshap["projected_shap"]["seq"]
		raw = deepshap["raw"]["seq"]
		avg_shap = deepshap["shap"]["seq"]

	else:
		assert(np.array_equal(raw, deepshap["raw"]["seq"]))
		avg_projected_shap += deepshap["projected_shap"]["seq"]
		avg_shap += deepshap["shap"]["seq"]

print("\taveraging...")
avg_projected_shap /= 5
avg_shap /= 5

print("\tsaving...")
deepshap_output = {"projected_shap": {"seq": avg_projected_shap},
                   "raw": {"seq": raw},
                   "shap": {"seq": avg_shap}}

deepdish.io.save(f"{base_dir}/{celltype}_average_counts_scores.h5", deepshap_output)