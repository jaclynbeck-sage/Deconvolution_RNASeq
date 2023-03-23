# Python script that runs AutogeneS to find marker genes for each cell type.
# Ideally this would have been done using reticulate in R to be consistent with
# other scripts, but Scanpy/AutogeneS functions sometimes crash when calling 
# them with reticulate in RStudio, so it was easier to just write a python 
# script and call it from here.
#
# This script will save 3 different versions of markers from one run of 
# AutogeneS:
#   1) The set of markers with the lowest expression correlation between cell types
#   2) The set of markers with the highest expression distance between cell types
#   3) The set of markers with the best score when lowest correlation and highest
#      distance are equally weighted
#
# This script assumes the following conda setup has already been run to install
# conda and set up a conda environment with all necessary packages:
# wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.1.0-1-Linux-x86_64.sh
# bash Miniconda3-py310_23.1.0-1-Linux-x86_64.sh
# conda create -n autogenes_env python=3.10 libffi=3.3
# conda activate autogenes_env
# pip install scanpy pandas numpy scipy anndata anndata2ri autogenes scikit-learn scikit-misc

import sys
import scanpy as sc
import pandas as pd
import numpy as np
import autogenes as ag
import anndata2ri
from rpy2.robjects import r
from rpy2.robjects import ListVector
anndata2ri.activate()

def find_ag_markers(dataset, granularity, output_filename_prefix):
  r('source(file.path("functions", "FileIO_HelperFunctions.R"))')
  
  # This is necessary to enable loading of seaRef, which is a Delayed Array
  cmd = "sce = Load_SingleCell(\"" + dataset + "\",\"" + granularity + "\");" \
             + "counts(sce) = as(counts(sce), \"CsparseMatrix\"); sce"
  adata = r(cmd)
  #adata = r['Load_SingleCell'](dataset, granularity, output_type = "counts")
  
  sc.pp.filter_cells(adata, min_genes=200)
  sc.pp.filter_genes(adata, min_cells=10)
  sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=4000)
  sc.pp.normalize_total(adata, target_sum=1e4)
  
  ag.init(adata, use_highly_variable=True, celltype_key='celltype')
  ag.optimize(ngen=5000, seed=0, mode='standard', offspring_size=100,
              verbose=False)
  
  wts = {"correlation": (-1,0), "combined": (-1,1), "distance": (0, 1)}
  
  for key, wt in wts.items():
    inds = ag.select(weights=wt)
    sc_means = ag.adata()[:,inds]
    X = pd.DataFrame(sc_means.X)
    markers = sc_means.var_names
    
    maxs = X.idxmax()
    cts = sc_means.obs_names[maxs]
    
    print("Markers for " + dataset + " / " + granularity + " cell types (" + key + "):")
    print(cts.value_counts())
  
    # This isn't exactly FC because it's the mean of means, but it's good enough 
    # for figuring out a relative ordering within each cell type
    def calc_fc(col, ct):
      val1 = X.iloc[ct, col]
      non_ct = [x for x in range(X.shape[0]) if x != ct]
      val2 = np.mean(X.iloc[non_ct, col])
      return val1 / max(val2, 0.001)
    
    fc = [calc_fc(m, maxs[m]) for m in range(len(maxs))]
    
    markers_df = pd.DataFrame({"gene": markers, "celltype": cts, "fc": fc})
    markers_df = markers_df.sort_values(by = "fc", ascending = False)
    tmp = markers_df.groupby("celltype")["gene"].apply(list).to_dict()
    
    r_list = ListVector(tmp)
    r_list = r['lapply'](r_list, "unlist")
    
    output_filename = output_filename_prefix + "_" + key + ".rds"
    r['saveRDS'](r_list, output_filename)
  
  # End of function


if __name__ == "__main__":
    dataset = sys.argv[1]
    granularity = sys.argv[2]
    output_filename_prefix = sys.argv[3]

    find_ag_markers(dataset, granularity, output_filename_prefix)
