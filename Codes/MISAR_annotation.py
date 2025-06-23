import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')



data_file = 'Data/E18_5-S1/'

adata = sc.read_h5ad(data_file+'adata_RNA.h5ad')
adata_results = sc.read_h5ad(data_file+'spatialDDM_results_noscale.h5ad')
adata_results.obs['SpatialDDM'] = adata_results.obs['SpatialDDM'].astype(str)

sc.pl.embedding(adata, basis='spatial', color='Combined_Clusters_annotation', s=70)
sc.pl.embedding(adata_results, basis='spatial', color='SpatialDDM', s=70)

labels = ['Cartilage_3', 'Cartilage_2', 'DPallv', 'Muscle', 'DPallm', 'Mesenchyme',
          'Subpallium_2', 'Cartilage_1', 'Basal_plate_of_hindbrain', 'Midbrain',
          'Diencephalon_and_hindbrain', 'Thalamus', 'Subpallium_1', 'Cartilage_4']



sc.pl.embedding(adata, basis='spatial', groups='Cartilage_4', color='Combined_Clusters_annotation', s=70)
sc.pl.embedding(adata_results, basis='spatial', groups='12', color='SpatialDDM', s=70)