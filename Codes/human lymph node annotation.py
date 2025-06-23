import scanpy as sc
import pandas as pd


data_file1 = 'Data/Dataset11_Human_Lymph_Node_A1/'
adata1 = sc.read_h5ad(data_file1 + 'adata_benchmarking_results.h5ad')
adata1.obs['SpatialDDM'] = adata1.obs['SpatialDDM'].astype(str)

annot1 = list(adata1.obs['SpatialDDM'])
mapping1 = {'1': 'medulla',
           '2': 'capsule',
           '3': 'cortex',
           '4': 'medulla',
           '5': 'follicle/subcapsular sinus',
           '6': 'pericapsular adipose tissue',
           }

annot1 = [mapping1[i] for i in annot1]
adata1.obs['annotation'] = pd.Categorical(annot1)
# sc.pl.embedding(adata1, basis='spatial', color='annotation', s=80)


import cosg
n_gene=30
groupby='annotation'
cosg.cosg(
   adata1,
   key_added='cosg',
   # use_raw=False,
   # layer='log1p',
   mu=100,
   expressed_pct=0.1,
   remove_lowly_expressed=True,
   n_genes_user=100,
   groupby=groupby
)

sc.tl.dendrogram(adata1, groupby=groupby, use_rep='SpatialDDM') ## Change use_rep to the cell embeddings key you'd like to use
df_tmp1=pd.DataFrame(adata1.uns['cosg']['names'][:5,]).T
df_tmp1=df_tmp1.reindex(adata1.uns['dendrogram_'+groupby]['categories_ordered'])
marker_genes_list1={idx: list(row.values) for idx, row in df_tmp1.iterrows()}
marker_genes_list1 = {k: v for k, v in marker_genes_list1.items() if not any(isinstance(x, float) for x in v)}

"""
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 18  # 全局字体大小
fig, ax = plt.subplots(figsize=(16, 6))
sc.pl.dotplot(
   adata1,
   marker_genes_list1,
   groupby=groupby,
   dendrogram=True,
   swap_axes=False,
   standard_scale='var',
   cmap='Spectral_r',
   ax=ax
 )
"""




data_file = 'Data/Dataset12_Human_Lymph_Node_D1/'
adata = sc.read_h5ad(data_file + 'adata_benchmarking_results.h5ad')
adata.obs['SpatialDDM'] = adata.obs['SpatialDDM'].astype(str)


annot = list(adata.obs['SpatialDDM'])
mapping = {'1': 'medulla',
           '6': 'capsule',
           '3': 'cortex',
           '2': 'medulla',
           '5': 'follicle/subcapsular sinus',
           '4': 'pericapsular adipose tissue',
           }

annot = [mapping[i] for i in annot]
adata.obs['annotation'] = pd.Categorical(annot)

sc.pl.embedding(adata, basis='spatial', color='annotation', s=80)



import cosg
n_gene=30
groupby='annotation'
cosg.cosg(
   adata,
   key_added='cosg',
   # use_raw=False,
   # layer='log1p',
   mu=100,
   expressed_pct=0.1,
   remove_lowly_expressed=True,
   n_genes_user=100,
   groupby=groupby
)

sc.tl.dendrogram(adata, groupby=groupby, use_rep='SpatialDDM') ## Change use_rep to the cell embeddings key you'd like to use
df_tmp=pd.DataFrame(adata.uns['cosg']['names'][:5,]).T
df_tmp=df_tmp.reindex(adata.uns['dendrogram_'+groupby]['categories_ordered'])
marker_genes_list={idx: list(row.values) for idx, row in df_tmp.iterrows()}
marker_genes_list = {k: v for k, v in marker_genes_list.items() if not any(isinstance(x, float) for x in v)}




import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 2, figsize=(11, 3),constrained_layout=True)
sc.pl.embedding(adata1, basis='spatial', color='Ground Truth', show=False, size=40, ax=axs[0],
              legend_fontsize=8)
# axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Ground Truth', fontsize=10)

sc.pl.embedding(adata1, basis='spatial', color='annotation', show=False, size=40, ax=axs[1],
              legend_fontsize=8)
# axs[0].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM annotation', fontsize=10)
plt.savefig('Results/Human_Lymph_Node/A1_LN_SpatialDDM_annotation.pdf', dpi=300)


fig, axs = plt.subplots(1, 2, figsize=(11, 3), constrained_layout=True)
sc.pl.embedding(adata, basis='spatial', color='Ground Truth', show=False, size=40, ax=axs[0],
              legend_fontsize=8)
# axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Ground Truth', fontsize=10)

sc.pl.embedding(adata, basis='spatial', color='annotation', show=False, size=40, ax=axs[1],
              legend_fontsize=8)
# axs[0].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM annotation', fontsize=10)
plt.savefig('Results/Human_Lymph_Node/D1_LN_SpatialDDM_annotation.pdf', dpi=300)




import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 1, figsize=(9, 8.5),constrained_layout=True)
sc.pl.dotplot(adata1, marker_genes_list1, groupby=groupby, dendrogram=True, swap_axes=False,
              standard_scale='var', cmap='Spectral_r', ax=axs[0])
# axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
# axs[0].set_title('Section 1 differential expression analysis', fontsize=12)

sc.pl.dotplot(adata, marker_genes_list, groupby=groupby, dendrogram=True, swap_axes=False,
              standard_scale='var', cmap='Spectral_r', ax=axs[1])
# axs[0].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
# axs[1].set_title('Section 1 differential expression analysis', fontsize=12)

plt.savefig('Results/Human_Lymph_Node/DEGs analysis.pdf', dpi=300)