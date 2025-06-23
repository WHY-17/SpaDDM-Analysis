import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from Codes.utils import mclust_R



os.environ['R_HOME'] = 'D:\software\R-4.3.3'
data_file = 'Data/Dataset12_Human_Lymph_Node_D1/'

# load data
adata = sc.read_h5ad(data_file+'adata_RNA.h5ad')
# Preprocessing the data
adata.var_names_make_unique()
adata.raw = adata
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]


adata_totalvi = sc.read_h5ad(data_file+'_totalVI_results.h5ad')
adata_multivi = sc.read_h5ad(data_file+'MultiVI_results.h5ad')
adata_ddm = sc.read_h5ad(data_file+'SpatialDDM_results.h5ad')
adata_seurat = sc.read_h5ad(data_file+'Seurat_results.h5ad')
adata_spatialglue = sc.read_h5ad(data_file+'SpatialGlue_results.h5ad')
stabmap_emb = pd.read_csv(data_file+'StabMap_embeddings.csv', index_col=0)
stabmap_cluster = pd.read_csv(data_file+'StabMap_clustering_results.csv', index_col=0)

adata.obsm['totalVI'] = adata_totalvi.obsm['X_totalVI']
adata.obsm['MultiVI'] = adata_multivi.obsm['X_MultiVI']
adata.obsm['SpatialDDM'] = adata_ddm.obsm['SpatialDDM']
# adata.obsm['SpatialDDM_pca'] = adata_ddm.obsm['SpatialDDM_pca']
adata.obsm["StabMap"] = stabmap_emb.values
adata.obsm['Seurat'] = adata_seurat.obsm['X_apca']
adata.obsm['SpatialGlue'] = adata_spatialglue.obsm['SpatialGlue']


emb_names = ['totalVI', 'MultiVI']
for emb_name in emb_names:
    adata = mclust_R(adata, num_cluster=6, used_obsm=emb_name)


adata.obs['Seurat'] = adata_seurat.obs['seurat_clusters'].astype('category')
adata.obs['SpatialDDM'] = adata_ddm.obs['SpatialDDM']
adata.obs['SpatialGlue'] = adata_spatialglue.obs['SpatialGlue']
adata.obs['StabMap'] = stabmap_cluster['seurat_clusters'].astype('category')
adata.obs['Ground Truth'] = adata_ddm.obs['Ground Truth']



plot_color = ['#F59B7B', '#ED8828','#A4DDD3', '#81B21F', '#8D73BA', '#ABD7EC']
adata.uns['SpatialDDM_colors'] = plot_color
adata.uns['SpatialGlue_colors'] = plot_color
adata.uns['Seurat_colors'] = plot_color
adata.uns['StabMap_colors'] = plot_color
adata.uns['MultiVI_colors'] = plot_color
adata.uns['totalVI_colors'] = plot_color



fig, axs = plt.subplots(1, 6, figsize=(8,1.33),constrained_layout=True)
sc.pl.embedding(adata, basis='spatial', color=['MultiVI'], ax=axs[0],
                show=False, size=10, legend_fontsize=0, legend_loc='on data')
# axs[0].invert_yaxis()
# axs[0].spines['right'].set_visible(False) # 去掉边框
# axs[0].spines['top'].set_visible(False)   # 去掉边框
# axs[0].spines['left'].set_visible(False) # 去掉边框
# axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('', fontsize=15)

sc.pl.embedding(adata, basis='spatial', color=['totalVI'], show=False, size=10, ax=axs[1],
              legend_fontsize=0, legend_loc='on data')
# axs[1].invert_yaxis()
# axs[1].spines['right'].set_visible(False) # 去掉边框
# axs[1].spines['top'].set_visible(False)   # 去掉边框
# axs[1].spines['left'].set_visible(False) # 去掉边框
# axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('', fontsize=15)

sc.pl.embedding(adata, basis='spatial', color=['StabMap'], show=False, size=10, ax=axs[2],
              legend_fontsize=0, legend_loc='on data')
# axs[2].invert_yaxis()
# axs[2].spines['right'].set_visible(False) # 去掉边框
# axs[2].spines['top'].set_visible(False)   # 去掉边框
# axs[2].spines['left'].set_visible(False) # 去掉边框
# axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('', fontsize=15)

sc.pl.embedding(adata, basis='spatial', color=['Seurat'], show=False, size=10, ax=axs[3],
              legend_fontsize=0, legend_loc='on data')
# axs[3].invert_yaxis()
# axs[3].spines['right'].set_visible(False) # 去掉边框
# axs[3].spines['top'].set_visible(False)   # 去掉边框
# axs[3].spines['left'].set_visible(False) # 去掉边框
# axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('', fontsize=15)

sc.pl.embedding(adata, basis='spatial', color=['SpatialGlue'], show=False, size=10, ax=axs[4],
              legend_fontsize=0, legend_loc='on data')
# axs[4].invert_yaxis()
# axs[4].spines['right'].set_visible(False) # 去掉边框
# axs[4].spines['top'].set_visible(False)   # 去掉边框
# axs[4].spines['left'].set_visible(False) # 去掉边框
# axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('', fontsize=15)

sc.pl.embedding(adata, basis='spatial', color=['SpatialDDM'], show=False, size=10, ax=axs[5],
              legend_fontsize=0, legend_loc='on data')
# axs[5].invert_yaxis()
# axs[5].spines['right'].set_visible(False) # 去掉边框
# axs[5].spines['top'].set_visible(False)   # 去掉边框
# axs[5].spines['left'].set_visible(False) # 去掉边框
# axs[5].spines['bottom'].set_visible(False)   # 去掉边框
axs[5].get_yaxis().set_visible(False)
axs[5].get_xaxis().set_visible(False)
axs[5].set_title('', fontsize=15)

plt.savefig('Results/Human_Lymph_Node/D1_LN_embedding-8.pdf', dpi=300)


Umap_MultiVI = sc.AnnData(adata.obsm['MultiVI'])
Umap_MultiVI.obs_names = adata.obs_names
Umap_MultiVI.obs['MultiVI'] = adata.obs['MultiVI']
sc.pp.neighbors(Umap_MultiVI)
sc.tl.umap(Umap_MultiVI)
# sc.pl.umap(Umap_MultiVI, color='MultiVI')

Umap_totalVI = sc.AnnData(adata.obsm['totalVI'])
Umap_totalVI.obs_names = adata.obs_names
Umap_totalVI.obs['totalVI'] = adata.obs['totalVI']
sc.pp.neighbors(Umap_totalVI)
sc.tl.umap(Umap_totalVI)

Umap_StabMap = sc.AnnData(adata.obsm['StabMap'])
Umap_StabMap.obs_names = adata.obs_names
Umap_StabMap.obs['StabMap'] = adata.obs['StabMap']
sc.pp.neighbors(Umap_StabMap)
sc.tl.umap(Umap_StabMap)

Umap_Seurat = sc.AnnData(adata_seurat.obsm['X_apca'])
Umap_Seurat.obs_names = adata.obs_names
Umap_Seurat.obs['Seurat'] = adata.obs['Seurat']
Umap_Seurat.obsm['X_umap'] = adata_seurat.obsm['X_wnn.umap']

Umap_SpatialGlue = sc.AnnData(adata.obsm['SpatialGlue'])
Umap_SpatialGlue.obs_names = adata.obs_names
Umap_SpatialGlue.obs['SpatialGlue'] = adata.obs['SpatialGlue']
sc.pp.neighbors(Umap_SpatialGlue)
sc.tl.umap(Umap_SpatialGlue)

Umap_SpatialDDM = sc.AnnData(adata.obsm['SpatialDDM'])
Umap_SpatialDDM.obs_names = adata.obs_names
Umap_SpatialDDM.obs['SpatialDDM'] = adata.obs['SpatialDDM']
sc.pp.neighbors(Umap_SpatialDDM)
sc.tl.umap(Umap_SpatialDDM)

sc.pp.neighbors(adata, use_rep='SpatialDDM')
sc.tl.umap(adata)

adata.obsm['MultiVI_umap'] = Umap_MultiVI.obsm['X_umap']
adata.obsm['totalVI_umap'] = Umap_totalVI.obsm['X_umap']
adata.obsm['StabMap_umap'] = Umap_StabMap.obsm['X_umap']
adata.obsm['Seurat_umap'] = Umap_Seurat.obsm['X_umap']
adata.obsm['SpatialGlue_umap'] = Umap_SpatialGlue.obsm['X_umap']
adata.obsm['SpatialDDM_umap'] = Umap_SpatialDDM.obsm['X_umap']
adata.write_h5ad(data_file+'adata_benchmarking_results.h5ad', compression='gzip')



Umap_SpatialDDM.uns['SpatialDDM_colors'] = plot_color
Umap_SpatialGlue.uns['SpatialGlue_colors'] = plot_color
Umap_Seurat.uns['Seurat_colors'] = plot_color
Umap_StabMap.uns['scAI_colors'] = plot_color
Umap_MultiVI.uns['MultiVI_colors'] = plot_color
Umap_totalVI.uns['totalVI_colors'] = plot_color



fig, axs = plt.subplots(1, 6, figsize=(8,1.33),constrained_layout=True)
sc.pl.umap(Umap_MultiVI, color='MultiVI', ax=axs[0], show=False,
           size=10, legend_loc='on data', legend_fontsize=0)
# axs[0].spines['right'].set_visible(False) # 去掉边框
# axs[0].spines['top'].set_visible(False)   # 去掉边框
# axs[0].spines['left'].set_visible(False) # 去掉边框
# axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('', fontsize=15)

sc.pl.umap(Umap_totalVI, color='totalVI', ax=axs[1], show=False,
           size=10, legend_loc='on data', legend_fontsize=0)
# axs[1].spines['right'].set_visible(False) # 去掉边框
# axs[1].spines['top'].set_visible(False)   # 去掉边框
# axs[1].spines['left'].set_visible(False) # 去掉边框
# axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('', fontsize=15)

sc.pl.umap(Umap_StabMap, color='StabMap', ax=axs[2], show=False,
           size=10, legend_loc='on data', legend_fontsize=0)
# axs[2].spines['right'].set_visible(False) # 去掉边框
# axs[2].spines['top'].set_visible(False)   # 去掉边框
# axs[2].spines['left'].set_visible(False) # 去掉边框
# axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('', fontsize=15)

sc.pl.umap(Umap_Seurat, color='Seurat', ax=axs[3], show=False,
           size=10, legend_loc='on data', legend_fontsize=0)
# axs[3].spines['right'].set_visible(False) # 去掉边框
# axs[3].spines['top'].set_visible(False)   # 去掉边框
# axs[3].spines['left'].set_visible(False) # 去掉边框
# axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('', fontsize=15)

sc.pl.umap(Umap_SpatialGlue, color='SpatialGlue', ax=axs[4], show=False,
           size=10, legend_loc='on data', legend_fontsize=0)
# axs[4].spines['right'].set_visible(False) # 去掉边框
# axs[4].spines['top'].set_visible(False)   # 去掉边框
# axs[4].spines['left'].set_visible(False) # 去掉边框
# axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('', fontsize=15)

sc.pl.umap(Umap_SpatialDDM, color='SpatialDDM', ax=axs[5], show=False,
           size=10, legend_loc='on data', legend_fontsize=0)
# axs[5].spines['right'].set_visible(False) # 去掉边框
# axs[5].spines['top'].set_visible(False)   # 去掉边框
# axs[5].spines['left'].set_visible(False) # 去掉边框
# axs[5].spines['bottom'].set_visible(False)   # 去掉边框
axs[5].get_yaxis().set_visible(False)
axs[5].get_xaxis().set_visible(False)
axs[5].set_title('', fontsize=15)

plt.savefig('Results/Human_Lymph_Node/D1_LN_Umap-8.pdf', dpi=300)


### 绘制图例###
# 1. 绘制 SpaDDM 的 spatial 图，用于获取 legend（不显示）
ax = sc.pl.embedding(adata, basis='spatial', color=['SpatialDDM'], s=10, show=False)
# 2. 兼容性处理：如果是列表，取第一个 Axes
if isinstance(ax, list):
    ax = ax[0]
# 3. 创建一个新的图例画布
fig_legend = plt.figure(figsize=(8, 1.1))  # 根据图例数量自定义尺寸
legend_ax = fig_legend.add_subplot(111)
legend_ax.axis("off")
# 4. 提取图例句柄与标签
handles, labels = ax.get_legend_handles_labels()
# 5. 横向排布图例
legend_ax.legend(
    handles,
    labels,
    loc="center",
    frameon=False,
    fontsize=8,
    ncol=6,           # 图例横向排列列数（根据需要调整）
    handletextpad=0.8,
    columnspacing=1.5
)
plt.show()
plt.savefig('Results/Human_Lymph_Node/A1_LN_legend.pdf', dpi=300)
