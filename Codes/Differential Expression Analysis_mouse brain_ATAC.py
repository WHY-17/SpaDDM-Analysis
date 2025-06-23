import scanpy as sc
import cosg
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')



data_file = "Data/Dataset7_Mouse_Brain_ATAC/"
adata = sc.read_h5ad(data_file+'adata_benchmarking_results.h5ad')
adata.obs['SpatialDDM'] = adata.obs['SpatialDDM'].astype(str)
adata.obs['SpatialDDM'] = adata.obs['SpatialDDM'].astype('category')

plot_color = ["#458A74", "#F5A216", "#57AF37", "#41B9c1", "#008B8B",
              "#4E5689", "#6A8EC9", "#652884", "#8A7355", "#CC5B45",
              "#EC3E31", "#B46DA9", "#F59B79", "#FBEA2E", "#A8D3A0", "#D65190", "#F898CB"]
adata.uns['SpatialDDM_colors'] = plot_color

pred_label = list(adata.obs['SpatialDDM'])
mapping = {'1': 'ACB',
           '2': 'CP',
           '3': 'CLA',
           '4': 'CTX-L1',
           '5': 'CTX-L5',
           '6': 'CTX-L5',
           '7': 'CCG',
           '8': 'CTX-L2/3',
           '9': 'CTX-L6b',
           '10': 'CTX-L2/3',
           '11': 'CP',
           '12': 'CTX-L4',
           '13': 'LS',
           '14': 'CTX-L6a',
           '15': 'VL',
           '16': 'ACB',
           '17': 'PIR',
           '18': 'LS',
           }

annot = [mapping[i] for i in pred_label]
adata.obs['annotation'] = pd.Categorical(annot)
adata.uns['annotation_colors'] = plot_color[:13]

# 应该将6和8设置为None
adata_new = adata[~adata.obs['SpatialDDM'].isin(['6', '8'])].copy()

n_gene=30
groupby='annotation'
cosg.cosg(
   adata_new,
   key_added='cosg',
   # use_raw=False,
   # layer='log1p',
   mu=100,
   expressed_pct=0.1,
   remove_lowly_expressed=True,
   n_genes_user=100,
   groupby=groupby
)

sc.tl.dendrogram(adata_new, groupby=groupby, use_rep='SpatialDDM') ## Change use_rep to the cell embeddings key you'd like to use
df_tmp=pd.DataFrame(adata_new.uns['cosg']['names'][:3,]).T
df_tmp=df_tmp.reindex(adata_new.uns['dendrogram_'+groupby]['categories_ordered'])
marker_genes_list={idx: list(row.values) for idx, row in df_tmp.iterrows()}
marker_genes_list = {k: v for k, v in marker_genes_list.items() if not any(isinstance(x, float) for x in v)}


import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 8  # 全局字体大小
fig, ax = plt.subplots(figsize=(16, 6))
sc.pl.dotplot(
   adata_new,
   marker_genes_list,
   groupby=groupby,
   dendrogram=True,
   swap_axes=False,
   standard_scale='var',
   cmap='Spectral_r',
   ax=ax
 )
plt.savefig('Results/Mouse Brain/ATAC_differential_expression_analysis.pdf', dpi=300)



fig, ax = plt.subplots(figsize=(1.6, 1.8))
# plt.rcParams['font.size'] = 16
sc.pl.embedding(
    adata,
    basis='spatial',
    color='annotation',
    s=10,
    title='SpaDDM annotation',
    legend_fontsize=6,
    ax=ax
)
# 调整标题字号
plt.xlabel('spatial 1', fontsize=7)
plt.ylabel('spatial 2', fontsize=7)
ax.title.set_size(7)
plt.show()
plt.savefig('Results/Mouse Brain/ATAC_annotation-1.6.pdf', dpi=300)


fig, axs = plt.subplots(2, 8, figsize=(8, 2.2),constrained_layout=True)
sc.pl.embedding(adata, basis='spatial', groups='CP', color='annotation', show=False, size=8, ax=axs[0,0],
              legend_fontsize=0, legend_loc='on data')
# axs[0].invert_yaxis()
# axs[0,0].spines['right'].set_visible(False) # 去掉边框
# axs[0,0].spines['top'].set_visible(False)   # 去掉边框
# axs[0,0].spines['left'].set_visible(False) # 去掉边框
# axs[0,0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,0].get_yaxis().set_visible(False)
axs[0,0].get_xaxis().set_visible(False)
axs[0,0].set_title('CP', fontsize=6)

sc.pl.embedding(adata, basis='spatial', groups='CLA', color='annotation', show=False, size=8, ax=axs[0,1],
              legend_fontsize=0, legend_loc='on data')
# axs[0].invert_yaxis()
# axs[0,1].spines['right'].set_visible(False) # 去掉边框
# axs[0,1].spines['top'].set_visible(False)   # 去掉边框
# axs[0,1].spines['left'].set_visible(False) # 去掉边框
# axs[0,1].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,1].get_yaxis().set_visible(False)
axs[0,1].get_xaxis().set_visible(False)
axs[0,1].set_title('CLA', fontsize=6)

sc.pl.embedding(adata, basis='spatial', groups='CCG', color='annotation', show=False, size=8, ax=axs[0,2],
              legend_fontsize=0, legend_loc='on data')
# axs[0].invert_yaxis()
# axs[0,2].spines['right'].set_visible(False) # 去掉边框
# axs[0,2].spines['top'].set_visible(False)   # 去掉边框
# axs[0,2].spines['left'].set_visible(False) # 去掉边框
# axs[0,2].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,2].get_yaxis().set_visible(False)
axs[0,2].get_xaxis().set_visible(False)
axs[0,2].set_title('CCG', fontsize=6)

sc.pl.embedding(adata, basis='spatial', groups='CTX-L6b', color='annotation', show=False, size=8, ax=axs[0,3],
              legend_fontsize=0, legend_loc='on data')
# axs[0].invert_yaxis()
# axs[0,3].spines['right'].set_visible(False) # 去掉边框
# axs[0,3].spines['top'].set_visible(False)   # 去掉边框
# axs[0,3].spines['left'].set_visible(False) # 去掉边框
# axs[0,3].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,3].get_yaxis().set_visible(False)
axs[0,3].get_xaxis().set_visible(False)
axs[0,3].set_title('CTX-L6b', fontsize=6)

sc.pl.embedding(adata, basis='spatial', groups='CTX-L2/3', color='annotation', show=False, size=8, ax=axs[0,4],
              legend_fontsize=0, legend_loc='on data')
# axs[0].invert_yaxis()
# axs[0,4].spines['right'].set_visible(False) # 去掉边框
# axs[0,4].spines['top'].set_visible(False)   # 去掉边框
# axs[0,4].spines['left'].set_visible(False) # 去掉边框
# axs[0,4].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,4].get_yaxis().set_visible(False)
axs[0,4].get_xaxis().set_visible(False)
axs[0,4].set_title('CTX-L2/3', fontsize=6)

sc.pl.embedding(adata, basis='spatial', groups='CTX-L4', color='annotation', show=False, size=8, ax=axs[0,5],
              legend_fontsize=0, legend_loc='on data')
# axs[0].invert_yaxis()
# axs[0,5].spines['right'].set_visible(False) # 去掉边框
# axs[0,5].spines['top'].set_visible(False)   # 去掉边框
# axs[0,5].spines['left'].set_visible(False) # 去掉边框
# axs[0,5].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,5].get_yaxis().set_visible(False)
axs[0,5].get_xaxis().set_visible(False)
axs[0,5].set_title('CTX-L4', fontsize=6)

sc.pl.embedding(adata, basis='spatial', groups='LS', color='annotation', show=False, size=8, ax=axs[0,6],
              legend_fontsize=0, legend_loc='on data')
# axs[0].invert_yaxis()
# axs[0,6].spines['right'].set_visible(False) # 去掉边框
# axs[0,6].spines['top'].set_visible(False)   # 去掉边框
# axs[0,6].spines['left'].set_visible(False) # 去掉边框
# axs[0,6].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,6].get_yaxis().set_visible(False)
axs[0,6].get_xaxis().set_visible(False)
axs[0,6].set_title('LS', fontsize=6)

sc.pl.embedding(adata, basis='spatial', groups='VL', color='annotation', show=False, size=8, ax=axs[0,7],
              legend_fontsize=0, legend_loc='on data')
# axs[0].invert_yaxis()
# axs[0,7].spines['right'].set_visible(False) # 去掉边框
# axs[0,7].spines['top'].set_visible(False)   # 去掉边框
# axs[0,7].spines['left'].set_visible(False) # 去掉边框
# axs[0,7].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,7].get_yaxis().set_visible(False)
axs[0,7].get_xaxis().set_visible(False)
axs[0,7].set_title('VL', fontsize=6)


sc.pl.embedding(adata, basis='spatial', color='Rgs9', cmap='Spectral_r', show=False, size=8, ax=axs[1,0],
              colorbar_loc=None)
# axs[1,0].spines['right'].set_visible(False) # 去掉边框
# axs[1,0].spines['top'].set_visible(False)   # 去掉边框
# axs[1,0].spines['left'].set_visible(False) # 去掉边框
# axs[1,0].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,0].get_yaxis().set_visible(False)
axs[1,0].get_xaxis().set_visible(False)
axs[1,0].set_title('Rgs9', fontsize=6)

sc.pl.embedding(adata, basis='spatial', color='Nr4a2', cmap='Spectral_r',  show=False, size=8, ax=axs[1,1],
              colorbar_loc=None)
# axs[0].invert_yaxis()
# axs[1,1].spines['right'].set_visible(False) # 去掉边框
# axs[1,1].spines['top'].set_visible(False)   # 去掉边框
# axs[1,1].spines['left'].set_visible(False) # 去掉边框
# axs[1,1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,1].get_yaxis().set_visible(False)
axs[1,1].get_xaxis().set_visible(False)
axs[1,1].set_title('Nr4a2', fontsize=6)

sc.pl.embedding(adata, basis='spatial', color='Gsn', cmap='Spectral_r', show=False, size=8, ax=axs[1,2],
              colorbar_loc=None)
# axs[0].invert_yaxis()
# axs[1,2].spines['right'].set_visible(False) # 去掉边框
# axs[1,2].spines['top'].set_visible(False)   # 去掉边框
# axs[1,2].spines['left'].set_visible(False) # 去掉边框
# axs[1,2].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,2].get_yaxis().set_visible(False)
axs[1,2].get_xaxis().set_visible(False)
axs[1,2].set_title('Gsn', fontsize=6)

sc.pl.embedding(adata, basis='spatial', color='Ctgf', cmap='Spectral_r', show=False, size=8, ax=axs[1,3],
              colorbar_loc=None)
# axs[0].invert_yaxis()
# axs[1,3].spines['right'].set_visible(False) # 去掉边框
# axs[1,3].spines['top'].set_visible(False)   # 去掉边框
# axs[1,3].spines['left'].set_visible(False) # 去掉边框
# axs[1,3].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,3].get_yaxis().set_visible(False)
axs[1,3].get_xaxis().set_visible(False)
axs[1,3].set_title('Ctgf', fontsize=6)

sc.pl.embedding(adata, basis='spatial', color='Cux2', cmap='Spectral_r', show=False, size=8, ax=axs[1,4],
              colorbar_loc=None)
# axs[0].invert_yaxis()
# axs[1,4].spines['right'].set_visible(False) # 去掉边框
# axs[1,4].spines['top'].set_visible(False)   # 去掉边框
# axs[1,4].spines['left'].set_visible(False) # 去掉边框
# axs[1,4].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,4].get_yaxis().set_visible(False)
axs[1,4].get_xaxis().set_visible(False)
axs[1,4].set_title('Cux2', fontsize=6)

sc.pl.embedding(adata, basis='spatial', color='Adam33', cmap='Spectral_r', show=False, size=8, ax=axs[1,5],
              colorbar_loc=None)
# axs[0].invert_yaxis()
# axs[1,5].spines['right'].set_visible(False) # 去掉边框
# axs[1,5].spines['top'].set_visible(False)   # 去掉边框
# axs[1,5].spines['left'].set_visible(False) # 去掉边框
# axs[1,5].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,5].get_yaxis().set_visible(False)
axs[1,5].get_xaxis().set_visible(False)
axs[1,5].set_title('Adam33', fontsize=6)

sc.pl.embedding(adata, basis='spatial', color='Trpc4', cmap='Spectral_r', show=False, size=8, ax=axs[1,6],
              colorbar_loc=None)
# axs[0].invert_yaxis()
# axs[1,6].spines['right'].set_visible(False) # 去掉边框
# axs[1,6].spines['top'].set_visible(False)   # 去掉边框
# axs[1,6].spines['left'].set_visible(False) # 去掉边框
# axs[1,6].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,6].get_yaxis().set_visible(False)
axs[1,6].get_xaxis().set_visible(False)
axs[1,6].set_title('Trpc4', fontsize=6)

sc.pl.embedding(adata, basis='spatial', color='Top2a', cmap='Spectral_r', show=False, size=8, ax=axs[1,7],
              colorbar_loc=None)
# axs[0].invert_yaxis()
# axs[1,7].spines['right'].set_visible(False) # 去掉边框
# axs[1,7].spines['top'].set_visible(False)   # 去掉边框
# axs[1,7].spines['left'].set_visible(False) # 去掉边框
# axs[1,7].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,7].get_yaxis().set_visible(False)
axs[1,7].get_xaxis().set_visible(False)
axs[1,7].set_title('Top2a', fontsize=6)

plt.savefig('Results/Mouse Brain/Marker_genes-8.pdf', dpi=300)







































sc.tl.rank_genes_groups(adata, groupby='annotation', method='wilcoxon')
sc.tl.dendrogram(adata_new, groupby='annotation', use_rep='SpatialDDM')
sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=3,
    groupby='annotation',
    standard_scale='var',
    dendrogram=True,
    cmap='Spectral_r'
)

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names  # 所有 cluster 的名字
# 将每个 cluster 的 marker gene 信息合并成 DataFrame
marker_genes = pd.DataFrame({
    group: result['names'][group] for group in groups
})
print(marker_genes.head())


group = 'Cortex-L5'  # 指定你感兴趣的 cluster
top_markers = result['names'][group][:5]  # 查看前10个基因
print(f"Cluster {group} 的前10个 marker genes:")
print(top_markers)

sc.pl.embedding(adata, basis='spatial', groups='Cortex-L5', color='annotation', s=40)
sc.pl.embedding(adata, basis='spatial', color='C1ql3', s=40)