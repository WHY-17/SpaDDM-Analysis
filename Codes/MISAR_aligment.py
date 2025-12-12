import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')


adata_11_s1 = sc.read_h5ad('./Data/E11_0-S1/spatialDDM_results.h5ad')
adata_11_s2 = sc.read_h5ad('./Data/E11_0-S2/spatialDDM_results.h5ad')
adata_13_s1 = sc.read_h5ad('./Data/E13_5-S1/spatialDDM_results.h5ad')
adata_13_s2 = sc.read_h5ad('./Data/E13_5-S2/spatialDDM_results_noscale.h5ad')
adata_15_s1 = sc.read_h5ad('./Data/E15_5-S1/spatialDDM_results.h5ad')
adata_15_s2 = sc.read_h5ad('./Data/E15_5-S2/spatialDDM_results_noscale.h5ad')
adata_18_s1 = sc.read_h5ad('./Data/E18_5-S1/spatialDDM_results_noscale.h5ad')
adata_18_s2 = sc.read_h5ad('./Data/E18_5-S2/spatialDDM_results_noscale.h5ad')



# 创建子图
fig, axs = plt.subplots(4, 2, figsize=(8, 11), constrained_layout=True)
sc.pl.embedding(adata_11_s1, basis='spatial', color='Combined_Clusters', s=40, legend_fontsize=8, ax=axs[0,0])
axs[0,0].invert_yaxis()  # 颠倒 y 轴方向
# axs[0,0].spines['right'].set_visible(False)
# axs[0,0].spines['top'].set_visible(False)
# axs[0,0].spines['left'].set_visible(False)
# axs[0,0].spines['bottom'].set_visible(False)
axs[0,0].get_yaxis().set_visible(False)
axs[0,0].get_xaxis().set_visible(False)
axs[0,0].set_title('Clusters', fontsize=10)

sc.pl.embedding(adata_11_s1, basis='spatial', color='Combined_Clusters_annotation', legend_fontsize=8, s=40, ax=axs[0,1])
axs[0,1].invert_yaxis()  # 颠倒 y 轴方向
# axs[0,1].spines['right'].set_visible(False)
# axs[0,1].spines['top'].set_visible(False)
# axs[0,1].spines['left'].set_visible(False)
# axs[0,1].spines['bottom'].set_visible(False)
axs[0,1].get_yaxis().set_visible(False)
axs[0,1].get_xaxis().set_visible(False)
axs[0,1].set_title('Annotation', fontsize=10)


sc.pl.embedding(adata_13_s1, basis='spatial', color='Combined_Clusters', legend_fontsize=8, s=40, ax=axs[1,0])
axs[1,0].invert_yaxis()  # 颠倒 y 轴方向
# axs[1,0].spines['right'].set_visible(False)
# axs[1,0].spines['top'].set_visible(False)
# axs[1,0].spines['left'].set_visible(False)
# axs[1,0].spines['bottom'].set_visible(False)
axs[1,0].get_yaxis().set_visible(False)
axs[1,0].get_xaxis().set_visible(False)
axs[1,0].set_title('Clusters', fontsize=10)

sc.pl.embedding(adata_13_s1, basis='spatial', color='Combined_Clusters_annotation', legend_fontsize=8, s=40, ax=axs[1,1])
axs[1,1].invert_yaxis()  # 颠倒 y 轴方向
# axs[1,1].spines['right'].set_visible(False)
# axs[1,1].spines['top'].set_visible(False)
# axs[1,1].spines['left'].set_visible(False)
# axs[1,1].spines['bottom'].set_visible(False)
axs[1,1].get_yaxis().set_visible(False)
axs[1,1].get_xaxis().set_visible(False)
axs[1,1].set_title('Annotation', fontsize=10)


sc.pl.embedding(adata_15_s1, basis='spatial', color='Combined_Clusters', s=40, legend_fontsize=8, ax=axs[2,0])
axs[2,0].invert_yaxis()  # 颠倒 y 轴方向
# axs[2,0].spines['right'].set_visible(False)
# axs[2,0].spines['top'].set_visible(False)
# axs[2,0].spines['left'].set_visible(False)
# axs[2,0].spines['bottom'].set_visible(False)
axs[2,0].get_yaxis().set_visible(False)
axs[2,0].get_xaxis().set_visible(False)
axs[2,0].set_title('Clusters', fontsize=10)

sc.pl.embedding(adata_15_s1, basis='spatial', color='Combined_Clusters_annotation', legend_fontsize=8, s=40, ax=axs[2,1])
axs[2,1].invert_yaxis()  # 颠倒 y 轴方向
# axs[2,1].spines['right'].set_visible(False)
# axs[2,1].spines['top'].set_visible(False)
# axs[2,1].spines['left'].set_visible(False)
# axs[2,1].spines['bottom'].set_visible(False)
axs[2,1].get_yaxis().set_visible(False)
axs[2,1].get_xaxis().set_visible(False)
axs[2,1].set_title('Annotation', fontsize=10)


sc.pl.embedding(adata_18_s1, basis='spatial', color='Combined_Clusters', s=40, legend_fontsize=8, ax=axs[3,0])
axs[3,0].invert_yaxis()  # 颠倒 y 轴方向
# axs[3,0].spines['right'].set_visible(False)
# axs[3,0].spines['top'].set_visible(False)
# axs[3,0].spines['left'].set_visible(False)
# axs[3,0].spines['bottom'].set_visible(False)
axs[3,0].get_yaxis().set_visible(False)
axs[3,0].get_xaxis().set_visible(False)
axs[3,0].set_title('Clusters', fontsize=10)

sc.pl.embedding(adata_18_s1, basis='spatial', color='Combined_Clusters_annotation', s=40, legend_fontsize=8, ax=axs[3,1])
axs[3,1].invert_yaxis()  # 颠倒 y 轴方向
# axs[3,1].spines['right'].set_visible(False)
# axs[3,1].spines['top'].set_visible(False)
# axs[3,1].spines['left'].set_visible(False)
# axs[3,1].spines['bottom'].set_visible(False)
axs[3,1].get_yaxis().set_visible(False)
axs[3,1].get_xaxis().set_visible(False)
axs[3,1].set_title('Annotation', fontsize=10)

plt.savefig('Results/MISAR/Combined_Clusters_annotation.pdf', dpi=300)





"""
sc.pl.embedding(adata_11_s1, basis='spatial', groups='Midbrain', color='Combined_Clusters_annotation')
plt.gca().invert_yaxis()
sc.pl.embedding(adata_11_s1, basis='spatial', color='SpatialDDM')
plt.gca().invert_yaxis()


sc.pl.embedding(adata_13_s1, basis='spatial', groups='Midbrain', color='Combined_Clusters_annotation')
plt.gca().invert_yaxis()
sc.pl.embedding(adata_13_s1, basis='spatial', color='SpatialDDM')
plt.gca().invert_yaxis()


sc.pl.embedding(adata_15_s1, basis='spatial', groups='Midbrain', color='Combined_Clusters_annotation')
plt.gca().invert_yaxis()
sc.pl.embedding(adata_15_s1, basis='spatial', color='SpatialDDM')
plt.gca().invert_yaxis()


sc.pl.embedding(adata_18_s1, basis='spatial', groups='Midbrain', color='Combined_Clusters_annotation')
plt.gca().invert_yaxis()
sc.pl.embedding(adata_18_s1, basis='spatial', color='SpatialDDM')
plt.gca().invert_yaxis()
"""


# 提取数据
df_11 = adata_11_s1.obs.copy()
df_11[['x', 'y']] = adata_11_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df_11['Combined_Clusters_annotation'] = df_11['Combined_Clusters_annotation'].astype(str)
df_11['SpatialDDM'] = df_11['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df_11, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df_11[df_11['Combined_Clusters_annotation'] == 'Muscle'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])  # 高亮 Muscle
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df_11, x='x', y='y', color=lowlight_color, s=10, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df_11[df_11['SpatialDDM'] == '6'],
                x='x', y='y', color=highlight_color, s=15,ax=axs[1])  # 高亮 6
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
"""
# 创建单独的图例画布
fig_legend = plt.figure(figsize=(0.1, 0.6))
handles = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=lowlight_color, markersize=6, label='Other'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=highlight_color, markersize=6, label='Dpallm'),
]
fig_legend.legend(handles=handles, loc='center', frameon=False, ncol=2, fontsize=7)
fig_legend.tight_layout()
# 显示图和图例
plt.show()
plt.savefig('Results/MISAR/Dpallm_region_legend.pdf', dpi=300)
"""
plt.savefig('Results/MISAR/E11_0-S1_Muscle.pdf', dpi=300)


# 提取数据
df_13 = adata_13_s1.obs.copy()
df_13[['x', 'y']] = adata_13_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df_13['Combined_Clusters_annotation'] = df_13['Combined_Clusters_annotation'].astype(str)
df_13['SpatialDDM'] = df_13['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df_13, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df_13[df_13['Combined_Clusters_annotation'] == 'Muscle'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])  # 高亮 Muscle
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df_13, x='x', y='y', color=lowlight_color, s=15, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df_13[df_13['SpatialDDM'] == '6'],
                x='x', y='y', color=highlight_color, s=30, ax=axs[1])  # 高亮 6
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E13_5-S1_Muscle.pdf', dpi=300)


# 提取数据
df_15 = adata_15_s1.obs.copy()
df_15[['x', 'y']] = adata_15_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df_15['Combined_Clusters_annotation'] = df_15['Combined_Clusters_annotation'].astype(str)
df_15['SpatialDDM'] = df_15['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df_15, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df_15[df_15['Combined_Clusters_annotation'] == 'Muscle'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])  # 高亮 Muscle
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df_15, x='x', y='y', color=lowlight_color, s=15, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df_15[df_15['SpatialDDM'] == '7'],
                x='x', y='y', color=highlight_color, s=30, ax=axs[1])  # 高亮 6
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E15_5-S1_Muscle.pdf', dpi=300)


# 提取数据
df_18 = adata_18_s1.obs.copy()
df_18[['x', 'y']] = adata_18_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df_18['Combined_Clusters_annotation'] = df_18['Combined_Clusters_annotation'].astype(str)
df_18['SpatialDDM'] = df_18['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df_18, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df_18[df_18['Combined_Clusters_annotation'] == 'Muscle'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])  # 高亮 Muscle
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df_18, x='x', y='y', color=lowlight_color, s=15, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df_18[df_18['SpatialDDM'] == '9'],
                x='x', y='y', color=highlight_color, s=30, ax=axs[1])  # 高亮 6
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E18_5-S1_Muscle.pdf', dpi=300)

"""
# 创建子图
fig, axs = plt.subplots(1, 8, figsize=(16, 2), constrained_layout=True)
sns.scatterplot(data=df_11, x='x', y='y', color=lowlight_color, s=15, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df_11[df_11['Combined_Clusters_annotation'] == 'Muscle'],
                x='x', y='y', color=highlight_color, s=30, ax=axs[0])  # 高亮 Muscle
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual', fontsize=18)

sns.scatterplot(data=df_11, x='x', y='y', color=lowlight_color, s=15, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df_11[df_11['SpatialDDM'] == '6'],
                x='x', y='y', color=highlight_color, s=30, ax=axs[1])  # 高亮 6
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpatialDDM', fontsize=18)
plt.show()


sns.scatterplot(data=df_13, x='x', y='y', color=lowlight_color, s=15, ax=axs[2])  # 画所有点（灰色）
sns.scatterplot(data=df_13[df_13['Combined_Clusters_annotation'] == 'Muscle'],
                x='x', y='y', color=highlight_color, s=30, ax=axs[2])  # 高亮 Muscle
axs[2].invert_yaxis()  # 颠倒 y 轴方向
axs[2].spines['right'].set_visible(False)
axs[2].spines['top'].set_visible(False)
axs[2].spines['left'].set_visible(False)
axs[2].spines['bottom'].set_visible(False)
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Manual', fontsize=18)

sns.scatterplot(data=df_13, x='x', y='y', color=lowlight_color, s=15, ax=axs[3])  # 画所有点（灰色）
sns.scatterplot(data=df_13[df_13['SpatialDDM'] == '6'],
                x='x', y='y', color=highlight_color, s=30, ax=axs[3])  # 高亮 6
axs[3].invert_yaxis()  # 颠倒 y 轴方向
axs[3].spines['right'].set_visible(False)
axs[3].spines['top'].set_visible(False)
axs[3].spines['left'].set_visible(False)
axs[3].spines['bottom'].set_visible(False)
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('SpatialDDM', fontsize=18)


sns.scatterplot(data=df_15, x='x', y='y', color=lowlight_color, s=15, ax=axs[4])  # 画所有点（灰色）
sns.scatterplot(data=df_15[df_15['Combined_Clusters_annotation'] == 'Muscle'],
                x='x', y='y', color=highlight_color, s=30, ax=axs[4])  # 高亮 Muscle
axs[4].invert_yaxis()  # 颠倒 y 轴方向
axs[4].spines['right'].set_visible(False)
axs[4].spines['top'].set_visible(False)
axs[4].spines['left'].set_visible(False)
axs[4].spines['bottom'].set_visible(False)
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('Manual', fontsize=18)

sns.scatterplot(data=df_15, x='x', y='y', color=lowlight_color, s=15, ax=axs[5])  # 画所有点（灰色）
sns.scatterplot(data=df_15[df_15['SpatialDDM'] == '7'],
                x='x', y='y', color=highlight_color, s=30, ax=axs[5])  # 高亮 6
axs[5].invert_yaxis()  # 颠倒 y 轴方向
axs[5].spines['right'].set_visible(False)
axs[5].spines['top'].set_visible(False)
axs[5].spines['left'].set_visible(False)
axs[5].spines['bottom'].set_visible(False)
axs[5].get_yaxis().set_visible(False)
axs[5].get_xaxis().set_visible(False)
axs[5].set_title('SpatialDDM', fontsize=18)


sns.scatterplot(data=df_18, x='x', y='y', color=lowlight_color, s=15, ax=axs[6])  # 画所有点（灰色）
sns.scatterplot(data=df_18[df_18['Combined_Clusters_annotation'] == 'Muscle'],
                x='x', y='y', color=highlight_color, s=30, ax=axs[6])  # 高亮 Muscle
axs[6].invert_yaxis()  # 颠倒 y 轴方向
axs[6].spines['right'].set_visible(False)
axs[6].spines['top'].set_visible(False)
axs[6].spines['left'].set_visible(False)
axs[6].spines['bottom'].set_visible(False)
axs[6].get_yaxis().set_visible(False)
axs[6].get_xaxis().set_visible(False)
axs[6].set_title('Manual', fontsize=18)

sns.scatterplot(data=df_18, x='x', y='y', color=lowlight_color, s=15, ax=axs[7])  # 画所有点（灰色）
sns.scatterplot(data=df_18[df_18['SpatialDDM'] == '9'],
                x='x', y='y', color=highlight_color, s=30, ax=axs[7])  # 高亮 6
axs[7].invert_yaxis()  # 颠倒 y 轴方向
axs[7].spines['right'].set_visible(False)
axs[7].spines['top'].set_visible(False)
axs[7].spines['left'].set_visible(False)
axs[7].spines['bottom'].set_visible(False)
axs[7].get_yaxis().set_visible(False)
axs[7].get_xaxis().set_visible(False)
axs[7].set_title('SpatialDDM', fontsize=18)
"""


# 提取数据
df = adata_13_s1.obs.copy()
df[['x', 'y']] = adata_13_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df['Combined_Clusters_annotation'] = df['Combined_Clusters_annotation'].astype(str)
df['SpatialDDM'] = df['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df[df['Combined_Clusters_annotation'] == 'DPallm'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])  # 高亮 Muscle
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df[df['SpatialDDM'] == '12'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])  # 高亮 6
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E13_5-S1_DPallm.pdf', dpi=300)


# 提取数据
df = adata_15_s1.obs.copy()
df[['x', 'y']] = adata_15_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df['Combined_Clusters_annotation'] = df['Combined_Clusters_annotation'].astype(str)
df['SpatialDDM'] = df['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df[df['Combined_Clusters_annotation'] == 'DPallm'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])  # 高亮 Muscle
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df[df['SpatialDDM'] == '6'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])  # 高亮 6
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E15_5-S1_DPallm.pdf', dpi=300)


# 提取数据
df = adata_18_s1.obs.copy()
df[['x', 'y']] = adata_18_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df['Combined_Clusters_annotation'] = df['Combined_Clusters_annotation'].astype(str)
df['SpatialDDM'] = df['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df[df['Combined_Clusters_annotation'] == 'DPallm'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])  # 高亮 Muscle
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df[df['SpatialDDM'] == '10'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])  # 高亮 6
sns.scatterplot(data=df[df['SpatialDDM'] == '13'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])  # 高亮 6
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E18_5-S1_DPallm.pdf', dpi=300)



# 提取数据
df = adata_13_s1.obs.copy()
df[['x', 'y']] = adata_13_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df['Combined_Clusters_annotation'] = df['Combined_Clusters_annotation'].astype(str)
df['SpatialDDM'] = df['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df[df['Combined_Clusters_annotation'] == 'Diencephalon_and_hindbrain'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df[df['SpatialDDM'] == '5'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])
sns.scatterplot(data=df[df['SpatialDDM'] == '9'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E13_5-S1_Diencephalon_and_hindbrain.pdf', dpi=300)


# 提取数据
df = adata_15_s1.obs.copy()
df[['x', 'y']] = adata_15_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df['Combined_Clusters_annotation'] = df['Combined_Clusters_annotation'].astype(str)
df['SpatialDDM'] = df['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df[df['Combined_Clusters_annotation'] == 'Diencephalon_and_hindbrain'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df[df['SpatialDDM'] == '4'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])
sns.scatterplot(data=df[df['SpatialDDM'] == '2'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E15_5-S1_Diencephalon_and_hindbrain.pdf', dpi=300)


# 提取数据
df = adata_18_s1.obs.copy()
df[['x', 'y']] = adata_18_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df['Combined_Clusters_annotation'] = df['Combined_Clusters_annotation'].astype(str)
df['SpatialDDM'] = df['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df[df['Combined_Clusters_annotation'] == 'Diencephalon_and_hindbrain'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=11, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df[df['SpatialDDM'] == '2'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E18_5-S1_Diencephalon_and_hindbrain.pdf', dpi=300)



# 提取数据
df = adata_11_s1.obs.copy()
df[['x', 'y']] = adata_11_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df['Combined_Clusters_annotation'] = df['Combined_Clusters_annotation'].astype(str)
df['SpatialDDM'] = df['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df[df['Combined_Clusters_annotation'] == 'Mesenchyme'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df[df['SpatialDDM'] == '7'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E11_0-S1_Mesenchyme-1.8.pdf', dpi=300)


# 提取数据
df = adata_13_s1.obs.copy()
df[['x', 'y']] = adata_13_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df['Combined_Clusters_annotation'] = df['Combined_Clusters_annotation'].astype(str)
df['SpatialDDM'] = df['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df[df['Combined_Clusters_annotation'] == 'Mesenchyme'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df[df['SpatialDDM'] == '11'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E13_5-S1_Mesenchyme-1.8.pdf', dpi=300)


# 提取数据
df = adata_15_s1.obs.copy()
df[['x', 'y']] = adata_15_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df['Combined_Clusters_annotation'] = df['Combined_Clusters_annotation'].astype(str)
df['SpatialDDM'] = df['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df[df['Combined_Clusters_annotation'] == 'Mesenchyme'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df[df['SpatialDDM'] == '11'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E15_5-S1_Mesenchyme-1.8.pdf', dpi=300)


# 提取数据
df = adata_18_s1.obs.copy()
df[['x', 'y']] = adata_18_s1.obsm['spatial']  # 获取空间坐标
# 确保类别变量是字符串类型
df['Combined_Clusters_annotation'] = df['Combined_Clusters_annotation'].astype(str)
df['SpatialDDM'] = df['SpatialDDM'].astype(str)
# 选择一个高亮颜色
lowlight_color = '#D3D3D3'
highlight_color = '#41B9C1'

# 创建子图
fig, axs = plt.subplots(1, 2, figsize=(1.8, 1.1), constrained_layout=True)
sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[0])  # 画所有点（灰色）
sns.scatterplot(data=df[df['Combined_Clusters_annotation'] == 'Mesenchyme'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[0])
axs[0].invert_yaxis()  # 颠倒 y 轴方向
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=7)

sns.scatterplot(data=df, x='x', y='y', color=lowlight_color, s=10, ax=axs[1])  # 画所有点（灰色）
sns.scatterplot(data=df[df['SpatialDDM'] == '1'],
                x='x', y='y', color=highlight_color, s=15, ax=axs[1])
axs[1].invert_yaxis()  # 颠倒 y 轴方向
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaDDM', fontsize=7)
plt.show()
plt.savefig('Results/MISAR/E18_5-S1_Mesenchyme-1.8.pdf', dpi=300)

