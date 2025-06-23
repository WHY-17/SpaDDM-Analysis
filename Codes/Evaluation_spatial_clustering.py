import squidpy as sq
import numpy as np
import pandas as pd
import scanpy as sc
# from pandas.conftest import compression
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import mantel
from scipy.stats import spearmanr
from sklearn.manifold import trustworthiness
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')



def Morans(ad, cols, coord_type='generic', **kwargs):
    col_data = []
    for col in cols:
        if pd.api.types.is_numeric_dtype(ad.obs[col]):
            col_data.append(ad.obs[col].to_list())
        else:
            col_data.append(ad.obs[col].astype('category').cat.codes)

    col_data = np.hstack(col_data).reshape(len(cols), -1).T
    ad_holder = sc.AnnData(col_data, obsm={'spatial': ad.obsm['spatial']})
    ad_holder.var_names = cols

    sq.gr.spatial_neighbors(ad_holder, coord_type=coord_type, **kwargs)
    sq.gr.spatial_autocorr(
        ad_holder,
        mode="moran",
        genes=cols,
        n_perms=100,
        n_jobs=6,
    )
    return ad_holder.uns["moranI"]


def Mantel_test(adata, col):
    if 'mental' not in adata.uns:
        adata.uns['mental'] = {}
    dist_low = squareform(pdist(adata.obsm[col], metric='euclidean'))
    dist_spatial = squareform(pdist(adata.obsm['spatial'], metric='euclidean'))
    # 曼特尔检验
    r, p_value, n = mantel(dist_low, dist_spatial, method='pearson', permutations=999)
    adata.uns['mental'][col] = (r, p_value, n)
    return adata

def Spearman_correlation(adata, col):
    if 'spearman' not in adata.uns:
        adata.uns['spearman'] = {}
    dist_low = squareform(pdist(adata.obsm[col], metric='euclidean'))
    dist_spatial = squareform(pdist(adata.obsm['spatial'], metric='euclidean'))
    r, p = spearmanr(dist_low.flatten(), dist_spatial.flatten())
    adata.uns['spearman'][col] = (r,p)
    return adata


def Trustworthiness(adata, col):
    if 'trustworthiness' not in adata.uns:
        adata.uns['trustworthiness'] = {}
    score = trustworthiness(adata.obsm['spatial'], adata.obsm[col], n_neighbors=10)
    adata.uns['trustworthiness'][col] = score
    return adata


def compute_cluster_spatial_corr(adata, emb_key='X_emb', cluster_key='cluster',
                                 spatial_key='spatial', corr_method='spearman'):
    """
        对每个聚类簇计算低维表示和空间坐标的距离矩阵之间的相关性，并存入 adata.uns['cluster_spatial_corr']

        Parameters:
            adata: AnnData 对象
            emb_key: obsm 中低维表示的 key
            cluster_key: obs 中聚类结果的 key
            spatial_key: obsm 中空间坐标的 key
            corr_method: 'spearman'（默认）或 'pearson'
        """

    X_emb = adata.obsm[emb_key]
    coords = adata.obsm[spatial_key]
    clusters = adata.obs[cluster_key]

    results = {}

    for cluster_id in clusters.unique():
        idx = clusters == cluster_id
        X_sub = X_emb[idx.values]
        coords_sub = coords[idx.values]

        if len(X_sub) < 5:
            continue

        dist_emb = pdist(X_sub, metric='euclidean')
        dist_spatial = pdist(coords_sub, metric='euclidean')

        if corr_method == 'spearman':
            r, p = spearmanr(dist_emb, dist_spatial)
        elif corr_method == 'pearson':
            from scipy.stats import pearsonr
            r, p = pearsonr(dist_emb, dist_spatial)
        else:
            raise ValueError("corr_method must be 'spearman' or 'pearson'")

        results[str(cluster_id)] = {
            'r': float(r),
            'p_value': float(p),
            'n_spots': int(len(X_sub))
        }

    # 初始化外部容器（只做一次）
    if 'cluster_spatial_corr' not in adata.uns:
        adata.uns['cluster_spatial_corr'] = {}

    adata.uns['cluster_spatial_corr'][cluster_key] = results
    print(f"✅ [{cluster_key}] 每个聚类簇的空间/表达相关性已计算完成！")
    return adata



data_file = "Data/E18_5-S1/"
adata = sc.read_h5ad(data_file+'adata_benchmarking_results.h5ad')

cols = ['totalVI', 'MultiVI', 'scMM', 'scAI', 'Seurat', 'SpatialDDM', 'SpatialGlue']

adata_morans = Morans(adata, cols)

for col in cols:
    adata = Trustworthiness(adata, col)

for col in cols:
    adata = compute_cluster_spatial_corr(adata, emb_key=col, cluster_key=col)




all_data = []
for method, clusters in adata.uns['cluster_spatial_corr'].items():
    for cluster_id, stats in clusters.items():
        all_data.append({
            'method': method,
            'cluster': cluster_id,
            'r': stats['r'],
            'p_value': stats['p_value'],
            'n_spots': stats['n_spots']
        })
df_corr = pd.DataFrame(all_data)
df = df_corr[['method', 'r']]
df.loc[df['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'
# algorithms = ["MultiVI", "totalVI", "scMM", "Seurat", "SpatialGlue", "SpaDDM"]
algorithms = ["MultiVI", "totalVI", "scMM", "scAI", "Seurat", "SpatialGlue", "SpaDDM"]
df["method"] = pd.Categorical(df["method"], categories=algorithms, ordered=True)

df.to_csv(data_file+'Trustworthiness_box_data.csv')

# 颜色优化方案
"""
algorithm_colors = {
    "SpatialDDM": "#F1AEA1",
    "Seurat": "#BACF90",
    "SpatialGlue": "#F5BDC8",
    "totalVI": "#BAD4EA",
    "MultiVI": "#D4EAE0",
    "scMM": "#DBC6E0",
    "StabMap": '#FEE281'
}

algorithm_colors = {
    "SpatialDDM": "#f9f871",
    "Seurat": "#ff9671",
    "SpatialGlue": "#ffc75f",
    "totalVI": "#d65db1",
    "MultiVI": "#845ec2",
    "scMM": "#ff6f91",
    "StabMap": '#FEE281'
}
"""

algorithm_colors = {
    "SpaDDM": "#EC3E31",
    "Seurat": "#A8D3A0",
    "SpatialGlue": "#A6D0E6",
    "totalVI": "#8582BD",
    "MultiVI": "#F8B072",
    "scMM": "#4F99C9",
    "scAI": '#FBEA2E',
}

plt.figure(figsize=(8, 3))
sns.boxplot(data=df, x="method", y="r", palette=algorithm_colors, width=0.7, linewidth=1.0, fliersize=0)
sns.stripplot(data=df, x="method", y="r", color="black", alpha=0.6, jitter=0.25, size=1.5)
plt.xlabel("")
plt.ylabel("Cluster correlation", fontsize=8)
plt.xticks(fontsize=8, rotation=0)  # X轴刻度字体
plt.yticks(fontsize=8)  # Y轴刻度字体
ax = plt.gca()  # 获取当前坐标轴
for spine in ['top', 'right', 'bottom', 'left']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(1.0)
from matplotlib.patches import Patch
handles = [Patch(color=color, label=method) for method, color in algorithm_colors.items()]
# plt.legend(handles=handles, title="", bbox_to_anchor=(0.25, 1.0), loc="upper right", fontsize=10, title_fontsize=10, frameon=False)
plt.tight_layout(pad=1.0)
plt.show()
plt.savefig(data_file+'Trustworthiness_box_plot.pdf', dpi=300)




# 准备数据
data = adata.uns['trustworthiness']
df_bar = pd.DataFrame({'method': list(data.keys()), 'value': list(data.values())})
df_bar.loc[df_bar['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'
df_bar["method"] = pd.Categorical(df_bar["method"], categories=algorithms, ordered=True)
df_bar.sort_values("method", inplace=True)
df_bar.to_csv(data_file+'Trustworthiness_data.csv')

# 画图
plt.figure(figsize=(4, 3))
sns.barplot(data=df_bar, x="method", y="value", width=0.6, palette=algorithm_colors, edgecolor="black")
# 美化
plt.ylabel("Global correlation (r)", fontsize=8)
plt.xlabel("")
plt.xticks(rotation=20, fontsize=8)
plt.yticks(fontsize=8)
# 去除刻度线
plt.tick_params(axis='both', which='both', length=0)
ax = plt.gca()  # 获取当前坐标轴
for spine in ['top', 'right', 'bottom', 'left']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(1.0)
# plt.ylim(0, 1.05)  # 可选：统一上限
plt.tight_layout(pad=1.0)
plt.show()
plt.savefig(data_file+'Trustworthiness_plot.pdf', dpi=300)






data_morans = pd.DataFrame(adata_morans['I'])
data_morans.insert(0, "method", data_morans.index)
data_morans.index = range(len(data_morans))
data_morans.loc[data_morans['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'
data_morans["method"] = pd.Categorical(data_morans["method"], categories=algorithms, ordered=True)
data_morans.sort_values("method", inplace=True)
data_morans.to_csv(data_file+'morans_data.csv')

# 画图
plt.figure(figsize=(4, 3))
sns.barplot(data=data_morans, x="method", y="I", width=0.6, palette=algorithm_colors, edgecolor="black")
# 美化
plt.ylabel("Moran's I", fontsize=8)
plt.xlabel("")
plt.xticks(rotation=20, fontsize=8)
plt.yticks(fontsize=8)
# 去除刻度线
plt.tick_params(axis='both', which='both', length=0)
ax = plt.gca()  # 获取当前坐标轴
for spine in ['top', 'right', 'bottom', 'left']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(1.0)
# plt.ylim(0, 1.05)  # 可选：统一上限
plt.tight_layout(pad=1.0)
plt.show()
plt.savefig(data_file+'Moran_plot.pdf', dpi=300)
