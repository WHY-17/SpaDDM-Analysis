import scanpy as sc
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score, adjusted_mutual_info_score, homogeneity_score




adata1 = sc.read_h5ad('Data/Dataset11_Human_Lymph_Node_A1/adata_benchmarking_results.h5ad')
adata2 = sc.read_h5ad('Data/Dataset12_Human_Lymph_Node_D1/adata_benchmarking_results.h5ad')


ground_truth_1 = adata1.obs['Ground Truth']
# 你想评估的各方法列名
methods = ['totalVI', 'MultiVI', 'Seurat', 'SpatialDDM', 'SpatialGlue', 'StabMap']
# 用于存储结果
results = []
for method in methods:
    preds = adata1.obs[method]
    ari = adjusted_rand_score(ground_truth_1, preds)
    nmi = normalized_mutual_info_score(ground_truth_1, preds)
    ami = adjusted_mutual_info_score(ground_truth_1, preds)
    hom = homogeneity_score(ground_truth_1, preds)
    results.append([ari, nmi, ami, hom])
# 构建 DataFrame
scores_df1 = pd.DataFrame(results, index=methods, columns=['ARI', 'NMI', 'AMI', 'Hom'])
# scores_df1.to_csv('Data/Dataset11_Human_Lymph_Node_A1/evaluation_metrics_scores.csv')



ground_truth_2 = adata2.obs['Ground Truth']
# 你想评估的各方法列名
methods = ['totalVI', 'MultiVI', 'Seurat', 'SpatialDDM', 'SpatialGlue', 'StabMap']
# 用于存储结果
results = []
for method in methods:
    preds = adata2.obs[method]
    ari = adjusted_rand_score(ground_truth_2, preds)
    nmi = normalized_mutual_info_score(ground_truth_2, preds)
    ami = adjusted_mutual_info_score(ground_truth_2, preds)
    hom = homogeneity_score(ground_truth_2, preds)
    results.append([ari, nmi, ami, hom])
# 构建 DataFrame
scores_df2 = pd.DataFrame(results, index=methods, columns=['ARI', 'NMI', 'AMI', 'Hom'])
# scores_df2.to_csv('Data/Dataset12_Human_Lymph_Node_D1/evaluation_metrics_scores.csv')



# 转换为竖长型数据
df_1 = (
    scores_df1
    .reset_index()                  # 将行索引转为列（默认列名为 'index'）
    .melt(
        id_vars=['index'],          # 保留的标识列（方法名称）
        var_name='Metric',          # 合并后的指标列名
        value_name='Score'          # 合并后的数值列名
    )
    .rename(columns={'index': 'Method'})  # 重命名列以增强可读性
)


df_2 = (
    scores_df2
    .reset_index()                  # 将行索引转为列（默认列名为 'index'）
    .melt(
        id_vars=['index'],          # 保留的标识列（方法名称）
        var_name='Metric',          # 合并后的指标列名
        value_name='Score'          # 合并后的数值列名
    )
    .rename(columns={'index': 'Method'})  # 重命名列以增强可读性
)

df_1.loc[df_1['Method'] == 'SpatialDDM', 'Method'] = 'SpaDDM'
df_2.loc[df_2['Method'] == 'SpatialDDM', 'Method'] = 'SpaDDM'

# 颜色优化方案
algorithm_colors = {
    "SpaDDM": "#EC3E31",   # 深蓝
    "Seurat": "#A8D3A0",       # 紫色
    "SpatialGlue": "#A6D0E6",  # 绿色
    "StabMap": "#4F99C9",         # 橙色
    "totalVI": "#F8B072",      # 金色
    "MultiVI": "#FBEA2E",      # 粉色
    "scMM": "#8582BD"          # 砖红色
}

metrics_order = ['ARI', 'NMI', 'AMI', 'Hom']
algorithms = ["MultiVI", "StabMap", "Seurat", "totalVI", "SpatialGlue", "SpaDDM"]

import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(4, 2))
sns.set_style("whitegrid")
ax = sns.barplot(
    x="Metric",
    y="Score",
    hue="Method",
    data=df_1,
    palette=algorithm_colors,
    order=metrics_order,
    hue_order=algorithms,
    width=0.6,
    dodge=True,
)

plt.xlabel('Section', fontsize=0)
plt.ylabel('Scores', fontsize=7)
plt.xticks(fontsize=7, rotation=0)  # X轴刻度字体
plt.yticks(fontsize=7)  # Y轴刻度字体

# 颜色覆盖增强（解决seaborn颜色透明度问题）‌:ml-citation{ref="2,4" data="citationList"}
for i, box in enumerate(ax.artists):
    hue_index = i % len(algorithms)
    alg = algorithms[hue_index]
    box.set_facecolor(algorithm_colors[alg])
    box.set_edgecolor(algorithm_colors[alg])
    box.set_alpha(0.9)  # 增强颜色饱和度
# 移除原图例
ax.legend_.remove()
plt.savefig('Results/human lymph node A1 plot-4.pdf', dpi=300)


plt.figure(figsize=(4, 2))
sns.set_style("whitegrid")
ax = sns.barplot(
    x="Metric",
    y="Score",
    hue="Method",
    data=df_2,
    palette=algorithm_colors,
    order=metrics_order,
    hue_order=algorithms,
    width=0.6,
    dodge=True,
)

plt.xlabel('Section', fontsize=0)
plt.ylabel('Scores', fontsize=7)
plt.xticks(fontsize=7, rotation=0)  # X轴刻度字体
plt.yticks(fontsize=7)  # Y轴刻度字体

# 颜色覆盖增强（解决seaborn颜色透明度问题）‌:ml-citation{ref="2,4" data="citationList"}
for i, box in enumerate(ax.artists):
    hue_index = i % len(algorithms)
    alg = algorithms[hue_index]
    box.set_facecolor(algorithm_colors[alg])
    box.set_edgecolor(algorithm_colors[alg])
    box.set_alpha(0.9)  # 增强颜色饱和度
# 移除原图例
ax.legend_.remove()
plt.savefig('Results/human lymph node D1 plot-4.pdf', dpi=300)


fig_legend = plt.figure(figsize=(8, 1.2))  # 宽一点，矮一点
legend_ax = fig_legend.add_subplot(111)
legend_ax.axis("off")  # 不显示轴
# 获取图例句柄与标签
handles, labels = ax.get_legend_handles_labels()
# 横向图例：设置 ncol 为算法数量
legend_ax.legend(
    handles,
    labels,
    loc="center",
    frameon=False,
    fontsize=8,
    ncol=len(algorithms),  # 横向排列
    handletextpad=0.8,
    columnspacing=1.5
)
plt.savefig('Results/human lymph node legend plot.pdf', dpi=300)
