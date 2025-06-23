import scanpy as sc
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score, adjusted_mutual_info_score, homogeneity_score

adata1 = sc.read_h5ad('Data/Dataset13_Simulation1/adata_benchmarking_results.h5ad')
adata2 = sc.read_h5ad('Data/Dataset14_Simulation2/adata_benchmarking_results.h5ad')
adata3 = sc.read_h5ad('Data/Dataset15_Simulation3/adata_benchmarking_results.h5ad')
adata4 = sc.read_h5ad('Data/Dataset16_Simulation4/adata_benchmarking_results.h5ad')
adata5 = sc.read_h5ad('Data/Dataset17_Simulation5/adata_benchmarking_results.h5ad')


# 你的数据列表
adata_list = [adata1, adata2, adata3, adata4, adata5]
# 需要计算的聚类标签列
methods = ['totalVI', 'MultiVI', 'scMM', 'Seurat', 'SpatialGlue', 'SpatialDDM']
# 指标函数映射
metrics = {
    'NMI': normalized_mutual_info_score,
    'ARI': adjusted_rand_score,
    'AMI': adjusted_mutual_info_score,
    'HOM': homogeneity_score
}

df = []

for i, adata in enumerate(adata_list, 1):
    gt = adata.obs['Ground Truth']
    for method in methods:
        pred = adata.obs[method]
        # 计算所有指标
        scores = {metric_name: metric_func(gt, pred) for metric_name, metric_func in metrics.items()}
        scores['Dataset'] = f'Simulation-{i}'
        scores['Method'] = method
        df.append(scores)

# 转成DataFrame
df_all = pd.DataFrame(df)
df_all.loc[df_all['Method'] == 'SpatialDDM', 'Method'] = 'SpaDDM'

# 重组数据结构
metrics_order = ['ARI', 'NMI', 'AMI', 'HOM']
algorithms = ["MultiVI", "totalVI", "scMM", "Seurat", "SpatialGlue", "SpaDDM"]

# 把结果变成长格式 (指标名、指标值两列)
df_long = df_all.melt(
    id_vars=['Dataset', 'Method'],
    value_vars=['ARI', 'NMI', 'AMI', 'HOM'],
    var_name='Metric',
    value_name='Score'
)


# 颜色优化方案
algorithm_colors = {
    "SpaDDM": "#EC3E31",   # 深蓝
    "Seurat": "#A8D3A0",       # 紫色
    "SpatialGlue": "#A6D0E6",  # 绿色
    "totalVI": "#F8B072",      # 金色
    "MultiVI": "#FBEA2E",      # 粉色
    "scMM": "#8582BD"          # 砖红色
}

import matplotlib.pyplot as plt
import seaborn as sns
plt.figure(figsize=(7, 3))
sns.set_style("whitegrid")

ax = sns.barplot(
    x="Metric",
    y="Score",
    hue="Method",
    data=df_long,
    palette=algorithm_colors,
    order=metrics_order,
    hue_order=algorithms,
    ci="sd",        # 显示标准差误差条，也可以ci=None关闭误差条
    capsize=0.1,    # 误差条帽宽度
)

plt.xlabel('', fontsize=12)
plt.ylabel('Score', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
ax.legend_.remove()
plt.savefig('Results/Simulation clustering evaluation.pdf', dpi=300)



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

plt.savefig('Results/Simulation clustering evaluation legend.pdf', dpi=300)




