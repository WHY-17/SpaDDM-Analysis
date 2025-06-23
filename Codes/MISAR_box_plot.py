import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


data_ari = {
    "totalVI": [0.1903, 0.1151, 0.0278, 0.1510, 0.1390, 0.1088, 0.2610, 0.1016],
    "MultiVI": [0.0597, 0.0197, 0.0407, 0.1108, 0.0812, 0.0460, 0.0415, 0.0469],
    "scAI": [0.3735, 0.3271, 0.1410, 0.2873, 0.4284, 0.3914, 0.3187, 0.2354],
    "Seurat": [0.4235, 0.1991, 0.1432, 0.3111, 0.5204, 0.4200, 0.3339, 0.2865],
    "SpatialGlue": [0.3576, 0.2512, 0.1239, 0.2319, 0.4634, 0.4465, 0.5537, 0.3730],
    "SpaDDM": [0.4645, 0.2923, 0.1426, 0.2804, 0.5363, 0.4103, 0.5927, 0.4051],
    "scMM": [0.3364, 0.2022, 0.2960, 0.2661, 0.4344, 0.2954, 0.1962, 0.1519]
}

data_nmi = {
    "totalVI": [0.2729, 0.2356, 0.1230, 0.3283, 0.2826, 0.2303, 0.3507, 0.1656],
    "MultiVI": [0.0681, 0.0688, 0.0875, 0.1511, 0.1188, 0.1084, 0.1020, 0.0955],
    "scAI": [0.4692, 0.3778, 0.2653, 0.4286, 0.4824, 0.5153, 0.4702, 0.4244],
    "Seurat": [0.5225, 0.3461, 0.3644, 0.5354, 0.5984, 0.5994, 0.5698, 0.5017],
    "SpatialGlue": [0.4771, 0.4219, 0.2895, 0.4939, 0.6198, 0.6294, 0.5934, 0.5367],
    "SpaDDM": [0.5094, 0.4146, 0.3730, 0.5096, 0.6014, 0.6056, 0.5939, 0.5381],
    "scMM": [0.4043, 0.2690, 0.3043, 0.3968, 0.4276, 0.4263, 0.2793, 0.2459]
}

data_ami = {
    "totalVI": [0.2647, 0.2171, 0.1084, 0.3168, 0.2722, 0.2123, 0.3378, 0.1452],
    "MultiVI": [0.0577, 0.0466, 0.0722, 0.1363, 0.1060, 0.0876, 0.0861, 0.0741],
    "scAI": [0.4630, 0.3626, 0.2530, 0.4185, 0.4745, 0.5034, 0.4608, 0.4108],
    "Seurat": [0.5172, 0.3306, 0.3540, 0.5274, 0.5924, 0.5899, 0.5621, 0.4901],
    "SpatialGlue": [0.4712, 0.4080, 0.2774, 0.4853, 0.6141, 0.6207, 0.5859, 0.5256],
    "SpaDDM": [0.5039, 0.4015, 0.3626, 0.5012, 0.5955, 0.5963, 0.5863, 0.5271],
    "scMM": [0.3972, 0.2506, 0.2921, 0.3859, 0.4181, 0.4125, 0.2652, 0.2269]
}

data_hom = {
    "totalVI": [0.3162, 0.2688, 0.1465, 0.3728, 0.3029, 0.2500, 0.3383, 0.1782],
    "MultiVI": [0.0787, 0.0800, 0.1023, 0.1696, 0.1268, 0.1182, 0.1124, 0.1067],
    "scAI": [0.5175, 0.4277, 0.3182, 0.6123, 0.6178, 0.6419, 0.6195, 0.5652],
    "Seurat": [0.6176, 0.4125, 0.4466, 0.6123, 0.6178, 0.6419, 0.6195, 0.5652],
    "SpatialGlue": [0.5569, 0.4887, 0.3388, 0.5677, 0.6447, 0.6779, 0.6213, 0.5898],
    "SpaDDM": [0.5844, 0.4788, 0.4528, 0.5845, 0.6329, 0.6554, 0.6150, 0.6035],
    "scMM": [0.4285, 0.2736, 0.3127, 0.4228, 0.3990, 0.4418, 0.2534, 0.2474]
}

"""
df = pd.DataFrame(data_hom)
fig, ax = plt.subplots(figsize=(13,5))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
sns.boxplot(data=combined_df, palette='viridis')
sns.stripplot(data=combined_df, color="orange", jitter=0.2, size=4)
plt.xticks(["totalVI", "MultiVI", "scAI", "Seurat", "SpatialGlue", "SpatialDDM", "scMM"],
           rotation=10, fontsize=10)
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=10)
plt.xlabel('Methods', fontsize=12)
plt.ylabel('ARI score', fontsize=12)
plt.savefig('C:/Users/haiyu/Desktop/true_label_box.pdf', dpi=300)
"""


import pandas as pd
# 重组数据结构
metrics_order = ['ARI', 'NMI', 'AMI', 'HOM']
algorithms = ["MultiVI", "totalVI", "scMM", "scAI", "Seurat", "SpatialGlue", "SpaDDM"]

# 创建横向长格式数据
combined_df = pd.DataFrame(
    [(alg, metric, value)
     for metric in metrics_order
     for alg in algorithms
     for value in eval(f"data_{metric.lower()}")[alg]],
    columns=['Algorithm', 'Metric', 'Value']
)



# # 定义算法颜色映射字典 ‌:ml-citation{ref="1,2" data="citationList"}
# algorithm_colors = {
#     "SpatialDDM": "#845ec2",  # 珊瑚红
#     "Seurat": "#d65db1",      # 紫晶色
#     "SpatialGlue": "#2c73d2", # 橄榄绿
#     "scAI": "#ff6f91",        # 海军蓝
#     "totalVI": "#ffc75f",     # 橙色
#     "MultiVI": "#f9f871",      # 棕红
#     "scMM": "#ff9671"          # 蓝绿色
# }

# 颜色优化方案
# algorithm_colors = {
#     "SpatialDDM": "#1F77B4",   # 深蓝
#     "Seurat": "#9467BD",       # 紫色
#     "SpatialGlue": "#2CA02C",  # 绿色
#     "scAI": "#FF7F0E",         # 橙色
#     "totalVI": "#FFD700",      # 金色
#     "MultiVI": "#E377C2",      # 粉色
#     "scMM": "#D62728"          # 砖红色
# }

# 颜色优化方案
algorithm_colors = {
    "SpaDDM": "#EC3E31",
    "Seurat": "#A8D3A0",
    "SpatialGlue": "#A6D0E6",
    "scAI": "#4F99C9",
    "totalVI": "#F8B072",
    "MultiVI": "#FBEA2E",
    "scMM": "#8582BD"
}



import matplotlib.pyplot as plt
import seaborn as sns
plt.figure(figsize=(8, 2.5))
sns.set_style("whitegrid")


ax = sns.boxplot(
    x="Metric",
    y="Value",
    hue="Algorithm",
    data=combined_df,
    palette=algorithm_colors,
    order=metrics_order,
    hue_order=algorithms,
    width=0.7,
    dodge=True,
    flierprops=dict(marker='o', markersize=2)
)

plt.xlabel('Metrics', fontsize=0)
plt.ylabel('Score', fontsize=7)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)

# 颜色覆盖增强（解决seaborn颜色透明度问题）‌:ml-citation{ref="2,4" data="citationList"}
for i, box in enumerate(ax.artists):
    hue_index = i % len(algorithms)
    alg = algorithms[hue_index]
    box.set_facecolor(algorithm_colors[alg])
    box.set_edgecolor(algorithm_colors[alg])
    box.set_alpha(0.9)  # 增强颜色饱和度

# 移除原图例
ax.legend_.remove()
plt.savefig('Results/MISAR/box_plot.pdf', dpi=300)


fig_legend = plt.figure(figsize=(0.5, 2))  # 宽一点，矮一点
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
    fontsize=7,
    ncol=1,  # 横向排列
    handletextpad=0.8,
    columnspacing=1.5
)

plt.savefig('Results/MISAR/box_plot_legend.pdf', dpi=300)