import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')


df1 = pd.read_csv('Data/Dataset7_Mouse_Brain_ATAC/Trustworthiness_box_data.csv', index_col=0)
df2 = pd.read_csv('Data/Dataset8_Mouse_Brain_H3K4me3/Trustworthiness_box_data.csv', index_col=0)
df3 = pd.read_csv('Data/Dataset9_Mouse_Brain_H3K27ac/Trustworthiness_box_data.csv', index_col=0)
df4 = pd.read_csv('Data/Dataset10_Mouse_Brain_H3K27me3/Trustworthiness_box_data.csv', index_col=0)



# 添加 dataset 标签
df1["dataset"] = "ATAC"
df2["dataset"] = "H3K4me3"
df3["dataset"] = "H3K27ac"
df4["dataset"] = "H3K27me3"

# 合并数据
df_all = pd.concat([df1, df2, df3, df4], ignore_index=True)
df_all.loc[df_all['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'
# 重组数据结构
section_order = ['ATAC', 'H3K4me3', 'H3K27ac', 'H3K27me3']
algorithms = ["MultiVI", "totalVI", "scMM", "Seurat", "SpatialGlue", "SpaDDM"]


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
plt.figure(figsize=(2.5, 2.0))
sns.set_style("whitegrid")

ax = sns.boxplot(
    x="dataset",
    y="r",
    hue="method",
    data=df_all,
    palette=algorithm_colors,
    order=section_order,
    hue_order=algorithms,
    width=0.6,
    dodge=True,
    flierprops=dict(marker='o', markersize=1.0)
)

plt.xlabel('Section', fontsize=0)
plt.ylabel('Cluster correlation', fontsize=7)
plt.xticks(fontsize=7, rotation=20)  # X轴刻度字体
plt.yticks(fontsize=7)  # Y轴刻度字体
ax = plt.gca()  # 获取当前坐标轴
for spine in ['top', 'right', 'bottom', 'left']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(1.0)
from matplotlib.patches import Patch
handles = [Patch(color=color, label=method) for method, color in algorithm_colors.items()]

# 颜色覆盖增强（解决seaborn颜色透明度问题）‌:ml-citation{ref="2,4" data="citationList"}
for i, box in enumerate(ax.artists):
    hue_index = i % len(algorithms)
    alg = algorithms[hue_index]
    box.set_facecolor(algorithm_colors[alg])
    box.set_edgecolor(algorithm_colors[alg])
    box.set_alpha(0.9)  # 增强颜色饱和度
# 移除原图例
ax.legend_.remove()
plt.tight_layout(pad=0.5)
plt.savefig('Results/mouse brain box plot-2.5.pdf', dpi=300)


data1 = pd.read_csv('Data/Dataset7_Mouse_Brain_ATAC/Trustworthiness_data.csv', index_col=0)
data2 = pd.read_csv('Data/Dataset8_Mouse_Brain_H3K4me3/Trustworthiness_data.csv', index_col=0)
data3 = pd.read_csv('Data/Dataset9_Mouse_Brain_H3K27ac/Trustworthiness_data.csv', index_col=0)
data4 = pd.read_csv('Data/Dataset10_Mouse_Brain_H3K27me3/Trustworthiness_data.csv', index_col=0)

# 添加 dataset 标签
data1["dataset"] = "ATAC"
data2["dataset"] = "H3K4me3"
data3["dataset"] = "H3K27ac"
data4["dataset"] = "H3K27me3"


# 合并数据
data_all = pd.concat([data1, data2, data3, data4], ignore_index=True)
data_all.loc[data_all['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'

import matplotlib.pyplot as plt
import seaborn as sns
plt.figure(figsize=(2.5, 2.0))
sns.set_style("whitegrid")

ax = sns.barplot(
    x="dataset",
    y="value",
    hue="method",
    data=data_all,
    palette=algorithm_colors,
    order=section_order,
    hue_order=algorithms,
    width=0.6,
    dodge=True,
)

plt.xlabel('Section', fontsize=0)
plt.ylabel('Global correlation', fontsize=7)
plt.xticks(fontsize=7, rotation=20)  # X轴刻度字体
plt.yticks(fontsize=7)  # Y轴刻度字体
ax = plt.gca()  # 获取当前坐标轴
for spine in ['top', 'right', 'bottom', 'left']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(1.0)
from matplotlib.patches import Patch
handles = [Patch(color=color, label=method) for method, color in algorithm_colors.items()]

# 颜色覆盖增强（解决seaborn颜色透明度问题）‌:ml-citation{ref="2,4" data="citationList"}
for i, box in enumerate(ax.artists):
    hue_index = i % len(algorithms)
    alg = algorithms[hue_index]
    box.set_facecolor(algorithm_colors[alg])
    box.set_edgecolor(algorithm_colors[alg])
    box.set_alpha(0.9)  # 增强颜色饱和度
# 移除原图例
ax.legend_.remove()
plt.tight_layout(pad=0.5)
plt.savefig('Results/mouse brain trustworthiness plot-2.5.pdf', dpi=300)


moran1 = pd.read_csv('Data/Dataset7_Mouse_Brain_ATAC/morans_data.csv', index_col=0)
moran2 = pd.read_csv('Data/Dataset8_Mouse_Brain_H3K4me3/morans_data.csv', index_col=0)
moran3 = pd.read_csv('Data/Dataset9_Mouse_Brain_H3K27ac/morans_data.csv', index_col=0)
moran4 = pd.read_csv('Data/Dataset10_Mouse_Brain_H3K27me3/morans_data.csv', index_col=0)

# 添加 dataset 标签
moran1["dataset"] = "ATAC"
moran2["dataset"] = "H3K4me3"
moran3["dataset"] = "H3K27ac"
moran4["dataset"] = "H3K27me3"

# 合并数据
moran_all = pd.concat([moran1, moran2, moran3, moran4], ignore_index=True)
moran_all.loc[moran_all['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'

import matplotlib.pyplot as plt
import seaborn as sns
plt.figure(figsize=(2.5, 2.0))
sns.set_style("whitegrid")

ax = sns.barplot(
    x="dataset",
    y="I",
    hue="method",
    data=moran_all,
    palette=algorithm_colors,
    order=section_order,
    hue_order=algorithms,
    width=0.6,
    dodge=True,
)

plt.xlabel('Section', fontsize=0)
plt.ylabel("Moran's I", fontsize=7)

plt.xticks(fontsize=7, rotation=20)  # X轴刻度字体
plt.yticks(fontsize=7)  # Y轴刻度字体
ax = plt.gca()  # 获取当前坐标轴
for spine in ['top', 'right', 'bottom', 'left']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(1.0)
from matplotlib.patches import Patch
handles = [Patch(color=color, label=method) for method, color in algorithm_colors.items()]

# 颜色覆盖增强（解决seaborn颜色透明度问题）‌:ml-citation{ref="2,4" data="citationList"}
for i, box in enumerate(ax.artists):
    hue_index = i % len(algorithms)
    alg = algorithms[hue_index]
    box.set_facecolor(algorithm_colors[alg])
    box.set_edgecolor(algorithm_colors[alg])
    box.set_alpha(0.9)  # 增强颜色饱和度
# 移除原图例
ax.legend_.remove()
plt.tight_layout(pad=0.5)
plt.savefig('Results/mouse brain morans plot-2.5.pdf', dpi=300)




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
    fontsize=8,
    ncol=1,  # 横向排列
    handletextpad=0.8,
    columnspacing=1.5
)
plt.savefig('Results/mouse brain legend plot.pdf', dpi=300)