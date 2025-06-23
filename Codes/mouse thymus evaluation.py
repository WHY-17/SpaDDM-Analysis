import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')



df1 = pd.read_csv('Data/Dataset3_Mouse_Thymus1/Trustworthiness_box_data.csv', index_col=0)
df2 = pd.read_csv('Data/Dataset4_Mouse_Thymus2/Trustworthiness_box_data.csv', index_col=0)
df3 = pd.read_csv('Data/Dataset5_Mouse_Thymus3/Trustworthiness_box_data.csv', index_col=0)
df4 = pd.read_csv('Data/Dataset6_Mouse_Thymus4/Trustworthiness_box_data.csv', index_col=0)


# 添加 dataset 标签
df1["dataset"] = "Section 1"
df2["dataset"] = "Section 2"
df3["dataset"] = "Section 3"
df4["dataset"] = "Section 4"

# 合并数据
df_all = pd.concat([df1, df2, df3, df4], ignore_index=True)
df_all.loc[df_all['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'
# 重组数据结构
section_order = ['Section 1', 'Section 2', 'Section 3', 'Section 4']
algorithms = ["MultiVI", "totalVI", "scMM", "StabMap", "Seurat", "SpatialGlue", "SpaDDM"]


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


import matplotlib.pyplot as plt
import seaborn as sns
plt.figure(figsize=(2.8, 2.0))
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
    flierprops=dict(marker='o', markersize=1.5)
)

plt.xlabel('Section', fontsize=0)
plt.ylabel('Cluster correlation', fontsize=7)
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
plt.tight_layout(pad=1.0)
plt.savefig('Results/mouse thymus box plot-2.pdf', dpi=300)


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
plt.savefig('Results/mouse thymus evaluation legend.pdf', dpi=300)


data1 = pd.read_csv('Data/Dataset3_Mouse_Thymus1/Trustworthiness_data.csv', index_col=0)
data2 = pd.read_csv('Data/Dataset4_Mouse_Thymus2/Trustworthiness_data.csv', index_col=0)
data3 = pd.read_csv('Data/Dataset5_Mouse_Thymus3/Trustworthiness_data.csv', index_col=0)
data4 = pd.read_csv('Data/Dataset6_Mouse_Thymus4/Trustworthiness_data.csv', index_col=0)

# 添加 dataset 标签
data1["dataset"] = "Section 1"
data2["dataset"] = "Section 2"
data3["dataset"] = "Section 3"
data4["dataset"] = "Section 4"

# 合并数据
data_all = pd.concat([data1, data2, data3, data4], ignore_index=True)
data_all.loc[data_all['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'

import matplotlib.pyplot as plt
import seaborn as sns
plt.figure(figsize=(2.8, 2.0))
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
plt.ylabel('Global correlation', fontsize=8)
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
plt.tight_layout(pad=1.0)
plt.savefig('Results/mouse thymus global correlation-2.8.pdf', dpi=300)



moran1 = pd.read_csv('Data/Dataset3_Mouse_Thymus1/morans_data.csv', index_col=0)
moran2 = pd.read_csv('Data/Dataset4_Mouse_Thymus2/morans_data.csv', index_col=0)
moran3 = pd.read_csv('Data/Dataset5_Mouse_Thymus3/morans_data.csv', index_col=0)
moran4 = pd.read_csv('Data/Dataset6_Mouse_Thymus4/morans_data.csv', index_col=0)

# 添加 dataset 标签
moran1["dataset"] = "Section 1"
moran2["dataset"] = "Section 2"
moran3["dataset"] = "Section 3"
moran4["dataset"] = "Section 4"

# 合并数据
moran_all = pd.concat([moran1, moran2, moran3, moran4], ignore_index=True)
moran_all.loc[moran_all['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'

import matplotlib.pyplot as plt
import seaborn as sns
plt.figure(figsize=(2.8, 2.0))
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
plt.tight_layout(pad=1.0)
plt.savefig('Results/mouse thymus morans plot-2.8.pdf', dpi=300)