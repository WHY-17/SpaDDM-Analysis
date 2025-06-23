import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')



df1 = pd.read_csv('Data/Dataset13_Simulation1/Trustworthiness_box_data.csv', index_col=0)
df2 = pd.read_csv('Data/Dataset14_Simulation2/Trustworthiness_box_data.csv', index_col=0)
df3 = pd.read_csv('Data/Dataset15_Simulation3/Trustworthiness_box_data.csv', index_col=0)
df4 = pd.read_csv('Data/Dataset16_Simulation4/Trustworthiness_box_data.csv', index_col=0)
df5 = pd.read_csv('Data/Dataset17_Simulation5/Trustworthiness_box_data.csv', index_col=0)

# 添加 dataset 标签
df1["dataset"] = "Simulation 1"
df2["dataset"] = "Simulation 2"
df3["dataset"] = "Simulation 3"
df4["dataset"] = "Simulation 4"
df5["dataset"] = "Simulation 5"

# 合并数据
df_all = pd.concat([df1, df2, df3, df4, df5], ignore_index=True)
df_all.loc[df_all['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'
# 重组数据结构
section_order = ['Simulation 1', 'Simulation 2', 'Simulation 3', 'Simulation 4']
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
plt.figure(figsize=(9, 3))
sns.set_style("whitegrid")

ax = sns.boxplot(
    x="dataset",
    y="r",
    hue="method",
    data=df_all,
    palette=algorithm_colors,
    order=section_order,
    hue_order=algorithms,
    width=0.7,
    dodge=True,
    flierprops=dict(marker='o', markersize=2)
)

plt.xlabel('', fontsize=0)
plt.ylabel('Cluster correlation', fontsize=9)

# 颜色覆盖增强（解决seaborn颜色透明度问题）‌:ml-citation{ref="2,4" data="citationList"}
for i, box in enumerate(ax.artists):
    hue_index = i % len(algorithms)
    alg = algorithms[hue_index]
    box.set_facecolor(algorithm_colors[alg])
    box.set_edgecolor(algorithm_colors[alg])
    box.set_alpha(0.9)  # 增强颜色饱和度
# 移除原图例
ax.legend_.remove()
plt.savefig('Results/Simulation spatial structure preservation box plot.pdf', dpi=300)


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
plt.savefig('Results/Simulation spatial structure preservation evaluation legend.pdf', dpi=300)


df1 = pd.read_csv('Data/Dataset13_Simulation1/Trustworthiness_data.csv', index_col=0)
df2 = pd.read_csv('Data/Dataset14_Simulation2/Trustworthiness_data.csv', index_col=0)
df3 = pd.read_csv('Data/Dataset15_Simulation3/Trustworthiness_data.csv', index_col=0)
df4 = pd.read_csv('Data/Dataset16_Simulation4/Trustworthiness_data.csv', index_col=0)
df5 = pd.read_csv('Data/Dataset17_Simulation5/Trustworthiness_data.csv', index_col=0)

# 添加 dataset 标签
df1["dataset"] = "Simulation 1"
df2["dataset"] = "Simulation 2"
df3["dataset"] = "Simulation 3"
df4["dataset"] = "Simulation 4"
df5["dataset"] = "Simulation 5"

# 合并数据
df_all = pd.concat([df1, df2, df3, df4, df5], ignore_index=True)
df_all.loc[df_all['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'


import matplotlib.pyplot as plt
import seaborn as sns
plt.figure(figsize=(9, 3))
sns.set_style("whitegrid")

ax = sns.barplot(
    x="dataset",
    y="value",
    hue="method",
    data=df_all,
    palette=algorithm_colors,
    order=section_order,
    hue_order=algorithms,
    width=0.7,
    dodge=True,
)

plt.xlabel('Section', fontsize=0)
plt.ylabel('Global correlation', fontsize=9)

# 颜色覆盖增强（解决seaborn颜色透明度问题）‌:ml-citation{ref="2,4" data="citationList"}
for i, box in enumerate(ax.artists):
    hue_index = i % len(algorithms)
    alg = algorithms[hue_index]
    box.set_facecolor(algorithm_colors[alg])
    box.set_edgecolor(algorithm_colors[alg])
    box.set_alpha(0.9)  # 增强颜色饱和度
# 移除原图例
ax.legend_.remove()
plt.savefig('Results/Simulation global correlation.pdf', dpi=300)



df1 = pd.read_csv('Data/Dataset13_Simulation1/morans_data.csv', index_col=0)
df2 = pd.read_csv('Data/Dataset14_Simulation2/morans_data.csv', index_col=0)
df3 = pd.read_csv('Data/Dataset15_Simulation3/morans_data.csv', index_col=0)
df4 = pd.read_csv('Data/Dataset16_Simulation4/morans_data.csv', index_col=0)
df5 = pd.read_csv('Data/Dataset17_Simulation5/morans_data.csv', index_col=0)

# 添加 dataset 标签
df1["dataset"] = "Simulation 1"
df2["dataset"] = "Simulation 2"
df3["dataset"] = "Simulation 3"
df4["dataset"] = "Simulation 4"
df5["dataset"] = "Simulation 5"

# 合并数据
df_all = pd.concat([df1, df2, df3, df4, df5], ignore_index=True)
df_all.loc[df_all['method'] == 'SpatialDDM', 'method'] = 'SpaDDM'

import matplotlib.pyplot as plt
import seaborn as sns
plt.figure(figsize=(9, 3))
sns.set_style("whitegrid")

ax = sns.barplot(
    x="dataset",
    y="I",
    hue="method",
    data=df_all,
    palette=algorithm_colors,
    order=section_order,
    hue_order=algorithms,
    width=0.7,
    dodge=True,
)

plt.xlabel('Section', fontsize=0)
plt.ylabel("Moran's I", fontsize=9)

# 颜色覆盖增强（解决seaborn颜色透明度问题）‌:ml-citation{ref="2,4" data="citationList"}
for i, box in enumerate(ax.artists):
    hue_index = i % len(algorithms)
    alg = algorithms[hue_index]
    box.set_facecolor(algorithm_colors[alg])
    box.set_edgecolor(algorithm_colors[alg])
    box.set_alpha(0.9)  # 增强颜色饱和度
# 移除原图例
ax.legend_.remove()
plt.savefig('Results/Simulation morans plot.pdf', dpi=300)