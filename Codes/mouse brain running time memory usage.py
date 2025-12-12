import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


data_time_spatialglue = {
    "ATAC": [63.51, 60.99, 62.10, 60.94, 60.85],
    "H3K4me3": [72.49, 71.98, 71.74, 72.07, 73.50],
    "H3K27ac": [67.34, 69.13, 66.68, 65.63, 65.72],
    "H3K27me3": [74.21, 75.50, 74.51, 73.07, 74.44],
}

data_time_spaddm = {
    "ATAC": [122.38, 123.04, 126.29, 124.48, 125.72],
    "H3K4me3": [125.85, 124.70, 124.64, 124.11, 125.53],
    "H3K27ac": [126.56, 127.18, 125.96, 130.61, 127.76],
    "H3K27me3": [130.95, 130.89, 132.89, 130.74, 132.19],
}


data_memory_spatialglue = {
    "ATAC": [4778.75, 4781.16, 4777.95, 4781.97, 4779.92],
    "H3K4me3": [4502.47, 4501.59, 4503.44, 4502.77, 4500.85],
    "H3K27ac": [4636.59, 4637.80, 4635.55, 4637.45, 4635.40],
    "H3K27me3": [4643.97, 4642.57, 4643.94, 4642.88, 4644.40],
}

data_memory_spaddm = {
    "ATAC": [4665.75, 4666.89, 4667.38, 4669.81, 4665.78],
    "H3K4me3": [3633.86, 3631.65, 3632.47, 3631.61, 3632.91],
    "H3K27ac": [4491.32, 4491.48, 4490.04, 4488.20, 4489.42],
    "H3K27me3": [4444.19, 3645.23, 3644.59, 3645.60, 3645.36],
}
datasets = ["ATAC", "H3K4me3", "H3K27ac", "H3K27me3"]


# 转为长格式 DataFrame
df_list = []

for dataset, values in data_time_spatialglue.items():
    df_list.append(pd.DataFrame({
        "Dataset": dataset,
        "Time": values,
        "Method": "SpatialGLUE"
    }))

for dataset, values in data_time_spaddm.items():
    df_list.append(pd.DataFrame({
        "Dataset": dataset,
        "Time": values,
        "Method": "SpaDDM"
    }))

df = pd.concat(df_list, ignore_index=True)

# ====== 开始绘制 ======
plt.figure(figsize=(5, 3))
# 箱线图（不同模型不同颜色）
ax = sns.boxplot(
    data=df,
    x="Dataset",
    y="Time",
    hue="Method",
    width=0.6,
    fliersize=0,  # 不显示默认异常点
    flierprops=dict(marker='o', markersize=0.5)
)

# 用 scatter / strip 点覆盖上去（颜色自动跟 hue 匹配）
ax = sns.stripplot(
    data=df,
    x="Dataset",
    y="Time",
    hue="Method",
    dodge=True,        # 分开两个模型的点
    jitter=True,       # 稍微散开
    alpha=0.1,
    legend=False
)
# ====== 在数据集之间添加竖直虚线 ======
for i in range(len(datasets) - 1):
    plt.axvline(x=i + 0.5, color="gray", linestyle="--", alpha=0.6)

# 避免重复图例
plt.legend(title="Method", bbox_to_anchor=(1.5, 1), loc="upper right", fontsize=8)
# ax.legend_.remove()
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.ylabel("Time (s)", fontsize=8)
plt.xlabel("Datasets", fontsize=8)
plt.title("Runtime comparison (SpatialGlue vs SpaDDM)", fontsize=8)
plt.tight_layout()
plt.show()
plt.savefig('Runtime comparison (SpatialGlue vs SpaDDM).pdf', dpi=300)






# 转为长格式 DataFrame
df_list = []

for dataset, values in data_memory_spatialglue.items():
    df_list.append(pd.DataFrame({
        "Dataset": dataset,
        "Time": values,
        "Method": "SpatialGLUE"
    }))

for dataset, values in data_memory_spaddm.items():
    df_list.append(pd.DataFrame({
        "Dataset": dataset,
        "Time": values,
        "Method": "SpaDDM"
    }))

df = pd.concat(df_list, ignore_index=True)

# ====== 开始绘制 ======
plt.figure(figsize=(5, 3))
# 箱线图（不同模型不同颜色）
ax = sns.boxplot(
    data=df,
    x="Dataset",
    y="Time",
    hue="Method",
    width=0.6,
    fliersize=0,  # 不显示默认异常点
    flierprops=dict(marker='o', markersize=0.5)
)

# 用 scatter / strip 点覆盖上去（颜色自动跟 hue 匹配）
ax = sns.stripplot(
    data=df,
    x="Dataset",
    y="Time",
    hue="Method",
    dodge=True,        # 分开两个模型的点
    jitter=True,       # 稍微散开
    alpha=0.1,
    legend=False
)
# ====== 在数据集之间添加竖直虚线 ======
for i in range(len(datasets) - 1):
    plt.axvline(x=i + 0.5, color="gray", linestyle="--", alpha=0.6)

# 避免重复图例
plt.legend(title="Method", bbox_to_anchor=(1.5, 1), loc="upper right", fontsize=8)
# ax.legend_.remove()
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.ylabel("Memory Usage (MB)", fontsize=8)
plt.xlabel("Datasets", fontsize=8)
plt.title("Memory usage comparison (SpatialGlue vs SpaDDM)", fontsize=8)
plt.tight_layout()
plt.show()
plt.savefig('Memory usage comparison (SpatialGlue vs SpaDDM).pdf', dpi=300)
