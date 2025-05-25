import typer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde

__version__ = "0.4.0"

def correlation(
    tss_table: str = typer.Option(..., "-i", "--input", help="Input TSS table (tab-delimited)"),
    output: str = typer.Option(None, "-o", "--output", help="Output correlation matrix file (csv, optional)"),
    plot: bool = typer.Option(False, "--plot", help="If set, also plot pairwise scatter/correlation matrix (like R pairs)", is_flag=True),
    plot_file: str = typer.Option("correlation_pairs.png", "--plot-file", help="Output plot file (png/pdf)"),
    version: bool = typer.Option(False, "--version", help="Show version and exit.")
):
    """
    Calculate and output sample correlation matrix. Optionally plot pairwise scatter/correlation matrix (R pairs style).
    """
    if version:
        typer.echo(f"TSSpy correlation version {__version__}")
        raise typer.Exit()
    typer.echo("[1/4] Reading input table...")
    df = pd.read_csv(tss_table, sep='\t')
    sample_cols = df.columns[3:]
    data = df[sample_cols]
    typer.echo(f"[2/4] Found {len(sample_cols)} samples, {len(df)} rows. Calculating correlation matrix...")
    # Remove rows where all samples are zero
    data = data.loc[~(data == 0).all(axis=1)]
    # log10(x+1) for plotting
    data_log = np.log10(data + 1)
    corr = data.corr(method='pearson')
    # 输出相关系数矩阵
    print("Correlation matrix (Pearson r):")
    print(corr.round(3))
    if output:
        corr.to_csv(output)
        typer.echo(f"Correlation matrix saved to {output}")
    if not plot:
        return
    typer.echo(f"[3/4] Plotting pairwise scatter/correlation matrix (R pairs style, blue points)...")
    n = len(sample_cols)
    fig, axes = plt.subplots(n, n, figsize=(2.5*n, 2.5*n), squeeze=False)
    for ax_row in axes:
        for ax in ax_row:
            ax.set_aspect('equal', adjustable='box')
    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if i == j:
                # 对角线：只显示样品名
                ax.annotate(sample_cols[i], (0.5, 0.5), xycoords="axes fraction", ha="center", va="center", fontsize=14, color='black', fontweight='bold')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_frame_on(False)
                ax.set_xlabel("")
                ax.set_ylabel("")
            elif i > j:
                # 下三角：散点图，全部用蓝色点
                xvals = data_log.iloc[:, j]
                yvals = data_log.iloc[:, i]
                ax.scatter(xvals, yvals, s=2, color='blue', alpha=0.7)
                ax.set_xlim(0, 5)
                ax.set_ylim(0, 5)
                if j == 0:
                    ax.set_ylabel(sample_cols[i])
                else:
                    ax.set_ylabel("")
                if i == n - 1:
                    ax.set_xlabel(sample_cols[j])
                else:
                    ax.set_xlabel("")
            else:
                # 右上角：相关系数文本，背景色映射相关系数，字体大小随|r|增大
                r = corr.iloc[i, j]
                from matplotlib import cm
                blue_red_cmap = LinearSegmentedColormap.from_list('blue_red', [(0,0,1), (1,1,1), (1,0,0)], N=100)
                color = blue_red_cmap(int(np.clip((r + 1) / 2 * 99, 0, 99)) / 99)
                ax.set_facecolor(color)
                font_size = 12 + abs(r)*16
                font_size = min(font_size, 28)
                ax.annotate(f"{r:.2f}", (0.5, 0.5), xycoords="axes fraction", ha="center", va="center", fontsize=font_size, color='black', fontweight='bold')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_frame_on(False)
                ax.set_xlabel("")
                ax.set_ylabel("")
    plt.tight_layout()
    # 不加suptitle，避免标题遮挡
    plt.savefig(plot_file, dpi=150)
    typer.echo(f"[4/4] Correlation plot saved to {plot_file}") 