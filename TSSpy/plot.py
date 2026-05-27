import typer
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

app = typer.Typer(help="Plotting utilities for TSSpy")

def correlation_color(r, cmap='coolwarm'):
    # r: -1~1, 映射到0~1
    norm = (r + 1) / 2
    return plt.get_cmap(cmap)(norm)

@app.command()
def correlation(
    tss_table: str = typer.Option(..., "-i", "--input", help="Input TSS table (tab-delimited)"),
    output: str = typer.Option("correlation_plot.png", "-o", "--output", help="Output image file (png/pdf)")
):
    """
    Plot pairwise correlation and scatter plots between samples in a TSS table.
    Each subplot's color reflects the correlation coefficient (red=high, blue=low).
    """
    typer.echo("[1/5] Reading input table...")
    df = pd.read_csv(tss_table, sep='\t')
    sample_cols = df.columns[3:]
    typer.echo(f"[2/5] Found {len(sample_cols)} samples, {len(df)} rows.")
    data = df[sample_cols]
    typer.echo("[3/5] Filtering rows where all samples are zero...")
    data = data.loc[~(data == 0).all(axis=1)]
    typer.echo(f"[4/5] Plotting pairwise scatter and correlation (color=correlation, red=high, blue=low)...")

    # 创建 PairGrid
    g = sns.PairGrid(data)
    for i, var1 in enumerate(sample_cols):
        for j, var2 in enumerate(sample_cols):
            if i > j:
                r = data[var1].corr(data[var2])
                color = correlation_color(r, cmap='coolwarm')
                ax = g.axes[j, i]
                ax.scatter(data[var1], data[var2], color=color, s=10, alpha=0.7)
                ax.annotate(f"r={r:.2f}", (0.5, 0.9), xycoords="axes fraction", ha="center", fontsize=10, color='black')
            elif i == j:
                ax = g.axes[j, i]
                ax.hist(data[var1], bins=20, color='gray', alpha=0.7)
            else:
                g.axes[j, i].axis('off')
    plt.suptitle("Pairwise Correlation of TSS Samples\n(color=correlation, red=high, blue=low)", y=1.02)
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    typer.echo(f"[5/5] Correlation plot saved to {output}")

if __name__ == "__main__":
    app() 