import typer
import pandas as pd
import numpy as np
import click

app = typer.Typer(help="TSS clustering utilities")

@app.command()
def main(
    version: bool = typer.Option(False, "--version", help="Show version and exit."),
    input_file: str = typer.Option(None, "-i", "--input", help="Input TSS raw signal table (tab-delimited, must include chr, pos, strand, and sample column)"),
    output_file: str = typer.Option(None, "-o", "--output", help="Output clustered result table"),
    sample_col: str = typer.Option(None, "-s", "--sample", help="Sample column name to use for clustering"),
    peak_distance: int = typer.Option(100, "--peak-distance", help="Minimum distance between peaks"),
    local_threshold: float = typer.Option(0.02, "--local-threshold", help="Local filtering threshold (TPM)"),
    extension_distance: int = typer.Option(30, "--extension-distance", help="Cluster boundary extension distance"),
    cluster_threshold: float = typer.Option(1, "--cluster-threshold", help="Cluster filtering threshold (TPM)")
):
    """
    Cluster TSSs to infer core promoters (Python version of cluster_by_peak)
    """
    __version__ = '0.2.0'
    if version:
        print(f'clustering.py version {__version__}')
        return
    if not (input_file and output_file and sample_col):
        print('Error: -i/--input, -o/--output, and -s/--sample are required unless --version is specified.')
        return
    df = pd.read_csv(input_file, sep='\t')
    if sample_col not in df.columns:
        raise ValueError(f"Sample column {sample_col} not found in input file")
    df[sample_col] = df[sample_col].astype(float)
    total_tags = df[sample_col].sum()
    if total_tags == 0:
        raise ValueError("Total tags for the selected sample column is zero.")
    df['TPM'] = df[sample_col] / total_tags * 1e6
    df = df[df['TPM'] > 0].copy()

    header_written = False
    for (chr_, strand), group in df.groupby(['chr', 'strand']):
        group = group.sort_values('pos').reset_index(drop=True)
        clusters = cluster_by_peak(group, peak_distance, local_threshold, extension_distance, sample_col, cluster_threshold)
        if clusters is not None and not clusters.empty:
            clusters['chr'] = chr_
            clusters['strand'] = strand
            out_df = clusters[['cluster', 'chr', 'start', 'end', 'strand', 'dominant_tss', 'tags', 'tags.dominant_tss', 'tags_TPM', 'tags_TPM.dominant_tss', 'q_0.1', 'q_0.9', 'interquantile_width']]
            out_df.to_csv(output_file, sep='\t', index=False, mode='a', header=not header_written)
            header_written = True

def cluster_by_peak(df, peak_distance, local_threshold, extension_distance, sample_col, cluster_threshold):
    peak_id = np.zeros(len(df), dtype=int)
    for i in range(len(df)):
        pos = df.iloc[i]['pos']
        tags = df.iloc[i]['TPM']
        mask = (df['pos'] > pos - peak_distance) & (df['pos'] < pos + peak_distance)
        local_max = df.loc[mask, 'TPM'].max() if mask.any() else 0
        if tags == local_max:
            peak_id[i] = i + 1
    df = df.copy()
    df['peak'] = peak_id
    df['ID'] = np.arange(1, len(df) + 1)
    keep_idx = set(df.index)
    for i in np.where(peak_id > 0)[0]:
        peak_tpm = df.iloc[i]['TPM']
        peak_pos = df.iloc[i]['pos']
        if df.iloc[i]['strand'] == '+':
            mask = (df['pos'] >= peak_pos) & (df['pos'] <= peak_pos + peak_distance)
        else:
            mask = (df['pos'] >= peak_pos - peak_distance) & (df['pos'] <= peak_pos)
        sub = df[mask]
        drop = sub.index[sub['TPM'] < peak_tpm * local_threshold]
        keep_idx -= set(drop)
    df = df.loc[list(keep_idx)].sort_values('pos').reset_index(drop=True)
    pos = df['pos'].values
    forward = np.zeros(len(df), dtype=int)
    reverse = np.zeros(len(df), dtype=int)
    forward[:-1] = (pos[1:] < pos[:-1] + extension_distance).astype(int)
    reverse[1:] = (pos[:-1] > pos[1:] - extension_distance).astype(int)
    df['forward'] = forward
    df['reverse'] = reverse
    rleid = (df['peak'].ne(0) | df['forward'].ne(0) | df['reverse'].ne(0)).cumsum()
    df['rleid'] = rleid
    clusters = []
    for _, sub in df.groupby('rleid'):
        if sub['peak'].max() == 0:
            continue
        start = sub['pos'].min()
        end = sub['pos'].max()
        cluster_data = sub.copy()
        tags_sum = cluster_data['TPM'].sum()
        tags_sum_raw = cluster_data[sample_col].sum()
        dominant_idx = cluster_data['TPM'].idxmax()
        dominant_tss = cluster_data.loc[dominant_idx, 'pos']
        tags_dom = cluster_data.loc[dominant_idx, sample_col]
        tags_dom_tpm = cluster_data.loc[dominant_idx, 'TPM']
        cumsum = cluster_data['TPM'].cumsum()
        q1 = cluster_data.loc[cumsum > 0.1 * tags_sum, 'pos'].min()
        q9 = cluster_data.loc[cumsum > 0.9 * tags_sum, 'pos'].min()
        interquantile_width = q9 - q1 + 1 if pd.notnull(q1) and pd.notnull(q9) else 0
        clusters.append({
            'cluster': len(clusters) + 1,
            'start': start,
            'end': end,
            'dominant_tss': dominant_tss,
            'tags': tags_sum_raw,
            'tags.dominant_tss': tags_dom,
            'tags_TPM': tags_sum,
            'tags_TPM.dominant_tss': tags_dom_tpm,
            'q_0.1': q1,
            'q_0.9': q9,
            'interquantile_width': interquantile_width
        })
    clusters_df = pd.DataFrame(clusters)
    clusters_df = clusters_df[clusters_df['tags_TPM'] > cluster_threshold].reset_index(drop=True)
    return clusters_df

if __name__ == '__main__':
    app() 