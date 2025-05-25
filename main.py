import typer

from TSSpy import tss_calling
from TSSpy import clustering
from TSSpy import gene_assign
from TSSpy import plot
from TSSpy.correlation import correlation
from TSSpy import bigwig

app = typer.Typer(help="TSSpy: Python CLI for TSS analysis (main command)")

app.add_typer(tss_calling.app, name="tssCalling")
app.add_typer(clustering.app, name="clustering")
app.add_typer(gene_assign.app, name="geneAssign")
app.add_typer(plot.app, name="plot")
app.add_typer(bigwig.app, name="bigwig")
app.command()(correlation)

if __name__ == "__main__":
    app() 