import typer

app = typer.Typer(help="Assign TSS clusters to genes (geneAssign)")

@app.command()
def main():
    """
    Assign TSS clusters to genes (to be implemented)
    """
    typer.echo("[geneAssign] This feature is under development.")

if __name__ == "__main__":
    app() 