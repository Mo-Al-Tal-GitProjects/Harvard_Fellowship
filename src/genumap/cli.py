from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from genumap.config import load_config
from genumap.pipeline import run_study

app = typer.Typer(
    name="genumap",
    help="Run the GENUMAP UMAP sensitivity study.",
    no_args_is_help=True,
)


@app.command()
def check(
    config: Annotated[
        Path,
        typer.Option("--config", "-c", exists=True, dir_okay=False, readable=True),
    ] = Path("configs/study.yaml"),
) -> None:
    """Validate a study configuration without downloading data or running analyses."""
    study = load_config(config)
    grid_size = (
        len(study.analysis.pca_components)
        * len(study.analysis.umap_neighbors)
        * len(study.analysis.umap_min_dist)
        * len(study.analysis.seeds)
    )
    typer.echo(f"Configuration valid: {study.study_name}")
    typer.echo(f"UMAP fits: {grid_size}")
    typer.echo(f"Maximum retained variants: {study.data.max_variants}")
    typer.echo(f"Configuration hash: {study.stable_hash()}")


@app.command()
def run(
    config: Annotated[
        Path,
        typer.Option("--config", "-c", exists=True, dir_okay=False, readable=True),
    ] = Path("configs/study.yaml"),
    synthetic: Annotated[
        bool,
        typer.Option(help="Use a deterministic synthetic fixture instead of human genomic data."),
    ] = False,
    quick: Annotated[
        bool,
        typer.Option(help="Use the reduced verification parameter grid."),
    ] = False,
) -> None:
    """Execute the study end to end and write a versioned artifact directory."""
    study = load_config(config)
    typer.echo("Running GENUMAP study...")
    output = run_study(study, synthetic=synthetic, quick=quick)
    typer.echo(f"Completed: {output}")
