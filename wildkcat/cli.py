import typer

from dotenv import load_dotenv
from wildkcat import run_extraction, run_retrieval, run_prediction_part1, run_prediction_part2, generate_summary_report


load_dotenv()

__version__ = "0.3.1"

app = typer.Typer(help="WILDkCAT CLI - Extract, Retrieve and Predict kcat values for a metabolic model.")

def version_callback(value: bool):
    """Show the application's version and exit."""
    if value:
        typer.echo(f"WILDkCAT version: {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    ctx: typer.Context,
    version: bool = typer.Option(
        False, "--version", callback=version_callback, is_eager=True,
        help="Show the application's version."
    )
):
    """Common CLI callback."""
    pass


@app.command()
def extraction(
    model_path: str, 
    output_folder: str, 
    report: bool = True
):
    """
    Extracts kcat-related data from a metabolic model and generates output files and an optional HTML report.

    Parameters:

        model_path (str): Path to the metabolic model file (JSON, MATLAB, or SBML format).

        output_folder (str): Path to the output folder where the files will be saved.

        report (bool, optional): Whether to generate an HTML report (default: True).
    """
    run_extraction(model_path=model_path, output_folder=output_folder, report=report)
    typer.echo(f"Extraction finished. Output saved at {output_folder}/kcat.tsv")


@app.command()
def retrieval(
    output_folder: str,
    organism: str,
    temperature_range: tuple[float, float],
    ph_range: tuple[float, float],
    database: str = 'both',
    report: bool = True
):
    """
    Retrieves closests kcat values from specified databases for entries in a kcat file, applies filtering criteria, 
    and saves the results to an output file.
    
    Parameters:

        output_folder (str): Path to the output folder where the files are and will be saved.

        organism (str): Organism scientific name (e.g. "Escherichia coli", "Homo sapiens").

        temperature_range (tuple): Acceptable temperature range for filtering (min, max).

        ph_range (tuple): Acceptable pH range for filtering (min, max).

        database (str, optional): Specifies which database(s) to query for kcat values. Options are 'both' (default), 'brenda', or 'sabio_rk'.

        report (bool, optional): Whether to generate an HTML report using the retrieved data (default: True).       
    """
    run_retrieval(
        output_folder=output_folder,
        organism=organism,
        temperature_range=temperature_range,
        pH_range=ph_range,
        database=database,
        report=report
    )
    typer.echo(f"Retrieval finished. Output saved at {output_folder}/kcat_retrieved.tsv")


@app.command()
def prediction_part1(
    output_folder: str, 
    limit_penalty_score: int,
    report: bool = True
):
    """
    Processes kcat data file to generate input files for CataPro prediction.
    Optionally, it can produce a summary report of the processed data.

    Parameters:
        output_folder (str): Path to the output folder where the files are and will be saved.

        limit_penalty_score (int): Threshold for filtering entries based on penalty score.

        report (bool, optional): Whether to generate a report using the retrieved data (default: True). 
    """
    run_prediction_part1(
        output_folder=output_folder,
        limit_penalty_score=limit_penalty_score, 
        report=report
    )
    typer.echo(f"Prediction Part 1 finished. Output saved at {output_folder}/machine_learning/catapro_input.csv")


@app.command()
def prediction_part2(
    output_folder: str,
    catapro_predictions_path: str,
    limit_penalty_score: int
):
    """
    Runs the second part of the kcat prediction pipeline by integrating Catapro predictions,
    mapping substrates to SMILES, formatting the output, and optionally generating a report.
    
    Parameters:
        output_folder (str): Path to the output folder where the files are and will be saved.

        catapro_predictions_path (str): Path to the Catapro predictions CSV file.

        limit_penalty_score (float): Threshold for taking predictions over retrieved values.

        report (bool, optional): If True, generates a report (default: True). 
    """
    run_prediction_part2(
        output_folder=output_folder,
        catapro_predictions_path=catapro_predictions_path,
        limit_penalty_score=limit_penalty_score
    )
    typer.echo(f"Prediction Part 2 finished. Output saved at {output_folder}/kcat_full.tsv")


@app.command()
def report(model_path: str, output_folder: str):
    """
    Generate a HTML report summarizing the kcat extraction, retrieval and prediction for a given model. 

    Parameters:

        model_path (str): Path to the metabolic model file (JSON, MATLAB, or SBML format).
        
        output_folder (str): Path to the output folder where the files are and will be saved.
    """
    generate_summary_report(
        model_path=model_path,
        output_folder=output_folder
    )
    typer.echo(f"Summary report generated. Output saved at {output_folder}/reports/general_report.html")


if __name__ == "__main__":
    app()
