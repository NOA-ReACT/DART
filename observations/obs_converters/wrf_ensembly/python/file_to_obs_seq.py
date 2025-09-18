"""
This script will convert the a set of input files to a CSV file, then call the converter
to convert the CSV file to an obs_seq file. The CSV file is deleted.
"""

import os
import subprocess
from pathlib import Path

import click

@click.command()
@click.argument("input_args", nargs=-1)
@click.argument("output_path", type=click.Path(path_type=Path))
@click.option("--conversion-script", type=click.Path(exists=True, path_type=Path))
def convert_to_csv(input_args: list[str], output_path: Path, conversion_script: Path) -> None:
    """
    Convert the input files to CSV using the given CONVERSION_SCRIPT, then call the Universal CSV converter to turn them into an obs_seq file, written at OUTPUT_PATH.

    The INPUT_ARGS are all the arguments needed to convert the input files to CSV.
    The CONVERSION_SCRIPT is the script that will convert the input files to CSV.
    The Universal CSV converter should reside in `../work` from the current directory or in
    the path indicated in the environment variable `UNIVERSAL_CSV_CONVERTER_PATH`.

    The CONVERSION_SCRIPT is called:
        python3 CONVERSION_SCRIPT INPUT_ARGS OUTPUT_PATH.csv

    The Universal CSV converter is called:
        ./universal_csv_converter OUTPUT_PATH.csv OUTPUT_PATH

    All output from the conversion script and the universal CSV converter is printed to the console.
    """

    # Check if we can locate the universal CSV converter
    universal_csv_converter_directory = Path(
        os.environ.get("UNIVERSAL_CSV_CONVERTER_PATH", "../work/")
    )
    universal_csv_converter_binary = universal_csv_converter_directory / "convert_universal_csv"
    if not universal_csv_converter_binary.exists() or not universal_csv_converter_binary.is_file():
        raise FileNotFoundError(
            f"Could not find the universal CSV converter at {universal_csv_converter_binary}"
        )

    intermediate_csv_path = output_path.with_suffix(".csv")

    command = ['python3', str(conversion_script), *input_args, str(intermediate_csv_path)]
    print(f"Running command: {' '.join(command)}")
    print('---')
    result = subprocess.run(
        command,
        cwd=conversion_script.parent,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    print(result.stdout)

    if result.returncode != 0:
        raise RuntimeError(f"Conversion script failed with return code {result.returncode}")

    command = [str(universal_csv_converter_binary), str(intermediate_csv_path.resolve()), str(output_path.resolve())]
    print(f"Running command: {' '.join(command)}")
    print('---')
    result = subprocess.run(
        command,
        cwd=universal_csv_converter_directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    print(result.stdout)

    if result.returncode != 0:
        raise RuntimeError(f"Universal CSV converter failed with return code {result.returncode}")

    # Clean up the CSV file
    # intermediate_csv_path.unlink()


if __name__ == "__main__":
    convert_to_csv()