"""
Convert a RemoTAP-SPEX AOD file to a 'universal CSV' for conversion to obs_sql.
The script is using the IMERG land/sea mask to mask out land retrievals over ocean and vice versa.
"""

import os
import re
from pathlib import Path

import click
import numpy as np
import pandas as pd
import xarray as xr


def parse_spex_date(utc_date: np.ndarray, fraction_of_day: np.ndarray) -> np.ndarray:
    """
    SPEX dates are provided in a UTC date as an integer and a float representing the fraction of day.
    This function converts them to a pandas timestamp.
    """

    utcdate_str = np.array([str(int(date)) for date in utc_date])
    utcdate_str[utcdate_str == '29990101'] = np.nan
    dates = pd.to_datetime(utcdate_str, format="%Y%m%d")
    dates += pd.to_timedelta(fraction_of_day, unit='D')

    return dates

def spex_aod_to_df(spex_file_path: Path) -> pd.DataFrame:
    """
    Read all SPEX AOD observations into a DataFrame.
    """

    # geolocation arrays & date parsing
    ds_geolocation = xr.open_dataset(spex_file_path, group="geolocation_data")

    latitude = ds_geolocation["latitude"].values.flatten()
    longitude = ds_geolocation["longitude"].values.flatten()
    timestamp = parse_spex_date(
        ds_geolocation["utc_date"].values.flatten(),
        ds_geolocation["fracday"].values.flatten(),
    )

    # AOD & uncertainty
    ds_geophysical = xr.open_dataset(spex_file_path, group="geophysical_data")
    aod550 = ds_geophysical["aot550"].values.flatten()
    aod550_uncertainty = ds_geophysical["aot550_uncertainty"].values.flatten()

    df = pd.DataFrame(
        {
            "latitude": latitude,
            "longitude": (longitude + 360) % 360, # ensure longitudes are in [0, 360) range
            "timestamp": timestamp,
            "aod550": aod550,
            "aod550_uncertainty": aod550_uncertainty,
        }
    )

    # Split date into components for easier parsing in FORTRAN
    df = df.sort_values(by="timestamp")
    df["year"] = df["timestamp"].dt.year
    df["month"] = df["timestamp"].dt.month
    df["day"] = df["timestamp"].dt.day
    df["hour"] = df["timestamp"].dt.hour
    df["minute"] = df["timestamp"].dt.minute
    df["second"] = df["timestamp"].dt.second
    df = df.drop(columns=["timestamp"])

    # Sort columns by putting geolocation first
    df = df[
        [
            "latitude",
            "longitude",
            "year",
            "month",
            "day",
            "hour",
            "minute",
            "second",
            "aod550",
            "aod550_uncertainty",
        ]
    ]

    # Drop nan values (e.g. no AOD retrieval or a bad date value)
    df = df.dropna()

    return df


@click.command()
@click.argument("spex_file_path", type=click.Path(exists=True, path_type=Path))
@click.argument("output_path", type=click.Path(path_type=Path))
def main(spex_file_path: Path, output_path: Path):
    """
    Convert a RemoTAP-SPEX AOD file to a 'universal CSV' for conversion to obs_sql.
    The script is using the IMERG land/sea mask to mask out land retrievals over ocean and vice versa.

    The IMERG land/sea mask must be either in the same directory called 'IMERG_land_sea_mask.nc'
    or defined in the environment variable 'IMERG_LAND_SEA_MASK'. The env. variable takes precedence.
    """

    # Locate the IMERG land/sea mask
    imerg_land_sea_mask_path = Path(
        os.environ.get("IMERG_LAND_SEA_MASK", "IMERG_land_sea_mask.nc")
    )
    if not imerg_land_sea_mask_path.exists():
        raise FileNotFoundError(
            f"Could not find IMERG land/sea mask at {imerg_land_sea_mask_path}"
        )
    imerg = xr.open_dataset(imerg_land_sea_mask_path)
    imerg['is_ocean'] = imerg['landseamask'] > 60

    # Check if the file is a valid SPEX file
    with xr.open_dataset(spex_file_path) as root_ds:
        product_name = root_ds.attrs.get("product_name", "unknown")
    if not re.match(r'PACE_SPEXONE.*AER_(?:LAND|OCEAN)_REMOTAP\.nc', product_name):
        raise ValueError(
            f"File {spex_file_path} is not a valid SPEX file. "
            r"Expected product name to match 'PACE_SPEXONE.*AER_(?:LAND|OCEAN)_REMOTAP\.nc', "
            f"but got '{product_name}'"
        )

    # Load each file into a big dataframe
    df = spex_aod_to_df(spex_file_path)
    df['is_ocean_retrieval'] = "OCEAN" in product_name
    length = len(df)
    print(f"Got {length} observations")

    if df.empty:
        print("No observations found. Exiting.")
        return

    # Filter out observations that don't match the land/sea mask (i.e. land retrievals over ocean and vice versa)
    df['imerg_is_ocean'] = df.apply(
        lambda row: imerg['is_ocean'].sel(lon=row['longitude'], lat=row['latitude'], method='nearest').item(),
        axis=1
    )
    df = df.loc[~(df["is_ocean_retrieval"]) ^ (df["imerg_is_ocean"])] # NOT XOR
    print(f"Filtered to {len(df)} observations matching the land/sea mask")

    # Fix columns to be as expected by the FORTRAN code
    df['obs_type'] = 'AIRSENSE_AOD'
    df['vert'] = 0.0
    df['obs_value'] = df['aod550']
    df['obs_error'] = df['aod550_uncertainty']
    df['obs_meta'] = ''
    df = df[[
        "obs_type",
        "longitude",
        "latitude",
        "vert",
        "year",
        "month",
        "day",
        "hour",
        "minute",
        "second",
        "obs_value",
        "obs_error",
        "obs_meta",
    ]]

    print(f"Writing observations to {output_path}")
    df.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
