# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polars",
#   "pandas",
#   "pyarrow",
#   "numpy",
# ]
# ///
"""Preprocess raw UKB covariates parquet to a processed parquet file.

Usage:
    uv run preprocess_covariates.py \
        --covariates-path    <path/to/covariates.parquet> \
        --mapping-path       <path/to/covariate_mapping.parquet> \
        --output-dir         <path/to/output_dir> \
        [--samples-path      <path/to/samples.csv>]

Outputs (written to --output-dir):
    covariates_preprocessed.parquet  – one row per sample, columns: sample + covariate features
"""

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Union, Protocol

import numpy as np
import pandas as pd
import polars as pl

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


class CovariateExtractor(Protocol):
    def __call__(self, cov_column: Union[int, tuple[int, ...]], cov_label: str, cov_df: pd.DataFrame) -> pd.DataFrame: ...


@dataclass
class Covariate:
    code: Union[int, tuple[int, ...]]
    label: str
    extract_function: CovariateExtractor


def _extract_categorical(cov_column: int, cov_label: str, cov_df: pd.DataFrame) -> pd.DataFrame:
    return pd.get_dummies(cov_df[str(cov_column)], prefix=cov_label)


def _extract_raw(cov_column: int, cov_label: str, cov_df: pd.DataFrame) -> pd.DataFrame:
    return cov_df[str(cov_column)].rename(cov_label)


def _extract_array(cov_column: int, cov_label: str, cov_df: pd.DataFrame) -> pd.DataFrame:
    array_size = len(cov_df[str(cov_column)].iloc[0])
    return pd.DataFrame(
        cov_df[str(cov_column)].tolist(),
        columns=[f"{cov_label}_{i}" for i in range(array_size)],
        index=cov_df.index,
    )


_CATEGORICAL = [
    (54, 'UK_Biobank_assessment_centre'),
    (1239, 'Current_tobacco_smoking'),
    (1249, 'Past_tobacco_smoking'),
    (1259, 'Smoking_or_smokers_in_household'),
    (1329, 'Oily_fish_intake'),
    (1339, 'Non-oily_fish_intake'),
    (1349, 'Processed_meat_intake'),
    (1359, 'Poultry_intake'),
    (1369, 'Beef_intake'),
    (1379, 'Lamb_or_mutton_intake'),
    (1389, 'Pork_intake'),
    (1408, 'Cheese_intake'),
    (1418, 'Milk_type_used'),
    (1428, 'Spread_type'),
    (1448, 'Bread_type'),
    (1468, 'Cereal_type'),
    (1478, 'Salt_added_to_food'),
    (1518, 'Hot_drink_temperature'),
    (1538, 'Major_dietary_changes_in_the_last_5_years'),
    (1548, 'Variation_in_diet'),
    (1558, 'Alcohol_intake_frequency'),
    (1628, 'Alcohol_intake_versus_10_years_previously'),
    (1677, 'Breastfed_as_a_baby'),
    (3089, 'Caffeine_drink_within_last_hour'),
    (21000, 'ethnicity'),
    (22021, 'Genetic_kinship_to_other_participants'),
    (24014, 'Close_to_major_road'),
    (2724, 'Had_menopause'),
    (2784, 'taken_contraceptive_pill'),
    (2814, 'used_hormone_replacement_therapy'),
    (2834, 'both_ovaries_removed'),
    (3140, 'pregnant'),
    (3591, 'had_womb_removed'),
    (53, 'Date_of_attending_assessment_centre'),
]

_RAW = [
    (22001, 'sex'),
    (21003, 'age'),
    (22189, 'Townsend_deprivation_index_at_recruitment'),
    (24003, 'Nitrogen_dioxide_air_pollution_2010'),
    (24004, 'Nitrogen_oxides_air_pollution_2010'),
    (24010, 'Inverse_distance_to_the_nearest_road'),
    (24012, 'Inverse_distance_to_the_nearest_major_road'),
    (24015, 'Sum_of_road_length_of_major_roads_within_100m'),
    (24016, 'Nitrogen_dioxide_air_pollution_2005'),
    (24017, 'Nitrogen_dioxide_air_pollution_2006'),
    (24018, 'Nitrogen_dioxide_air_pollution_2007'),
    (24019, 'Particulate_matter_air_pollution_pm10_2007'),
    (24020, 'Average_daytime_sound_level_of_noise_pollution'),
    (24021, 'Average_evening_sound_level_of_noise_pollution'),
    (24022, 'Average_night_time_sound_level_of_noise_pollution'),
    (24023, 'Average_16_hour_sound_level_of_noise_pollution'),
    (24024, 'Average_24_hour_sound_level_of_noise_pollution'),
    (74, 'Fasting_time'),
    (1289, 'Cooked_vegetable_intake'),
    (1299, 'Salad_or_raw_vegetable_intake'),
    (1309, 'Fresh_fruit_intake'),
    (1319, 'Dried_fruit_intake'),
    (1438, 'Bread_intake'),
    (1458, 'Cereal_intake'),
    (1488, 'Tea_intake'),
    (1498, 'Coffee_intake'),
    (1528, 'Water_intake'),
    (24009, 'Traffic_intensity_on_the_nearest_road'),
    (24011, 'Traffic_intensity_on_the_nearest_major_road'),
    (24013, 'Total_traffic_load_on_major_roads'),
]

_ARRAY = [
    (22009, 'Genetic_PCs'),
]

_CUSTOM = [
    (21003, 'age^2',
     lambda col, label, d: pd.Series(d[str(col)] ** 2).to_frame(name=label)),
    ((21003, 22001), 'ageXsex',
     lambda cols, label, d: pd.Series(d[str(cols[0])] * d[str(cols[1])]).to_frame(name=label)),
    ((21003, 22001), 'age^2Xsex',
     lambda cols, label, d: pd.Series(d[str(cols[0])] * 2 * d[str(cols[1])]).to_frame(name=label)),
]

COVARIATES = [
    *[Covariate(code, label, _extract_categorical) for code, label in _CATEGORICAL],
    *[Covariate(code, label, _extract_raw) for code, label in _RAW],
    *[Covariate(code, label, _extract_array) for code, label in _ARRAY],
    *[Covariate(code, label, fn) for code, label, fn in _CUSTOM],
]


def transform_covariates(df: pd.DataFrame, covariates: list[Covariate]) -> pd.DataFrame:
    assert 'eid' not in df.columns
    return pd.concat(
        [cov.extract_function(cov.code, cov.label, df).astype(np.float32) for cov in covariates],
        axis=1,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--covariates-path", required=True, type=Path)
    parser.add_argument("--mapping-path", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--samples-path", type=Path, default=None)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    mapping_df = pl.read_parquet(args.mapping_path)
    mapping_df = (
        mapping_df
        .with_columns(instance=pl.col('title').str.extract(r' \| Instance (\d+)').cast(pl.Int32))
        .with_columns(title=pl.col('title').str.split(' | ').list.get(0))
        .group_by('ukb-code', 'title', 'instance')
        .agg(pl.col('field_name'))
        .filter((pl.col('instance') == 0) | (pl.col('instance').is_null()))
    )

    fields_to_keep = [y for x in mapping_df['field_name'].to_list() for y in x]
    cov_df = pl.read_csv(args.covariates_path).select(['eid'] + fields_to_keep)

    for ukb_code, field_names in mapping_df.select(['ukb-code', 'field_name']).rows():
        if len(field_names) == 1:
            cov_df = cov_df.with_columns(**{str(ukb_code): pl.col(field_names[0])}).drop(field_names)
        else:
            cov_df = cov_df.with_columns(**{str(ukb_code): pl.concat_list(field_names)}).drop(field_names)

    cov_df = cov_df.rename({'eid': 'sample'})

    cov_pd = transform_covariates(cov_df.to_pandas().set_index('sample'), covariates=COVARIATES)
    cov_df = pl.from_pandas(cov_pd.reset_index())
    cov_samples = cov_df['sample'].to_list()

    if args.samples_path:
        logger.info("Subsetting to samples in %s", args.samples_path)
        samples = list(
            pl.read_csv(args.samples_path, has_header=True)
            .select(pl.nth(0).cast(pl.String))
            .to_series()
        )
        missing = set(samples) - set(cov_samples)
        if missing:
            logger.warning("%d samples not found and will be skipped: %s", len(missing), missing)
            samples = [s for s in samples if s not in missing]
        cov_df = pl.DataFrame({'sample': samples}).join(cov_df, on='sample', how='left')

    cov_out = args.output_dir / "covariates_preprocessed.parquet"
    cov_df.write_parquet(cov_out)
    logger.info("Saved: %s", cov_out)


if __name__ == "__main__":
    main()
