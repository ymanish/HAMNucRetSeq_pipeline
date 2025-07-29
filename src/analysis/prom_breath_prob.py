## **5. Promoter Region Nucleosome Breathing Propensity**


from src.config.path import DATA_DIR, RESULTS_DIR
import polars as pl
from src.config.path import EXAMPLES_DATA, EXAMPLES_OUT
from pathlib import Path
import csv
from src.utils.logger_util import get_logger
import logging

def get_promoter_breathing_data(input_file: Path, id_lookup_path: Path) -> pl.LazyFrame:
    """
    Load promoter breathing data from the example dataset.
    """
    schema = {
        "id": pl.Int64,
        "subid": pl.Int64,
        "left_index": pl.Int64,
        "right_index": pl.Int64,
        "dF": pl.Float64,
        "F_enthalpy": pl.Float64,
        "F_freedna": pl.Float64
    }

    schema_id = {
        "id_original": pl.Utf8,
        "id_code": pl.Int64
    }

    lz_promdt = pl.scan_csv(input_file,
                                separator="\t",
                                  infer_schema_length=0,
                                    dtypes=schema)

    lz_promid = pl.scan_csv(id_lookup_path,
                                    separator="\t",
                                      infer_schema_length=0,
                                        dtypes=schema_id)

    return lz_promdt, lz_promid


def refine_lookup_tb(lz_prombound_id: pl.LazyFrame)-> pl.LazyFrame:
  NUM_PARTS = 5

  df_id_split = (
              lz_prombound_id
              .with_columns(
                            pl.col("id_original")
                      
                              .str.split_exact(by="|", n=NUM_PARTS-1).struct.rename_fields(
                                [f"part_{i+1}" for i in range(NUM_PARTS)]
                              ).alias("parts")
                            )
              .unnest("parts")
              .drop("id_original", "part_1", "part_3", "part_4")
              .with_columns(
                            pl.col("part_5").str.split_exact(by=":", n=1).struct.rename_fields(
                                ["pos_", "subid"]
                              ).alias("positions")
                            )
              .unnest("positions")
              .drop("part_5", "pos_")
              .cast({"subid": pl.Int64})
              .rename({"part_2": "id_"})
          )
  return df_id_split



def promoter_breath_summary(lz_prombound_joined: pl.LazyFrame) -> pl.LazyFrame:

    min_dF_lf = (
        lz_prombound_joined
          .group_by(["id","subid"])
          .agg(
            pl.col("dF").min().alias("min_dF")
          )
    )

    df_with_min_dF = lz_prombound_joined.join(
        min_dF_lf,
        on=["id", "subid"],
        how="left"
    ).with_columns(
        (pl.col("dF") - pl.col("min_dF")).alias("dF_reduced")
    )

    df_with_exp = df_with_min_dF.with_columns(
        (-pl.col('dF_reduced')).exp().alias('exp_neg_beta_dF')
    )

    # Compute Z per group (id, subid): sum(exp_neg_beta_dF)
    z_df = df_with_exp.group_by(['id', 'subid']).agg(
        pl.col('exp_neg_beta_dF').sum().alias('Z')
    )

    # Join Z back to df to compute P = exp_neg_beta_dF / Z
    df_with_p = df_with_exp.join(z_df, on=['id', 'subid'], how='left').with_columns(
        (pl.col('exp_neg_beta_dF') / pl.col('Z')).alias('P')
    )

    print(df_with_p.limit(105).collect())



    position_lf = df_with_p.group_by(['id', 'subid']).agg(
        (pl.col('P') * pl.col('left_open')).sum().alias('avg_left_unwrap'),
        (pl.col('P') * pl.col('right_open')).sum().alias('avg_right_unwrap'),
        ((pl.col('left_open') + pl.col('right_open')) * pl.col('P')).sum().alias('avg_total_unwrap')
    )

    return position_lf


def main(breath_data: Path, id_lookup: Path, out_file: Path) -> None:
    lz_promdt, lz_promid = get_promoter_breathing_data(breath_data, id_lookup)

    lz_lookup_refined = refine_lookup_tb(lz_promid)

    lz_prom_joined = (
        lz_promdt.drop("subid")
          .join(
              lz_lookup_refined,
              left_on="id",
              right_on="id_code",
              how="left"   
          ).rename({"F_enthalpy": "F_freedna", ##### this is a hack to rename the column because of the mistake in 06_multiharmonic_fe.job. Remove it later
                    "F_freedna": "F_enthalpy"})
            .with_columns((13-pl.col("right_index")).alias("right_open"), 
                                pl.col("left_index").alias("left_open")
                        )
            .drop("right_index", "left_index", "id")
            .rename({"id_": "id"})
    )

    position_lf = promoter_breath_summary(lz_prom_joined)
    position_df = position_lf.collect().to_pandas()
    
    position_df.to_csv(out_file, sep="\t", index=False)
    parent_logger.info(f"Results saved to {out_file}")




def arg_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Process promoter breathing data.")
    parser.add_argument("--input_file", type=Path, required=True, help="Path to the input promoter breathing data file.")
    parser.add_argument("--id_lookup_path", type=Path, required=True, help="Path to the ID lookup file.")
    parser.add_argument("--outfile", type=Path, help="Output file path for results.")
    return parser.parse_args()

if __name__ == "__main__":
  import time 
  start = time.perf_counter()

  BREATH_DT = EXAMPLES_DATA / f'minpoint_boundpromoter_regions_breath/dF/001.tsv'
  ID_LOOKUP = EXAMPLES_DATA / f'minpoint_boundpromoter_regions_breath/dF/id_lookup.tsv'
  OUTFILE = EXAMPLES_OUT / f'minpoint_boundpromoter_regions_breath/breath_summary/001.tsv'


  parent_logger = get_logger(__name__, level=logging.INFO)

  # args = arg_parser()
  # if args.input_file:
  #   BREATH_DT = args.input_file
  # if args.id_lookup_path:
  #   ID_LOOKUP = args.id_lookup_path
  # if args.outfile:
  #   OUTFILE = args.outfile

  if not BREATH_DT.exists():
    raise FileNotFoundError(f"Input file {BREATH_DT} does not exist.")

  if not ID_LOOKUP.exists():
    raise FileNotFoundError(f"ID lookup file {ID_LOOKUP} does not exist.")

  if not OUTFILE.parent.exists():
    OUTFILE.parent.mkdir(parents=True, exist_ok=True)

  parent_logger.info(f"Processing promoter breathing data from {BREATH_DT} with ID lookup {ID_LOOKUP}")
  main(breath_data=BREATH_DT, id_lookup=ID_LOOKUP, out_file=OUTFILE)

  end = time.perf_counter()
  parent_logger.info(f"Execution time: {end - start:.2f} seconds")