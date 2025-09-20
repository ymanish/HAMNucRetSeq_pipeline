## **5. Promoter Region Nucleosome Breathing Propensity**


from src.config.path import DATA_DIR, RESULTS_DIR, NUCMONTE_DIR
import polars as pl
from src.config.path import EXAMPLES_DATA, EXAMPLES_OUT
from pathlib import Path
import csv
from src.utils.logger_util import get_logger
import logging
from py_analysis.adsorbpipe.extract_eads import main as extract_eads_func
from py_analysis.config.eads_var import Parameters
import numpy as np



def add_fully_open_states_lazy(lzdf: pl.LazyFrame, length: int = 14) -> pl.LazyFrame:
    """
    Ensure that each id has all combinations where left_open + right_open == 'length'.
    Works fully lazily inside Polars.
    """

    # 1. Do we already have at least one fully-open state somewhere?
    fully_open_present = (
        lzdf
        .filter((pl.col("left_open") + pl.col("right_open")) == length)
        .select(pl.len())
        .collect(engine="streaming")            
        .item() > 0
    )
    if fully_open_present:
        print("Fully open states already present, skipping addition.")
        return lzdf                      

    fully_open_states = (
        pl.DataFrame({
            "left_open": list(range(length + 1)),
        })
        .with_columns((pl.lit(length) - pl.col("left_open")).alias("right_open"))
        .lazy()
    )

    per_id = (
        lzdf
        .select("id", "subid", "F_freedna")
        .unique()
    )

    # 4. Cross-join -> 15 rows per id, then add the constant columns
    additions = (
        per_id
        .join(fully_open_states, how="cross")
        .with_columns([
            pl.col("F_freedna").alias("F"),
            pl.col("F_freedna").alias("F_entropy"),
            pl.lit(0.0).alias("F_enthalpy"),
            pl.col("F_freedna").alias("F_freedna"),
            pl.lit(0.0).alias("dF"),
        ])
    )

    original_schema = dict(lzdf.collect_schema())
    original_cols   = list(original_schema.keys())

    # for any column in the original that additions is missing, add it as null of the right dtype
    additions_schema = dict(additions.collect_schema())
    additions = additions.with_columns([
       pl.lit(None).cast(original_schema[col]).alias(col)
       for col in original_cols
       if col not in additions_schema
        ]).select(original_cols)



    # additions = additions.with_columns([
    #     pl.lit(None).cast(original_schema[col]).alias(col)
    #     for col in original_cols
    #     if col not in additions.collect_schema().to_dict().keys()
    # ]).select(original_cols)

    # reorder additions to match the original column order
    additions = additions.select(original_cols)

    return (
        pl.concat([lzdf, additions])
          .unique(subset=["id", "subid", "left_open", "right_open"])
    )



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
                                   schema_overrides=schema)

    lz_promid = pl.scan_csv(id_lookup_path,
                                    separator="\t",
                                      infer_schema_length=0,
                                        schema_overrides=schema_id)
    


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



def add_adsorption_energy(lz_prombound_joined: pl.LazyFrame, eads_array: np.ndarray) -> pl.LazyFrame:
  # E_ads_list = eads_array.tolist()
  # L = len(E_ads_list)

  # lz_prombound_joined = lz_prombound_joined.with_columns(
  #     pl.struct(["left_open", "right_open"])
  #       .map_elements(lambda row: -sum(E_ads_list[row["left_open"] : L - row["right_open"]]), return_dtype=pl.Float64)
  #       .alias("Adsorp_F")
  # )

  # lz_prombound_joined = lz_prombound_joined.with_columns((pl.col("dF") + pl.col("Adsorp_F")).alias("dF_total"))

  L = len(eads_array)
  cum_E = np.cumsum(np.insert(eads_array, 0, 0.0)).tolist()
  # 3) compute Adsorp_F = -(cum_E[L-right_open] - cum_E[left_open])
  #    and then dF_total = dF + Adsorp_F
  lz_prombound_joined_ = lz_prombound_joined.with_columns(
                                                          Adsorp_F = -(
                                                                          pl.lit(cum_E, dtype=pl.List(pl.Float64)).list.get(pl.lit(L) - pl.col("right_open")) 
                                                                          - pl.lit(cum_E, dtype=pl.List(pl.Float64)).list.get(pl.col("left_open"))
                                                                        )
                                                        ).with_columns(
                                                                    dF_total = pl.col("dF") + pl.col("Adsorp_F")
                                                                      )

  return lz_prombound_joined_




def add_constant_adsorption_energy(lz_prombound_joined: pl.LazyFrame, Eads_val: float) -> pl.LazyFrame:

  lz_prombound_joined_ = lz_prombound_joined.with_columns(Adsorp_F = -Eads_val*(14-(pl.col("left_open") + pl.col("right_open")))).with_columns(
                                dF_total = pl.col("dF") + pl.col("Adsorp_F")
                                )
  return lz_prombound_joined_





def promoter_breath_summary(lz_prombound_joined: pl.LazyFrame) -> pl.LazyFrame:

    min_dF_lf = (
        lz_prombound_joined
          .group_by(["id","subid"])
          .agg(
            pl.col("dF_total").min().alias("min_dF")
          )
    )

    df_with_min_dF = lz_prombound_joined.join(
        min_dF_lf,
        on=["id", "subid"],
        how="left"
    ).with_columns(
        (pl.col("dF_total") - pl.col("min_dF")).alias("dF_reduced")
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


    position_lf = df_with_p.group_by(['id', 'subid']).agg(
        (pl.col('P') * pl.col('left_open')).sum().alias('avg_left_unwrap'),
        (pl.col('P') * pl.col('right_open')).sum().alias('avg_right_unwrap'),
        ((pl.col('left_open') + pl.col('right_open')) * pl.col('P')).sum().alias('avg_total_unwrap')
    )

    return position_lf




def main(breath_data: Path, 
         id_lookup: Path, 
         out_file: Path, 
         outfile_energy: Path, 
         params: Parameters | None = None, 
         homo_adsorption_val: float | None = None) -> None:
    
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
    lz_prom_joined = add_fully_open_states_lazy(lz_prom_joined, length=14)

    if homo_adsorption_val is not None:
        # homogeneous adsorption with provided constant
        lz_full = add_constant_adsorption_energy(lz_prom_joined, homo_adsorption_val)
    else:      
        eads_array = extract_eads_func(eads_params=params, 
                                  sequence="601",
                                  leftdwell_path=NUCMONTE_DIR/"left.csv",
                                  rightdwell_path=NUCMONTE_DIR/"right.csv",
                                  dna_histone_dir=NUCMONTE_DIR/"dna_histone")

        lz_full = add_adsorption_energy(lz_prom_joined, eads_array)



    position_lf = promoter_breath_summary(lz_full)
    position_df = position_lf.collect(engine="streaming").to_pandas()
    
    position_df.to_csv(out_file, sep="\t", index=False)
    parent_logger.info(f"Results saved to {out_file}")

    # Save energy landscape data
    #energy_df = lz_full.drop("Adsorp_F").collect(engine="streaming").to_pandas()
    #energy_df.to_csv(outfile_energy, sep="\t", index=False)
    lz_full.drop("Adsorp_F").sink_csv(outfile_energy, separator="\t")
    parent_logger.info(f"Energy landscape saved to {outfile_energy}")




def arg_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Process promoter breathing data.")
    parser.add_argument("--input_file", type=Path, required=True, help="Path to the input promoter breathing data file.")
    parser.add_argument("--id_lookup_path", type=Path, required=True, help="Path to the ID lookup file.")
    parser.add_argument("--outfile", type=Path, help="Output file path for results.")
    parser.add_argument("--outfile_energy", type=Path, help="Output file path for energy landscape results.")
    parser.add_argument("--eads_param_nucmethod", type=str, default="crystal", help="Nucleosome method for Eadsorption extraction.")
    parser.add_argument("--eads_param_freednamethod", type=str, default="md", help="Free DNA method for Eadsorption extraction.")
    parser.add_argument("--eads_param_krescfactor", type=float, default=1.0, help="Krescfactor for Eadsorption extraction.")
    parser.add_argument("--eads_param_Eoffset", type=float, default=0.0, help="Energy offset for Eadsorption extraction.")
    parser.add_argument(
        "--eads_homo_adsorption",
        nargs='?',            # flag may take an optional float
        const=16.32,          # default constant if no value is provided
        type=float,
        default=None,         # None means "use heterogeneous model"
        help=(
            "Use homogeneous adsorption model. "
            "Optionally specify constant adsorption energy (float). "
            "If flag is given without a value, uses default const=16.32."
        )
    )
    return parser.parse_args()

if __name__ == "__main__":
  import time 
  start = time.perf_counter()

  # BREATH_DT = EXAMPLES_DATA / f'minpoint_boundpromoter_regions_breath/dF/001.tsv'
  # ID_LOOKUP = EXAMPLES_DATA / f'minpoint_boundpromoter_regions_breath/dF/id_lookup.tsv'
  # OUTFILE = EXAMPLES_OUT / f'minpoint_boundpromoter_regions_breath/breath_summary/001.tsv'
  # OUTFILE_ENERGY = EXAMPLES_OUT / f'minpoint_boundpromoter_regions_breath/STATE_ENERGY/001.tsv'


  parent_logger = get_logger(__name__, level=logging.INFO)

  args = arg_parser()
  if args.input_file:
    BREATH_DT = args.input_file
  if args.id_lookup_path:
    ID_LOOKUP = args.id_lookup_path
  if args.outfile:
    OUTFILE = args.outfile
  if args.outfile_energy:
    OUTFILE_ENERGY = args.outfile_energy


  # decide between heterogeneous vs homogeneous adsorption
  if args.eads_homo_adsorption is None:
      EADS_PARAMS = Parameters(
          NUCMETHOD       = args.eads_param_nucmethod,
          FREEDNA_METHOD  = args.eads_param_freednamethod,
          KRESCFACTOR     = args.eads_param_krescfactor,
          E_OFFSET        = args.eads_param_Eoffset,
      )
      homo_val = None
      parent_logger.info(f"Using heterogeneous adsorption with parameters: {EADS_PARAMS}")
  else:
      EADS_PARAMS = None
      homo_val = args.eads_homo_adsorption
      parent_logger.info(f"Using homogeneous adsorption model with constant Eads = {homo_val}")



  if not BREATH_DT.exists():
    raise FileNotFoundError(f"Input file {BREATH_DT} does not exist.")

  if not ID_LOOKUP.exists():
    raise FileNotFoundError(f"ID lookup file {ID_LOOKUP} does not exist.")

  if not OUTFILE.parent.exists():
    OUTFILE.parent.mkdir(parents=True, exist_ok=True)

  if not OUTFILE_ENERGY.parent.exists():
    OUTFILE_ENERGY.parent.mkdir(parents=True, exist_ok=True)

  parent_logger.info(f"Processing promoter breathing data from {BREATH_DT} with ID lookup {ID_LOOKUP}")
  parent_logger.info(f"Output will be written to {OUTFILE}")
  parent_logger.info("Starting promoter breathing analysis...")
  parent_logger.info(f"Using homogeneous adsorption model: {homo_val}")
  parent_logger.info(f"Using Eadsorption parameters: {EADS_PARAMS}")


  main(
    breath_data=BREATH_DT,
    id_lookup=ID_LOOKUP,
    out_file=OUTFILE,
    outfile_energy=OUTFILE_ENERGY,
    params=EADS_PARAMS,
    homo_adsorption_val=homo_val,
)

  end = time.perf_counter()
  parent_logger.info(f"Execution time: {end - start:.2f} seconds")
