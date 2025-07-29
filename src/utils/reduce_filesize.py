
import os
import glob
import pandas as pd
from src.config.path import DATA_DIR, OUT_DIR
from src.utils.logger_util import get_logger


def main():
    LOOKUP_FILE = os.path.join(OUTPUT_DIR, "id_lookup.tsv")


    parent_logger.info("Scanning for unique IDs ....")
    unique_ids = set()
    for fn in glob.glob(os.path.join(LARGE_DIR, "*.tsv")):
        for chunk in pd.read_csv(fn, sep="\t", usecols=["id"], chunksize=CHUNKSIZE):
            unique_ids |= set(chunk["id"].unique())

    parent_logger.info(f"Found {len(unique_ids):,} unique IDs. Building lookup table ...")
    # map each long ID -> small int code
    id_list = sorted(unique_ids)
    id_to_code = {sid: idx for idx, sid in enumerate(id_list)}

    # write the lookup table 
    pd.DataFrame({
        "id_original": id_list,
        "id_code":      list(range(len(id_list)))
    }).to_csv(LOOKUP_FILE, sep="\t", index=False)
    parent_logger.info(f"Lookup saved to {LOOKUP_FILE!r}")

    parent_logger.info("Re-writing and rounding float....")
    for fn in glob.glob(os.path.join(LARGE_DIR, "*.tsv")):
        basename = os.path.basename(fn)
        out_fn = os.path.join(OUTPUT_DIR, basename)
        parent_logger.info(f" -> {basename}")
        first = True

        for chunk in pd.read_csv(fn, sep="\t", chunksize=CHUNKSIZE, dtype={"subid":str}):
            # a) replace long IDs with small int codes
            chunk["id"] = chunk["id"].map(id_to_code).astype("int32")

            # b) round all floats
            float_cols = chunk.select_dtypes("float64").columns
            chunk[float_cols] = chunk[float_cols].round(FLOAT_DP)

            # c) write out
            chunk.to_csv(
                out_fn,
                sep="\t",
                index=False,
                header=first,
                mode="w" if first else "a",
                float_format=f"%.{FLOAT_DP}f"
            )
            first = False


def arg_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Reduce file size of large TSV files by rounding floats and replacing long IDs with small integers.")
    parser.add_argument("--input_dir", type=str,  help="Directory containing input TSV files.")
    parser.add_argument("--output_dir", type=str,  help="Directory to save reduced TSV files.")
    parser.add_argument("--chunksize", type=int, default=1_000_000, help="Number of rows to read in each chunk.")
    parser.add_argument("--float_dp", type=int  , default=4, help="Number of decimal places to round floats.")
    return parser.parse_args()


if __name__ == "__main__":
    import logging
    parent_logger = get_logger(__name__, log_file=None, level=logging.INFO)

    LARGE_DIR = f"{OUT_DIR}/minpoint_boundpromoter_regions_breath/crystal_freedna_md_merged"
    OUTPUT_DIR = f"{LARGE_DIR}/reduced"


    args = arg_parser()
    if args.input_dir:
        LARGE_DIR = args.input_dir
    if args.output_dir:
        OUTPUT_DIR = args.output_dir
    if args.chunksize:
        CHUNKSIZE = args.chunksize
    if args.float_dp:
        FLOAT_DP = args.float_dp 

    if not os.path.exists(LARGE_DIR):
        raise FileNotFoundError(f"Input directory {LARGE_DIR} does not exist. Please check the path.")
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    parent_logger.info(f"Processing files in: {LARGE_DIR}")
    parent_logger.info(f"Output will be written to: {OUTPUT_DIR}")
    parent_logger.info(f"Using chunk size: {CHUNKSIZE}")
    parent_logger.info(f"Using float decimal places: {FLOAT_DP}")

    main()
