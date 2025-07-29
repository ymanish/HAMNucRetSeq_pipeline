
import argparse
from pathlib import Path

import pandas as pd

CHUNKSIZE_DEFAULT = 500_000 
FLOAT_FMT = "%.6f"

def main(infile: Path, outfile: Path, chunksize: int) -> None:
    # only read the columns you actually need
    usecols = [
        "id", "subid", "left_index", "right_index",
        "Ffree_bound", "Ffree_unbound", "F", "F_enthalpy"
    ]

    reader = pd.read_csv(
        infile,
        sep="\t",
        usecols=usecols,
        dtype={"id": str, "subid": str},
        chunksize=chunksize
    )

    first = True
    for chunk in reader:
        # compute new cols
        chunk["F_new"]          = chunk["F"] + chunk["Ffree_unbound"]
        chunk["F_freedna_new"]  = chunk["Ffree_unbound"] + chunk["Ffree_bound"]
        # chunk["F_enthalpy_new"] = chunk["F_enthalpy"]
        chunk["F_entropy_new"]   = chunk["F_new"] - chunk["F_enthalpy"]
        chunk.drop(columns=["F", "Ffree_bound", "Ffree_unbound"], inplace=True)


        chunk = chunk.rename(columns={
            "F_new": "F",
            "F_freedna_new": "F_freedna",
            "F_entropy_new": "F_entropy"
        })

        out_cols = [
            "id", "subid", "left_index", "right_index",
            "F", "F_freedna", "F_enthalpy", "F_entropy"
        ]

        chunk.to_csv(
            outfile,
            sep="\t",
            columns=out_cols,
            index=False,
            header=first,
            mode="w" if first else "a",
            float_format=FLOAT_FMT
        )
        first = False

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Compute multiharmonic energies in streaming mode"
    )
    p.add_argument("--infile",     type=str, help="Input TSV file")
    p.add_argument("--outfile",    type=str, help="Output TSV file")
    p.add_argument(
        "--chunksize",
        type=int,
        default=CHUNKSIZE_DEFAULT,
        help="Rows per pandas chunk"
    )
    args = p.parse_args()
    main(Path(args.infile), Path(args.outfile), args.chunksize)