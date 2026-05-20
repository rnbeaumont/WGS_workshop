#!/usr/bin/env python3
"""Convert REGENIE association output to manc_cojo summary-stat format.

Default output columns (tab-delimited, GCTA-COJO input):
SNP A1 A2 freq b se p N

Assumptions for REGENIE single-variant output:
- A1 is the tested allele (REGENIE ALLELE1)
- A2 is the other allele (REGENIE ALLELE0)
- p is derived from REGENIE LOG10P as 10^(-LOG10P)
"""

import argparse
import gzip
import math
import sys
from pathlib import Path
from typing import Iterable, Optional, TextIO, Tuple, List


REQUIRED_COLUMNS = {
    "CHROM",
    "GENPOS",
    "ID",
    "ALLELE0",
    "ALLELE1",
    "A1FREQ",
    "N",
    "BETA",
    "SE",
    "LOG10P",
}


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert REGENIE output to manc_cojo/GCTA columns: "
            "SNP A1 A2 freq b se p N"
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input REGENIE file (.regenie or .regenie.gz)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output file path",
    )
    parser.add_argument(
        "--p_thresh",
        default=0,
        type=float,
        help="-log10 P-value threshold for inclusion",
    )
    parser.add_argument(
        "--test",
        default="ADD",
        help=(
            "Only include rows where TEST matches this value. "
            "Use --test all to keep all tests. Default: ADD"
        ),
    )
    parser.add_argument(
        "--n-override",
        type=int,
        default=None,
        help="Override N with a fixed sample size for all rows",
    )
    parser.add_argument(
        "--with-chr-bp",
        action="store_true",
        help=(
            "Write an extended 10-column file with leading Chr and bp columns "
            "(Chr SNP bp A1 A2 freq b se p n). "
            "By default, writes strict 8-column GCTA-COJO input."
        ),
    )
    return parser.parse_args()


def safe_p_from_log10p(log10p_str: str) -> str:
    value = log10p_str.strip()
    if value in {"", "NA", "NaN", "nan"}:
        return "NA"

    log10p = float(value)
    if math.isnan(log10p):
        return "NA"
    if math.isinf(log10p):
        return "0" if log10p > 0 else "1"
    if log10p < 0:
        return "NA"

    p = 10.0 ** (-log10p)

    # manc_cojo can be sensitive to scientific-notation p-values,
    # so write fixed-decimal strings instead of e-notation.
    # Keep enough decimals to preserve very small p-values.
    decimals = min(max(int(math.ceil(log10p)) + 8, 8), 400)
    p_str = f"{p:.{decimals}f}".rstrip("0").rstrip(".")
    if p_str == "":
        return "0"
    return p_str


def validate_header(columns: List[str]) -> None:
    missing = sorted(REQUIRED_COLUMNS.difference(columns))
    if missing:
        raise ValueError(
            "Input is missing required columns: " + ", ".join(missing)
        )


def iter_rows(handle: TextIO) -> Iterable[List[str]]:
    for line in handle:
        line = line.strip()
        if not line:
            continue
        yield line.split()


def convert(
    input_path: Path,
    output_path: Path,
    test_filter: str,
    n_override: Optional[int],
    with_chr_bp: bool,
    pval_thresh: float,
) -> Tuple[int, int]:
    written = 0
    skipped = 0

    with open_text(input_path) as fin:
        rows = iter_rows(fin)
        try:
            header = next(rows)
        except StopIteration as exc:
            raise ValueError("Input file is empty") from exc

        validate_header(header)
        idx = {name: i for i, name in enumerate(header)}

        with output_path.open("w", encoding="utf-8", newline="") as fout:
            if with_chr_bp:
                fout.write("Chr\tSNP\tbp\tA1\tA2\tfreq\tb\tse\tp\tn\n")
            else:
                fout.write("SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n")

            for fields in rows:
                if len(fields) < len(header):
                    skipped += 1
                    continue

                if test_filter.lower() != "all" and "TEST" in idx:
                    if fields[idx["TEST"]] != test_filter:
                        continue

                p_value = safe_p_from_log10p(fields[idx["LOG10P"]])
                n_value = str(n_override) if n_override is not None else fields[idx["N"]]

                if with_chr_bp:
                    out = [
                        fields[idx["CHROM"]],
                        fields[idx["ID"]],
                        fields[idx["GENPOS"]],
                        fields[idx["ALLELE1"]],  # tested allele
                        fields[idx["ALLELE0"]],  # non-tested allele
                        fields[idx["A1FREQ"]],
                        fields[idx["BETA"]],
                        fields[idx["SE"]],
                        p_value,
                        n_value,
                    ]
                else:
                    out = [
                        fields[idx["ID"]],
                        fields[idx["ALLELE1"]],  # tested allele
                        fields[idx["ALLELE0"]],  # non-tested allele
                        fields[idx["A1FREQ"]],
                        fields[idx["BETA"]],
                        fields[idx["SE"]],
                        p_value,
                        n_value,
                    ]
                if float(fields[idx["LOG10P"]].strip()) >= pval_thresh:
                    fout.write("\t".join(out) + "\n")
                    written += 1
                else:
                    skipped += 1

    return written, skipped


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        print(f"ERROR: input file not found: {input_path}", file=sys.stderr)
        return 2

    try:
        written, skipped = convert(
            input_path=input_path,
            output_path=output_path,
            test_filter=args.test,
            n_override=args.n_override,
            with_chr_bp=args.with_chr_bp,
            pval_thresh=args.p_thresh,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"Wrote {written} rows to {output_path} "
        f"(skipped malformed rows: {skipped})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
