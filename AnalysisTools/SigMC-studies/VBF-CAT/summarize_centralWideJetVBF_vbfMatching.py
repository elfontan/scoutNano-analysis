#!/usr/bin/env python3

"""
Summarize VBF-pair truth matching for the Central-to-VBFTag strategy.

For each input tree, it reports:
 - number of truth-available events with an LHE VBF pair
 - number of reco VBF-pair candidates
 - number of tagged VBF events
 - matched efficiency with respect to truth-available events
 - matched purity among tagged VBF events

It can also write a CSV for later comparisons across masses/configurations.
"""

import os
import re
import csv
import glob
import argparse

import ROOT

ROOT.gROOT.SetBatch(True)

DEFAULT_INPUT_DIR = os.path.dirname(os.path.abspath(__file__))


parser = argparse.ArgumentParser(description="Summarize VBF reco matching for central-widejet study outputs")
parser.add_argument("--inputs", nargs="+", default=None, help="Input ROOT files or glob patterns.")
parser.add_argument(
    "--selection",
    type=str,
    default=None,
    help="Selection tag used in the ROOT filename, e.g. centPt30-Deta1p3_fwdPt24_VBF-Deta4-Mjj500.",
)
parser.add_argument(
    "--inputDir",
    type=str,
    default=DEFAULT_INPUT_DIR,
    help="Directory to search when --selection is used.",
)
parser.add_argument(
    "--samplePrefix",
    type=str,
    default="CentralVBFHTo2B",
    help="Filename prefix used when --selection is used.",
)
parser.add_argument("--tree", type=str, default="Events", help="Tree name inside the input ROOT files.")
parser.add_argument("--writeCSV", action="store_true", help="Write a CSV summary.")
parser.add_argument("--csv", type=str, default="summary_centralWideJetVBF_vbfMatching.csv", help="CSV output path.")
args = parser.parse_args()


def input_patterns():
    if args.inputs:
        return args.inputs
    if args.selection:
        return [
            os.path.join(
                args.inputDir,
                f"{args.samplePrefix}_M*_centralWideJetVBF_scoutNano_{args.selection}.root",
            )
        ]
    raise RuntimeError("Please provide either --inputs or --selection.")


def expand_inputs(patterns):
    files = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            files.extend(matches)
        elif os.path.exists(pattern):
            files.append(pattern)
    if not files:
        raise RuntimeError("No input files found.")
    return sorted(set(files))


def open_tree(path):
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open input file: {path}")
    tree = root_file.Get(args.tree)
    if tree is None:
        root_file.Close()
        raise RuntimeError(f"Could not find tree '{args.tree}' in {path}")
    return root_file, tree


def count(tree, selection):
    return int(tree.GetEntries(selection))


def infer_mass(path):
    match = re.search(r"(?:^|[_-])M(\d+)(?:[_.-]|$)", os.path.basename(path), re.IGNORECASE)
    return int(match.group(1)) if match else -1


def infer_selection_tag(path):
    base = os.path.basename(path)
    match = re.search(
        r"(centPt\d+(?:-Deta[0-9p]+)?_fwdPt\d+_(?:VBF-)?Deta[0-9p]+(?:-Mjj\d+)?)",
        base,
        re.IGNORECASE,
    )
    return match.group(1) if match else "unknown"


def infer_sample_label(path):
    base = os.path.basename(path)
    sample = re.sub(r"_M\d+.*$", "", base)
    sample = re.sub(r"_centralWideJetVBF.*$", "", sample)
    return sample


def summarize_strategy(tree):
    baseline = "pass_trigger_baseline == 1"
    truth = f"{baseline} && has_lhe_vbf_pair == 1"
    pair_sel = f"{baseline} && has_vbf_pair == 1"
    tagged_sel = f"{baseline} && isVBF == 1"
    pair_matched_sel = f"{pair_sel} && vbf_pair_lhematched == 1"
    tagged_matched_sel = f"{tagged_sel} && vbf_pair_lhematched == 1"

    n_truth = count(tree, truth)
    n_pair = count(tree, pair_sel)
    n_pair_matched = count(tree, pair_matched_sel)
    n_tagged = count(tree, tagged_sel)
    n_tagged_matched = count(tree, tagged_matched_sel)

    return {
        "n_truth_lhe": n_truth,
        "n_pair": n_pair,
        "n_pair_matched": n_pair_matched,
        "n_tagged": n_tagged,
        "n_tagged_matched": n_tagged_matched,
        "pair_match_purity_pct": 100.0 * n_pair_matched / n_pair if n_pair > 0 else 0.0,
        "tag_match_eff_pct": 100.0 * n_tagged_matched / n_truth if n_truth > 0 else 0.0,
        "tag_match_purity_pct": 100.0 * n_tagged_matched / n_tagged if n_tagged > 0 else 0.0,
    }


def summarize_file(path):
    root_file, tree = open_tree(path)
    row = {
        "file": path,
        "mass": infer_mass(path),
        "sample": infer_sample_label(path),
        "selection": infer_selection_tag(path),
    }

    central = summarize_strategy(tree)
    for key, value in central.items():
        row[f"central_{key}"] = value

    root_file.Close()
    return row


def print_rows(rows):
    selections = sorted({row["selection"] for row in rows})
    if len(selections) == 1:
        print(f"Selection: {selections[0]}")
    else:
        print(f"Selections: {', '.join(selections)}")
    print()
    print(
        f"{'Mass':>6s} {'Selection':>34s} "
        f"{'CtrEff%':>8s} {'CtrPur%':>8s} "
        f"{'CtrTag':>8s}"
    )
    for row in rows:
        print(
            f"{row['mass']:6d} "
            f"{row['selection'][:34]:>34s} "
            f"{row['central_tag_match_eff_pct']:8.2f} "
            f"{row['central_tag_match_purity_pct']:8.2f} "
            f"{row['central_n_tagged']:8d}"
        )


def write_csv(path, rows):
    fieldnames = [
        "file",
        "sample",
        "mass",
        "selection",
        "central_n_truth_lhe",
        "central_n_pair",
        "central_n_pair_matched",
        "central_n_tagged",
        "central_n_tagged_matched",
        "central_pair_match_purity_pct",
        "central_tag_match_eff_pct",
        "central_tag_match_purity_pct",
    ]
    csv_dir = os.path.dirname(path)
    if csv_dir:
        os.makedirs(csv_dir, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    files = expand_inputs(input_patterns())
    rows = [summarize_file(path) for path in files]
    rows.sort(key=lambda row: (row["mass"], row["file"]))

    print_rows(rows)

    if args.writeCSV:
        write_csv(args.csv, rows)
        print(f"\n[INFO] Wrote CSV summary to: {args.csv}")


if __name__ == "__main__":
    main()
