#!/usr/bin/env python3

"""
Summarize VBF-tag purity vs signal mass for the Central-to-VBFTag strategy.

The headline quantity is:
 - Central-to-VBFTag purity:
   N(pass_trigger_baseline && isVBF && vbf_pair_lhematched) /
   N(pass_trigger_baseline && isVBF)
"""

import os
import re
import csv
import glob
import argparse

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

DEFAULT_INPUT_DIR = os.path.dirname(os.path.abspath(__file__))


def hex_to_rootcolor(hexcode):
    hexcode = hexcode.lstrip("#")
    r = int(hexcode[0:2], 16) / 255.0
    g = int(hexcode[2:4], 16) / 255.0
    b = int(hexcode[4:6], 16) / 255.0
    return ROOT.TColor.GetColor(r, g, b)


STRATEGY_STYLES = {
    "central": (hex_to_rootcolor("#3f90da"), "Central-to-VBFTag"),
}


parser = argparse.ArgumentParser(description="Plot VBF-tag purity vs mass for central-widejet study outputs")
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
parser.add_argument("--outdir", type=str, default=".", help="Output directory.")
parser.add_argument("--label", type=str, default="VBFHTo2B", help="Label drawn on the canvas.")
parser.add_argument("--outname", type=str, default="summary_vbfPurity_vs_mass", help="Output basename.")
parser.add_argument("--csv", type=str, default=None, help="Optional CSV output path.")
parser.add_argument(
    "--selectionLine",
    action="append",
    default=[],
    help="Selection line drawn in the annotation box. Can be passed multiple times.",
)
args = parser.parse_args()


def cms_style():
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.14)
    ROOT.gStyle.SetPadTopMargin(0.12)
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleOffset(1.1, "X")
    ROOT.gStyle.SetTitleOffset(1.2, "Y")


def draw_label(extra_text):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.12, 0.895, "CMS")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextFont(52)
    latex2.SetTextSize(0.036)
    latex2.DrawLatex(0.20, 0.895, "Simulation Preliminary")

    latex3 = ROOT.TLatex()
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.036)
    latex3.SetTextAlign(31)
    latex3.DrawLatex(0.95, 0.895, extra_text)


def draw_note(lines, x1=0.16, y1=0.69, x2=0.82, y2=0.86):
    box = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    box.SetFillColor(ROOT.kGray)
    box.SetFillColorAlpha(ROOT.kGray, 0.12)
    box.SetLineColor(ROOT.kGray + 1)
    box.SetLineWidth(2)
    box.SetTextAlign(12)
    box.SetTextFont(42)
    box.SetTextSize(0.030)
    for line in lines:
        box.AddText(line)
    box.Draw()
    return box


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


def describe_selection(tag):
    if tag == "unknown":
        return ["VBF Selection applied"]
    central_match = re.search(r"centPt(\d+)(?:-Deta([0-9p]+))?_fwdPt", tag)
    forward_match = re.search(r"_fwdPt(\d+)_", tag)
    vbf_match = re.search(r"_fwdPt\d+_(?:VBF-)?Deta([0-9p]+)", tag)
    mjj_match = re.search(r"-Mjj(\d+)", tag)

    lines = ["VBF Selection applied"]
    if central_match:
        line = f"Central pair: p_{{T}} > {central_match.group(1)} GeV"
        if central_match.group(2):
            line += f", |#Delta#eta| < {central_match.group(2).replace('p', '.')}"
        lines.append(line)
    if forward_match and vbf_match:
        line = (
            f"Forward pair: p_{{T}} > {forward_match.group(1)} GeV, "
            f"|#Delta#eta| > {vbf_match.group(1).replace('p', '.')}"
        )
        if mjj_match:
            line += f", #it{{m}}_{{jj}} > {mjj_match.group(1)} GeV"
        lines.append(line)
    return lines


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


def summarize_file(path):
    root_file, tree = open_tree(path)
    baseline = "pass_trigger_baseline == 1"
    central_tagged = f"{baseline} && isVBF == 1"
    central_matched = f"{central_tagged} && vbf_pair_lhematched == 1"

    n_central_tagged = count(tree, central_tagged)
    n_central_matched = count(tree, central_matched)

    row = {
        "file": path,
        "mass": infer_mass(path),
        "selection": infer_selection_tag(path),
        "central_tagged": n_central_tagged,
        "central_matched": n_central_matched,
        "central_purity_pct": 100.0 * n_central_matched / n_central_tagged if n_central_tagged > 0 else 0.0,
    }
    root_file.Close()
    return row


def write_csv(path, rows):
    fieldnames = [
        "file",
        "mass",
        "selection",
        "central_tagged",
        "central_matched",
        "central_purity_pct",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def style_hist(hist, strategy):
    color, _ = STRATEGY_STYLES[strategy]
    hist.SetFillColorAlpha(color, 0.45)
    hist.SetLineColor(color)
    hist.SetLineWidth(3)


def draw_bar_labels(hist, values, x_shift, y_shift, color):
    labels = []
    for idx, value in enumerate(values, start=1):
        latex = ROOT.TLatex()
        latex.SetTextFont(42)
        latex.SetTextSize(0.028)
        latex.SetTextAlign(21)
        latex.SetTextColor(color)
        latex.DrawLatex(idx + x_shift, value + y_shift, f"{value:.1f}%")
        labels.append(latex)
    return labels


def save_summary_plot(rows, outdir):
    if not rows:
        raise RuntimeError("No rows to plot.")
    masses = sorted(row["mass"] for row in rows)
    canvas = ROOT.TCanvas("c_vbf_purity_summary", "", 1000, 800)

    h_central = ROOT.TH1F("h_vbfpur_central", ";Signal mass [GeV];VBF-tag purity [%]", len(masses), 0.5, len(masses) + 0.5)

    values_central = []
    for idx, row in enumerate(rows, start=1):
        h_central.GetXaxis().SetBinLabel(idx, f"{row['mass']}")
        h_central.SetBinContent(idx, row["central_purity_pct"])
        values_central.append(row["central_purity_pct"])

    style_hist(h_central, "central")

    h_central.SetMaximum(100.0)
    h_central.SetMinimum(0.0)
    h_central.SetBarWidth(0.58)
    h_central.SetBarOffset(0.21)
    h_central.GetXaxis().SetTitleOffset(1.45)

    h_central.Draw("b")

    labels = []
    labels.extend(draw_bar_labels(h_central, values_central, 0.0, 1.0, STRATEGY_STYLES["central"][0]))

    selection_lines = args.selectionLine[:]
    if not selection_lines:
        selections = sorted({row["selection"] for row in rows})
        if len(selections) == 1:
            selection_lines = describe_selection(selections[0])
        else:
            selection_lines = ["VBF Selection applied", f"{len(selections)} working points combined"]
    note = draw_note(selection_lines)

    draw_label(args.label)
    os.makedirs(outdir, exist_ok=True)
    outbase = os.path.join(outdir, args.outname)
    canvas.SaveAs(outbase + ".png")
    canvas.SaveAs(outbase + ".pdf")
    return outbase + ".png", outbase + ".pdf"


def main():
    cms_style()
    files = expand_inputs(input_patterns())
    rows = [summarize_file(path) for path in files]
    rows = sorted(rows, key=lambda row: row["mass"])

    if args.csv:
        csv_dir = os.path.dirname(args.csv)
        if csv_dir:
            os.makedirs(csv_dir, exist_ok=True)
        write_csv(args.csv, rows)
        print(f"[INFO] Wrote CSV summary: {args.csv}")

    png_path, pdf_path = save_summary_plot(rows, args.outdir)
    print(f"[INFO] Wrote purity-summary plot: {png_path}")
    print(f"[INFO] Wrote purity-summary plot: {pdf_path}")


if __name__ == "__main__":
    main()
