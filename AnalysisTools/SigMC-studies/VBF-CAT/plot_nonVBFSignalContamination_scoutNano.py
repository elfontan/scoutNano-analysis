#!/usr/bin/env python3

"""
Plot the rate at which non-VBF ScoutNano signal samples are selected by the
current VBF-category strategies.

Input:
 - CSV produced by summarize_nonVBFSignalContamination_scoutNano.py

Outputs:
 - one summary plot per working point showing efficiencies vs mass for all signals
 - one summary plot per signal family showing efficiencies vs mass for all selections
"""

import os
import re
import glob
import csv
import argparse
from collections import defaultdict

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


parser = argparse.ArgumentParser(
    description="Plot non-VBF signal contamination vs mass from the ScoutNano VBF-category summary CSV."
)
parser.add_argument(
    "--csv",
    action="append",
    required=True,
    help="Input CSV from summarize_nonVBFSignalContamination_scoutNano.py. Repeat to combine files.",
)
parser.add_argument("--outdir", default=".", help="Output directory.")
parser.add_argument("--label", default="2024 (13.6 TeV)", help="Text label drawn on plots.")
parser.add_argument(
    "--families",
    nargs="+",
    default=None,
    help="Optional subset of sample families to plot, e.g. zprime_qq rsgraviton_2q.",
)
parser.add_argument(
    "--compareKey",
    default="eff_central_pct",
    choices=["eff_central_pct"],
    help="Efficiency column to plot on the y axis.",
)
parser.add_argument(
    "--workingPoint",
    default=None,
    help="Optional working_point value to restrict the inputs to one selection.",
)
args = parser.parse_args()


FAMILY_LABELS = {
    "zprime_qq": "Z' #rightarrow qq",
    "glugluspin0_2b": "GluGluSpin0 #rightarrow bb",
    "rsgraviton_2glu": "RS Graviton #rightarrow gg",
    "rsgraviton_2q": "RS Graviton #rightarrow qq",
}

_DRAWN_OBJECTS = []

CASE_STYLES = [
    (ROOT.TColor.GetColor("#3f90da"), 20),
    (ROOT.TColor.GetColor("#ffa90e"), 21),
    (ROOT.TColor.GetColor("#bd1f01"), 22),
    (ROOT.TColor.GetColor("#94a4a2"), 23),
    (ROOT.TColor.GetColor("#832db6"), 33),
    (ROOT.TColor.GetColor("#2ca02c"), 34),
]

FAMILY_STYLES = {
    "zprime_qq": (ROOT.TColor.GetColor("#3f90da"), 20),
    "glugluspin0_2b": (ROOT.TColor.GetColor("#832db6"), 21),
    "rsgraviton_2glu": (ROOT.TColor.GetColor("#ffa90e"), 22),
    "rsgraviton_2q": (ROOT.TColor.GetColor("#bd1f01"), 33),
}


def infer_mass_from_filename(path):
    base = os.path.basename(path)
    patterns = [
        r"_M(\d+)\.root$",
        r"Par-M-(\d+)\.root$",
        r"Par-[^-]+-M-(\d+)\.root$",
    ]
    for pattern in patterns:
        match = re.search(pattern, base)
        if match:
            return int(match.group(1))
    return -1


def infer_case_label(path):
    base = os.path.splitext(os.path.basename(path))[0]
    suffix = re.sub(r"^summary_nonVBFSignalContamination_", "", base)
    preset_match = re.search(r"(fwdPt\d+_VBF-Deta[0-9p]+(?:-Mjj\d+)?)", suffix)
    if preset_match:
        return preset_match.group(1)
    tokens = suffix.split("_")
    if len(tokens) <= 1:
        return "Default"
    case = "_".join(tokens[1:])
    deta_match = re.fullmatch(r"DEta(\d+)", case)
    if deta_match:
        return f"|#Delta#eta| > {deta_match.group(1)}"
    mjj_match = re.fullmatch(r"Mjj(\d+)", case)
    if mjj_match:
        return f"#it{{m}}_{{jj}} > {mjj_match.group(1)} GeV"
    return case.replace("_", " ")


def sanitize_token(text):
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(text)).strip("_")


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
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleOffset(1.10, "X")
    ROOT.gStyle.SetTitleOffset(1.15, "Y")


def read_rows(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row_case_label = row.get("working_point", "").strip() or infer_case_label(path)
            if args.workingPoint and row_case_label != args.workingPoint:
                continue

            row["mass"] = int(row["mass"])
            if row["mass"] < 0:
                row["mass"] = infer_mass_from_filename(row.get("file", ""))

            for key in [
                "n_total",
                "n_after_trigger",
                "n_central_selected",
                "eff_central_pct",
            ]:
                row[key] = float(row[key]) if "eff_" in key else int(row[key])

            row["case_label"] = row_case_label
            row["source_csv"] = path
            rows.append(row)
    return rows


def expand_csv_inputs(items):
    paths = []
    for item in items:
        matches = sorted(glob.glob(item))
        if matches:
            paths.extend(matches)
        elif os.path.exists(item):
            paths.append(item)
    if not paths:
        raise RuntimeError("No CSV inputs found.")
    return sorted(set(paths))


def deduplicate_rows(rows):
    deduped = []
    seen = set()
    for row in rows:
        key = (
            row.get("sample", ""),
            row.get("mass", -1),
            row.get("file", ""),
            row.get("case_label", ""),
            row.get("n_central_selected", -1),
        )
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    return deduped


def draw_label(extra_text):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.045)
    cms = latex.DrawLatex(0.12, 0.905, "CMS")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextFont(52)
    latex2.SetTextSize(0.036)
    prelim = latex2.DrawLatex(0.20, 0.905, "Simulation Preliminary")

    latex3 = ROOT.TLatex()
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.036)
    latex3.SetTextAlign(31)
    label = latex3.DrawLatex(0.95, 0.905, extra_text)
    _DRAWN_OBJECTS.extend([latex, latex2, latex3, cms, prelim, label])


def draw_note(lines, x1=0.15, y1=0.70, x2=0.45, y2=0.84):
    note = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    note.SetFillColor(ROOT.kGray)
    note.SetFillColorAlpha(ROOT.kGray, 0.12)
    note.SetLineColor(ROOT.kGray + 1)
    note.SetLineWidth(2)
    note.SetTextAlign(12)
    note.SetTextFont(42)
    note.SetTextSize(0.030)
    for line in lines:
        note.AddText(line)
    note.Draw()
    _DRAWN_OBJECTS.append(note)


def make_graph(points, key, color, marker, title):
    graph = ROOT.TGraph(len(points))
    for idx, row in enumerate(points):
        graph.SetPoint(idx, float(row["mass"]), float(row[key]))
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(marker)
    graph.SetMarkerSize(1.2)
    graph.SetLineWidth(3)
    graph.SetTitle(title)
    return graph


def make_y_frame(name, masses, ymax):
    xmin = 0.95 * min(masses)
    xmax = 1.05 * max(masses)
    if min(masses) == max(masses):
        xmin = min(masses) - 1.0
        xmax = max(masses) + 1.0
    frame = ROOT.TH1F(
        name,
        ";Signal mass [GeV];Selected fraction [%]",
        100,
        xmin,
        xmax,
    )
    frame.SetMinimum(0.0)
    frame.SetMaximum(40.0)
    return frame


def save_working_point_signal_overlay_plot(case_label, grouped_rows, outdir):
    _DRAWN_OBJECTS.clear()
    canvas = ROOT.TCanvas(f"c_{sanitize_token(case_label)}_signals", "", 1000, 800)
    masses = sorted({row["mass"] for rows in grouped_rows.values() for row in rows})
    if not masses:
        return

    ymax = max(row[args.compareKey] for rows in grouped_rows.values() for row in rows)
    frame = make_y_frame(f"frame_{sanitize_token(case_label)}_signals", masses, ymax)
    frame.Draw()
    _DRAWN_OBJECTS.append(frame)

    legend = ROOT.TLegend(0.54, 0.62, 0.90, 0.86)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.032)

    for family, rows in sorted(grouped_rows.items()):
        rows = sorted(rows, key=lambda row: row["mass"])
        color, marker = FAMILY_STYLES.get(family, (ROOT.kBlack, 20))
        graph = make_graph(rows, args.compareKey, color, marker, FAMILY_LABELS.get(family, family))
        graph.Draw("LP same")
        legend.AddEntry(graph, FAMILY_LABELS.get(family, family), "lp")
        _DRAWN_OBJECTS.append(graph)

    legend.Draw()
    _DRAWN_OBJECTS.append(legend)
    draw_note(
        [
            "Selection:",
            case_label,
            "VBF Selection applied",
        ]
    )
    draw_label(args.label)

    outbase = os.path.join(
        outdir,
        f"VBF_vs_NoVBF_{sanitize_token(case_label)}",
    )
    canvas.SaveAs(outbase + ".png")
    canvas.SaveAs(outbase + ".pdf")


def save_family_selection_overlay_plot(family, case_rows_map, outdir):
    _DRAWN_OBJECTS.clear()
    canvas = ROOT.TCanvas(f"c_{family}_cases", "", 1000, 800)
    masses = sorted({row["mass"] for rows in case_rows_map.values() for row in rows})
    if not masses:
        return

    ymax = max(row[args.compareKey] for rows in case_rows_map.values() for row in rows)
    frame = make_y_frame(f"frame_{family}_cases", masses, ymax)
    frame.Draw()
    _DRAWN_OBJECTS.append(frame)

    legend = ROOT.TLegend(0.52, 0.60, 0.90, 0.86)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.032)

    for idx, (case_label, rows) in enumerate(sorted(case_rows_map.items())):
        rows = sorted(rows, key=lambda row: row["mass"])
        color, marker = CASE_STYLES[idx % len(CASE_STYLES)]
        graph = make_graph(rows, args.compareKey, color, marker, case_label)
        graph.Draw("LP same")
        legend.AddEntry(graph, case_label, "lp")
        _DRAWN_OBJECTS.append(graph)

    legend.Draw()
    _DRAWN_OBJECTS.append(legend)
    draw_note(
        [
            "Signal model:",
            FAMILY_LABELS.get(family, family),
            "VBF Selection applied",
        ]
    )
    draw_label(args.label)

    outbase = os.path.join(outdir, f"VBF_vs_NoVBF_{sanitize_token(family)}")
    canvas.SaveAs(outbase + ".png")
    canvas.SaveAs(outbase + ".pdf")


def main():
    cms_style()
    os.makedirs(args.outdir, exist_ok=True)

    csv_paths = expand_csv_inputs(args.csv)
    all_rows = []
    for path in csv_paths:
        all_rows.extend(read_rows(path))

    all_rows = deduplicate_rows(all_rows)
    if args.families:
        all_rows = [row for row in all_rows if row["sample"] in args.families]

    if not all_rows:
        raise RuntimeError("No rows left after applying the family and/or working-point selection.")

    grouped_by_case = defaultdict(lambda: defaultdict(list))
    grouped_by_family = defaultdict(lambda: defaultdict(list))
    for row in all_rows:
        grouped_by_case[row["case_label"]][row["sample"]].append(row)
        grouped_by_family[row["sample"]][row["case_label"]].append(row)

    for case_label, grouped_rows in sorted(grouped_by_case.items()):
        save_working_point_signal_overlay_plot(case_label, grouped_rows, args.outdir)

    for family, case_rows_map in sorted(grouped_by_family.items()):
        save_family_selection_overlay_plot(family, case_rows_map, args.outdir)

    print(f"[INFO] Saved plots in: {args.outdir}")


if __name__ == "__main__":
    main()
