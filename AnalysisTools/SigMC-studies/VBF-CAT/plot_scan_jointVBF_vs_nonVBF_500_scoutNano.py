#!/usr/bin/env python3

"""
Plot the joint 500 GeV VBF vs non-VBF forward-cut scan.

Input:
 - CSV from scan_jointVBF_vs_nonVBF_500_scoutNano.py

Outputs:
 - tradeoff scatter: VBF efficiency vs non-VBF efficiency
 - heatmaps per forward-jet pT threshold
"""

import os
import csv
import math
import argparse
from collections import defaultdict

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


parser = argparse.ArgumentParser(description="Plot the joint 500 GeV VBF vs non-VBF forward-cut scan")
parser.add_argument("--csv", required=True, help="Input CSV from scan_jointVBF_vs_nonVBF_500_scoutNano.py")
parser.add_argument("--outdir", default="plots_scan_jointVBF_vs_nonVBF_500", help="Output directory")
parser.add_argument("--label", default="500 GeV ScoutNano optimization", help="Text label")
parser.add_argument("--topN", type=int, default=12, help="Number of best points to annotate in the scatter")
parser.add_argument("--paretoCsv", default=None, help="Optional Pareto CSV from scan_jointVBF_vs_nonVBF_500_scoutNano.py")
args = parser.parse_args()


def cms_style():
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.12)
    ROOT.gStyle.SetPadBottomMargin(0.14)
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleOffset(1.10, "X")
    ROOT.gStyle.SetTitleOffset(1.15, "Y")


def draw_label(extra_text):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.12, 0.905, "CMS")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextFont(52)
    latex2.SetTextSize(0.036)
    latex2.DrawLatex(0.185, 0.905, "Simulation Preliminary")

    latex3 = ROOT.TLatex()
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.036)
    latex3.SetTextAlign(31)
    latex3.DrawLatex(0.95, 0.905, extra_text)


def draw_note(lines, x1=0.15, y1=0.72, x2=0.47, y2=0.86):
    note = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    note.SetFillColor(ROOT.kGray)
    note.SetFillColorAlpha(ROOT.kGray, 0.12)
    note.SetLineColor(ROOT.kGray + 1)
    note.SetLineWidth(2)
    note.SetTextAlign(12)
    note.SetTextFont(42)
    note.SetTextSize(0.028)
    for line in lines:
        note.AddText(line)
    note.Draw()
    return note


def read_rows(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            for key in [
                "forwardJetPtMin",
                "forwardPairDEtaMin",
                "forwardPairExtraMjjMin",
                "vbf_eff_selected_pct",
                "nonvbf_eff_selected_pct",
                "eff_diff_pct",
                "eff_ratio",
                "score",
                "significance_like",
                "purity_like",
            ]:
                if key in row and row[key] != "":
                    row[key] = float(row[key])
            rows.append(row)
    return rows


def pareto_frontier(rows):
    frontier = []
    for row in rows:
        dominated = False
        for other in rows:
            if other is row:
                continue
            if (
                other["vbf_eff_selected_pct"] >= row["vbf_eff_selected_pct"]
                and other["nonvbf_eff_selected_pct"] <= row["nonvbf_eff_selected_pct"]
                and (
                    other["vbf_eff_selected_pct"] > row["vbf_eff_selected_pct"]
                    or other["nonvbf_eff_selected_pct"] < row["nonvbf_eff_selected_pct"]
                )
            ):
                dominated = True
                break
        if not dominated:
            frontier.append(row)
    return sorted(frontier, key=lambda r: r["vbf_eff_selected_pct"])


def save_scatter(rows, outdir):
    canvas = ROOT.TCanvas("c_joint_tradeoff", "", 1450, 900)
    canvas.SetLeftMargin(0.11)
    canvas.SetRightMargin(0.34)
    canvas.SetBottomMargin(0.14)
    canvas.SetTopMargin(0.10)
    frame = ROOT.TH1F(
        "frame_joint_tradeoff",
        ";VBF selected efficiency [%];non-VBF selected efficiency [%]",
        100,
        0.0,
        max(1.0, 1.10 * max(row["vbf_eff_selected_pct"] for row in rows)),
    )
    frame.SetMinimum(0.0)
    frame.SetMaximum(max(1.0, 1.15 * max(row["nonvbf_eff_selected_pct"] for row in rows)))
    frame.Draw()

    graph = ROOT.TGraph(len(rows))
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.1)
    graph.SetMarkerColor(ROOT.TColor.GetColor("#3f90da"))
    graph.SetLineColor(ROOT.TColor.GetColor("#3f90da"))
    for idx, row in enumerate(rows):
        graph.SetPoint(idx, row["vbf_eff_selected_pct"], row["nonvbf_eff_selected_pct"])
    graph.Draw("P same")

    frontier_rows = pareto_frontier(rows)
    frontier = ROOT.TGraph(len(frontier_rows))
    frontier.SetLineColor(ROOT.TColor.GetColor("#bd1f01"))
    frontier.SetLineWidth(3)
    frontier.SetMarkerStyle(24)
    frontier.SetMarkerSize(1.2)
    frontier.SetMarkerColor(ROOT.TColor.GetColor("#bd1f01"))
    for idx, row in enumerate(frontier_rows):
        frontier.SetPoint(idx, row["vbf_eff_selected_pct"], row["nonvbf_eff_selected_pct"])
    frontier.Draw("LP same")

    top_rows = sorted(rows, key=lambda r: r["score"], reverse=True)[:args.topN]
    top_graph = ROOT.TGraph(len(top_rows))
    top_graph.SetMarkerStyle(29)
    top_graph.SetMarkerSize(1.6)
    top_graph.SetMarkerColor(ROOT.kBlack)
    top_graph.SetLineColor(ROOT.kBlack)
    for idx, row in enumerate(top_rows):
        top_graph.SetPoint(idx, row["vbf_eff_selected_pct"], row["nonvbf_eff_selected_pct"])
    top_graph.Draw("P same")

    legend = ROOT.TLegend(0.69, 0.70, 0.96, 0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.030)
    legend.AddEntry(graph, "All working points", "p")
    legend.AddEntry(frontier, "Best trade-off", "lp")
    legend.AddEntry(top_graph, f"Top {len(top_rows)} score points", "p")
    legend.Draw()

    shortlist = ROOT.TPaveText(0.69, 0.20, 0.97, 0.66, "NDC")
    shortlist.SetFillColor(ROOT.kGray)
    shortlist.SetFillColorAlpha(ROOT.kGray, 0.10)
    shortlist.SetLineColor(ROOT.kGray + 1)
    shortlist.SetLineWidth(2)
    shortlist.SetTextAlign(12)
    shortlist.SetTextFont(42)
    shortlist.SetTextSize(0.024)
    header = shortlist.AddText(f"Top {len(top_rows)} score points")
    header.SetTextFont(62)
    shortlist.AddText("(fwdPt, dEta, extraMjj)")
    for idx, row in enumerate(top_rows, start=1):
        shortlist.AddText(
            f"{idx}. ({row['forwardJetPtMin']:.0f}, {row['forwardPairDEtaMin']:.2f}, {row['forwardPairExtraMjjMin']:.0f})"
        )
    shortlist.Draw()

    draw_label(args.label)
    canvas.SaveAs(os.path.join(outdir, "joint_tradeoff_scatter.png"))
    canvas.SaveAs(os.path.join(outdir, "joint_tradeoff_scatter.pdf"))


def build_heatmap(rows, zkey, title):
    xvals = sorted({row["forwardPairDEtaMin"] for row in rows})
    yvals = sorted({row["forwardPairExtraMjjMin"] for row in rows})
    if len(xvals) < 1 or len(yvals) < 1:
        return None

    xstep = xvals[1] - xvals[0] if len(xvals) > 1 else 0.25
    ystep = yvals[1] - yvals[0] if len(yvals) > 1 else 100.0
    hist = ROOT.TH2F(
        f"h_{zkey}_{int(rows[0]['forwardJetPtMin'])}",
        title,
        len(xvals),
        min(xvals) - 0.5 * xstep,
        max(xvals) + 0.5 * xstep,
        len(yvals),
        min(yvals) - 0.5 * ystep,
        max(yvals) + 0.5 * ystep,
    )
    for row in rows:
        hist.Fill(row["forwardPairDEtaMin"], row["forwardPairExtraMjjMin"], row[zkey])
    return hist


def save_heatmaps(rows, outdir):
    rows_by_fwdpt = defaultdict(list)
    for row in rows:
        rows_by_fwdpt[row["forwardJetPtMin"]].append(row)

    metrics = [
        ("vbf_eff_selected_pct", "VBF efficiency [%]"),
        ("nonvbf_eff_selected_pct", "non-VBF selected efficiency [%]"),
        ("score", "Score = eff_{VBF} - eff_{nonVBF}"),
        ("purity_like", "Purity-like = eff_{VBF} / (eff_{VBF} + eff_{nonVBF})"),
    ]

    for fwdpt, pt_rows in sorted(rows_by_fwdpt.items()):
        for key, ztitle in metrics:
            canvas = ROOT.TCanvas(f"c_{key}_{int(fwdpt)}", "", 1150, 850)
            canvas.SetLeftMargin(0.12)
            canvas.SetRightMargin(0.20)
            canvas.SetBottomMargin(0.14)
            canvas.SetTopMargin(0.10)
            hist = build_heatmap(
                pt_rows,
                key,
                f";Forward pair |#Delta#eta| min;Extra forward-pair mjj min [GeV];{ztitle}",
            )
            if hist is None:
                continue
            hist.SetContour(50)
            hist.GetZaxis().SetTitleOffset(1.45)
            hist.GetZaxis().SetTitleSize(0.042)
            hist.GetZaxis().SetLabelSize(0.033)
            hist.Draw("COLZ TEXT")

            best_row = max(pt_rows, key=lambda row: row["score"])
            marker = ROOT.TMarker(best_row["forwardPairDEtaMin"], best_row["forwardPairExtraMjjMin"], 29)
            marker.SetMarkerColor(ROOT.kBlack)
            marker.SetMarkerSize(2.0)
            marker.Draw()

            draw_note(
                [
                    f"forwardJetPtMin > {fwdpt:.0f} GeV",
                    f"Best score: {best_row['score']:.2f}",
                    f"Best point: dEta>{best_row['forwardPairDEtaMin']:.2f}, mjj>{best_row['forwardPairExtraMjjMin']:.0f}",
                ]
            )
            draw_label(args.label)
            outbase = os.path.join(outdir, f"heatmap_{key}_fwdPt{int(round(fwdpt))}")
            canvas.SaveAs(outbase + ".png")
            canvas.SaveAs(outbase + ".pdf")


def write_recommendation_outputs(rows, pareto_rows, outdir):
    top_score = sorted(rows, key=lambda row: row["score"], reverse=True)[:10]
    top_purity = sorted(rows, key=lambda row: row["purity_like"], reverse=True)[:10]
    top_rejection = sorted(rows, key=lambda row: row["nonvbf_eff_selected_pct"])[:10]

    csv_path = os.path.join(outdir, "recommended_working_points.csv")
    fieldnames = [
        "category",
        "working_point",
        "forwardJetPtMin",
        "forwardPairDEtaMin",
        "forwardPairExtraMjjMin",
        "vbf_eff_selected_pct",
        "nonvbf_eff_selected_pct",
        "eff_diff_pct",
        "score",
        "significance_like",
        "purity_like",
        "motivation",
    ]
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for category, subset, motivation in [
            ("best_score", top_score, "Best balance of high VBF efficiency and low non-VBF contamination for the chosen linear score."),
            ("best_purity", top_purity, "Cleanest selected sample according to VBF% / (VBF% + NoVBF%)."),
            ("best_rejection", top_rejection, "Lowest non-VBF mistag efficiency, typically with a stronger cost in VBF efficiency."),
            ("pareto", pareto_rows[:20], "Non-dominated shortlist: no other point is simultaneously better in both VBF efficiency and non-VBF rejection."),
        ]:
            for row in subset:
                out_row = {field: row.get(field, "") for field in fieldnames}
                out_row["category"] = category
                out_row["motivation"] = motivation
                writer.writerow(out_row)

    txt_path = os.path.join(outdir, "selection_recommendation.txt")
    with open(txt_path, "w") as handle:
        handle.write("Joint VBF vs non-VBF scan summary\n")
        handle.write("=================================\n\n")
        handle.write(f"Label: {args.label}\n")
        handle.write("Score definitions:\n")
        handle.write(" - score = VBF% - NoVBF% by default, or more generally VBF% - penalty * NoVBF%\n")
        handle.write(" - significance_like = VBF% / sqrt(NoVBF%)\n")
        handle.write(" - purity_like = VBF% / (VBF% + NoVBF%)\n\n")
        handle.write("Interpretation:\n")
        handle.write(" - score is the most direct compromise metric if you want to keep VBF efficiency high while explicitly penalizing contamination.\n")
        handle.write(" - significance_like is useful when you care more strongly about suppressing non-VBF leakage.\n")
        handle.write(" - purity_like is useful when you want to argue that the selected category is cleaner, even if efficiency drops.\n\n")
        handle.write("Recommended points:\n")
        for title, row, note in [
            ("Best by score", top_score[0], "Balanced default choice."),
            ("Best by purity_like", top_purity[0], "Best if category cleanliness is the main priority."),
            ("Best by non-VBF rejection", top_rejection[0], "Best if minimizing mistag is the main priority."),
        ]:
            handle.write(
                f" - {title}: fwdPt>{row['forwardJetPtMin']:.0f} GeV, dEta>{row['forwardPairDEtaMin']:.2f}, extraMjj>{row['forwardPairExtraMjjMin']:.0f} GeV"
                f" | VBF%={row['vbf_eff_selected_pct']:.3f}, NoVBF%={row['nonvbf_eff_selected_pct']:.3f}, Diff%={row['eff_diff_pct']:.3f},"
                f" Score={row['score']:.3f}, PurityLike={row['purity_like']:.4f}. {note}\n"
            )
        handle.write("\nPareto guidance:\n")
        handle.write(" - Use the Pareto shortlist when you want to show there is no strictly better point in both signal efficiency and contamination.\n")
        handle.write(" - A practical motivation is to choose the Pareto point nearest your preferred operating region, then cite the scatter and heatmaps.\n")


def print_top(rows):
    print(
        f"{'Rank':>4s} {'fwdPt':>7s} {'dEta':>7s} {'xMjj':>7s} "
        f"{'VBF%':>8s} {'NoVBF%':>8s} {'Score':>8s}"
    )
    for idx, row in enumerate(sorted(rows, key=lambda r: r["score"], reverse=True)[:15], start=1):
        print(
            f"{idx:4d} {row['forwardJetPtMin']:7.1f} {row['forwardPairDEtaMin']:7.2f} "
            f"{row['forwardPairExtraMjjMin']:7.1f} {row['vbf_eff_selected_pct']:8.3f} "
            f"{row['nonvbf_eff_selected_pct']:8.3f} {row['score']:8.3f}"
        )


def main():
    cms_style()
    os.makedirs(args.outdir, exist_ok=True)
    rows = read_rows(args.csv)
    if not rows:
        raise RuntimeError("No rows found in the input CSV.")
    pareto_rows = read_rows(args.paretoCsv) if args.paretoCsv else pareto_frontier(rows)
    save_scatter(rows, args.outdir)
    save_heatmaps(rows, args.outdir)
    write_recommendation_outputs(rows, pareto_rows, args.outdir)
    print_top(rows)
    print(f"[INFO] Saved recommendation table: {os.path.join(args.outdir, 'recommended_working_points.csv')}")
    print(f"[INFO] Saved recommendation note:  {os.path.join(args.outdir, 'selection_recommendation.txt')}")
    print(f"[INFO] Saved plots in: {args.outdir}")


if __name__ == "__main__":
    main()
