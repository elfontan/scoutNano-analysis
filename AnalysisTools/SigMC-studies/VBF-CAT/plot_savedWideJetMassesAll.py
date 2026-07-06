#!/usr/bin/env python3

"""
Overlay wide-jet mass distributions across signal mass points.

Supported inputs:
 - legacy `signal_peak_histograms.root` files produced by
   plot_centralWideJetVBF_scoutNano.py
 - current `CentralVBFHTo2B_M*_centralWideJetVBF_scoutNano_<WP>.root` files,
   auto-discovered for a chosen working point

In tree-input mode, the script overlays the VBF-tagged wide-jet mass
(`pass_trigger_baseline == 1 && isVBF == 1`) and optionally the matched subset
(`central_pair_genmatched == 1`) using the same multi-mass style.
"""

import argparse
import glob
import os
import re

import ROOT

ROOT.gROOT.SetBatch(True)


parser = argparse.ArgumentParser(description="Overlay saved wide-jet mass histograms by strategy and signal mass")
parser.add_argument("--input", "-i", type=str, action="append", default=[], help="One signal_peak_histograms.root file for a signal mass point. Repeat for multiple masses.")
parser.add_argument("--massLabel", type=str, action="append", default=[], help="Optional mass label for each repeated --input entry, e.g. M300. If omitted, the mass label is inferred from the filename.")
parser.add_argument("--input300", type=str, help="Legacy signal_peak_histograms.root for M=300")
parser.add_argument("--input500", type=str, help="Legacy signal_peak_histograms.root for M=500")
parser.add_argument("--input750", type=str, help="Legacy signal_peak_histograms.root for M=750")
parser.add_argument("--input1000", type=str, help="Legacy signal_peak_histograms.root for M=1000")
parser.add_argument("--input2000", type=str, help="Legacy signal_peak_histograms.root for M=2000")
parser.add_argument("--input3000", type=str, help="Legacy signal_peak_histograms.root for M=3000")
parser.add_argument("--outdir", type=str, default="plots_savedWideJetMasses_byStrategy", help="Output directory")
parser.add_argument("--label", type=str, default="ScoutNano VBF study", help="Text label drawn on plots")
parser.add_argument("--plotMatched", action="store_true", help="Also make matched-only overlays")
parser.add_argument("--norm", action="store_true", help="Normalize each histogram to unit area before plotting.")
parser.add_argument("--xmax", type=float, default=None, help="Maximum x-axis value [GeV]. If omitted, choose it adaptively from the largest mass point.")
parser.add_argument("--workingPoint", type=str, default=None, help="Current tree-output mode: working-point suffix to auto-discover, e.g. fwdPt24_VBF-Deta4-Mjj500")
parser.add_argument("--treeDir", type=str, default=os.path.dirname(os.path.abspath(__file__)), help="Directory containing CentralVBFHTo2B_M*_centralWideJetVBF_scoutNano_*.root files")
parser.add_argument("--tree", type=str, default="Events", help="Tree name for --workingPoint mode")
args = parser.parse_args()


_DRAWN_OBJECTS = []


def cms_style():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleOffset(1.05, "X")
    ROOT.gStyle.SetTitleOffset(1.35, "Y")
    ROOT.gStyle.SetEndErrorSize(0)
    ROOT.gStyle.SetLegendFont(42)


def hex_to_rootcolor(hexcode):
    hexcode = hexcode.lstrip("#")
    r = int(hexcode[0:2], 16) / 255.0
    g = int(hexcode[2:4], 16) / 255.0
    b = int(hexcode[4:6], 16) / 255.0
    return ROOT.TColor.GetColor(r, g, b)


def draw_label(extra_text):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.05)
    cms_text = latex.DrawLatex(0.132, 0.932, "CMS")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextFont(52)
    latex2.SetTextSize(0.04)
    sim_text = latex2.DrawLatex(0.235, 0.932, "Simulation Preliminary")

    latex3 = ROOT.TLatex()
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.04)
    latex3.SetTextAlign(31)
    extra_label = latex3.DrawLatex(0.95, 0.932, extra_text)
    _DRAWN_OBJECTS.extend([latex, latex2, latex3, cms_text, sim_text, extra_label])


def draw_note(selection_text, requirements_text, x=0.16, y=0.84):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(13)
    drawn = [latex]

    latex.SetTextSize(0.032)
    selection = latex.DrawLatex(x, y, f"#bf{{#it{{Selection: VBF ({selection_text}}}}})")
    drawn.append(selection)

    latex.SetTextSize(0.027)
    requirements = latex.DrawLatex(
        x,
        y - 0.050,
        f"#it{{{requirements_text}}}",
    )
    drawn.append(requirements)
    _DRAWN_OBJECTS.extend(drawn)


def load_hist(root_file, hist_name):
    hist = root_file.Get(hist_name)
    if not hist:
        raise RuntimeError(f"Could not find histogram '{hist_name}' in {root_file.GetName()}")
    clone = hist.Clone(f"{hist_name}_{os.path.basename(root_file.GetName()).replace('.', '_')}")
    clone.SetDirectory(0)
    return clone


def infer_mass_label(path):
    base = os.path.basename(path)
    tokens = base.replace('.', '_').replace('-', '_').split('_')
    for token in tokens:
        if token.upper().startswith('M') and token[1:].isdigit():
            return token.upper()
        if token.isdigit():
            return f"M{token}"
    import re

    match = re.search(r'(\d{3,4})', base)
    if match:
        return f"M{match.group(1)}"
    raise RuntimeError(f"Could not infer mass label from input path: {path}")


def make_hist_from_tree(tree, name, expr, selection, nbins, xmin, xmax, xtitle, ytitle="Events / 5 GeV"):
    hist = ROOT.TH1F(name, f";{xtitle};{ytitle}", nbins, xmin, xmax)
    hist.Sumw2()
    tree.Draw(f"{expr}>>{name}", selection, "goff")
    hist.SetDirectory(0)
    return hist


def build_input_files(args):
    file_inputs = {}
    for mass in [300, 500, 750, 1000, 2000, 3000]:
        arg_name = f"input{mass}"
        path = getattr(args, arg_name, None)
        if path:
            file_inputs[f"M{mass}"] = path

    for idx, path in enumerate(args.input):
        if idx < len(args.massLabel) and args.massLabel[idx]:
            label = args.massLabel[idx].upper()
            if not label.startswith('M'):
                label = f"M{label}"
        else:
            label = infer_mass_label(path)
        if label in file_inputs:
            raise RuntimeError(f"Duplicate input provided for {label}.")
        file_inputs[label] = path

    if not file_inputs:
        raise RuntimeError("No input files provided. Use --input/-i repeatedly or legacy --input300/500/750/1000/2000/3000 options.")

    def mass_value(label):
        return int(''.join(filter(str.isdigit, label)))

    return sorted(file_inputs.items(), key=lambda item: mass_value(item[0]))


def discover_tree_working_point_files(args):
    pattern = os.path.join(
        args.treeDir,
        f"CentralVBFHTo2B_M*_centralWideJetVBF_scoutNano_*{args.workingPoint}.root",
    )
    file_inputs = {}
    for path in sorted(glob.glob(pattern)):
        base = os.path.basename(path)
        match = re.search(r"_M(\d+)_", base)
        if not match:
            continue
        selection_suffix = base.split("centralWideJetVBF_scoutNano_", 1)[1].rsplit(".root", 1)[0]
        if not selection_suffix.endswith(args.workingPoint):
            continue
        label = f"M{int(match.group(1))}"
        if label in file_inputs:
            raise RuntimeError(f"Duplicate tree input provided for {label}.")
        file_inputs[label] = path

    if not file_inputs:
        raise RuntimeError(f"No tree files found for working point '{args.workingPoint}' in {args.treeDir}")

    return sorted(file_inputs.items(), key=lambda item: int(item[0].replace("M", "")))


def style_hist(hist, color, line_style=1):
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetLineWidth(3)
    hist.SetLineStyle(line_style)


def normalize(hist):
    if hist.Integral() > 0:
        hist.Scale(1.0 / hist.Integral())


def compute_xmax(mass_values):
    if args.xmax is not None:
        return args.xmax
    if not mass_values:
        return 2000.0
    return max(1600.0, 2.0 * max(mass_values))


def save_mass_overlay(hists, labels, outbase, selection_text, legend_title=None):
    _DRAWN_OBJECTS.clear()
    canvas = ROOT.TCanvas(f"c_{os.path.basename(outbase)}", "", 1000, 800)
    legend = ROOT.TLegend(0.6, 0.40, 0.90, 0.7)
    legend.SetTextSize(0.035)
    legend.SetFillStyle(0)
    if legend_title:
        legend.SetHeader(legend_title, "C")

    for hist in hists:
        if args.norm:
            normalize(hist)

    ymax = max(hist.GetMaximum() for hist in hists) if hists else 0.0
    ytitle = "Normalized events" if args.norm else "Events / 5 GeV"
    frame = ROOT.TH1F(f"frame_{os.path.basename(outbase)}", f";Wide-jet dijet mass [GeV];{ytitle}", 100, 0.0, args.xmax)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.35 * ymax)
    frame.Draw()
    _DRAWN_OBJECTS.append(frame)

    for hist in hists:
        hist.SetMaximum(1.35 * ymax)
        hist.SetMinimum(0.0)

    for idx, hist in enumerate(hists):
        hist.Draw("hist same")
        legend.AddEntry(hist, labels[idx], "l")

    legend.Draw()
    draw_note(
        selection_text=selection_text,
        requirements_text="Central wide-jets (#DeltaR = 1.1; |#eta| < 2.5; p_{T} > 30 GeV)",
    )
    draw_label(args.label)
    canvas.SaveAs(outbase + ".png")
    canvas.SaveAs(outbase + ".pdf")


def open_root(path):
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {path}")
    return root_file


def main():
    cms_style()
    os.makedirs(args.outdir, exist_ok=True)

    tree_mode = bool(args.workingPoint)
    file_entries = discover_tree_working_point_files(args) if tree_mode else build_input_files(args)
    colors_palette = [
        "#3f90da",
        "#ffa90e",
        "#bd1f01",
        "#2ca02c",
        "#9467bd",
        "#d62728",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]

    files = {label: open_root(path) for label, path in file_entries}
    colors = {label: hex_to_rootcolor(colors_palette[idx % len(colors_palette)]) for idx, (label, _) in enumerate(file_entries)}

    selected_hists = []
    matched_hists = []
    central_final_hists = []
    inverted_final_hists = []
    central_matched_hists = []
    inverted_matched_hists = []
    labels = []
    mass_values = []
    mass_keys = [label for label, _ in file_entries]

    for mass_label in mass_keys:
        root_file = files[mass_label]
        mass_value = float(mass_label.replace("M", ""))
        mass_values.append(mass_value)
        labels.append(mass_label.replace("M", "m_{H} = ") + " GeV")

        if tree_mode:
            tree = root_file.Get(args.tree)
            if tree is None:
                raise RuntimeError(f"Could not find tree '{args.tree}' in {root_file.GetName()}")
            xmax_guess = max(1600.0, 2.0 * mass_value)
            nbins = int(round(xmax_guess / 5.0))

            h_selected = make_hist_from_tree(
                tree,
                f"h_widejet_peak_isVBF_{mass_label}",
                "widejet_pair_mass",
                "pass_trigger_baseline == 1 && isVBF == 1",
                nbins,
                0.0,
                xmax_guess,
                "Wide-jet dijet mass [GeV]",
            )
            style_hist(h_selected, colors[mass_label])
            selected_hists.append(h_selected)

            if args.plotMatched:
                h_matched = make_hist_from_tree(
                    tree,
                    f"h_widejet_peak_isVBF_matched_{mass_label}",
                    "widejet_pair_mass",
                    "pass_trigger_baseline == 1 && isVBF == 1 && central_pair_genmatched == 1",
                    nbins,
                    0.0,
                    xmax_guess,
                    "Wide-jet dijet mass [GeV]",
                )
                style_hist(h_matched, colors[mass_label])
                matched_hists.append(h_matched)
        else:
            h_central_final = load_hist(root_file, "h_widejet_peak_centralFirst_final_5GeV")
            style_hist(h_central_final, colors[mass_label])
            central_final_hists.append(h_central_final)

            h_inverted_final = load_hist(root_file, "h_widejet_peak_inverted_final_5GeV")
            style_hist(h_inverted_final, colors[mass_label])
            inverted_final_hists.append(h_inverted_final)

            if args.plotMatched:
                h_central_matched = load_hist(root_file, "h_widejet_peak_centralFirst_matched_5GeV")
                style_hist(h_central_matched, colors[mass_label])
                central_matched_hists.append(h_central_matched)

                h_inverted_matched = load_hist(root_file, "h_widejet_peak_inverted_matched_5GeV")
                style_hist(h_inverted_matched, colors[mass_label])
                inverted_matched_hists.append(h_inverted_matched)

    args.xmax = compute_xmax(mass_values)
    mass_key = "_".join(mass_keys)

    if tree_mode:
        wp_tag = args.workingPoint.replace("/", "_")
        save_mass_overlay(
            selected_hists,
            labels,
            os.path.join(args.outdir, f"widejetMass_isVBF_{wp_tag}_{mass_key}"),
            args.workingPoint,
            legend_title="VBFHTo2B",
        )
        if args.plotMatched:
            save_mass_overlay(
                matched_hists,
                labels,
                os.path.join(args.outdir, f"widejetMass_isVBF_matched_{wp_tag}_{mass_key}"),
                args.workingPoint,
                legend_title="VBFHTo2B (matched)",
            )
    else:
        save_mass_overlay(
            central_final_hists,
            labels,
            os.path.join(args.outdir, f"widejetMass_centralFirst_final_{mass_key}"),
            "Central-to-VBFTag",
            legend_title="Selected events",
        )
        save_mass_overlay(
            inverted_final_hists,
            labels,
            os.path.join(args.outdir, f"widejetMass_inverted_final_{mass_key}"),
            "VBFTag-to-Central",
            legend_title="Selected events",
        )

        if args.plotMatched:
            save_mass_overlay(
                central_matched_hists,
                labels,
                os.path.join(args.outdir, f"widejetMass_centralFirst_matched_{mass_key}"),
                "Central-to-VBFTag",
                legend_title="Selected events (matched)",
            )
            save_mass_overlay(
                inverted_matched_hists,
                labels,
                os.path.join(args.outdir, f"widejetMass_inverted_matched_{mass_key}"),
                "VBFTag-to-Central",
                legend_title="Selected events (matched)",
            )

    for root_file in files.values():
        root_file.Close()

    print(f"[INFO] Saved overlays in: {args.outdir}")


if __name__ == "__main__":
    main()
