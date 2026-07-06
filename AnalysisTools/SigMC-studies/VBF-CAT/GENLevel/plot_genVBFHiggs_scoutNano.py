#!/usr/bin/env python3

"""
Gen-level summary plotter for the ScoutNano Higgs study tree.

It produces:
 - generated Higgs kinematics
 - decay-product comparisons before/after minimal acceptance cuts
 - gen-jet multiplicity comparisons
 - forward gen-jet pair properties
 - a simple VBF acceptance cutflow based on forward gen-jets
"""

import os
import argparse

import ROOT

ROOT.gROOT.SetBatch(True)


parser = argparse.ArgumentParser(description="Plot gen-level ScoutNano Higgs summary")
parser.add_argument(
    "--input",
    type=str,
    default="VBFHTo2B_M300_genHiggs_scoutNano.root",
    help="Input ROOT file containing the Events tree",
)
parser.add_argument(
    "--outdir",
    type=str,
    default="plots_genVBFHiggs_scoutNano",
    help="Output directory for plots",
)
parser.add_argument(
    "--label",
    type=str,
    default="VBFHTo2B M300 (2024)",
    help="Label shown on plots",
)
parser.add_argument(
    "--forwardEtaMin",
    type=float,
    default=3.0,
    help="Minimum |eta| defining a forward gen jet for the cutflow labels.",
)
parser.add_argument(
    "--forwardEtaMax",
    type=float,
    default=5.0,
    help="Maximum |eta| defining a forward gen jet for the cutflow labels.",
)
parser.add_argument(
    "--forwardPairDEtaMin",
    type=float,
    default=5.0,
    help="Minimum |DeltaEta| threshold used by the dEta-based VBF tag modes.",
)
parser.add_argument(
    "--decayPtAcceptanceMin",
    type=float,
    default=30.0,
    help="Minimum pT [GeV] used for the Gen-level Higgs-daughter acceptance selection.",
)
parser.add_argument(
    "--genJetLheMatchDRMax",
    type=float,
    default=0.4,
    help="Maximum DeltaR used to define a successful GenJet-to-LHE VBF pair match.",
)
parser.add_argument(
    "--vbfTagMode",
    type=str,
    default="twoForward_maxDEta_dEtaGtMin",
    choices=[
        "twoForward_maxDEta",
        "twoForward_maxDEta_dEtaGtMin",
        "oneForward_maxDEta",
        "oneForward_maxDEta_dEtaGtMin",
        "oneForward_maxMjj",
        "maxDEta",
        "maxMjj",
    ],
    help="VBF tag mode used when producing the input tree.",
)
args = parser.parse_args()


def split_output_dirs(outdir, vbf_tag_mode):
    expected_leaf = f"VBFTagMode-{vbf_tag_mode}"
    normalized = os.path.normpath(outdir)
    if os.path.basename(normalized) == expected_leaf:
        return os.path.dirname(normalized), normalized
    return normalized, os.path.join(normalized, expected_leaf)


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
    ROOT.gStyle.SetTitleOffset(1.1, "X")
    ROOT.gStyle.SetTitleOffset(1.3, "Y")


def draw_label(extra_text, cms_x=0.14, cms_y=0.93, extra_x=0.95, extra_y=0.93):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.05)
    latex.DrawLatex(cms_x, cms_y, "CMS")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextFont(52)
    latex2.SetTextSize(0.04)
    latex2.DrawLatex(cms_x + 0.10, cms_y, "Simulation Preliminary")

    latex3 = ROOT.TLatex()
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.04)
    latex3.SetTextAlign(31)
    latex3.DrawLatex(extra_x, extra_y, extra_text)


def make_hist(tree, name, expr, selection, nbins, xmin, xmax, xtitle, ytitle="Events"):
    hist = ROOT.TH1F(name, f";{xtitle};{ytitle}", nbins, xmin, xmax)
    hist.Sumw2()
    draw_expr = f"{expr}>>{name}"
    tree.Draw(draw_expr, selection, "goff")
    return hist


def make_hist_from_expressions(tree, name, expressions, selection, nbins, xmin, xmax, xtitle, ytitle="Events"):
    hist = ROOT.TH1F(name, f";{xtitle};{ytitle}", nbins, xmin, xmax)
    hist.Sumw2()
    for idx, expr in enumerate(expressions):
        draw_expr = f"{expr}>>{name}" if idx == 0 else f"{expr}>>+{name}"
        tree.Draw(draw_expr, selection, "goff")
    return hist


def style_hist(hist, color):
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetLineWidth(3)


def style_single_filled(hist, line_color, fill_color, fill_alpha=0.35):
    hist.SetLineColor(line_color)
    hist.SetMarkerColor(line_color)
    hist.SetLineWidth(3)
    hist.SetFillColorAlpha(fill_color, fill_alpha)


def style_overlay_line(hist, color, line_style=1, line_width=3):
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetLineStyle(line_style)
    hist.SetLineWidth(line_width)
    hist.SetFillStyle(0)


def style_overlay_filled(hist, line_color, fill_color, fill_alpha=0.30, line_width=3):
    hist.SetLineColor(line_color)
    hist.SetMarkerColor(line_color)
    hist.SetLineWidth(line_width)
    hist.SetFillColorAlpha(fill_color, fill_alpha)


def save_single(hist, outdir, outname, label, logy=False):
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1000, 800)
    if logy:
        canvas.SetLogy()
        hist.SetMinimum(0.5)
    draw_opt = "hist"
    if hist.GetFillColor() and hist.GetFillStyle() != 0:
        draw_opt = "hist"
    hist.Draw(draw_opt)
    draw_label(label)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def save_overlay(hists, labels, outdir, outname, label, logy=False, normalize=False):
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1000, 800)
    if logy:
        canvas.SetLogy()

    legend = ROOT.TLegend(0.5, 0.72, 0.88, 0.89)
    legend.SetTextSize(0.035)
    legend.SetFillStyle(0)

    max_y = 0.0
    for hist in hists:
        if normalize and hist.Integral() > 0:
            hist.Scale(1.0 / hist.Integral())
        max_y = max(max_y, hist.GetMaximum())

    for idx, hist in enumerate(hists):
        hist.SetMaximum(max_y * (20.0 if logy else 1.35))
        if logy:
            hist.SetMinimum(1e-4 if normalize else 0.5)
        option = "hist" if idx == 0 else "hist same"
        hist.Draw(option)
        legend.AddEntry(hist, labels[idx], "l")

    legend.Draw()
    draw_label(label)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def build_cutflow_hist(counts):
    hist = ROOT.TH1F("h_vbf_cutflow", ";Selection;Events", len(counts), 0.5, len(counts) + 0.5)
    for idx, (label, value) in enumerate(counts, start=1):
        hist.SetBinContent(idx, value)
        hist.GetXaxis().SetBinLabel(idx, label)
    hist.SetFillColor(ROOT.kAzure - 9)
    hist.SetLineColor(ROOT.kAzure + 2)
    hist.SetLineWidth(2)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetXaxis().SetLabelOffset(0.02)
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleOffset(2.1)
    return hist


def save_cutflow(hist, counts, outdir, outname, label):
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1200, 1000)
    canvas.SetBottomMargin(0.24)
    canvas.SetTopMargin(0.10)
    canvas.SetRightMargin(0.06)
    hist.SetBarWidth(0.8)
    hist.SetBarOffset(0.1)
    hist.SetMaximum(hist.GetMaximum() * 1.65)
    hist.Draw("hist text0")
    hist.LabelsOption("u", "X")

    total = float(counts[0][1]) if counts else 0.0
    box = ROOT.TPaveText(0.36, 0.5, 0.91, 0.8, "NDC")
    box.SetFillColor(ROOT.kGray)
    box.SetFillColorAlpha(ROOT.kGray, 0.15)
    box.SetFillStyle(1001)
    box.SetLineColor(ROOT.kGray + 2)
    box.SetLineWidth(2)
    box.SetShadowColor(0)
    box.SetTextAlign(12)
    box.SetTextFont(42)
    box.SetTextSize(0.026)
    title = box.AddText("Gen-level study: summary")
    title.SetTextFont(62)
    for idx, (step_label, value) in enumerate(counts):
        frac = 100.0 * value / total if total > 0 else 0.0
        line = f"{step_label}: {value} ({frac:.2f}%)"
        box.AddText(line)

    box.Draw()
    draw_label(label, cms_x=0.14, cms_y=0.845, extra_x=0.945, extra_y=0.912)
    #draw_label(label, cms_x=0.14, cms_y=0.845, extra_x=0.955, extra_y=0.905)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def pair_exists_selection():
    return "forward_genjet_pair_mass > -900"


def tag_selection():
    return "has_forward_genjet_pair == 1"


def vbf_tag_step_label():
    if args.vbfTagMode == "twoForward_maxDEta":
        return "VBFTag passed: >=2 forward GenJets, max |#Delta#eta|"
    if args.vbfTagMode == "twoForward_maxDEta_dEtaGtMin":
        return (
            "VBFTag passed: >=2 forward GenJets, max |#Delta#eta|, "
            f"|#Delta#eta| > {args.forwardPairDEtaMin:g}"
        )
    if args.vbfTagMode == "oneForward_maxDEta":
        return "VBFTag passed: >=1 forward GenJet, max |#Delta#eta|"
    if args.vbfTagMode == "oneForward_maxDEta_dEtaGtMin":
        return (
            "VBFTag passed: >=1 forward GenJet, max |#Delta#eta|, "
            f"|#Delta#eta| > {args.forwardPairDEtaMin:g}"
        )
    if args.vbfTagMode == "oneForward_maxMjj":
        return "VBFTag passed: >=1 forward GenJet, max m_{jj}"
    if args.vbfTagMode == "maxDEta":
        return "VBFTag passed: max |#Delta#eta|"
    if args.vbfTagMode == "maxMjj":
        return "VBFTag passed: max m_{jj}"
    return "VBFTag passed"


cms_style()
common_outdir, mode_outdir = split_output_dirs(args.outdir, args.vbfTagMode)
os.makedirs(common_outdir, exist_ok=True)
os.makedirs(mode_outdir, exist_ok=True)

root_file = ROOT.TFile.Open(args.input)
if not root_file or root_file.IsZombie():
    raise RuntimeError(f"Could not open input file: {args.input}")

tree = root_file.Get("Events")
if tree is None:
    raise RuntimeError(f"Could not find tree 'Events' in {args.input}")

base_higgs = "higgs_mass > -900"
base_decay = "decay_mass > -900"
accepted_decay = (
    "decay_mass > -900"
    f" && lead_decay_pt > {args.decayPtAcceptanceMin:g} && sublead_decay_pt > {args.decayPtAcceptanceMin:g}"
    " && abs(lead_decay_eta) < 3 && abs(sublead_decay_eta) < 3"
)
has_forward_pair = "has_forward_genjet_pair == 1"
selected_pair = pair_exists_selection()
has_vbf_tag = tag_selection()
has_truth_lhe_vbf_pair = "has_lhe_vbf_pair == 1"
has_genjet_lhe_comparison = "has_lhe_vbf_pair == 1 && forward_genjet_pair_mass > -900"
has_genjet_lhe_match_pair = "has_genjet_lhe_match_pair == 1"
has_genjet_lhe_mismatch_pair = "has_lhe_vbf_pair == 1 && forward_genjet_pair_mass > -900 && has_genjet_lhe_match_pair == 0"
vbf_tag_and_decayacc = (
    f"{has_vbf_tag}"
    f" && lead_decay_pt > {args.decayPtAcceptanceMin:g} && sublead_decay_pt > {args.decayPtAcceptanceMin:g}"
    " && abs(lead_decay_eta) < 3 && abs(sublead_decay_eta) < 3"
)

# Higgs properties
h_higgs_pt = make_hist(tree, "h_higgs_pt", "higgs_pt", base_higgs, 80, 0, 800, "Generated Higgs p_{T} [GeV]")
h_higgs_eta = make_hist(tree, "h_higgs_eta", "higgs_eta", base_higgs, 60, -6, 6, "Generated Higgs #eta")
h_higgs_phi = make_hist(tree, "h_higgs_phi", "higgs_phi", base_higgs, 64, -3.2, 3.2, "Generated Higgs #phi")
h_higgs_mass = make_hist(tree, "h_higgs_mass", "higgs_mass", base_higgs, 800, 0, 3200, "Generated Higgs mass [GeV]")

for hist, outname in [
    (h_higgs_pt, "higgs_pt"),
    (h_higgs_eta, "higgs_eta"),
    (h_higgs_phi, "higgs_phi"),
    (h_higgs_mass, "higgs_mass"),
]:
    style_single_filled(hist, ROOT.kAzure + 2, ROOT.kAzure - 9, 0.45)
    save_single(hist, common_outdir, outname, args.label)

# Decay comparisons
decay_specs = [
    ("lead_decay_pt", 200, 0, 1000, "Leading decay-product p_{T} [GeV]"),
    ("sublead_decay_pt", 200, 0, 1000, "Subleading decay-product p_{T} [GeV]"),
    ("lead_decay_eta", 60, -6, 6, "Leading decay-product #eta"),
    ("sublead_decay_eta", 60, -6, 6, "Subleading decay-product #eta"),
    ("decay_dR", 60, 0, 6, "Decay-product #DeltaR"),
    ("decay_dEta", 60, 0, 6, "Decay-product |#Delta#eta|"),
    ("decay_dPhi", 64, 0, 3.2, "Decay-product |#Delta#phi|"),
    ("decay_mass", 800, 0, 3200, "Decay-product invariant mass [GeV]"),
]

for idx, (expr, nbins, xmin, xmax, xtitle) in enumerate(decay_specs):
    h_all = make_hist(tree, f"h_decay_all_{idx}", expr, base_decay, nbins, xmin, xmax, xtitle, "Normalized events")
    h_acc = make_hist(tree, f"h_decay_acc_{idx}", expr, accepted_decay, nbins, xmin, xmax, xtitle, "Normalized events")
    style_overlay_line(h_all, ROOT.kGray + 2, line_style=2, line_width=3)
    style_overlay_filled(h_acc, ROOT.kAzure + 2, ROOT.kAzure - 9, 0.35, line_width=3)
    save_overlay(
        [h_all, h_acc],
        ["No acceptance cuts", f"gen p_{{T}} > {args.decayPtAcceptanceMin:g} GeV, |#eta| < 3"],
        common_outdir,
        f"compare_{expr}",
        args.label,
        normalize=True,
    )

# LHE VBF-pair properties
base_lhe_vbf = "has_lhe_vbf_pair == 1"
lhe_count_specs = [
    ("n_lheparts_final_state", 20, -0.5, 19.5, "Final-state LHE multiplicity"),
    ("n_lhe_final_state_partons", 12, -0.5, 11.5, "Final-state LHE parton multiplicity"),
    ("n_lhe_final_state_quarks", 12, -0.5, 11.5, "Final-state LHE quark multiplicity"),
    ("n_lhe_final_state_gluons", 12, -0.5, 11.5, "Final-state LHE gluon multiplicity"),
    ("n_lhe_vbf_quarks", 12, -0.5, 11.5, "Mother-tagged LHE VBF-quark multiplicity"),
]
for expr, nbins, xmin, xmax, xtitle in lhe_count_specs:
    hist = make_hist(tree, f"h_{expr}", expr, "n_lheparts_total > 0", nbins, xmin, xmax, xtitle)
    style_single_filled(hist, ROOT.kBlack, ROOT.kGray + 1, 0.25)
    save_single(hist, common_outdir, expr, args.label)

lhe_vbf_specs = [
    ("lhe_vbf1_pt", 300, 0, 1500, "LHE VBF quark 1 p_{T} [GeV]"),
    ("lhe_vbf2_pt", 300, 0, 1500, "LHE VBF quark 2 p_{T} [GeV]"),
    ("lhe_vbf1_eta", 30, -6, 6, "LHE VBF quark 1 #eta"),
    ("lhe_vbf2_eta", 30, -6, 6, "LHE VBF quark 2 #eta"),
    ("lhe_vbf_pair_dEta", 35, 0, 14, "LHE VBF pair |#Delta#eta|"),
    ("lhe_vbf_pair_dPhi", 32, 0, 3.2, "LHE VBF pair |#Delta#phi|"),
    ("lhe_vbf_pair_dR", 60, 0, 12, "LHE VBF pair #DeltaR"),
    ("lhe_vbf_pair_mass", 800, 0, 8000, "LHE VBF pair mass [GeV]"),
]

for expr, nbins, xmin, xmax, xtitle in lhe_vbf_specs:
    hist = make_hist(tree, f"h_{expr}", expr, base_lhe_vbf, nbins, xmin, xmax, xtitle)
    style_single_filled(hist, ROOT.kOrange + 7, ROOT.kOrange - 2, 0.35)
    save_single(hist, common_outdir, expr, args.label, logy=("mass" in expr or "_pt" in expr))

# LHE VBF vs Higgs-decay comparisons
comparison_specs = [
    ("lhe_vbf_pair_dEta", "decay_dEta", 35, 0, 14, "Pair |#Delta#eta|"),
    ("lhe_vbf_pair_dPhi", "decay_dPhi", 32, 0, 3.2, "Pair |#Delta#phi|"),
    ("lhe_vbf_pair_dR", "decay_dR", 60, 0, 12, "Pair #DeltaR"),
    ("lhe_vbf_pair_mass", "decay_mass", 800, 0, 8000, "Pair invariant mass [GeV]"),
]
comparison_selection = "has_lhe_vbf_pair == 1 && decay_mass > -900"

for idx, (vbf_expr, decay_expr, nbins, xmin, xmax, xtitle) in enumerate(comparison_specs):
    h_vbf = make_hist(
        tree,
        f"h_compare_vbf_{idx}",
        vbf_expr,
        comparison_selection,
        nbins,
        xmin,
        xmax,
        xtitle,
        "Normalized events",
    )
    h_decay = make_hist(
        tree,
        f"h_compare_decay_{idx}",
        decay_expr,
        comparison_selection,
        nbins,
        xmin,
        xmax,
        xtitle,
        "Normalized events",
    )
    style_overlay_filled(h_vbf, ROOT.kOrange + 7, ROOT.kOrange - 2, 0.28, line_width=3)
    style_overlay_filled(h_decay, ROOT.kAzure + 2, ROOT.kAzure - 9, 0.28, line_width=3)
    save_overlay(
        [h_vbf, h_decay],
        ["LHE VBF pair", "Gen Higgs daughters"],
        common_outdir,
        f"compare_vbfpair_vs_higgs_{decay_expr}",
        args.label,
        normalize=True,
    )

# Inclusive LHE VBF-jet vs Higgs-daughter kinematics
inclusive_object_specs = [
    (
        ["lhe_vbf1_pt", "lhe_vbf2_pt"],
        ["lead_decay_pt", "sublead_decay_pt"],
        150,
        0,
        600,
        "Object p_{T} [GeV]",
        "compare_lheVBF_vs_higgsDaughters_ptInclusive",
        True,
    ),
    (
        ["lhe_vbf1_eta", "lhe_vbf2_eta"],
        ["lead_decay_eta", "sublead_decay_eta"],
        60,
        -6,
        6,
        "Object #eta",
        "compare_lheVBF_vs_higgsDaughters_etaInclusive",
        False,
    ),
]
inclusive_selection = "has_lhe_vbf_pair == 1 && decay_mass > -900"

for idx, (lhe_exprs, decay_exprs, nbins, xmin, xmax, xtitle, outname, logy) in enumerate(inclusive_object_specs):
    h_lhe = make_hist_from_expressions(
        tree,
        f"h_compare_lhe_obj_{idx}",
        lhe_exprs,
        inclusive_selection,
        nbins,
        xmin,
        xmax,
        xtitle,
        "Normalized objects",
    )
    h_decay = make_hist_from_expressions(
        tree,
        f"h_compare_decay_obj_{idx}",
        decay_exprs,
        inclusive_selection,
        nbins,
        xmin,
        xmax,
        xtitle,
        "Normalized objects",
    )
    style_overlay_filled(h_lhe, ROOT.kOrange + 7, ROOT.kOrange - 2, 0.28, line_width=3)
    style_overlay_filled(h_decay, ROOT.kAzure + 2, ROOT.kAzure - 9, 0.28, line_width=3)
    save_overlay(
        [h_lhe, h_decay],
        ["LHE VBF jets (inclusive)", "Gen Higgs daughters (inclusive)"],
        common_outdir,
        outname,
        args.label,
        normalize=True,
        logy=logy,
    )

# Gen-jet multiplicities
h_n_genjets_total = make_hist(tree, "h_n_genjets_total", "n_genjets_total", "", 16, -0.5, 15.5, "GenJet multiplicity", "Normalized events")
h_n_genjets_central = make_hist(tree, "h_n_genjets_central", "n_genjets_central", "", 16, -0.5, 15.5, "GenJet multiplicity", "Normalized events")
h_n_genjets_forward = make_hist(tree, "h_n_genjets_forward", "n_genjets_forward", "", 16, -0.5, 15.5, "GenJet multiplicity", "Normalized events")
style_hist(h_n_genjets_total, ROOT.kBlack)
style_hist(h_n_genjets_central, ROOT.kAzure + 2)
style_hist(h_n_genjets_forward, ROOT.kOrange + 7)
save_overlay(
    [h_n_genjets_total, h_n_genjets_central, h_n_genjets_forward],
    ["All GenJets", "|#eta| < 3", "3 #leq |#eta| #leq 5"],
    mode_outdir,
    "genjet_multiplicities",
    args.label,
    normalize=True,
)

# Forward-pair properties
forward_specs = [
    ("forward_genjet1_pt", 200, 0, 1000, "Forward GenJet 1 p_{T} [GeV]"),
    ("forward_genjet2_pt", 200, 0, 1000, "Forward GenJet 2 p_{T} [GeV]"),
    ("forward_genjet1_eta", 30, -6, 6, "Forward GenJet 1 #eta"),
    ("forward_genjet2_eta", 30, -6, 6, "Forward GenJet 2 #eta"),
    ("forward_genjet_pair_dEta", 35, 0, 14, "Forward GenJet pair |#Delta#eta|"),
    ("forward_genjet_pair_dPhi", 15, 0, 6, "Forward GenJet pair |#Delta#phi|"),
    ("forward_genjet_pair_dR", 30, 0, 12, "Forward GenJet pair #DeltaR"),
    ("forward_genjet_pair_mass", 500, 0, 5000, "Forward GenJet pair mass [GeV]"),
]

for expr, nbins, xmin, xmax, xtitle in forward_specs:
    hist = make_hist(tree, f"h_{expr}", expr, selected_pair, nbins, xmin, xmax, xtitle)
    style_single_filled(hist, ROOT.kOrange + 7, ROOT.kOrange - 2, 0.35)
    save_single(hist, mode_outdir, expr, args.label, logy=("mass" in expr or "_pt" in expr))

# GenJet-to-LHE VBF matching diagnostics
matching_specs = [
    ("genjet_lhe_match_n", 4, -0.5, 3.5, "Matched selected GenJets"),
    ("genjet_lhe_match_dr1", 50, 0, 2.5, "Best-assignment #DeltaR_{1}"),
    ("genjet_lhe_match_dr2", 50, 0, 2.5, "Best-assignment #DeltaR_{2}"),
    ("genjet_lhe_match_maxdr", 50, 0, 2.5, "Best-assignment max #DeltaR"),
    ("genjet_lhe_match_sumdr", 60, 0, 4.0, "Best-assignment #Sigma#DeltaR"),
]
for expr, nbins, xmin, xmax, xtitle in matching_specs:
    selection = has_genjet_lhe_comparison if expr != "genjet_lhe_match_n" else has_truth_lhe_vbf_pair
    hist = make_hist(tree, f"h_{expr}", expr, selection, nbins, xmin, xmax, xtitle)
    style_single_filled(hist, ROOT.kOrange + 7, ROOT.kOrange - 2, 0.35)
    save_single(hist, mode_outdir, expr, args.label)

matching_bias_specs = [
    ("forward_genjet_pair_dEta", 35, 0, 14, "Selected GenJet pair |#Delta#eta|"),
    ("forward_genjet_pair_dPhi", 32, 0, 3.2, "Selected GenJet pair |#Delta#phi|"),
    ("forward_genjet_pair_dR", 50, 0, 10, "Selected GenJet pair #DeltaR"),
    ("forward_genjet_pair_mass", 500, 0, 5000, "Selected GenJet pair mass [GeV]"),
]
for idx, (expr, nbins, xmin, xmax, xtitle) in enumerate(matching_bias_specs):
    h_match = make_hist(
        tree,
        f"h_match_bias_pass_{idx}",
        expr,
        has_genjet_lhe_match_pair,
        nbins,
        xmin,
        xmax,
        xtitle,
        "Normalized events",
    )
    h_fail = make_hist(
        tree,
        f"h_match_bias_fail_{idx}",
        expr,
        has_genjet_lhe_mismatch_pair,
        nbins,
        xmin,
        xmax,
        xtitle,
        "Normalized events",
    )
    style_overlay_filled(h_match, ROOT.kAzure + 2, ROOT.kAzure - 9, 0.30, line_width=3)
    style_overlay_line(h_fail, ROOT.kGray + 2, line_style=2, line_width=3)
    save_overlay(
        [h_match, h_fail],
        [f"Matched pair (both #DeltaR < {args.genJetLheMatchDRMax:g})", "Unmatched comparable pair"],
        mode_outdir,
        f"compare_matchBias_{expr}",
        args.label,
        normalize=True,
        logy=("mass" in expr),
    )

# VBF acceptance summary
n_all = tree.GetEntries()
n_tag = tree.GetEntries(has_vbf_tag)
n_truth_lhe_vbf = tree.GetEntries(has_truth_lhe_vbf_pair)
n_genjet_lhe_cmp = tree.GetEntries(has_genjet_lhe_comparison)
n_genjet_lhe_match = tree.GetEntries(has_genjet_lhe_match_pair)
n_tagged_genjet_lhe_cmp = tree.GetEntries(f"{has_genjet_lhe_comparison} && {has_vbf_tag}")
n_tagged_genjet_lhe_match = tree.GetEntries(f"{has_genjet_lhe_match_pair} && {has_vbf_tag}")
n_decayacc_after_tag = tree.GetEntries(vbf_tag_and_decayacc)

cutflow_counts = [
    ("Generated events", n_all),
    (vbf_tag_step_label(), n_tag),
    (f"Central pair acc. (p_{{T}} > {args.decayPtAcceptanceMin:g}, |#eta| < 3)", n_decayacc_after_tag),
]
h_cutflow = build_cutflow_hist(cutflow_counts)
save_cutflow(h_cutflow, cutflow_counts, mode_outdir, "vbf_genjet_cutflow", args.label)

print("[INFO] Summary counts")
for label, value in cutflow_counts:
    frac = 100.0 * value / n_all if n_all > 0 else 0.0
    print(f"  {label:20s}: {value:8d} ({frac:6.2f}%)")
if n_genjet_lhe_cmp > 0:
    eff = 100.0 * n_genjet_lhe_match / n_genjet_lhe_cmp
    print(f"  truth LHE VBF pair     : {n_truth_lhe_vbf:8d}")
    print(f"  GenJet/LHE comparable  : {n_genjet_lhe_cmp:8d}")
    print(f"  pair-match eff.        : {n_genjet_lhe_match:8d} ({eff:6.2f}%) for #DeltaR < {args.genJetLheMatchDRMax:g}")
if n_tagged_genjet_lhe_cmp > 0:
    eff_tag = 100.0 * n_tagged_genjet_lhe_match / n_tagged_genjet_lhe_cmp
    print(f"  tagged pair-match eff. : {n_tagged_genjet_lhe_match:8d} ({eff_tag:6.2f}%) for #DeltaR < {args.genJetLheMatchDRMax:g}")

root_file.Close()
print(f"[INFO] Saved common plots in: {common_outdir}")
print(f"[INFO] Saved VBF-tag-mode plots in: {mode_outdir}")
