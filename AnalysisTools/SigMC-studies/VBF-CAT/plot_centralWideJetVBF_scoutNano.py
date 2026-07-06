#!/usr/bin/env python3

"""
Plotter for the central-widejet VBF ScoutNano study tree.

It produces:
 - central AK4-pair kinematics split by isVBF
 - central widejet-pair kinematics split by isVBF
 - forward/VBF-pair kinematics split by isVBF
 - cutflows and selection summaries
 - AK4 vs widejet signal-peak comparisons for the isVBF-selected sample
"""

import argparse
import os
import re

import ROOT

ROOT.gROOT.SetBatch(True)


_DRAWN_OBJECTS = []


parser = argparse.ArgumentParser(description="Plot the central widejet + VBF ScoutNano study tree")
parser.add_argument(
    "--input",
    type=str,
    default="VBFHTo2B_M300_centralWideJetVBF_scoutNano_centPt30_fwdPt20.root",
    help="Input ROOT file containing the Events tree.",
)
parser.add_argument(
    "--outdir",
    type=str,
    default="plots_centralWideJetVBF_scoutNano",
    help="Directory where plots will be written.",
)
parser.add_argument(
    "--label",
    type=str,
    default="VBFHTo2B M300 (2024)",
    help="Text label drawn on plots.",
)
parser.add_argument(
    "--mass",
    type=float,
    default=None,
    help="Signal mass hypothesis used to place a guide line. If omitted, try to infer it from the input filename.",
)
parser.add_argument(
    "--fitOutput",
    type=str,
    default=None,
    help="Optional fit-ready ROOT output. Defaults to VBF-dijetMass-Histos_ForFIT/vbf-m<MASS>-leading.root.",
)
args = parser.parse_args()


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
    cms_text = latex.DrawLatex(cms_x, cms_y, "CMS")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextFont(52)
    latex2.SetTextSize(0.04)
    sim_text = latex2.DrawLatex(cms_x + 0.10, cms_y, "Simulation Preliminary")

    latex3 = ROOT.TLatex()
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.04)
    latex3.SetTextAlign(31)
    extra_label = latex3.DrawLatex(extra_x, extra_y, extra_text)
    _DRAWN_OBJECTS.extend([latex, latex2, latex3, cms_text, sim_text, extra_label])


def draw_note(lines, x1=0.55, y1=0.68, x2=0.91, y2=0.87, text_size=0.030):
    box = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    box.SetFillColor(ROOT.kGray)
    box.SetFillColorAlpha(ROOT.kGray, 0.15)
    box.SetFillStyle(1001)
    box.SetLineColor(ROOT.kGray + 1)
    box.SetLineWidth(2)
    box.SetShadowColor(0)
    box.SetTextAlign(12)
    box.SetTextFont(42)
    box.SetTextSize(text_size)
    for line in lines:
        box.AddText(line)
    box.Draw()
    _DRAWN_OBJECTS.append(box)
    return box


def draw_metric_box(title, lines, x1, y1, x2, y2, text_size=0.028):
    box = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    box.SetFillColor(ROOT.kGray)
    box.SetFillColorAlpha(ROOT.kGray, 0.10)
    box.SetFillStyle(1001)
    box.SetLineColor(ROOT.kGray + 1)
    box.SetLineWidth(2)
    box.SetShadowColor(0)
    box.SetTextAlign(12)
    box.SetTextFont(42)
    box.SetTextSize(text_size)
    title_text = box.AddText(title)
    title_text.SetTextFont(62)
    for line in lines:
        box.AddText(line)
    box.Draw()
    _DRAWN_OBJECTS.append(box)
    return box


def infer_mass(input_path):
    if args.mass is not None:
        return args.mass
    match = re.search(r"_M(\d+)", os.path.basename(input_path))
    return float(match.group(1)) if match else None


def infer_forward_deta_threshold(input_path):
    base = os.path.basename(input_path)
    match = re.search(r"VBF-?Deta(\d+(?:p\d+)?)", base)
    if not match:
        return None
    return float(match.group(1).replace("p", "."))


def infer_forward_mjj_threshold(input_path):
    base = os.path.basename(input_path)
    match = re.search(r"Mjj(\d+(?:p\d+)?)", base)
    if not match:
        return None
    return float(match.group(1).replace("p", "."))


def format_threshold(value):
    if value is None:
        return None
    if abs(value - round(value)) < 1e-6:
        return f"{int(round(value))}"
    return f"{value:.1f}"


def get_mass_plot_max(mass_guess):
    if mass_guess is None:
        return 1200.0
    return max(700.0, 2.2 * float(mass_guess))


def get_mass_plot_nbins(mass_guess, bin_width=5.0):
    xmax = get_mass_plot_max(mass_guess)
    return int(round(xmax / bin_width)), xmax


def default_fit_output(mass_guess):
    if args.fitOutput:
        return args.fitOutput
    fit_dir = os.path.join(os.getcwd(), "VBF-dijetMass-Histos_ForFIT")
    mass_token = f"m{int(round(mass_guess))}" if mass_guess is not None else "mUnknown"
    return os.path.join(fit_dir, f"vbf-{mass_token}-leading.root")


def make_hist(tree, name, expr, selection, nbins, xmin, xmax, xtitle, ytitle="Events"):
    hist = ROOT.TH1F(name, f";{xtitle};{ytitle}", nbins, xmin, xmax)
    hist.Sumw2()
    tree.Draw(f"{expr}>>{name}", selection, "goff")
    return hist


def style_hist(hist, color, line_style=1, line_width=3):
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetLineStyle(line_style)
    hist.SetLineWidth(line_width)


def set_hist_max(hists, logy=False, max_scale=1.35):
    max_y = max(hist.GetMaximum() for hist in hists) if hists else 0.0
    for hist in hists:
        hist.SetMaximum(max_y * (20.0 if logy else max_scale))
        if logy:
            hist.SetMinimum(1e-4 if hist.Integral() > 0 else 1e-5)


def save_overlay(
    hists,
    labels,
    outdir,
    outname,
    label,
    normalize=False,
    logy=False,
    note_lines=None,
    legend_coords=None,
    note_coords=None,
    note_text_size=0.030,
    vline_x=None,
    vline_label=None,
    legend_title=None,
):
    _DRAWN_OBJECTS.clear()
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1000, 800)
    if logy:
        canvas.SetLogy()

    if legend_coords is None:
        legend_coords = (0.58, 0.70, 0.90, 0.88)
    legend = ROOT.TLegend(*legend_coords)
    legend.SetTextSize(0.035)
    legend.SetFillStyle(0)
    if legend_title:
        legend.SetHeader(legend_title, "C")

    for hist in hists:
        if normalize and hist.Integral() > 0:
            hist.Scale(1.0 / hist.Integral())

    set_hist_max(hists, logy=logy, max_scale=1.40)

    for idx, hist in enumerate(hists):
        hist.Draw("hist" if idx == 0 else "hist same")
        legend.AddEntry(hist, labels[idx], "l")

    line = None
    if vline_x is not None and hists:
        y1 = hists[0].GetMinimum() if logy else 0.0
        y2 = hists[0].GetMaximum()
        line = ROOT.TLine(vline_x, y1, vline_x, y2)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(7)
        line.SetLineWidth(3)
        line.Draw()
        if vline_label:
            legend.AddEntry(line, vline_label, "l")

    legend.Draw()
    note_box = None
    if note_lines:
        if note_coords is None:
            note_box = draw_note(note_lines, text_size=note_text_size)
        else:
            note_box = draw_note(note_lines, *note_coords, text_size=note_text_size)
    draw_label(label)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def build_cutflow_hist(name, counts):
    hist = ROOT.TH1F(name, ";Selection;Events", len(counts), 0.5, len(counts) + 0.5)
    for idx, (text, value) in enumerate(counts, start=1):
        hist.SetBinContent(idx, value)
        hist.GetXaxis().SetBinLabel(idx, text)
    hist.SetFillColor(ROOT.kAzure - 9)
    hist.SetLineColor(ROOT.kAzure + 2)
    hist.SetLineWidth(2)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetXaxis().SetTitleOffset(2.0)
    return hist


def effective_sigma(hist, frac=0.683):
    nbins = hist.GetNbinsX()
    total = hist.Integral(1, nbins)
    if total <= 0:
        return 0.0, None, None

    target = frac * total
    best_width = None
    best_low = None
    best_high = None

    for ibin_low in range(1, nbins + 1):
        integral = 0.0
        for ibin_high in range(ibin_low, nbins + 1):
            integral += hist.GetBinContent(ibin_high)
            if integral >= target:
                x_low = hist.GetBinLowEdge(ibin_low)
                x_high = hist.GetBinLowEdge(ibin_high + 1)
                width = x_high - x_low
                if best_width is None or width < best_width:
                    best_width = width
                    best_low = x_low
                    best_high = x_high
                break

    if best_width is None:
        return 0.0, None, None

    return 0.5 * best_width, best_low, best_high


def fit_peak_gaussian(hist, mass_guess):
    if hist.Integral() <= 0:
        return None

    peak_bin = hist.GetMaximumBin()
    peak_x = hist.GetBinCenter(peak_bin)
    peak_y = hist.GetBinContent(peak_bin)
    if peak_y <= 0:
        return None

    center_guess = peak_x if peak_x > 0 else mass_guess
    sigma_guess = max(15.0, 0.12 * mass_guess)
    fit_min = max(hist.GetXaxis().GetXmin(), center_guess - 3.0 * sigma_guess)
    fit_max = min(hist.GetXaxis().GetXmax(), center_guess + 3.0 * sigma_guess)
    if fit_max <= fit_min:
        return None

    func_name = f"f_gaus_{hist.GetName()}"
    func = ROOT.TF1(func_name, "gaus", fit_min, fit_max)
    func.SetParameters(peak_y, center_guess, sigma_guess)
    func.SetParLimits(1, max(hist.GetXaxis().GetXmin(), 0.6 * mass_guess), min(hist.GetXaxis().GetXmax(), 1.4 * mass_guess))
    func.SetParLimits(2, 5.0, max(30.0, 0.50 * mass_guess))
    fit_result = hist.Fit(func, "RQ0S")
    if int(fit_result.Status()) != 0:
        return None

    sigma = abs(func.GetParameter(2))
    mean = func.GetParameter(1)
    return {
        "func": func,
        "mean": mean,
        "mean_err": func.GetParError(1),
        "sigma": sigma,
        "sigma_err": func.GetParError(2),
        "fwhm": 2.355 * sigma,
        "yield_window": hist.Integral(hist.FindBin(mean - sigma), hist.FindBin(mean + sigma)),
    }


def save_peak_comparison(
    h_ak4,
    h_wide,
    outdir,
    outname,
    label,
    mass_guess,
    note_lines=None,
    fit_selection_label=None,
):
    _DRAWN_OBJECTS.clear()
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1100, 850)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetBottomMargin(0.12)

    h_ak4.SetLineWidth(4)
    h_wide.SetLineWidth(4)
    h_ak4.SetMarkerStyle(20)
    h_ak4.SetMarkerSize(0.0)
    h_wide.SetMarkerStyle(20)
    h_wide.SetMarkerSize(0.0)

    draw_min = 0.0
    draw_max = min(h_ak4.GetXaxis().GetXmax(), get_mass_plot_max(mass_guess))
    h_ak4.GetXaxis().SetRangeUser(draw_min, draw_max)
    h_wide.GetXaxis().SetRangeUser(draw_min, draw_max)

    max_y = max(h_ak4.GetMaximum(), h_wide.GetMaximum())
    h_ak4.SetMaximum(1.35 * max_y)
    h_ak4.SetMinimum(0.0)
    h_ak4.Draw("hist")
    h_wide.Draw("hist same")

    fit_ak4 = fit_peak_gaussian(h_ak4, mass_guess)
    fit_wide = fit_peak_gaussian(h_wide, mass_guess)

    mass_line = ROOT.TLine(mass_guess, 0.0, mass_guess, 1.35 * max_y)
    mass_line.SetLineColor(ROOT.kBlack)
    mass_line.SetLineStyle(7)
    mass_line.SetLineWidth(3)
    mass_line.Draw()

    legend = ROOT.TLegend(0.58, 0.69, 0.90, 0.88)
    legend.SetTextSize(0.033)
    legend.SetFillStyle(0)
    legend.AddEntry(h_ak4, "AK4 central pair", "l")
    legend.AddEntry(h_wide, "Widejet central pair", "l")
    legend.AddEntry(mass_line, "Nominal mass", "l")
    legend.Draw()

    if note_lines:
        draw_note(note_lines, 0.17, 0.76, 0.41, 0.88, 0.028)

    sigma_eff_ak4, _, _ = effective_sigma(h_ak4)
    sigma_eff_wide, _, _ = effective_sigma(h_wide)
    ak4_mass_ref = fit_ak4["mean"] if fit_ak4 is not None and fit_ak4["mean"] != 0 else mass_guess
    wide_mass_ref = fit_wide["mean"] if fit_wide is not None and fit_wide["mean"] != 0 else mass_guess
    ak4_sigma_eff_over_m = 100.0 * sigma_eff_ak4 / ak4_mass_ref if ak4_mass_ref else 0.0
    wide_sigma_eff_over_m = 100.0 * sigma_eff_wide / wide_mass_ref if wide_mass_ref else 0.0
    draw_metric_box(
        fit_selection_label if fit_selection_label else "Peak comparison",
        [],
        0.53,
        0.62,
        0.91,
        0.69,
        text_size=0.028,
    )
    draw_metric_box(
        "AK4 jets",
        [
            f"Selected: {int(h_ak4.Integral())}",
            f"FWHM: {fit_ak4['fwhm']:.1f} GeV" if fit_ak4 is not None else "FWHM: n/a",
            f"sigma_eff/m: {ak4_sigma_eff_over_m:.2f}%",
        ],
        0.53,
        0.43,
        0.71,
        0.61,
        text_size=0.026,
    )
    draw_metric_box(
        "Wide-jets",
        [
            f"Selected: {int(h_wide.Integral())}",
            f"FWHM: {fit_wide['fwhm']:.1f} GeV" if fit_wide is not None else "FWHM: n/a",
            f"sigma_eff/m: {wide_sigma_eff_over_m:.2f}%",
        ],
        0.73,
        0.43,
        0.91,
        0.61,
        text_size=0.026,
    )

    draw_label(label)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def save_cutflow(hist, counts, outdir, outname, label, title):
    _DRAWN_OBJECTS.clear()
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1200, 950)
    canvas.SetBottomMargin(0.24)
    hist.SetBarWidth(0.8)
    hist.SetBarOffset(0.1)
    hist.SetMaximum(hist.GetMaximum() * 1.35)
    hist.Draw("hist text0")
    hist.LabelsOption("u", "X")

    box = ROOT.TPaveText(0.54, 0.50, 0.93, 0.86, "NDC")
    box.SetFillColor(ROOT.kGray)
    box.SetFillColorAlpha(ROOT.kGray, 0.15)
    box.SetLineColor(ROOT.kGray + 1)
    box.SetLineWidth(2)
    box.SetTextAlign(12)
    box.SetTextFont(42)
    box.SetTextSize(0.026)
    head = box.AddText(title)
    head.SetTextFont(62)
    total = float(counts[0][1]) if counts else 0.0
    for step, value in counts:
        frac = 100.0 * value / total if total > 0 else 0.0
        box.AddText(f"{step}: {value} ({frac:.2f}%)")
    box.Draw()

    draw_label(label, cms_x=0.14, cms_y=0.865, extra_x=0.955, extra_y=0.925)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def branch_exists(tree, branch_name):
    return tree.GetBranch(branch_name) is not None


def build_multiplicity_hist(tree, name, branch_name, selection, xtitle):
    values = []
    n_values = int(tree.Draw(branch_name, selection, "goff"))
    vector = tree.GetV1()
    max_value = 0
    for idx in range(n_values):
        value = int(round(vector[idx]))
        values.append(value)
        if value > max_value:
            max_value = value

    xmax = max(12, max_value + 2)
    hist = ROOT.TH1F(name, f";{xtitle};Events", xmax, -0.5, xmax - 0.5)
    hist.Sumw2()
    for value in values:
        hist.Fill(value)
    return hist


def save_multiplicity_comparison(tree, branch_names, branch_labels, branch_colors, selection, outdir, outname, label, xtitle, note_lines):
    hists = []
    labels = []
    for idx, branch_name in enumerate(branch_names):
        if not branch_exists(tree, branch_name):
            continue
        hist = build_multiplicity_hist(
            tree,
            f"h_{outname}_{idx}",
            branch_name,
            selection,
            xtitle,
        )
        style_hist(hist, branch_colors[idx])
        hists.append(hist)
        labels.append(branch_labels[idx])

    if not hists:
        return

    save_overlay(
        hists,
        labels,
        outdir,
        outname,
        label,
        normalize=False,
        logy=False,
        note_lines=note_lines,
        legend_coords=(0.55, 0.72, 0.90, 0.88),
        note_coords=(0.18, 0.70, 0.50, 0.88),
        note_text_size=0.028,
        legend_title="Reco-jet multiplicity",
    )


cms_style()
os.makedirs(args.outdir, exist_ok=True)

root_file = ROOT.TFile.Open(args.input)
if not root_file or root_file.IsZombie():
    raise RuntimeError(f"Could not open input file: {args.input}")

tree = root_file.Get("Events")
if tree is None:
    raise RuntimeError(f"Could not find tree 'Events' in {args.input}")

signal_mass = infer_mass(args.input)
forward_deta_threshold = infer_forward_deta_threshold(args.input)
forward_mjj_threshold = infer_forward_mjj_threshold(args.input)
mass_plot_nbins, mass_plot_max = get_mass_plot_nbins(signal_mass, 5.0)

baseline = "pass_trigger_baseline == 1"
has_central = f"{baseline} && has_central_pair == 1"
step1_pre = f"{baseline} && pass_step1_kinematics == 1"
step1 = f"{baseline} && pass_central_pair_selection == 1"
has_forward = f"{step1} && has_vbf_pair == 1"
isvbf = f"{baseline} && isVBF == 1"
ispp = f"{baseline} && isPP == 1"
step1_matched = f"{step1} && central_pair_genmatched == 1"
step1_unmatched = f"{step1} && central_pair_genmatched == 0"
isvbf_matched = f"{isvbf} && central_pair_genmatched == 1"
isvbf_unmatched = f"{isvbf} && central_pair_genmatched == 0"
ispp_matched = f"{ispp} && central_pair_genmatched == 1"
ispp_unmatched = f"{ispp} && central_pair_genmatched == 0"
vbf_lhe_matched = f"{isvbf} && vbf_pair_lhematched == 1"
vbf_lhe_unmatched = f"{isvbf} && vbf_pair_lhematched == 0"

n_all = tree.GetEntries()
n_baseline = tree.GetEntries(baseline)
n_has_central = tree.GetEntries(has_central)
n_step1_pre = tree.GetEntries(step1_pre)
n_step1 = tree.GetEntries(step1)
n_has_forward = tree.GetEntries(has_forward)
n_isvbf = tree.GetEntries(isvbf)
n_ispp = tree.GetEntries(ispp)
n_step1_matched = tree.GetEntries(step1_matched)
n_step1_unmatched = tree.GetEntries(step1_unmatched)
n_isvbf_matched = tree.GetEntries(isvbf_matched)
n_isvbf_unmatched = tree.GetEntries(isvbf_unmatched)
n_ispp_matched = tree.GetEntries(ispp_matched)
n_ispp_unmatched = tree.GetEntries(ispp_unmatched)
n_vbf_truth = tree.GetEntries(f"{baseline} && has_lhe_vbf_pair == 1")
n_vbf_match = tree.GetEntries(vbf_lhe_matched)

step1_match_purity = 100.0 * n_step1_matched / n_step1 if n_step1 else 0.0
isvbf_match_purity = 100.0 * n_isvbf_matched / n_isvbf if n_isvbf else 0.0
ispp_match_purity = 100.0 * n_ispp_matched / n_ispp if n_ispp else 0.0
vbf_match_eff = 100.0 * n_vbf_match / n_vbf_truth if n_vbf_truth else 0.0
vbf_match_purity = 100.0 * n_vbf_match / n_isvbf if n_isvbf else 0.0

summary_note = [
    f"Events: {n_all}",
    f"Trigger baseline: {n_baseline}",
    f"CentralSel: {n_step1}",
    f"isVBF / isPP: {n_isvbf} / {n_ispp}",
]


def build_mass_hist(name, expr, selection):
    return make_hist(
        tree,
        name,
        expr,
        selection,
        mass_plot_nbins,
        0,
        mass_plot_max,
        "Dijet mass [GeV]",
        "Events / 5 GeV",
    )


central_specs = [
    ("central_pair_dEta", 60, 0, 6, "Central AK4 seed pair |#Delta#eta|"),
    ("central_pair_dPhi", 64, 0, 3.2, "Central AK4 seed pair |#Delta#phi|"),
    ("central_pair_mass", mass_plot_nbins, 0, mass_plot_max, "Central AK4 seed pair mass [GeV]"),
    ("widejet_pair_dEta", 60, 0, 6, "Central wide-jet pair |#Delta#eta|"),
    ("widejet_pair_dPhi", 64, 0, 3.2, "Central wide-jet pair |#Delta#phi|"),
    ("widejet_pair_mass", mass_plot_nbins, 0, mass_plot_max, "Central wide-jet pair mass [GeV]"),
]

for expr, nbins, xmin, xmax, xtitle in central_specs:
    h_step1 = make_hist(tree, f"h_step1_{expr}", expr, step1, nbins, xmin, xmax, xtitle, "Events")
    h_isvbf = make_hist(tree, f"h_isvbf_{expr}", expr, isvbf, nbins, xmin, xmax, xtitle, "Events")
    h_ispp = make_hist(tree, f"h_ispp_{expr}", expr, ispp, nbins, xmin, xmax, xtitle, "Events")
    style_hist(h_step1, ROOT.kBlack)
    style_hist(h_isvbf, ROOT.kRed + 1)
    style_hist(h_ispp, ROOT.kBlue + 1)
    save_overlay(
        [h_step1, h_isvbf, h_ispp],
        ["CentralSel", "isVBF", "isPP"],
        args.outdir,
        f"compare_step1_isVBF_isPP_{expr}",
        args.label,
        normalize=False,
        logy=False if "mass" in expr else False,
        note_lines=summary_note,
        legend_coords=(0.62, 0.72, 0.90, 0.88),
        note_coords=(0.18, 0.72, 0.48, 0.88),
        note_text_size=0.028,
    )

forward_specs = [
    ("forward_pair_dEta", "forward_pair_dEta", 80, 0, 10, "Forward leftover pair |#Delta#eta|"),
    ("forward_pair_dPhi", "forward_pair_dPhi", 64, 0, 3.2, "Forward leftover pair |#Delta#phi|"),
    ("forward_pair_dR", "forward_pair_dR", 80, 0, 10, "Forward leftover pair #DeltaR"),
    ("forward_pair_mass", "forward_pair_mass", 120, 0, 4000, "Forward leftover pair mass [GeV]"),
    ("forward_jet1_pt", "forward_jet1_pt", 100, 0, 800, "Leading forward leftover jet p_{T} [GeV]"),
    ("forward_jet2_pt", "forward_jet2_pt", 100, 0, 600, "Subleading forward leftover jet p_{T} [GeV]"),
    ("abs(forward_jet1_eta)", "forward_jet1_abseta", 60, 0, 6, "Leading forward leftover jet |#eta|"),
    ("abs(forward_jet2_eta)", "forward_jet2_abseta", 60, 0, 6, "Subleading forward leftover jet |#eta|"),
]

forward_deta_label = (
    f"VBF |#Delta#eta| thr. = {format_threshold(forward_deta_threshold)}"
    if forward_deta_threshold is not None
    else None
)
forward_mjj_label = (
    f"VBF m_{{jj}} thr. = {format_threshold(forward_mjj_threshold)} GeV"
    if forward_mjj_threshold is not None
    else None
)

for expr, suffix, nbins, xmin, xmax, xtitle in forward_specs:
    h_isvbf = make_hist(tree, f"h_forward_isvbf_{suffix}", expr, isvbf, nbins, xmin, xmax, xtitle, "Events")
    h_ispp = make_hist(tree, f"h_forward_ispp_{suffix}", expr, ispp, nbins, xmin, xmax, xtitle, "Events")
    style_hist(h_isvbf, ROOT.kOrange + 7)
    style_hist(h_ispp, ROOT.kAzure + 2)
    save_overlay(
        [h_isvbf, h_ispp],
        ["isVBF", "isPP"],
        args.outdir,
        f"compare_forward_isVBF_isPP_{suffix}",
        args.label,
        normalize=False,
        logy=False if "mass" in expr else False,
        note_lines=[
            f"VBF categorization: {n_has_forward}",
            f"isVBF: {n_isvbf}",
            f"isPP: {n_ispp}",
        ],
        legend_coords=(0.67, 0.76, 0.90, 0.88),
        note_coords=(0.18, 0.72, 0.48, 0.88),
        note_text_size=0.028,
        vline_x=(
            forward_deta_threshold if expr == "forward_pair_dEta"
            else forward_mjj_threshold if expr == "forward_pair_mass"
            else None
        ),
        vline_label=(
            forward_deta_label if expr == "forward_pair_dEta"
            else forward_mjj_label if expr == "forward_pair_mass"
            else None
        ),
    )

multiplicity_specs = [
    ("total", "Inclusive reco-jet multiplicity", "nRecoJets_inclusive"),
    ("central_eta2p5", "Central reco-jet multiplicity, |#eta| < 2.5", "nRecoJets_central_eta2p5"),
    ("forward_eta2p5", "Forward reco-jet multiplicity, 2.5 #leq |#eta| < 5.0", "nRecoJets_forward_eta2p5"),
]

for region, xtitle, suffix in multiplicity_specs:
    for pt_threshold in (20, 30):
        branch_names = [
            f"n_recojets_raw_{region}_pt{pt_threshold}",
            f"n_recojets_jec_{region}_pt{pt_threshold}",
            f"n_recojets_postwide_{region}_pt{pt_threshold}",
        ]
        branch_labels = [
            f"Raw jets, p_{{T}} > {pt_threshold} GeV",
            f"JEC jets, p_{{T}} > {pt_threshold} GeV",
            f"Post wide-jet, p_{{T}} > {pt_threshold} GeV",
        ]
        branch_colors = [ROOT.kBlack, ROOT.kAzure + 2, ROOT.kRed + 1]
        save_multiplicity_comparison(
            tree,
            branch_names,
            branch_labels,
            branch_colors,
            baseline,
            args.outdir,
            f"compare_multiplicity_{suffix}_pt{pt_threshold}",
            args.label,
            xtitle,
            None,
        )

matched_compare_specs = [
    ("central_pair_mass", mass_plot_nbins, 0, mass_plot_max, "Central AK4 seed pair mass [GeV]"),
    ("widejet_pair_mass", mass_plot_nbins, 0, mass_plot_max, "Central wide-jet pair mass [GeV]"),
    ("widejet_pair_dEta", 60, 0, 6, "Central wide-jet pair |#Delta#eta|"),
]

matched_compare_cases = [
    ("step1", step1_matched, step1_unmatched, n_step1_matched, n_step1_unmatched, "CentralSel"),
    ("isvbf", isvbf_matched, isvbf_unmatched, n_isvbf_matched, n_isvbf_unmatched, "isVBF"),
    ("ispp", ispp_matched, ispp_unmatched, n_ispp_matched, n_ispp_unmatched, "isPP"),
]

for suffix, matched_sel, unmatched_sel, n_match, n_unmatch, label_text in matched_compare_cases:
    for expr, nbins, xmin, xmax, xtitle in matched_compare_specs:
        h_match = make_hist(tree, f"h_match_{suffix}_{expr}", expr, matched_sel, nbins, xmin, xmax, xtitle, "Events")
        h_unmatch = make_hist(tree, f"h_unmatch_{suffix}_{expr}", expr, unmatched_sel, nbins, xmin, xmax, xtitle, "Events")
        style_hist(h_match, ROOT.kRed + 1)
        style_hist(h_unmatch, ROOT.kBlack)
        save_overlay(
            [h_match, h_unmatch],
            ["Matched", "Not matched"],
            args.outdir,
            f"compare_matched_unmatched_{suffix}_{expr}",
            args.label,
            normalize=False,
            logy=False if "mass" in expr else False,
            note_lines=[
                f"Selection: {label_text}",
                f"Matched: {n_match}",
                f"Not matched: {n_unmatch}",
            ],
            legend_coords=(0.67, 0.76, 0.90, 0.88),
            note_coords=(0.18, 0.72, 0.48, 0.88),
            note_text_size=0.030,
        )

h_step1_wide = build_mass_hist("h_step1_wide", "widejet_pair_mass", step1)
h_step1_match = build_mass_hist("h_step1_match", "widejet_pair_mass", step1_matched)
h_step1_unmatch = build_mass_hist("h_step1_unmatch", "widejet_pair_mass", step1_unmatched)
style_hist(h_step1_wide, ROOT.kBlue + 1)
style_hist(h_step1_match, ROOT.kRed + 1)
style_hist(h_step1_unmatch, ROOT.kBlack)
save_overlay(
    [h_step1_wide, h_step1_match, h_step1_unmatch],
    ["CentralSel all", "Gen-matched", "Not gen-matched"],
    args.outdir,
    "compare_mass_step1_matchedUnmatched",
    args.label,
    normalize=False,
    logy=False,
    note_lines=[
        f"CentralSel: {n_step1}",
        f"Matched purity: {step1_match_purity:.1f}%",
    ],
    legend_coords=(0.56, 0.72, 0.90, 0.88),
    note_coords=(0.18, 0.74, 0.46, 0.88),
    vline_x=signal_mass,
    vline_label="Nominal signal mass",
)

h_isvbf_wide = build_mass_hist("h_isvbf_wide", "widejet_pair_mass", isvbf)
h_ispp_wide = build_mass_hist("h_ispp_wide", "widejet_pair_mass", ispp)
style_hist(h_isvbf_wide, ROOT.kRed + 1)
style_hist(h_ispp_wide, ROOT.kBlue + 1)
save_overlay(
    [h_isvbf_wide, h_ispp_wide],
    ["isVBF", "isPP"],
    args.outdir,
    "compare_mass_isVBF_isPP",
    args.label,
    normalize=False,
    logy=False,
    note_lines=[
        f"isVBF: {n_isvbf}",
        f"isPP: {n_ispp}",
    ],
    legend_coords=(0.68, 0.76, 0.90, 0.88),
    note_coords=(0.18, 0.70, 0.48, 0.88),
    note_text_size=0.028,
    vline_x=signal_mass,
    vline_label="Nominal signal mass",
)

h_ak4_isvbf = build_mass_hist("h_ak4_isvbf", "central_pair_mass", isvbf)
h_wide_isvbf = build_mass_hist("h_wide_isvbf", "widejet_pair_mass", isvbf)
style_hist(h_ak4_isvbf, ROOT.kAzure - 3)
style_hist(h_wide_isvbf, ROOT.kMagenta - 7)
save_peak_comparison(
    h_ak4_isvbf,
    h_wide_isvbf,
    args.outdir,
    "compare_signal_peak_isVBF",
    args.label,
    signal_mass,
    note_lines=[
        f"isVBF events: {n_isvbf}",
        f"Gen-matched: {n_isvbf_matched}",
        f"VBF LHE matched: {n_vbf_match}",
    ],
    fit_selection_label="Selection: isVBF",
)

h_ak4_step1_match = build_mass_hist("h_ak4_step1_match", "central_pair_mass", step1_matched)
h_wide_step1_match = build_mass_hist("h_wide_step1_match", "widejet_pair_mass", step1_matched)
style_hist(h_ak4_step1_match, ROOT.kAzure - 3)
style_hist(h_wide_step1_match, ROOT.kMagenta - 7)
save_peak_comparison(
    h_ak4_step1_match,
    h_wide_step1_match,
    args.outdir,
    "compare_signal_peak_step1_higgsMatched",
    args.label,
    signal_mass,
    note_lines=[
        f"Gen-matched: {n_step1_matched}",
        f"CentralSel purity: {step1_match_purity:.1f}%",
    ],
    fit_selection_label="Selection: CentralSel gen-matched",
)

h_mass_shift_step1 = make_hist(
    tree,
    "h_mass_shift_step1",
    "widejet_pair_mass-central_pair_mass",
    step1,
    120,
    -300,
    300,
    "#it{m}_{jj}^{wide} - #it{m}_{jj}^{AK4} [GeV]",
    "Events",
)
h_mass_shift_isvbf = make_hist(
    tree,
    "h_mass_shift_isvbf",
    "widejet_pair_mass-central_pair_mass",
    isvbf,
    120,
    -300,
    300,
    "#it{m}_{jj}^{wide} - #it{m}_{jj}^{AK4} [GeV]",
    "Events",
)
style_hist(h_mass_shift_step1, ROOT.kBlack)
style_hist(h_mass_shift_isvbf, ROOT.kRed + 1)
save_overlay(
    [h_mass_shift_step1, h_mass_shift_isvbf],
    ["CentralSel", "isVBF"],
    args.outdir,
    "compare_signal_peak_massShift_step1_vs_isVBF",
    args.label,
    normalize=False,
    logy=False,
    note_lines=[
        "Positive values mean the wide-jet mass is larger",
        f"CentralSel: {n_step1}",
        f"isVBF: {n_isvbf}",
    ],
    legend_coords=(0.69, 0.78, 0.90, 0.88),
    note_coords=(0.18, 0.72, 0.50, 0.88),
    note_text_size=0.028,
)

h_vbf_match = build_mass_hist("h_vbf_match", "widejet_pair_mass", vbf_lhe_matched)
h_vbf_unmatch = build_mass_hist("h_vbf_unmatch", "widejet_pair_mass", vbf_lhe_unmatched)
style_hist(h_vbf_match, ROOT.kRed + 1)
style_hist(h_vbf_unmatch, ROOT.kBlack)
save_overlay(
    [h_vbf_match, h_vbf_unmatch],
    ["VBF LHE matched", "VBF not LHE matched"],
    args.outdir,
    "compare_mass_isVBF_vbfMatch",
    args.label,
    normalize=False,
    logy=False,
    note_lines=[
        f"isVBF: {n_isvbf}",
        f"VBF-tag match purity: {vbf_match_purity:.1f}%",
    ],
    legend_coords=(0.56, 0.76, 0.90, 0.88),
    note_coords=(0.18, 0.72, 0.50, 0.88),
    note_text_size=0.028,
)

cutflow = [
    ("All events", n_all),
    ("Trigger", n_baseline),
    ("Central pair candidate", n_step1_pre),
    ("CentralSel (#Delta#eta < 1.3)", n_step1),
    ("VBF cand", n_has_forward),
    ("isVBF", n_isvbf),
    ("Gen-matched", n_step1_matched),
]
h_cutflow = build_cutflow_hist("h_cutflow_centralWideJetVBF", cutflow)
save_cutflow(h_cutflow, cutflow, args.outdir, "cutflow_centralWideJetVBF", args.label, "Central wide-jet + VBF study")

fit_root_path = os.path.join(args.outdir, "signal_peak_histograms.root")
fit_file = ROOT.TFile(fit_root_path, "RECREATE")
fit_hists = [
    build_mass_hist("h_widejet_peak_step1_5GeV", "widejet_pair_mass", step1),
    build_mass_hist("h_widejet_peak_step1_matched_5GeV", "widejet_pair_mass", step1_matched),
    build_mass_hist("h_widejet_peak_step1_unmatched_5GeV", "widejet_pair_mass", step1_unmatched),
    build_mass_hist("h_widejet_peak_isVBF_5GeV", "widejet_pair_mass", isvbf),
    build_mass_hist("h_widejet_peak_isPP_5GeV", "widejet_pair_mass", ispp),
    build_mass_hist("h_widejet_peak_isVBF_matched_5GeV", "widejet_pair_mass", isvbf_matched),
    build_mass_hist("h_widejet_peak_isVBF_unmatched_5GeV", "widejet_pair_mass", isvbf_unmatched),
    build_mass_hist("h_ak4_peak_isVBF_5GeV", "central_pair_mass", isvbf),
]
for hist in fit_hists:
    hist.Write()
fit_file.Close()
root_file.Close()

print("[INFO] Summary")
print(f"  Total events                 : {n_all}")
print(f"  Trigger baseline             : {n_baseline}")
print(f"  Central pair candidate       : {n_step1_pre}")
print(f"  CentralSel (dEta < 1.3)      : {n_step1}")
print(f"  VBF candidate found          : {n_has_forward}")
print(f"  isVBF                        : {n_isvbf}")
print(f"  isPP                         : {n_ispp}")
print(f"  Gen-matched                  : {n_step1_matched}")
print(f"  Not gen-matched              : {n_step1_unmatched}")
print(f"  isVBF matched                : {n_isvbf_matched}")
print(f"  isVBF unmatched              : {n_isvbf_unmatched}")
print(f"  VBF LHE matched              : {n_vbf_match}")
print(f"  CentralSel match purity [%]  : {step1_match_purity:.2f}")
print(f"  isVBF match purity [%]       : {isvbf_match_purity:.2f}")
print(f"  isPP match purity [%]        : {ispp_match_purity:.2f}")
print(f"  VBF-tag match efficiency [%] : {vbf_match_eff:.2f}")
print(f"  VBF-tag match purity [%]     : {vbf_match_purity:.2f}")
print(f"[INFO] Saved plots in: {args.outdir}")
print(f"[INFO] Saved signal-peak histos in: {fit_root_path}")
