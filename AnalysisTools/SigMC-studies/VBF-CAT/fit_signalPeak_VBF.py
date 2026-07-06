#!/usr/bin/env python3

import argparse
import os
import re

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fit VBF signal peaks with a double-sided Crystal Ball model."
    )
    parser.add_argument(
        "--indir",
        type=str,
        default=".",
        help="Directory containing ROOT files to fit.",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="plots_VBFSignalFits",
        help="Directory where fit plots and ROOT outputs will be saved.",
    )
    parser.add_argument(
        "--hist",
        type=str,
        default="h_signal_peak_selected_5GeV",
        help="Histogram name to fit when reading prebuilt histograms.",
    )
    parser.add_argument(
        "--algo",
        type=str,
        default=None,
        help="Only process files matching this algorithm/selection token.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Scan subdirectories of --indir as well.",
    )
    parser.add_argument(
        "--tree",
        type=str,
        default="Events",
        help="Tree name to use for current central-widejet study outputs.",
    )
    parser.add_argument(
        "--expr",
        type=str,
        default="widejet_pair_mass",
        help="Branch/expression to histogram when fitting from a tree.",
    )
    parser.add_argument(
        "--selection",
        type=str,
        default="pass_trigger_baseline == 1 && isVBF == 1",
        help="Tree selection used when fitting from a tree.",
    )
    parser.add_argument(
        "--binWidth",
        type=float,
        default=5.0,
        help="Histogram bin width in GeV when building the fit histogram from a tree.",
    )
    parser.add_argument(
        "--xmaxFactor",
        type=float,
        default=2.1,
        help="Maximum x-range relative to the signal mass when building the fit histogram from a tree.",
    )
    return parser.parse_args()


def hex_to_rootcolor(hexcode):
    hexcode = hexcode.lstrip("#")
    r = int(hexcode[0:2], 16) / 255.0
    g = int(hexcode[2:4], 16) / 255.0
    b = int(hexcode[4:6], 16) / 255.0
    return ROOT.TColor.GetColor(r, g, b)


MASS_COLORS = {
    300: hex_to_rootcolor("#3f90da"),
    500: hex_to_rootcolor("#ffa90e"),
    750: hex_to_rootcolor("#bd1f01"),
    1000: hex_to_rootcolor("#2ca02c"),
    2000: hex_to_rootcolor("#9467bd"),
    3000: hex_to_rootcolor("#8c564b"),
}


def set_style():
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetPadLeftMargin(0.14)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleOffset(1.05, "X")
    ROOT.gStyle.SetTitleOffset(1.35, "Y")


def draw_cms_label():
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.050)
    latex.DrawLatex(0.145, 0.935, "CMS")

    latex.SetTextFont(52)
    latex.SetTextSize(0.038)
    latex.DrawLatex(0.235, 0.935, "Simulation Preliminary")

    latex.SetTextFont(42)
    latex.SetTextSize(0.038)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.95, 0.935, "2024 (13.6 TeV)")


def parse_signal_file(path, relpath):
    match = re.search(r"vbf-m(?P<mass>\d+)-(?P<algo>[^/]+)\.root$", relpath)
    if match:
        parent_dir = os.path.dirname(relpath)
        return {
            "path": path,
            "mass": int(match.group("mass")),
            "algo": match.group("algo"),
            "tag": os.path.basename(parent_dir) if parent_dir not in ("", ".") else "base",
            "format": "hist",
        }

    base = os.path.basename(relpath)
    match = re.search(
        r"CentralVBFHTo2B_M(?P<mass>\d+)_centralWideJetVBF_scoutNano_(?P<selection>.+)\.root$",
        base,
    )
    if match:
        selection = match.group("selection")
        compact_match = re.search(r"(fwdPt.*)$", selection)
        algo = compact_match.group(1) if compact_match else selection
        return {
            "path": path,
            "mass": int(match.group("mass")),
            "algo": algo,
            "tag": "centralWideJetVBF",
            "full_selection": selection,
            "format": "tree",
        }

    return None


def discover_files(indir, recursive, algo_filter):
    file_infos = []
    if recursive:
        for root_dir, _, filenames in os.walk(indir):
            for filename in filenames:
                if not filename.endswith(".root"):
                    continue
                full_path = os.path.join(root_dir, filename)
                relpath = os.path.relpath(full_path, indir)
                info = parse_signal_file(full_path, relpath)
                if info is None:
                    continue
                if algo_filter and info["algo"] != algo_filter:
                    continue
                file_infos.append(info)
    else:
        for filename in os.listdir(indir):
            if not filename.endswith(".root"):
                continue
            full_path = os.path.join(indir, filename)
            info = parse_signal_file(full_path, filename)
            if info is None:
                continue
            if algo_filter and info["algo"] != algo_filter:
                continue
            file_infos.append(info)
    return sorted(file_infos, key=lambda item: (item["tag"], item["algo"], item["mass"]))


def sanitize_token(token):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", token)


def color_for_mass(mass):
    if mass in MASS_COLORS:
        return MASS_COLORS[mass]
    palette = [
        "#3f90da",
        "#ffa90e",
        "#bd1f01",
        "#2ca02c",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
    return hex_to_rootcolor(palette[mass % len(palette)])


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


def compute_fwhm_from_roofit_pdf(pdf, xvar, xmin, xmax, nscan=5000):
    xs = []
    ys = []
    step = (xmax - xmin) / float(nscan - 1)
    ymax = -1.0
    xpeak = xmin
    norm_set = ROOT.RooArgSet(xvar)

    for idx in range(nscan):
        xval = xmin + idx * step
        xvar.setVal(xval)
        yval = pdf.getVal(norm_set)
        xs.append(xval)
        ys.append(yval)
        if yval > ymax:
            ymax = yval
            xpeak = xval

    if ymax <= 0:
        return None

    yhalf = 0.5 * ymax

    x_left = None
    for idx in range(1, len(xs)):
        if xs[idx] > xpeak:
            break
        if ys[idx - 1] < yhalf <= ys[idx]:
            x1, x2 = xs[idx - 1], xs[idx]
            y1, y2 = ys[idx - 1], ys[idx]
            x_left = x1 if y2 == y1 else x1 + (yhalf - y1) * (x2 - x1) / (y2 - y1)
            break

    x_right = None
    passed_peak = False
    for idx in range(1, len(xs)):
        if xs[idx] >= xpeak:
            passed_peak = True
        if not passed_peak:
            continue
        if ys[idx - 1] >= yhalf > ys[idx]:
            x1, x2 = xs[idx - 1], xs[idx]
            y1, y2 = ys[idx - 1], ys[idx]
            x_right = x2 if y2 == y1 else x1 + (yhalf - y1) * (x2 - x1) / (y2 - y1)
            break

    if x_left is None or x_right is None or x_right <= x_left:
        return None

    fwhm = x_right - x_left
    return {
        "x_left": x_left,
        "x_right": x_right,
        "fwhm": fwhm,
        "sigma_eq": fwhm / 2.355,
    }


def get_fit_window(mass, hist_xmin, hist_xmax):
    return max(hist_xmin, 0.30 * mass), min(hist_xmax, 2.1 * mass)


def get_draw_range(mass, hist_xmin, hist_xmax):
    return max(hist_xmin, 0.20 * mass), min(hist_xmax, 2.1 * mass)


def initialize_dscb_parameters(hist, info, sigma_eff, mean, sigma_left, sigma_right, alpha_left, n_left, alpha_right, n_right):
    peak_bin = hist.GetMaximumBin()
    peak_x = hist.GetBinCenter(peak_bin)
    hist_rms = hist.GetRMS()
    core_sigma = sigma_eff if sigma_eff > 0 else max(10.0, 0.10 * info["mass"])
    core_sigma = min(core_sigma, max(30.0, 0.35 * info["mass"]))
    peak_guess = peak_x if peak_x > 0 else info["mass"]

    mean.setVal(peak_guess)
    sigma_left.setVal(max(10.0, min(1.10 * core_sigma, max(25.0, 0.30 * info["mass"]))))
    sigma_right.setVal(max(8.0, min(0.70 * core_sigma, max(20.0, 0.20 * info["mass"]))))
    alpha_left.setVal(2.0)
    n_left.setVal(3.0)
    alpha_right.setVal(1.5)
    n_right.setVal(3.0)

    mean.setRange(max(0.70 * info["mass"], peak_guess - max(80.0, 1.5 * hist_rms)),
                  min(1.20 * info["mass"], peak_guess + max(80.0, 1.5 * hist_rms)))


def refit_if_needed(model, datahist, fitres, mean, sigma_left, sigma_right, alpha_left, n_left, alpha_right, n_right, sigma_eff, info):
    if fitres.status() == 0 and fitres.covQual() >= 3:
        return fitres

    core_sigma = sigma_eff if sigma_eff > 0 else max(10.0, 0.10 * info["mass"])
    mean.setVal(min(max(mean.getVal(), 0.75 * info["mass"]), 1.15 * info["mass"]))
    sigma_left.setVal(max(10.0, min(core_sigma, max(25.0, 0.25 * info["mass"]))))
    sigma_right.setVal(max(8.0, min(0.85 * core_sigma, max(18.0, 0.18 * info["mass"]))))
    alpha_left.setVal(1.5)
    n_left.setVal(5.0)
    alpha_right.setVal(2.0)
    n_right.setVal(4.0)

    return model.fitTo(
        datahist,
        ROOT.RooFit.Save(True),
        ROOT.RooFit.PrintLevel(-1),
        ROOT.RooFit.Strategy(2),
    )


def build_histogram_from_tree(root_file, info, args):
    tree = root_file.Get(args.tree)
    if tree is None:
        print(f"[WARNING] Missing tree {args.tree} in {info['path']}")
        return None

    xmax = max(700.0, args.xmaxFactor * float(info["mass"]))
    nbins = int(round(xmax / args.binWidth))
    hist_name = f"h_fitTree_m{info['mass']}_{sanitize_token(info['algo'])}_{sanitize_token(info['tag'])}"
    hist = ROOT.TH1F(hist_name, f";m_{{jj}} [GeV];Events / {args.binWidth:g} GeV", nbins, 0.0, xmax)
    hist.Sumw2()
    draw_expr = f"{args.expr}>>{hist_name}"
    n_drawn = int(tree.Draw(draw_expr, args.selection, "goff"))
    if n_drawn <= 0 or hist.Integral() <= 0:
        print(f"[WARNING] Empty tree selection for {info['path']} with selection: {args.selection}")
        return None
    hist.SetDirectory(0)
    return hist


def load_histogram(info, args):
    root_file = ROOT.TFile.Open(info["path"])
    if not root_file or root_file.IsZombie():
        print(f"[WARNING] Could not open file: {info['path']}")
        return None

    hist = None
    if info.get("format") == "tree":
        try:
            hist = build_histogram_from_tree(root_file, info, args)
        finally:
            root_file.Close()
        return hist

    hist = root_file.Get(args.hist)
    if not hist:
        print(f"[WARNING] Missing histogram {args.hist} in {info['path']}")
        root_file.Close()
        return None
    if not hist.InheritsFrom("TH1"):
        print(f"[WARNING] Object {args.hist} is not a TH1 in {info['path']}")
        root_file.Close()
        return None

    hist = hist.Clone(f"{hist.GetName()}_{info['mass']}_{info['algo']}_{info['tag']}")
    hist.SetDirectory(0)
    root_file.Close()
    return hist


def fit_histogram(info, args, outdir):
    hist = load_histogram(info, args)
    if hist is None:
        return

    if hist.Integral() <= 0:
        print(f"[WARNING] Empty histogram in {info['path']}")
        return

    sigma_eff, x_eff_low, x_eff_high = effective_sigma(hist, frac=0.683)
    hist_xmin = hist.GetXaxis().GetXmin()
    hist_xmax = hist.GetXaxis().GetXmax()
    fit_xmin, fit_xmax = get_fit_window(info["mass"], hist_xmin, hist_xmax)
    draw_xmin, draw_xmax = get_draw_range(info["mass"], hist_xmin, hist_xmax)

    xvar = ROOT.RooRealVar(f"x_{info['mass']}_{info['algo']}", "m_{jj} [GeV]", fit_xmin, fit_xmax)
    datahist = ROOT.RooDataHist(
        f"datahist_{info['mass']}_{info['algo']}",
        f"datahist_{info['mass']}_{info['algo']}",
        ROOT.RooArgList(xvar),
        ROOT.RooFit.Import(hist),
    )

    mean = ROOT.RooRealVar(
        f"mean_{info['mass']}_{info['algo']}",
        "mean",
        info["mass"],
        0.75 * info["mass"],
        1.25 * info["mass"],
    )
    sigma_left = ROOT.RooRealVar(
        f"sigmaL_{info['mass']}_{info['algo']}",
        "sigmaL",
        max(10.0, 0.08 * info["mass"]),
        10.0,
        600.0,
    )
    sigma_right = ROOT.RooRealVar(
        f"sigmaR_{info['mass']}_{info['algo']}",
        "sigmaR",
        max(12.0, 0.12 * info["mass"]),
        1.0,
        500.0,
    )
    #alpha_left = ROOT.RooRealVar(f"alphaL_{info['mass']}_{info['algo']}", "alphaL", 1.0, 0.1, 10.0)
    #n_left = ROOT.RooRealVar(f"nL_{info['mass']}_{info['algo']}", "nL", 2.0, 0.1, 50.0)
    alpha_left = ROOT.RooRealVar(f"alphaL_{info['mass']}_{info['algo']}", "alphaL", 5.0, 0.05, 30.0)
    n_left = ROOT.RooRealVar(f"nL_{info['mass']}_{info['algo']}", "nL", 1.0, 0.5, 50.0)
    alpha_right = ROOT.RooRealVar(f"alphaR_{info['mass']}_{info['algo']}", "alphaR", 1.0, 0.01, 5.0)
    n_right = ROOT.RooRealVar(f"nR_{info['mass']}_{info['algo']}", "nR", 2.0, 0.1, 50.0)

    initialize_dscb_parameters(
        hist,
        info,
        sigma_eff,
        mean,
        sigma_left,
        sigma_right,
        alpha_left,
        n_left,
        alpha_right,
        n_right,
    )

    model = ROOT.RooCrystalBall(
        f"dscb_{info['mass']}_{info['algo']}",
        f"dscb_{info['mass']}_{info['algo']}",
        xvar,
        mean,
        sigma_left,
        sigma_right,
        alpha_left,
        n_left,
        alpha_right,
        n_right,
    )

    fitres = model.fitTo(
        datahist,
        ROOT.RooFit.Save(True),
        ROOT.RooFit.PrintLevel(-1),
        ROOT.RooFit.Strategy(1),
    )
    fitres = refit_if_needed(
        model,
        datahist,
        fitres,
        mean,
        sigma_left,
        sigma_right,
        alpha_left,
        n_left,
        alpha_right,
        n_right,
        sigma_eff,
        info,
    )

    fwhm_info = compute_fwhm_from_roofit_pdf(model, xvar, fit_xmin, fit_xmax)
    rel_width_pct = -1.0
    if fwhm_info is not None and mean.getVal() != 0:
        rel_width_pct = 100.0 * fwhm_info["sigma_eq"] / mean.getVal()

    print(f"\n[INFO] {os.path.basename(info['path'])}")
    print(f"       mean       = {mean.getVal():.3f} +/- {mean.getError():.3f}")
    print(f"       sigmaL     = {sigma_left.getVal():.3f} +/- {sigma_left.getError():.3f}")
    print(f"       sigmaR     = {sigma_right.getVal():.3f} +/- {sigma_right.getError():.3f}")
    print(f"       fit status = {fitres.status()}, covQual = {fitres.covQual()}")
    print(f"       sigma_eff  = {sigma_eff:.3f} GeV")
    if fwhm_info is not None:
        print(f"       FWHM       = {fwhm_info['fwhm']:.3f} GeV")
        print(f"       sigma/m    = {rel_width_pct:.3f} %")

    frame = xvar.frame()
    datahist.plotOn(
        frame,
        ROOT.RooFit.Name("data"),
        ROOT.RooFit.MarkerStyle(20),
        ROOT.RooFit.MarkerSize(0.9),
        ROOT.RooFit.LineColor(ROOT.kBlack),
    )
    model.plotOn(
        frame,
        ROOT.RooFit.Name("model"),
        ROOT.RooFit.LineColor(color_for_mass(info["mass"])),
        ROOT.RooFit.LineWidth(3),
    )

    canvas = ROOT.TCanvas(
        f"c_{sanitize_token(info['tag'])}_{info['mass']}_{sanitize_token(info['algo'])}",
        "",
        1200,
        1000,
    )
    canvas.SetLeftMargin(0.14)
    canvas.SetBottomMargin(0.13)

    frame.SetTitle("")
    frame.GetXaxis().SetTitle("m_{jj} [GeV]")
    frame.GetYaxis().SetTitle("Events / 5 GeV")
    frame.GetXaxis().SetRangeUser(draw_xmin, draw_xmax)
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.GetXaxis().SetTitleOffset(1.05)
    frame.Draw()

    ymax_frame = frame.GetMaximum()
    fwhm_arrow = None
    if fwhm_info is not None:
        y_arrow = 0.55 * ymax_frame
        y_tick_low = 0.51 * ymax_frame
        y_tick_high = 0.59 * ymax_frame

        line_half = ROOT.TLine(fwhm_info["x_left"], y_arrow, fwhm_info["x_right"], y_arrow)
        line_half.SetLineColor(ROOT.kGray + 2)
        line_half.SetLineStyle(2)
        line_half.SetLineWidth(2)
        line_half.Draw()

        fwhm_arrow = ROOT.TArrow(
            fwhm_info["x_left"], y_arrow,
            fwhm_info["x_right"], y_arrow,
            0.015, "<|>",
        )
        fwhm_arrow.SetLineColor(ROOT.kCyan - 6)
        fwhm_arrow.SetFillColor(ROOT.kCyan - 6)
        fwhm_arrow.SetLineWidth(3)
        fwhm_arrow.Draw()

        tick_left = ROOT.TLine(fwhm_info["x_left"], y_tick_low, fwhm_info["x_left"], y_tick_high)
        tick_right = ROOT.TLine(fwhm_info["x_right"], y_tick_low, fwhm_info["x_right"], y_tick_high)
        tick_left.SetLineColor(ROOT.kGray + 2)
        tick_right.SetLineColor(ROOT.kGray + 2)
        tick_left.SetLineWidth(2)
        tick_right.SetLineWidth(2)
        tick_left.Draw()
        tick_right.Draw()

    legend = ROOT.TLegend(0.62, 0.68, 0.90, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.030)
    legend.SetHeader(f"VBFHTo2B M={info['mass']} GeV", "C")
    #legend.SetHeader(f"VBFHTo2B M={info['mass']} GeV ({info['algo']})", "C")
    legend.AddEntry(frame.findObject("data"), "Data", "lep")
    legend.AddEntry(frame.findObject("model"), "Double-CB fit", "l")
    if fwhm_arrow:
        legend.AddEntry(fwhm_arrow, "FWHM interval", "l")
    legend.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.030)
    x0 = 0.6
    y0 = 0.6
    dy = 0.045
    text.DrawLatex(x0, y0, f"#mu = {mean.getVal():.2f} #pm {mean.getError():.2f} GeV")
    text.DrawLatex(x0, y0 - dy, f"#sigma_{{L}} = {sigma_left.getVal():.2f} GeV")
    text.DrawLatex(x0, y0 - 2 * dy, f"#sigma_{{R}} = {sigma_right.getVal():.2f} GeV")
    text.DrawLatex(x0, y0 - 3 * dy, f"#sigma_{{eff}} = {sigma_eff:.1f} GeV")
    if x_eff_low is not None and x_eff_high is not None:
        text.DrawLatex(0.6, 0.40, f"68.3% interval = [{x_eff_low:.0f}, {x_eff_high:.0f}] GeV")
    if fwhm_info is not None:
        text.DrawLatex(0.6, 0.35, f"FWHM = {fwhm_info['fwhm']:.1f} GeV")
        text.DrawLatex(0.6, 0.30, f"#sigma/m = {rel_width_pct:.2f}%")
   # if info.get("format") == "tree":
        #text.DrawLatex(0.6, 0.25, "Selection: pass trigger baseline, isVBF")

    draw_cms_label()

    stem = f"vbf_m{info['mass']}_{sanitize_token(info['algo'])}_{sanitize_token(info['tag'])}"
    outbase = os.path.join(outdir, stem + "_fit")
    canvas.SaveAs(outbase + ".png")
    canvas.SaveAs(outbase + ".pdf")
    #canvas.SaveAs(outbase + "_VBFTag-to-Central.png")
    #canvas.SaveAs(outbase + "_VBFTag-to-Central.pdf")

    output_root = ROOT.TFile(outbase + ".root", "RECREATE")
    hist.Write("hist")
    datahist.Write("datahist")
    fitres.Write("fitResult")
    canvas.Write("canvas")
    output_root.Close()


def main():
    args = parse_args()
    set_style()
    os.makedirs(args.outdir, exist_ok=True)

    file_infos = discover_files(args.indir, args.recursive, args.algo)
    if not file_infos:
        raise RuntimeError(f"No matching ROOT files found in {args.indir}")

    for info in file_infos:
        fit_histogram(info, args, args.outdir)

    print(f"\n[INFO] Done. Outputs saved in: {args.outdir}")


if __name__ == "__main__":
    main()
