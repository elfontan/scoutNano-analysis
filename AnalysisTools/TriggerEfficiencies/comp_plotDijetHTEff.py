#!/usr/bin/env python3
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import os
from argparse import ArgumentParser

ROOT.gROOT.SetBatch(True)

# CMS style
hep.style.use("CMS")

outdir = "/eos/user/e/elfontan/www/dijetAnaRun3/TRIGGER_EFF/INCLUSIVE_TrgEff/"

def get_efficiency(num, den):
    """Compute efficiency points with asymmetric errors using ROOT.TEfficiency."""
    eff_obj = ROOT.TEfficiency(num, den)
    centers, eff, err_low, err_up = [], [], [], []
    nbins = num.GetNbinsX()

    for b in range(1, nbins + 1):
        total = den.GetBinContent(b)
        if total <= 0:
            continue
        centers.append(num.GetBinCenter(b))
        eff.append(eff_obj.GetEfficiency(b))
        err_low.append(eff_obj.GetEfficiencyErrorLow(b))
        err_up.append(eff_obj.GetEfficiencyErrorUp(b))

    return np.array(centers), np.array(eff), np.array(err_low), np.array(err_up)

    #"""Compute efficiency and bin centers with uncertainties."""
    #values_num = np.array([num.GetBinContent(i+1) for i in range(num.GetNbinsX())])
    #values_den = np.array([den.GetBinContent(i+1) for i in range(den.GetNbinsX())])
    #edges = np.array([den.GetBinLowEdge(i+1) for i in range(den.GetNbinsX())] + [den.GetXaxis().GetXmax()])
    #
    #with np.errstate(divide='ignore', invalid='ignore'):
    #    eff = np.divide(values_num, values_den, out=np.zeros_like(values_num), where=values_den>0)
    #
    #err = np.sqrt(eff * (1 - eff) / np.maximum(values_den, 1))
    #centers = 0.5 * (edges[1:] + edges[:-1])
    #return centers, eff, err, edges

def main(args):
    #variables = ["pt_leading"]
    #variables = ["ht_inclusive", "pt_leading"]
    variables = ["ht_inclusive"]
    #variables = ["minv_1", "minv_2", "minv_3"]
    triggers = ["passed"]

    # Map file name --> label
    file_labels = []
    for fpath in args.rfiles:
        label = os.path.basename(fpath)
        label = label.replace("histos_Data2024_inclTrigEff_", "").replace(".root", "")
        #label = label.replace("histos_DijetHTTrigNanoAOD_", "").replace(".root", "")
        file_labels.append(label)

    colors = ["royalblue", "mediumseagreen", "tomato", "orchid", "orange", "deepskyblue", "goldenrod"]

    os.makedirs(f"{outdir}", exist_ok=True)

    for var in variables:
        fig, ax = plt.subplots(figsize=(10, 9))

        for i, (fpath, label) in enumerate(zip(args.rfiles, file_labels)):
            f = ROOT.TFile.Open(fpath)
            if not f or f.IsZombie():
                print(f"[Warning] Could not open file: {fpath}")
                continue

            fdir = f.Get("InclusiveTrigNanoAOD")
            #fdir = f.Get("DijethtTrigAnalyzerNanoAOD")
            if not fdir:
                print(f"[Warning] Missing directory in {fpath}")
                continue

            den = fdir.Get(f"h_{var}_all")
            num = fdir.Get(f"h_{var}_passed")
            if not num or not den:
                print(f"[Warning] Missing histograms for {var} in {fpath}")
                continue

            centers, eff, err_low, err_up = get_efficiency(num, den)

            ax.errorbar(
                centers, eff,
                yerr=[err_low, err_up],
                fmt='o', ms=4, lw=2.5,
                color=colors[i % len(colors)],
                label=label,
                capsize=2,
            )
            #centers, eff, err, edges = get_efficiency(num, den)
            # 
            #ax.errorbar(
            #    centers, eff, yerr=err,
            #    fmt='o', ms=4, lw=2.0,
            #    color=colors[i % len(colors)],
            #    label=label,
            #    capsize=2
            #)

        # CMS label + lumi text
        hep.cms.label("Preliminary", data=True, year=args.year, com="13.6", ax=ax)

        xlabel_map = {
            #"ht_inclusive": r"$H_{T}$ [GeV]",
            "pt_leading": r"AK4 leading jet $p_{T}$ [GeV]",
        }
        name_map = {
            #"ht_inclusive": "ht_Inclusive",
            "pt_leading": "pt_leading",
        }
        
        #xlabel_map = {
        #    "minv_1": r"$m_{jj}$ VBF [GeV]",
        #    "minv_2": r"$m_{jj}$ ISR [GeV]",
        #    "minv_3": r"$m_{jj}$ pp [GeV]",
        #}
        #name_map = {
        #    "minv_1": "mjjVBF",
        #    "minv_2": "mjjISR",
        #    "minv_3": "mjjPP",
        #}
        
        #ax.set_xlabel(xlabel_map.get(var, r"$m_{jj}$ [GeV]"))
        #ax.set_xlabel(r"$m_{jj}$ [GeV]")
        #ax.set_xlabel(r"AK4 leading jet $p_{T}$ [GeV]")
        ax.set_xlabel(r"$H_{T}$ [GeV]")
        ax.set_ylabel("Trigger efficiency")        
        ax.set_ylim(0, 1.1)
        ax.legend(fontsize=24, loc="lower right", title="Scenarios")
        plt.tight_layout()

        ax.axvline(280, color='black', linestyle='--', linewidth=1.4)

        var_tag = name_map.get(var, var)
        for fmt in args.formats:
            outname = f"{outdir}/TEffTrigEff_{var_tag}_comparison{fmt}"
            fig.savefig(outname, dpi=300)
            print(f"Saved: {outname}")

        plt.close(fig)


if __name__ == "__main__":
    parser = ArgumentParser(description="CMS-style trigger efficiency comparison (mplhep)")
    parser.add_argument(
        "--rfiles", nargs="+", required=True,
        help="List of ROOT files to compare"
    )
    parser.add_argument("--year", default="2024I", help="Year label for CMS text")
    parser.add_argument("--formats", nargs="+", default=[".png", ".pdf"], help="Output formats")
    args = parser.parse_args()
    main(args)
