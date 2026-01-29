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

# -----------------------
# Input / output
# -----------------------
DATA_FILE = "/eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff/histos_Data2024_inclTrigEff_wMuCleaning_GoldenJSON.root"
MC_FILE   = "/eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff/histos_QCD_inclTrigEff_noL1BitRequirement.root"
outdir = "/eos/user/e/elfontan/www/dijetAnaRun3/TRIGGER_EFF/INCLUSIVE_TrgEff/"
os.makedirs(f"{outdir}", exist_ok=True)

# -----------------------
# Helpers
# -----------------------
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


def ratio_and_error(data, data_err, mc, mc_err):
    """Compute Data/MC ratio with symmetric errors (approx)."""
    ratio = data / mc
    err = ratio * np.sqrt(
        (data_err / data) ** 2 + (mc_err / mc) ** 2
    )
    return ratio, err

def match_efficiency_bins(x_d, eff_d, err_d, x_m, eff_m, err_m, tol=1e-6):
    """Match Data and MC efficiency points by bin center."""
    xd, xm = [], []
    ed, em = [], []
    ederr, emerr = [], []

    for i, x in enumerate(x_d):
        diff = np.abs(x_m - x)
        j = np.argmin(diff)
        if diff[j] < tol:
            xd.append(x)
            ed.append(eff_d[i])
            em.append(eff_m[j])
            ederr.append(err_d[i])
            emerr.append(err_m[j])

    return (
        np.array(xd),
        np.array(ed),
        np.array(ederr),
        np.array(em),
        np.array(emerr),
    )

# -----------------------
# Main
# -----------------------
variables = ["ht_inclusive", "pt_leading"]
    
xlabel_map = {
    "ht_inclusive": r"$H_{T}$ [GeV]",
    "pt_leading": r"AK4 leading jet $p_{T}$ [GeV]",
}
name_map = {
    "ht_inclusive": "ht_Inclusive",
    "pt_leading": "pt_leading",
}        

#colors = ["royalblue", "mediumseagreen", "tomato", "orchid", "orange", "deepskyblue", "goldenrod"]

for var in variables:
    
    # --- Figure with ratio pad ---
    fig, (ax, rax) = plt.subplots(
        2, 1,
        figsize=(10, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )
    
    # -----------------------
    # Load ROOT files
    # -----------------------
    f_data = ROOT.TFile.Open(DATA_FILE)
    f_mc   = ROOT.TFile.Open(MC_FILE)
    ddir = f_data.Get("InclusiveTrigNanoAOD")
    mdir = f_mc.Get("InclusiveTrigNanoAOD")
    
    den_d = ddir.Get(f"h_{var}_all")
    num_d = ddir.Get(f"h_{var}_passed")
    
    den_m = mdir.Get(f"h_{var}_all")
    num_m = mdir.Get(f"h_{var}_passed")
    
    # -----------------------
    # Efficiencies
    # -----------------------
    x_d, eff_d, elo_d, ehi_d = get_efficiency(num_d, den_d)
    x_m, eff_m, elo_m, ehi_m = get_efficiency(num_m, den_m)
    
    # -----------------------
    # Upper pad: efficiencies
    # -----------------------
    ax.errorbar(
        x_d, eff_d,
        yerr=[elo_d, ehi_d],
        fmt="o", ms=4, lw=2,
        color="royalblue",
        label="Data",
        capsize=2,
    )
    
    ax.errorbar(
        x_m, eff_m,
        yerr=[elo_m, ehi_m],
        fmt="s", ms=4, lw=2,
        color="orchid",
        label="QCD MC",
        capsize=2,
    )
    
    ax.set_ylabel("Trigger efficiency")
    ax.set_ylim(0, 1.1)
    ax.legend(fontsize=22, loc="lower right")
    
    hep.cms.label(
        "Preliminary",
        data=True,
        year="2024I",
        com="13.6",
        ax=ax,
    )

    if (var == "ht_inclusive"):
        ax.axvline(280, color="gray", linestyle="--", lw=1.4)
    
    # -----------------------
    # Ratio pad: Data / MC
    # -----------------------
    # Symmetric errors
    err_d = 0.5 * (elo_d + ehi_d)
    err_m = 0.5 * (elo_m + ehi_m)
    
    # Match bins
    x_common, eff_d_m, err_d_m, eff_m_m, err_m_m = match_efficiency_bins(
        x_d, eff_d, err_d,
        x_m, eff_m, err_m,
    )
    
    # Compute ratio
    ratio = eff_d_m / eff_m_m
    ratio_err = ratio * np.sqrt(
        (err_d_m / eff_d_m) ** 2 + (err_m_m / eff_m_m) ** 2
    )
    
    rax.errorbar(
        x_common, ratio,
        yerr=ratio_err,
        fmt="o", ms=4,
        color="black",
        capsize=2,
    )
    
    rax.axhline(1.0, color="gray", linestyle="--", lw=1.5)
    rax.set_ylabel("Data / MC")
    rax.set_ylim(0.8, 1.2)
    rax.set_xlabel(xlabel_map[var])
    
    # -----------------------
    # Save
    # -----------------------
    var_tag = name_map.get(var, var)
    for ext in [".png", ".pdf"]:
        outname = f"{outdir}/dataMC_trigEff_{var_tag}_ratio{ext}"
        fig.savefig(outname, dpi=300)
        print(f"Saved: {outname}")
        
    plt.close(fig)
        
