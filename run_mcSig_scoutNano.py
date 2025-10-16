import ROOT
from CMGRDF import *
import os                                                                                                                                                                                                          
import math                                                                                                             
import glob

mHyp = 200

# -----------------------------------------------
# Define the dataset: studying signal MC events
# -----------------------------------------------
P = "/eos/cms/store/cmst3/user/elfontan/scoutAna/samples/mc/2024/Signal/"
signal_name = f"RSGraviton_M{mHyp}"
signal_file = P + signal_name +  "/hadd_nano.root" 

# Collect all ROOT files in the signal folder
# -------------------------------------------
#signal_files = glob.glob(f"{P}/{signal_name}/nano*.root")
#if not signal_files:
#    raise FileNotFoundError(f"No files found for pattern {P}/{signal_name}/nano*.root")

#print(f"Found {len(signal_files)} files for {signal_name}:")
#for f in signal_files:
#    print("  ", f)

signal_process = [
    Process("Signal", 
            [MCSample(signal_name, signal_file, xsec=1.0)],
#            [MCSample(signal_name, f, xsec=1.0) for f in signal_files],
            label=f"GravToGG ({mHyp} GeV)", 
            fillColor=ROOT.kAzure + 7,
            signal=True)
]


# -----------------------------------------------
# Trigger selection
# -----------------------------------------------
from CMGRDF import Flow, Cut, Define, Plot

def get_trigger_cut(year: int, run: int = None) -> Cut:
    return Cut("trigger", "DST_PFScouting_JetHT")


# -----------------------------------------------
# Build flows
# -----------------------------------------------
def build_flows(year: int, run: int = None):

    # -----------------------------------------
    # 1) AK4 PF jets flow
    # -----------------------------------------
    # NOTES: Nonzero(...) is a CMGRDF helper that returns the indices (integers) where the condition is true. Equivalent to NumPy’s np.nonzero().
    #        Sum(...) is conceptually very similar to the Nonzero(...) line:
    #        instead of returning the indices of the good jets, it returns the count of how many pass the selection.

    cuts_PFJets = Flow("PFJets",
        get_trigger_cut(year, run),
        Define("nPFJets", "nScoutingPFJet"),
        Define("goodPFJetsIdx", "Nonzero(ScoutingPFJet_pt > 70 && abs(ScoutingPFJet_eta) < 5)"),
        Define("nGoodPFJets", "Sum(ScoutingPFJet_pt > 70 && abs(ScoutingPFJet_eta) < 5)"),
        Cut("nGoodPFJets", "nGoodPFJets >= 1"),
        Define("iPFJet1", "goodPFJetsIdx[0]"),
        Define("pfj1pt",  "ScoutingPFJet_pt[iPFJet1]"),
        Define("pfj1eta", "ScoutingPFJet_eta[iPFJet1]"),
        Define("pfj1phi", "ScoutingPFJet_phi[iPFJet1]"),
    )

    plots_PFJets = [
        Plot("nPFJets", "nPFJets", (70, -0.5, 69.5), xTitle="PF AK4 jet multiplicity", legend="TR", integer=True),
        Plot("pfj1pt", "pfj1pt", (80, 0, 800), xTitle="Leading AK4 jet p_{T} (GeV)", legend="TR"),
        Plot("pfj1eta", "pfj1eta", (100, -5, 5), xTitle="Leading AK4 jet #eta", legend="TR"),
        Plot("pfj1phi", "pfj1phi", (70, -3.5, 3.5), xTitle="Leading AK4 jet #phi", legend="TR"),
    ]

    # -----------------------------------------
    # 2) Fat jets (AK8 reclustered) + GloParT
    # -----------------------------------------
    cuts_fatPFJetsRecluster = Flow("fatPFJetsRecluster",
        get_trigger_cut(year, run),
        Define("nFatJets",  "nScoutingFatPFJetRecluster"),
        Cut("nFatJets_cut", "nFatJets >= 1"),
        Define("fat_jpt",   "ScoutingFatPFJetRecluster_pt"),
        Define("fat_jeta",  "ScoutingFatPFJetRecluster_eta"),
        Define("fat_jphi",  "ScoutingFatPFJetRecluster_phi"),
        Define("fat_jmass", "ScoutingFatPFJetRecluster_mass"),
        Define("fat_jmassSD_mask", "ScoutingFatPFJetRecluster_msoftdrop > 20"),
        Define("fat_jmassSD_good", "ScoutingFatPFJetRecluster_msoftdrop[fat_jmassSD_mask]"),
        Cut("fat_jmassSD_cut", "Sum(fat_jmassSD_mask) > 0"),
        Define("fat_prob_xqq", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xqq"),
        Define("fat_prob_xgg", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xgg"),
        Define("fat_prob_qcd", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_QCD"),
        Define("gloParT_score", "fat_prob_xgg / (fat_prob_xgg + fat_prob_qcd)")
    )

    plots_fatPFJetsRecluster = [
        Plot("nFatJets", "nFatJets", (10, -0.5, 9.5), xTitle="Fat jet multiplicity", legend="TR", integer=True),
        Plot("fat_jpt", "fat_jpt", (80, 0, 800), xTitle="Fat jet p_{T} (GeV)", legend="TR"),
        Plot("fat_jeta", "fat_jeta", (100, -5, 5), xTitle="Fat jet #eta", legend="TR"),
        Plot("fat_jphi", "fat_jphi", (70, -3.5, 3.5), xTitle="Fat jet #phi", legend="TR"),
        Plot("fat_jmass", "fat_jmass", (50, 0, 500), xTitle="Fat jet mass (GeV)", legend="TR"),
        Plot("fat_jmassSD_good", "fat_jmassSD_good", (50, 0, 500), xTitle="Fat jet SD mass (GeV)", legend="TR"),
        Plot("fat_prob_xqq", "fat_prob_xqq", (50, 0, 1), xTitle="GloParT P(X_{qq})", legend="TR"),
        Plot("fat_prob_xgg", "fat_prob_xgg", (50, 0, 1), xTitle="GloParT P(X_{gg})", legend="TR"),
        Plot("fat_prob_qcd", "fat_prob_qcd", (50, 0, 1), xTitle="GloParT P(QCD)", legend="TR"),
        Plot("gloParT_score", "gloParT_score", (50, 0, 1), xTitle="GloParT score", legend="TR"),
    ]

    # -----------------------------------------
    # 3) ISR-like subset (AK4&AK8 back-to-back) 
    # -----------------------------------------    
    cuts_ISR_matched = Flow("fatAK8_with_AK4_ISR",
        get_trigger_cut(year, run),
        Define("nFatJets",  "nScoutingFatPFJetRecluster"),
        Define("nPFJets",   "nScoutingPFJet"),
        Cut("nFatJets_cut", "nFatJets >= 1"),
        Cut("nPFJets_cut",  "nPFJets >= 1"),
        Define("fat_jphi",  "ScoutingFatPFJetRecluster_phi"),
        Define("fat_jeta",  "ScoutingFatPFJetRecluster_eta"),
        Define("fat_jpt",   "ScoutingFatPFJetRecluster_pt"),
        Define("fat_jmass", "ScoutingFatPFJetRecluster_mass"),
        Define("dphi_fat0_vs_pf", "abs(acos(cos(fat_jphi[0] - ScoutingPFJet_phi)))"),
        Define("nBackToBackAK4", "Sum(dphi_fat0_vs_pf > 2.7)"),
        Cut("hasBackToBackAK4", "nBackToBackAK4 >= 1"),
        Define("backToBackIdx", "Nonzero(dphi_fat0_vs_pf > 2.7)"),
        Define("iAK4_ISR", "backToBackIdx[0]"),
        Define("isr_ak4_pt",  "ScoutingPFJet_pt[iAK4_ISR]"),
        Define("isr_ak4_eta", "ScoutingPFJet_eta[iAK4_ISR]"),
        Define("isr_ak4_phi", "ScoutingPFJet_phi[iAK4_ISR]"),
        Define("ak8_lead_pt",  "ScoutingFatPFJetRecluster_pt[0]"),
        Define("ak8_lead_eta", "ScoutingFatPFJetRecluster_eta[0]"),
        Define("ak8_lead_phi", "ScoutingFatPFJetRecluster_phi[0]"),
        Define("ak8_lead_mass","ScoutingFatPFJetRecluster_mass[0]"),
        Define("ak8_lead_prob_xqq", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xqq[0]"),
        Define("ak8_lead_prob_xgg", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xgg[0]"),
        Define("ak8_lead_prob_qcd", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_QCD[0]"),
        Define("ak8_lead_gloParT", "ak8_lead_prob_xgg / (ak8_lead_prob_xgg + ak8_lead_prob_qcd)")
    )

    plots_ISR = [
        Plot("ak8_lead_pt", "ak8_lead_pt", (80, 0, 800), xTitle="Leading AK8 p_{T} (GeV) - ISR subset", legend="TR"),
        Plot("ak8_lead_eta", "ak8_lead_eta", (100, -5, 5), xTitle="Leading AK8 #eta - ISR subset", legend="TR"), #logy=True
        Plot("ak8_lead_phi", "ak8_lead_phi", (70, -3.5, 3.5), xTitle="Leading AK8 #phi - ISR subset", legend="TR"),
        Plot("ak8_lead_mass", "ak8_lead_mass", (50, 0, 500), xTitle="Leading AK8 mass (GeV) - ISR subset", legend="TR"),
        Plot("isr_ak4_pt", "isr_ak4_pt", (80, 0, 800), xTitle="ISR AK4 p_{T} (GeV)", legend="TR"),
        Plot("isr_ak4_eta", "isr_ak4_eta", (100, -5, 5), xTitle="ISR AK4 #eta", legend="TR"),
        Plot("isr_ak4_phi", "isr_ak4_phi", (70, -3.5, 3.5), xTitle="ISR AK4 #phi", legend="TR"),
        Plot("ak8_lead_gloParT", "ak8_lead_gloParT", (50, 0, 1), xTitle="Leading AK8 GloParT - ISR subset", legend="TR"),
    ]

    # -----------------------------------------
    # 4) Subset using GloParT Score selection
    # -----------------------------------------
    cuts_gloParT_high = Flow("fatGloParT_high",
        get_trigger_cut(year, run),
        Define("nFatJets",  "nScoutingFatPFJetRecluster"),
        Cut("nFatJets_cut", "nFatJets >= 1"),
        Define("lead_prob_xqq", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xqq[0]"),
        Define("lead_prob_xgg", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xgg[0]"),
        Define("lead_prob_qcd", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_QCD[0]"),
        Define("lead_gloParT", "lead_prob_xgg / (lead_prob_xgg + lead_prob_qcd)"),
        Cut("gloParT_cut", "lead_gloParT > 0.4"),
        Define("ak8_hpt", "ScoutingFatPFJetRecluster_pt[0]"),
        Define("ak8_heta","ScoutingFatPFJetRecluster_eta[0]"),
        Define("ak8_hphi","ScoutingFatPFJetRecluster_phi[0]"),
        Define("ak8_hmass","ScoutingFatPFJetRecluster_mass[0]"),
    )

    plots_gloParT_high = [
        Plot("lead_gloParT", "lead_gloParT", (50, 0, 1), xTitle="Leading AK8 GloParT score (>0.4)", legend="TR"),
        Plot("ak8_hpt", "ak8_hpt", (100, 0, 1000), xTitle="Leading AK8 p_{T} (GeV) - GloParT score>0.4", legend="TR"),
        Plot("ak8_heta", "ak8_heta", (100, -5, 5), xTitle="Leading AK8 #eta - GloParT score>0.4", legend="TR"),
        Plot("ak8_hphi", "ak8_hphi", (100, -3.2, 3.2), xTitle="Leading AK8 #phi - GloParT score>0.4", legend="TR"),
        Plot("ak8_hmass", "ak8_hmass", (50, 0, 500), xTitle="Leading AK8 mass (GeV) - GloParT score>0.4", legend="TR"),
    ]

    return [
        (cuts_PFJets, plots_PFJets),
        (cuts_fatPFJetsRecluster, plots_fatPFJetsRecluster),
        (cuts_ISR_matched, plots_ISR),
        (cuts_gloParT_high, plots_gloParT_high)
    ]

# ----------------------------------------------------------------------------
# Run all flows
# ----------------------------------------------------------------------------
lumi = 1.0  # arbitrary for MC studies
ROOT.EnableImplicitMT(8)

flows = build_flows(2024)
maker = Processor()

for cuts, plots in flows:
    maker.book(signal_process, lumi, cuts, plots)

result_plots = maker.runPlots()

printer = PlotSetPrinter(
    topRightText="2024 (13.6 TeV)",
    topLeftText="#bf{CMS} #it{Simulation Preliminary}",
    legendTextSize = 0.04,
    showRatio=False,
    widePlot=True
)

printer.printSet(result_plots, "/eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/scoutNanoVal/M"+str(mHyp)+"/{flow}/")

print("----------CMGRDF signal flow finished successfully. ----------")
