import ROOT
from CMGRDF import *
import os                                                                                                                                                                                                          
import math                                                                                                             
import glob
import argparse

# ------------------------
# --- Argument parsing ---
# ------------------------
parser = argparse.ArgumentParser(description="MiniAOD dijet analysis")
parser.add_argument("--mHyp", type=int, default=200,
                    help="Mass hypothesis (default: 200 GeV)")
args = parser.parse_args()

# --------------------------
# --- User configuration ---
# --------------------------
mHyp = args.mHyp
print(f"[INFO] Running with mHyp = {mHyp} GeV")

# -----------------------------------------------
# Define the dataset: studying signal MC events
# -----------------------------------------------
P = "/eos/cms/store/cmst3/user/elfontan/scoutAna/samples/mc/2024/Signal/"
signal_name = f"RSGravitonToQuarkQuark_M{mHyp}_wGEN"
#signal_name = f"RSGraviton_M{mHyp}_wGEN"
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
            #label=f"GravToGG ({mHyp} GeV)", 
            label=f"GravToQQ ({mHyp} GeV)", 
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
        Define("goodPFJetsIdx", "Nonzero(ScoutingPFJet_pt > 30 && abs(ScoutingPFJet_eta) < 5)"),
        Define("nGoodPFJets", "Sum(ScoutingPFJet_pt > 30 && abs(ScoutingPFJet_eta) < 5)"),
        Cut("nGoodPFJets", "nGoodPFJets >= 2"),
        Define("iPFJet1", "goodPFJetsIdx[0]"),
        Define("iPFJet2", "goodPFJetsIdx[1]"),
        Define("pfj1pt",  "ScoutingPFJet_pt[iPFJet1]"),
        Define("pfj1eta", "ScoutingPFJet_eta[iPFJet1]"),
        Define("pfj1phi", "ScoutingPFJet_phi[iPFJet1]"),
        Define("pfj2pt",  "ScoutingPFJet_pt[iPFJet2]"),
        Define("pfj2eta", "ScoutingPFJet_eta[iPFJet2]"),
        Define("pfj2phi", "ScoutingPFJet_phi[iPFJet2]"),

        # --------------------------------------------------------
        # Build HT from AK4 jets (sum of AK4 pt with quality cuts)
        # --------------------------------------------------------
        # Use ScoutingPFJet_pt for AK4; require pt > 30 and |eta| < 5 in HT sum
        Define("ht_ak4_all", "Sum(ScoutingPFJet_pt * (ScoutingPFJet_pt > 30 && abs(ScoutingPFJet_eta) < 2.5))"),
        Define("HT", "ht_ak4_all"),
        # Apply HT selection, if needed
        Cut("HT_cut", "HT > 0"),

        # --- VBF-like intermediate category ---
        # forward jets: abs(eta) > 3 and pt > 30
        Define("forwardJetMask", "ScoutingPFJet_pt > 30 && abs(ScoutingPFJet_eta) > 3"),
        Define("forwardIdx", "Nonzero(forwardJetMask)"),
        Define("nForwardJets", "Sum(forwardJetMask)"),
        Cut("nForwardJets", "nForwardJets >= 2"),
        # calculate deltaEta between the two leading forward jets (indices from forwardIdx)
        Define("fjet1_eta", "ScoutingPFJet_eta[forwardIdx[0]]"),
        Define("fjet2_eta", "ScoutingPFJet_eta[forwardIdx[1]]"),
        Define("deltaEta_forward", "abs(fjet1_eta - fjet2_eta)"),
        Cut("deltaEta_forward_cut", "deltaEta_forward > 5"),

        Define("centralJetMask", "ScoutingPFJet_pt > 30 && abs(ScoutingPFJet_eta) < 2.5"),
        Define("centralIdx", "Nonzero(centralJetMask)"),
        Define("nCentralJets", "Sum(centralJetMask)"),
        Cut("nCentralJets", "nCentralJets >= 2"),

        # ------------------------------------------------------------------------
        # Choice of the pairing of the jets:
        # take indices for two highest-pt central jets 
        # ------------------------------------------------------------------------
        Define("iCent1", "centralIdx[0]"),
        Define("iCent2", "centralIdx[1]"),

        # ------------------------------------------------------------------------
        # Choice of the pairing of the jets:
        # select the pair of central jets with the *smallest* deltaPhi separation
        # ------------------------------------------------------------------------
        Define("central_dphi_matrix", R"""
        ROOT::VecOps::RVec<float> dphis;
        for (size_t i = 0; i < ScoutingPFJet_phi.size(); ++i) {
           if (!(ScoutingPFJet_pt[i] > 30 && fabs(ScoutingPFJet_eta[i]) < 2.5)) continue;
           for (size_t j = i + 1; j < ScoutingPFJet_phi.size(); ++j) {
              if (!(ScoutingPFJet_pt[j] > 30 && fabs(ScoutingPFJet_eta[j]) < 2.5)) continue;
              float dphi = fabs(TVector2::Phi_mpi_pi(ScoutingPFJet_phi[i] - ScoutingPFJet_phi[j]));
              dphis.push_back(dphi);
           }
        }
        return dphis;
        """),

        # Get indices of the two central jets with smallest deltaPhi
        Define("closest_central_pair", R"""
        std::array<int,2> pair = {-1, -1};
        float min_dphi = 999;
        for (size_t i = 0; i < ScoutingPFJet_phi.size(); ++i) {
          if (!(ScoutingPFJet_pt[i] > 30 && fabs(ScoutingPFJet_eta[i]) < 2.5)) continue;
           for (size_t j = i + 1; j < ScoutingPFJet_phi.size(); ++j) {
              if (!(ScoutingPFJet_pt[j] > 30 && fabs(ScoutingPFJet_eta[j]) < 2.5)) continue;
              float dphi = fabs(TVector2::Phi_mpi_pi(ScoutingPFJet_phi[i] - ScoutingPFJet_phi[j]));
              if (dphi < min_dphi) {
                 min_dphi = dphi;
                 pair = {static_cast<int>(i), static_cast<int>(j)};
              }
           }
        }
        return pair;
        """),

        # Jet indices for the chosen pair
        # -------------------------------
        #Define("iCent1", "closest_central_pair[0]"),
        #Define("iCent2", "closest_central_pair[1]"),

        # Kinematics of the selected central jets
        # ---------------------------------------
        Define("cent1_pt", "ScoutingPFJet_pt[iCent1]"),
        Define("cent1_eta","ScoutingPFJet_eta[iCent1]"),
        Define("cent1_phi","ScoutingPFJet_phi[iCent1]"),
        Define("cent2_pt", "ScoutingPFJet_pt[iCent2]"),
        Define("cent2_eta","ScoutingPFJet_eta[iCent2]"),
        Define("cent2_phi","ScoutingPFJet_phi[iCent2]"),
        # Dijet mass for central two jets
        # -------------------------------
        Define("mjj_central",
               "sqrt(2*cent1_pt*cent2_pt*(cosh(cent1_eta - cent2_eta) - cos(cent1_phi - cent2_phi)))")
    )
    

    plots_PFJets = [
        Plot("nPFJets", "nPFJets", (70, -0.5, 69.5), xTitle="PF AK4 jet multiplicity", legend="TR", moreY=1.5, integer=True),
        Plot("pfj1pt", "pfj1pt", (100, 0, 1000), xTitle="Leading AK4 jet p_{T} (GeV)", legend="TR", moreY=1.5),
        Plot("pfj1eta", "pfj1eta", (20, -5, 5), xTitle="Leading AK4 jet #eta", legend="TR", moreY=1.5),
        Plot("pfj1phi", "pfj1phi", (14, -3.5, 3.5), xTitle="Leading AK4 jet #phi", legend="TR", moreY=1.5),
        Plot("pfj2pt", "pfj2pt", (100, 0, 1000), xTitle="Subleading AK4 jet p_{T} (GeV)", legend="TR", moreY=1.5),
        Plot("pfj2eta", "pfj2eta", (20, -5, 5), xTitle="Subleading AK4 jet #eta", legend="TR", moreY=1.5),
        Plot("pfj2phi", "pfj2phi", (14, -3.5, 3.5), xTitle="Subleading AK4 jet #phi", legend="TR", moreY=1.5),
        Plot("HT", "HT", (70, 0, 3500), xTitle="HT from AK4 (GeV)", legend="TR", moreY=1.5),

        # ----- inclusive (array-level) distributions for all AK4 jets -----
        Plot("pfj_pt_inclusive",  "ScoutingPFJet_pt",  (100, 0, 300), xTitle="Inclusive AK4 jet p_{T} (GeV)", moreY=1.5),
        Plot("pfj_eta_inclusive", "ScoutingPFJet_eta", (20, -5, 5), xTitle="Inclusive AK4 jet #eta", moreY=1.5),
        Plot("pfj_phi_inclusive", "ScoutingPFJet_phi", (14, -3.5, 3.5), xTitle="Inclusive AK4 jet #phi", moreY=1.5),

        # ----- VBF-like category plots  -----
        Plot("nForwardJets", "nForwardJets", (20, -0.5, 19.5), xTitle="Number of forward jets (|#eta|>3, p_{T}>30)", moreY=1.5, integer=True),
        Plot("deltaEta_forward", "deltaEta_forward", (80, 0, 20), xTitle="|#Delta#eta| (two leading forward jets)", moreY=1.5),
        Plot("fjet1_pt", "ScoutingPFJet_pt[forwardIdx[0]]", (100, 0, 200), xTitle="Forward jet 1 p_{T}", moreY=1.5),
        Plot("fjet1_eta","ScoutingPFJet_eta[forwardIdx[0]]", (50, -5, 5), xTitle="Forward jet 1 #eta", moreY=1.5),
        Plot("fjet2_pt", "ScoutingPFJet_pt[forwardIdx[1]]", (100, 0, 200), xTitle="Forward jet 2 p_{T}", moreY=1.5),
        Plot("fjet2_eta","ScoutingPFJet_eta[forwardIdx[1]]", (50, -5, 5), xTitle="Forward jet 2 #eta", moreY=1.5),

        # Central jets (two highest-pt central jets) and central dijet mass
        Plot("cent1_pt", "cent1_pt", (150, 0, 1500), xTitle="Central jet 1 p_{T}", moreY=1.5),
        Plot("cent1_eta","cent1_eta", (50, -5, 5), xTitle="Central jet 1 #eta", moreY=1.5),
        Plot("cent2_pt", "cent2_pt", (150, 0, 1500), xTitle="Central jet 2 p_{T}", moreY=1.5),
        Plot("cent2_eta","cent2_eta", (50, -5, 5), xTitle="Central jet 2 #eta", moreY=1.5),
        Plot("mjj_central", "mjj_central", (100, 0, 2000), xTitle="Central dijet mass (GeV)", moreY=1.5),
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
        Define("fat_jmassSD_mask", "ScoutingFatPFJetRecluster_msoftdrop > 25"),
        Define("fat_jmassSD_good", "ScoutingFatPFJetRecluster_msoftdrop[fat_jmassSD_mask]"),
        Cut("fat_jmassSD_cut", "Sum(fat_jmassSD_mask) > 25"),
        Define("fat_prob_xqq", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xqq"),
        Define("fat_prob_xbb", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xbb"),
        Define("fat_prob_xcc", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xcc"),
        Define("fat_prob_xgg", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xgg"),
        Define("fat_prob_qcd", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_QCD"),
        Define("gloParT_score", "(fat_prob_xqq + fat_prob_xbb + fat_prob_xcc + fat_prob_xgg) / (fat_prob_xqq + fat_prob_xbb + fat_prob_xcc + fat_prob_xgg + fat_prob_qcd)")
    )

    plots_fatPFJetsRecluster = [
        Plot("nFatJets", "nFatJets", (10, -0.5, 9.5), xTitle="Fat jet multiplicity", legend="TR", moreY=1.5, integer=True),
        Plot("fat_jpt", "fat_jpt", (150, 0, 1500), xTitle="Fat jet p_{T} (GeV)", legend="TR", moreY=1.5),
        Plot("fat_jeta", "fat_jeta", (20, -5, 5), xTitle="Fat jet #eta", legend="TR", moreY=1.5),
        Plot("fat_jphi", "fat_jphi", (14, -3.5, 3.5), xTitle="Fat jet #phi", legend="TR", moreY=1.5),
        Plot("fat_jmass", "fat_jmass", (50, 0, 500), xTitle="Fat jet mass (GeV)", legend="TR", moreY=1.5),
        Plot("fat_jmassSD_good", "fat_jmassSD_good", (50, 0, 500), xTitle="Fat jet SD mass (GeV)", legend="TR", moreY=1.5),
        Plot("fat_prob_xqq", "fat_prob_xqq", (50, 0, 1), xTitle="GloParT P(X_{qq})", legend="TR", moreY=1.5, logy=True),
        Plot("fat_prob_xbb", "fat_prob_xbb", (50, 0, 1), xTitle="GloParT P(X_{bb})", legend="TR", moreY=1.5, logy=True),
        Plot("fat_prob_xcc", "fat_prob_xcc", (50, 0, 1), xTitle="GloParT P(X_{cc})", legend="TR", moreY=1.5, logy=True),
        Plot("fat_prob_xgg", "fat_prob_xgg", (50, 0, 1), xTitle="GloParT P(X_{gg})", legend="TR", moreY=1.5, logy=True),
        Plot("fat_prob_qcd", "fat_prob_qcd", (50, 0, 1), xTitle="GloParT P(QCD)", legend="TR", moreY=1.5, logy=True),
        Plot("gloParT_score", "gloParT_score", (50, 0, 1), xTitle="GloParT score", legend="TR", moreY=1.5, logy=True),
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

        # --- Jet-level masks for pt and SD mass ---
        #Define("goodFatJetsIdx", "Nonzero(ScoutingFatPFJetRecluster_pt > 200 && ScoutingFatPFJetRecluster_msoftdrop > 25)"),
        #Define("nGoodFatJets", "Sum(ScoutingFatPFJetRecluster_pt > 200 && ScoutingFatPFJetRecluster_msoftdrop > 25)"),
        #Cut("nGoodFatJets", "nGoodFatJets >= 1"),
        #Define("bestFatIdx", "goodFatJetsIdx[0]"),
        Define("nGoodFatJets", "Sum(ScoutingFatPFJetRecluster_pt > 200 && ScoutingFatPFJetRecluster_msoftdrop > 25)"),
        Cut("hasGoodFatJet", "nGoodFatJets >= 1"),
        Define("bestFatIdx", "ArgMax( (ScoutingFatPFJetRecluster_pt > 200 && ScoutingFatPFJetRecluster_msoftdrop > 25) * ScoutingFatPFJetRecluster_msoftdrop )"),
                            
        # --- Kinematics of the selected AK8 jet ---
        Define("ak8_pt",     "ScoutingFatPFJetRecluster_pt[bestFatIdx]"),
        Define("ak8_eta",    "ScoutingFatPFJetRecluster_eta[bestFatIdx]"),
        Define("ak8_phi",    "ScoutingFatPFJetRecluster_phi[bestFatIdx]"),
        Define("ak8_mass",   "ScoutingFatPFJetRecluster_mass[bestFatIdx]"),
        Define("ak8_massSD", "ScoutingFatPFJetRecluster_msoftdrop[bestFatIdx]"),

        # --- GloParT ---
        Define("ak8_probxqq", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xqq[bestFatIdx]"),
        Define("ak8_probxbb", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xbb[bestFatIdx]"),
        Define("ak8_probxcc", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xcc[bestFatIdx]"),
        Define("ak8_probxgg", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xgg[bestFatIdx]"),
        Define("ak8_probqcd", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_QCD[bestFatIdx]"),
        Define("ak8_lead_gloParT", "(ak8_probxqq + ak8_probxbb + ak8_probxcc + ak8_probxgg)/ (ak8_probxqq + ak8_probxbb + ak8_probxcc + ak8_probxgg + ak8_probqcd)"),

        # --- Match to AK4 jets for ISR ---
        #Define("dphi_fat_vs_pf", "abs(TVector2::Phi_mpi_pi(ak8_phi - ScoutingPFJet_phi))"),
        #Define("nBackToBackAK4", "Sum(dphi_fat_vs_pf > 2.5)"),
        #Cut("hasBackToBackAK4", "nBackToBackAK4 >= 1"),

        Define("dphi_fat_vs_pf", "abs(acos(cos(ak8_phi - ScoutingPFJet_phi)))"),
        Define("nBackToBackAK4", "Sum(dphi_fat_vs_pf > 2.5)"),
        Cut("hasBackToBackAK4", "nBackToBackAK4 >= 1"),
        Define("backToBackIdx", "Nonzero(dphi_fat_vs_pf > 2.5)"),
        Define("iAK4_ISR", "backToBackIdx[0]"),
        Define("isr_ak4_pt",  "ScoutingPFJet_pt[iAK4_ISR]"),
        Define("isr_ak4_eta", "ScoutingPFJet_eta[iAK4_ISR]"),
        Define("isr_ak4_phi", "ScoutingPFJet_phi[iAK4_ISR]"),    
    )

    plots_ISR = [
        Plot("ak8_pt", "ak8_pt", (150, 0, 1500), xTitle="Leading AK8 p_{T} (GeV) - ISR subset", legend="TR", moreY=1.5),
        Plot("ak8_eta", "ak8_eta", (20, -5, 5), xTitle="Leading AK8 #eta - ISR subset", legend="TR", moreY=1.5), #logy=True
        Plot("ak8_phi", "ak8_phi", (14, -3.5, 3.5), xTitle="Leading AK8 #phi - ISR subset", legend="TR", moreY=1.5),
        Plot("ak8_mass", "ak8_mass", (50, 0, 500), xTitle="Leading AK8 mass (GeV) - ISR subset", legend="TR", moreY=1.5),
        Plot("ak8_massSD", "ak8_massSD", (50, 0, 500), xTitle="Leading AK8 SD mass (GeV) - ISR subset", legend="TR", moreY=1.5),
        Plot("ak8_lead_gloParT", "ak8_lead_gloParT", (50, 0, 1), xTitle="Lead AK8 GloParT - ISR subset", legend="TR", moreY=1.5, logy=True),
        Plot("isr_ak4_pt", "isr_ak4_pt", (80, 0, 1600), xTitle="ISR AK4 p_{T} (GeV)", legend="TR", moreY=1.5),
        Plot("isr_ak4_eta", "isr_ak4_eta", (20, -5, 5), xTitle="ISR AK4 #eta", legend="TR", moreY=1.5),
        Plot("isr_ak4_phi", "isr_ak4_phi", (14, -3.5, 3.5), xTitle="ISR AK4 #phi", legend="TR", moreY=1.5),
        #Plot("debug_goodFatIdx_size", "debug_goodFatIdx_size", (10, 0, 10), xTitle="Size of goodFatIdx", legend="TR")
    ]

    # -----------------------------------------
    # 4) Subset using GloParT Score selection
    # -----------------------------------------
    cuts_gloParT_high = Flow("fatGloParT_high",
        get_trigger_cut(year, run),

        Define("nFatJets",  "nScoutingFatPFJetRecluster"),
        Define("nPFJets",   "nScoutingPFJet"),
        Cut("nFatJets_cut", "nFatJets >= 1"),
        Cut("nPFJets_cut",  "nPFJets >= 1"),                            

        # --- Jet-level masks for pt and SD mass ---
        Define("nGoodFatJets", "Sum(ScoutingFatPFJetRecluster_pt > 200 && ScoutingFatPFJetRecluster_msoftdrop > 25)"),
        Cut("hasGoodFatJet", "nGoodFatJets >= 1"),
        Define("bestFatIdx", "ArgMax( (ScoutingFatPFJetRecluster_pt > 200 && ScoutingFatPFJetRecluster_msoftdrop > 25) * ScoutingFatPFJetRecluster_msoftdrop )"),

        # --- GloParT ---
        Define("ak8_prob_xqq", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xqq[bestFatIdx]"),
        Define("ak8_prob_xbb", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xbb[bestFatIdx]"),
        Define("ak8_prob_xcc", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xcc[bestFatIdx]"),
        Define("ak8_prob_xgg", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_Xgg[bestFatIdx]"),
        Define("ak8_prob_qcd", "ScoutingFatPFJetRecluster_scoutGlobalParT_prob_QCD[bestFatIdx]"),
        Define("ak8_gloParT", "(ak8_prob_xqq + ak8_prob_xbb + ak8_prob_xcc + ak8_prob_xgg) / "
                          "(ak8_prob_xqq + ak8_prob_xbb + ak8_prob_xcc + ak8_prob_xgg + ak8_prob_qcd)"),
        Cut("gloParT_cut", "ak8_gloParT > 0.4"),
        Define("ak8Tagged_hpt",     "ScoutingFatPFJetRecluster_pt[bestFatIdx]"),
        Define("ak8Tagged_heta",    "ScoutingFatPFJetRecluster_eta[bestFatIdx]"),
        Define("ak8Tagged_hphi",    "ScoutingFatPFJetRecluster_phi[bestFatIdx]"),
        Define("ak8Tagged_hmass",   "ScoutingFatPFJetRecluster_mass[bestFatIdx]"),
        Define("ak8Tagged_hmassSD", "ScoutingFatPFJetRecluster_msoftdrop[bestFatIdx]"),
    )

    plots_gloParT_high = [
        Plot("ak8_gloParT", "ak8_gloParT", (50, 0, 1), xTitle="Leading AK8 GloParT score (>0.4)", legend="TR", moreY=1.5, logy=True),
        Plot("ak8Tagged_hpt", "ak8Tagged_hpt", (150, 0, 1500), xTitle="Leading AK8 p_{T} (GeV) - GloParT score>0.4", legend="TR", moreY=1.5),
        Plot("ak8Tagged_heta", "ak8Tagged_heta", (50, -5, 5), xTitle="Leading AK8 #eta - GloParT score>0.4", legend="TR", moreY=1.5),
        Plot("ak8Tagged_hphi", "ak8Tagged_hphi", (14, -3.5, 3.5), xTitle="Leading AK8 #phi - GloParT score>0.4", legend="TR", moreY=1.5),
        Plot("ak8Tagged_hmass", "ak8Tagged_hmass", (50, 0, 500), xTitle="Leading AK8 mass (GeV) - GloParT score>0.4", legend="TR", moreY=1.5),
        Plot("ak8Tagged_hmassSD", "ak8Tagged_hmassSD", (50, 0, 500), xTitle="Leading AK8 SDmass (GeV) - GloParT score>0.4", legend="TR", moreY=1.5),
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

#printer.printSet(result_plots, "/eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/scoutNanoVal/GravToGG_M"+str(mHyp)+"/{flow}/")
printer.printSet(result_plots, "/eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/scoutNanoVal/GravToQQ_M"+str(mHyp)+"/{flow}/")

print("----------CMGRDF signal flow finished successfully. ----------")
