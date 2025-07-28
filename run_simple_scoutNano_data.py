from CMGRDF import *
import ROOT
from CMGRDF.cms.DASSource import DASSource
import os

# ================================
#  CONFIGURATION
# ================================
DEFAULT_OUTPUT_DIR = "/eos/user/e/elfontan/www/dijetAnaRun3/CMGRDF"

YEARS = {
    "2025": { # Dataset /ScoutingPFRun3/Run2025*-ScoutNano-v1/NANOAOD
        "lumi": 1.0,  # placeholder
        "eras": {
            "2025C": [
                "root://xrootd-cms.infn.it//store/data/Run2025C/ScoutingPFRun3/NANOAOD/PromptReco-v1/000/393/087/00000/704e629a-ad74-442e-a108-319da7289aa5.root"
            ],
            # "2025D": [ ... ],
        }
    },
    "2024": { # Dataset /ScoutingPFRun3/Run2024*-ScoutNano-v1/NANOAOD
        "lumi": 1.0,  # placeholder
        "eras": {
            # "2024C": [ ... ],
            # "2024D": [ ... ],
            # "2024E": [ ... ],
            # "2024F": [ ... ],
            # "2024G": [ ... ],
            "2024H": [
                "root://xrootd-cms.infn.it//store/data/Run2024H/ScoutingPFRun3/NANOAOD/ScoutNano-v1/2520000/008642b9-07cb-43ad-8d99-4aad9b24f56d.root",
            ],
            # "2024I": [ ... ],
        }
    },
    "2023": { # Dataset /ScoutingPFRun3/Run2023*-ScoutNano-v1/NANOAOD
        "lumi": 1.0,  # placeholder
        "eras": {
            "2023C": [
                "root://xrootd-cms.infn.it//store/data/Run2023C/ScoutingPFRun3/NANOAOD/ScoutNano-v1/110000/00039154-59ff-4403-aa5d-581d27925cee.root"
            ],
            "2023D": [
                "root://xrootd-cms.infn.it//store/data/Run2023D/ScoutingPFRun3/NANOAOD/ScoutNano-v1/90000/cbe5939d-a464-4789-b260-748c3da1bec4.root"
            ],
        }
    },
    "2022": { # Dataset /ScoutingPFRun3/Run2022*-ScoutNano-v1/NANOAOD
        "lumi": 1.0,  # placeholder
        "eras": {
            "2022C": [
                "root://xrootd-cms.infn.it//store/data/Run2022C/ScoutingPFRun3/NANOAOD/ScoutNano-v1/2520000/00b84aa6-8ed5-420f-9b9e-581a3fceda28.root"
            ],
            # "2022D": [],
            # "2022E": [],
            # "2022F": [],
            # "2022G": [],
        }
    }
}

# ================================
#  STANDARD FUNCTIONS
# ================================

def ensure_dir(directory: str):
    """Ensure output directory exists."""
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

def build_data_samples(year: str, era: str, paths: list):
    """Return list of Data objects for given year/era."""
    samples = [DataSample(f"ScoutingPF_{era}", p) for p in paths]
    print("SAMPLES = ", samples)
    print("🔍 Sample debug:", samples[0])
    #print("🔍 Sample debug:", samples[0], dir(samples[0]))
    return [Data(samples=samples, era=era)]

def get_run_number(file_path: str) -> int:
    """Return the first run number found in the ScoutingNano tree."""
    f = ROOT.TFile.Open(file_path)
    tree = f.Get("Events")  # or "scoutingNanoAOD" depending on your structure
    run_val = None
    if tree:
        tree.GetEntry(0)
        run_val = int(tree.run)
        print("********* run_val = ", run_val)
    f.Close()
    return run_val


# ================================
#  APPLYING TRIGGER SELECTION
# ================================

from CMGRDF import Flow, Cut, Define

def get_trigger_cut(year: int, run: int = None) -> Cut:
    """
    Return the appropriate trigger Cut depending on the year (and run if needed).

    Parameters
    ----------
    year : int
        Data-taking year (2022, 2023, 2024, 2025).
    run : int, optional
        Run number (only used for 2023 split before/after 367622).
    """

    print(f"[DEBUG] get_trigger_cut({year}, run={run})")

    if year in [2024, 2025]:
        return Cut("trigger", "DST_PFScouting_JetHT")

    elif year == 2023: 
        if (run == None):         
            raise ValueError("For 2023, you must specify a run number!")
        elif (run > 367622):                                                                                                             
            return (                                                                                                                                                                         
                Cut("trigger", "DST_Run3_JetHT_PFScoutingPixelTracking")                                                 
                )                                                                                                                                                           
        else:                                                                  
            return(                                                                                                                                                                     
                Cut(                                                                                                                                             
                    "trigger", "DST_Run3_PFScoutingPixelTracking && "                                                                                            
                    "(L1_HTT280er || L1_HTT360er || L1_ETT2000 || L1_SingleJet180 || "    
                    "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5 || "                     
                    "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5)"                                                                                                                               
                )                                                                                                                                                        
            )                                                                                  

    elif year == 2022:
        return Cut(
            "trigger", "DST_Run3_PFScoutingPixelTracking && "
            "(L1_HTT280er || L1_HTT360er || L1_ETT2000 || L1_SingleJet180 || "
            "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5 || "
            "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5)"
        )    

    else:
        raise ValueError(f"Unsupported year: {year}")


# ============================================
#  BUILDING FLOWS: Selection and set of plots
# ============================================

def build_flows(year: int, run: int = None):
    """Define all analysis flows (PFJets_AK4 & FatJets_AK8)."""
    cuts_pfjets = Flow("pfjets",
        get_trigger_cut(year, run), #Cut("trigger", "DST_PFScouting_JetHT"),
        Define("nPFJets", "nScoutingPFJetRecluster"),
        Define("goodJetsIdx", "Nonzero(ScoutingPFJetRecluster_pt > 30 && abs(ScoutingPFJetRecluster_eta) < 5)"),
        Define("sortedJets", "Take(ScoutingPFJetRecluster_pt, goodJetsIdx)"),
        Define("nGoodJets", "Sum(ScoutingPFJetRecluster_pt > 30 && abs(ScoutingPFJetRecluster_eta) < 5)"),
        Cut("nGoodJets", "nGoodJets >= 2"),
        Define("iJet1", "goodJetsIdx[0]"),
        Define("iJet2", "goodJetsIdx[1]"),
        Define("j1pt",  "ScoutingPFJetRecluster_pt[iJet1]"),
        Define("j2pt",  "ScoutingPFJetRecluster_pt[iJet2]"),
        Define("j1eta", "ScoutingPFJetRecluster_eta[iJet1]"),
        Define("j2eta", "ScoutingPFJetRecluster_eta[iJet2]"),
        Define("j1phi", "ScoutingPFJetRecluster_phi[iJet1]"),
        Define("j2phi", "ScoutingPFJetRecluster_phi[iJet2]"),
        Define("mjj", "sqrt(2*j1pt*j2pt*(cosh(j1eta - j2eta) - cos(j1phi - j2phi)))"),
        Define("detajj", "abs(j1eta-j2eta)")
    )

    plots_pfjets = [
        Plot("j1pt", "j1pt", (50, 0, 500), xTitle="Leading jet p_{T} (GeV)", legend="TR"),
        Plot("j2pt", "j2pt", (50, 0, 500), xTitle="Subleading jet p_{T} (GeV)", legend="TR"),
        Plot("mjj", "mjj", (800, 0, 1600), xTitle="m_{jj} (GeV)", legend="TR", logy=True),
        Plot("detajj", "detajj", (20, 0, 10.0), xTitle="DeltaEta_{jj}", legend="TR"),
    ]

    cuts_fatjets = Flow("fatjets",
        get_trigger_cut(year, run), #Cut("trigger", "DST_PFScouting_JetHT"),
        Define("nFatJets",  "nScoutingFatPFJetRecluster"),
        Define("fat_jpt",   "ScoutingFatPFJetRecluster_pt"),
        Define("fat_jeta",  "ScoutingFatPFJetRecluster_eta"),
        Define("fat_jphi",  "ScoutingFatPFJetRecluster_phi"),
        Define("fat_jmass", "ScoutingFatPFJetRecluster_mass")
    )

    plots_fatjets = [
        Plot("fat_jpt", "fat_jpt", (100, 0, 500), xTitle="Fat jet p_{T} (GeV)", legend="TR"),
        Plot("fat_jmass", "fat_jmass", (50, 0, 200), xTitle="Fat jet mass (GeV)", legend="TR"),
    ]

    return [(cuts_pfjets, plots_pfjets), (cuts_fatjets, plots_fatjets)]

def run_analysis(year: str, era: str, data, lumi: float, output_dir: str, run=None, max_events=None):
    """Run the CMGRDF processor."""
    ROOT.EnableImplicitMT(8)
    maker = Processor()

    # If year is 2023, extract the run number from the first file in the sample
    #if year == "2023":
    #    for sample in data.samples:
    #        for file_path in sample.source().files:
    #            print("File path:", file_path)
    #            run = get_run_number(file_path) 
    #             print(f"🔸 Extracted run: {run}")
                    
    #run = None
    #print(f"Booking {era} ({year}) with lumi={lumi}, run={run}")
    #print("---------------------------------------------------")

    for cuts, plots in build_flows(year=year, run=run):
        maker.book(data_objects, lumi, cuts, plots)
            
    # Run plots
    result_plots = maker.runPlots()

    # Output organization
    save_dir = os.path.join(output_dir, f"{year}/{era}")
    ensure_dir(save_dir)
    printer = PlotSetPrinter(
        topRightText=f"Run{era} ({year}, 13.6 TeV)",
        showRatio=False,
        plotFormats="png,pdf"
    )
    printer.printSet(result_plots, save_dir, widePlot=True)
    print(f"✅ Plots saved to: {save_dir}")


# ================================
#  MAIN EXECUTION
# ================================
if __name__ == "__main__":
    for year, yinfo in YEARS.items():
        lumi = yinfo["lumi"]
        for era, paths in yinfo["eras"].items():
            print(f"--------------------------------------------------")
            print(f"🔹 Processing {year} - {era} with {len(paths)} files")

            # Build data samples for this year/era
            data_objects = build_data_samples(year, era, paths)

            for data in data_objects:
                if str(year) == "2023":
                    for sample in data.samples:
                        for file_path in sample.source().files:
                            run = get_run_number(file_path)
                            print(f"🔸 Extracted run: {run} from file {file_path}")
                            
                            # Run analysis for this specific file/run
                            run_analysis(
                                year=int(year),
                                era=era,
                                data=data,
                                lumi=lumi,
                                run=run,
                                output_dir=DEFAULT_OUTPUT_DIR
                            )
                else:
                    # For non-2023, run once per full data object (no run number needed)
                    run_analysis(
                        year=int(year),
                        era=era,
                        data=data,
                        lumi=lumi,
                        run=None,
                        output_dir=DEFAULT_OUTPUT_DIR
                    )

