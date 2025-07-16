from CMGRDF import *
import ROOT
import math
from CMGRDF.cms.DASSource import DASSource

# --------------------------------------------------- #
# Reading a dataset from EOS or a single file
# --------------------------------------------------- #
#P = localOrEOS("2018", "/data/gpetrucc/NanoTrees_TTH_v6", "/eos/cms/store/cmst3/group/tthlep/peruzzi/NanoTrees_TTH_091019_v6pre")
P = "root://xrootd-cms.infn.it//store/data/Run2025C/ScoutingPFRun3/NANOAOD/PromptReco-v1/000/393/087/00000/704e629a-ad74-442e-a108-319da7289aa5.root"

data_dijet = [
    Data([DataSample("ScoutingPF_Run2025C", P)]),
]
print(f"Reading from {P}")

# --------------------------------------------------- #
# Reading a dataset from DAS through DASSource
# --------------------------------------------------- #
# DataSample implements .source(era), needed when the book() method is called.
# Internally, DataSample delegates to DASSource to resolve files.

#data_dijet = [
#    Data(
#        samples=[
#            DataSample(
#                "ScoutingPFRun3_Run2025C", # Sample name 
#                DASSource(
#                    "ScoutingPFRun3_Run2025C", # Sample name
#                    "/ScoutingPFRun3/Run2025C-PromptReco-v1/NANOAOD" # Sample dataset
#                )
#            )
#        ],
#        era="2025C"
#    )
#]

#max_events = 100

cuts_pfjets = Flow("pfjets",
    #Cut("event_limit", "Entry$ < {}".format(max_events)),
    Cut("trigger", "DST_PFScouting_JetHT"),
    Define("nPFJets", "nScoutingPFJetRecluster"),
    Define("jpt",  "ScoutingPFJetRecluster_pt"),
)

plots_pfjets = [
    Plot("jpt", "jpt", (50, 0, 300), xTitle="PF jet p_{T} (GeV)", legend="pfjets"),
    Plot("nPFJets", "nPFJets", (10, -0.5, 9.5), xTitle="PF jet multiplicity", legend="pfjets", integer=True),
]

cuts_fatjets = Flow("fatjets",
    #Cut("event_limit", "Entry$ < {}".format(max_events)),
    #Cut("trigger", "DST_PFScouting_JetHT"),
    Cut("nFatJets", "nScoutingFatPFJetRecluster >= 1"),
    #Cut("nFatJets", "Sum(ScoutingFatPFJetRecluster_pt > 30 && abs(ScoutingFatPFJetRecluster_eta) < 5) >= 2"),
    Define("iJet1",  "ScoutingFatPFJetRecluster_pt.at(0)"),
    #Define("iJet1",  "(ScoutingFatPFJetRecluster_pt > 30 && abs(ScoutingFatPFJetRecluster_eta) < 5).at(0)"),
    #Define("iJet2",  "(ScoutingFatPFJetRecluster_pt > 30 && abs(ScoutingFatPFJetRecluster_eta) < 5).at(1)"),
    #Define("fj1pt",  "ScoutingFatPFJetRecluster_pt[iJet1]"),
    #Define("fj2pt",  "ScoutingFatPFJetRecluster_pt[iJet2]"),
    #Define("fj1eta", "ScoutingFatPFJetRecluster_eta[iJet1]"),
    #Define("fj2eta", "ScoutingFatPFJetRecluster_eta[iJet2]"),
    #Define("fj1phi", "ScoutingFatPFJetRecluster_phi[iJet1]"),
    #Define("fj2phi", "ScoutingFatPFJetRecluster_phi[iJet2]"),
    #Define("fj1m",   "ScoutingFatPFJetRecluster_mass[iJet1]"),
    #Define("fj2m",   "ScoutingFatPFJetRecluster_mass[iJet2]"),
    #Define("mjj", "sqrt(2*fj1pt*fj2pt*(cosh(fj1eta - fj2eta) - cos(fj1phi - fj2phi)))"),
    #Define("detajj", "abs(fj1eta - fj2eta)"),
)

plots_fatjet = [
    #Plot("fj1pt", "fj1pt", (50, 0, 300), xTitle="Leading fat jet p_{T} (GeV)", legend="Dijet"),
    #Plot("fj2pt", "fj2pt", (50, 0, 300), xTitle="Subleading fat jet p_{T} (GeV)", legend="Dijet"),
    #Plot("fj1eta", "fj1eta", (50, -3, 3), xTitle="Leading fat jet eta", legend="Dijet"),
    #Plot("fj2eta", "fj2eta", (50, -3, 3), xTitle="Subleading fat jet eta", legend="Dijet"),
    #Plot("fj1m", "fj1m", (50, 0, 200), xTitle="Leading fat jet mass (GeV)", legend="Dijet"),
    #Plot("fj2m", "fj2m", (50, 0, 200), xTitle="Subleading fat jet mass (GeV)", legend="Dijet"),
    #Plot("mjj", "mjj", (80, 0, 2000), xTitle="m_{jj} (GeV)", legend="Dijet", logy=True, moreY=10),
    #Plot("detajj", "detajj", (50, 0, 6), xTitle="DeltaEta_{jj}", legend="Dijet"),
    Plot("nJets", "Sum(ScoutingFatPFJetRecluster_pt > 0)", (10, -0.5, 9.5),
         xTitle="Fat jet multiplicity (p_{T} > 30, |eta| < 3)", legend="Dijet", integer=True),
]

lumi = 1.0  # placeholder

ROOT.EnableImplicitMT(8)
maker = Processor()
maker.book(data_dijet, lumi, cuts_pfjets, plots_pfjets)
#maker.book(data_dijet, lumi, cuts_fatjets, plots_fatjets)

result_plots = maker.runPlots()
printer = PlotSetPrinter(topRightText="Run2025C (13.6 TeV)", showRatio=False)
printer.printSet(result_plots, "/eos/user/e/elfontan/www/dijetAnaRun3/CMGRDF/")
