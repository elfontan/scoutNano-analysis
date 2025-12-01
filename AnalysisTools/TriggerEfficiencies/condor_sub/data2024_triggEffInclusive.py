# **********************************************
# Run 3 dijet scouting analysis:               #
# trigger efficiency studies                   #
# **********************************************

#!/usr/bin/env python3

import os, sys, math, copy
import ROOT
from ROOT import TLorentzVector, TMath

ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

import array
import numpy as np
import fastjet as fj
#print(fj.__version__) 

jetdef = fj.JetDefinition(fj.antikt_algorithm, 0.4)

# Importing tools from nanoAOD processing set up to store the ratio histograms in a root file
# -------------------------------------------------------------------------------------------
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

def PseudoJ(list_jets):
    pseudojets = []
    for j in list_jets:

        # Convert Collection Objects (pt, eta, phi, mass) to (px, py, pz, E)
        px = j.pt * TMath.Cos(j.phi)
        py = j.pt * TMath.Sin(j.phi)
        pz = j.pt * TMath.SinH(j.eta)
        E = TMath.Sqrt(px**2 + py**2 + pz**2 + j.m**2) #j.mass**2)
        
        # Input jets (px, py, pz, E)
        pj = fj.PseudoJet(px, py, pz, E)
        pseudojets.append(pj)

    return pseudojets

def Lorentz(pseudojets):    
    # Convert PseudoJets into TLorentz vectors (able to calculate Eta)
    return [TLorentzVector(pj.px(), pj.py(), pj.pz(), pj.E()) for pj in pseudojets]

def deltaR(jet, mu):
    deta = jet.eta - mu.eta
    dphi = math.fabs(math.atan2(math.sin(jet.phi - mu.phi), math.cos(jet.phi - mu.phi)))
    return math.sqrt(deta*deta + dphi*dphi)

# ---------------------------------------------------------------------------------------------------
def JetID(jet):
    eta = jet.eta
    aeta = abs(eta)
    
    totalE = (
        jet.neutralHadronEnergy + jet.HFHadronEnergy + 
        jet.photonEnergy + jet.HFEMEnergy +
        jet.muonEnergy + jet.electronEnergy +
        jet.chargedHadronEnergy
    )
    if totalE <= 0:
        return False

    NHF = (jet.neutralHadronEnergy + jet.HFHadronEnergy) / float(totalE)
    NEMF = (jet.photonEnergy + jet.HFEMEnergy) / float(totalE)
    muFrac = jet.muonEnergy / float(totalE)

    chargedMult = jet.chargedHadronMultiplicity + jet.HFHadronMultiplicity
    neutralMult = jet.neutralHadronMultiplicity + jet.HFEMMultiplicity
    nconst = jet.chargedHadronMultiplicity + jet.neutralHadronMultiplicity + jet.muonMultiplicity + jet.electronMultiplicity + jet.photonMultiplicity 

    # -------- |Eta| < 2.6 --------
    if aeta < 2.6:
        if NHF >= 0.99: return False
        if NEMF >= 0.90: return False
        if nconst <= 1: return False
        if chargedMult <= 0: return False
        if muFrac >= 0.80: return False
        return True

    # --- |Eta| = [2.6,2.7] --------
    if (aeta >= 2.6 and aeta < 2.7):
        if NEMF >= 0.99: return False
        if muFrac >= 0.80: return False
        return True
    
    # --- |Eta| = [2.7,3.0] --------
    if (aeta >= 2.7 and aeta < 3.0):

        if NEMF >= 0.99: return False
        if neutralMult <= 1: return False
        return True

    # --- |Eta| = [3.0,5.0] --------
    if (aeta >= 3.0 and aeta < 5.0):

        if NEMF >= 0.10: return False
        return True

    return False
# ---------------------------------------------------------------------------------------------------
def MuonID(mu):
    if mu.pt <= 30: return False
    if abs(mu.eta) >= 0.8: return False

    if abs(mu.trk_dxy) >= 0.2: return False
    if abs(mu.trk_dz) >= 0.5: return False
    if mu.normchi2 >= 3: return False

    if mu.nValidRecoMuonHits <= 0: return False
    if mu.nRecoMuonMatchedStations <= 3: return False
    if mu.nValidPixelHits <= 1: return False
    if mu.nTrackerLayersWithMeasurement <= 7: return False

    return True


class TrigDijetHTAnalysis(Module):
    def __init__(self):
        self.writeHistFile=True
        self.reference_paths=reference_paths
        print(reference_paths)
        self.signal_paths=signal_paths
        
    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)

        # Counters
        # --------
        self.n_totEvents_refTrig = 0 
        self.n_totEvents_refTrigJetId = 0 

        # Histos
        # --------
        self.h_passreftrig = ROOT.TH1F("h_passreftrig" , "; passed ref trigger", 2, 0. , 2.)
        self.h_ht_inclusive_all = ROOT.TH1F("h_ht_inclusive_all", "; H_{T} Inclusive [GeV]; Efficiency;", 150, 0., 1500.)
        self.h_ht_inclusive_passed = ROOT.TH1F("h_ht_inclusive_passed", "; H_{T} Inclusive [GeV]; Efficiency;", 150, 0., 1500.)
        self.h_pt_leading_all = ROOT.TH1F("h_pt_leading_all", "; AK4 Leading jet p_{T} [GeV]; Efficiency;", 50, 0., 500.)
        self.h_pt_leading_passed = ROOT.TH1F("h_pt_leading_passed", "; AK4 Leading jet p_{T} [GeV]; Efficiency;", 50, 0., 500.)

        for h in [
                self.h_passreftrig,
                self.h_ht_inclusive_all, self.h_ht_inclusive_passed,
                self.h_pt_leading_all, self.h_pt_leading_passed,
        ]:
            self.addObject(h)
            
    def analyze(self, event):
        HT = 0.0
        
        dst = Object(event, "DST")

        # -------------------------
        # --- Reference trigger ---
        # -------------------------

        refAccept = any(getattr(dst, path) for path in self.reference_paths)
        self.h_passreftrig.Fill(refAccept)

        if not refAccept:
            return False

        self.n_totEvents_refTrig += 1

        #print("*** ----------- ***")
        #print("*** Event = ", event)
        #print("*** ----------- ***")
        #print("*** n_total_events_refTrigger: ", self.n_totEvents_refTrig)
        #print("*** ----------- ***")

        jets = Collection(event, "ScoutingPFJet")
        njets = getattr(event, "n" + jets._prefix)

        # ---------------------------------
        # --- Preselection requirements ---
        # ---------------------------------

        # --- Jet preselection ---
        njetAcc = [
            j for j in jets
            if j.pt > 30
            and abs(j.eta) < 5
            and JetID(j)          
        ]

        # --- Muon veto (ScoutingMuonNoVtx with dR(mu,jet) < 0.4 ---
        muonsVtx = Collection(event, "ScoutingMuonVtx")
        muons = Collection(event, "ScoutingMuonNoVtx")
        #print("Number of muons = ", len(muons))
        ngood_muonsVtx = [mu for mu in muonsVtx if MuonID(mu)]

        if (len(ngood_muonsVtx) < 1):
            return False
            
        # ----------------------------
        # --- Global HT definition ---
        # ----------------------------
        njetHt = [j for j in njetAcc if j.pt > 30 and abs(j.eta) < 2.5]
        HT = sum(j.pt for j in njetHt)

        # **********************************
        # Check trigger as a function of HT:
        # **********************************
        self.h_ht_inclusive_all.Fill(HT) ## Inclusive

        # --- Unprescaled L1 bits in JetHT scouting trigger
        #unprescaled_l1Triggers = [
        #    "L1_HTT280er",
        #    "L1_SingleJet180",
        #    "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5",
        #    "L1_ETT2000",
        #]
        # --- Prescaled L1 bits in JetHT scouting trigger
        #prescaled_l1Triggers = [
        #    "L1_HTT200er",
        #    "L1_HTT250er",
        #]

        # --- Keep only events in which one of the unprescaled triggers in the JetHT scouting path fired ---          
        if not (                                                                                              
                getattr(event, "L1_HTT280er", False) or                                                        
                getattr(event, "L1_SingleJet180", False) or                                                        
                getattr(event, "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", False) or                                                        
                getattr(event, "L1_ETT2000", False)                                                       
        ):                                                                                    
            return False                                                                                                          
        
        if dst.PFScouting_JetHT == 1:
            self.h_ht_inclusive_passed.Fill(HT) ## Inclusive

        self.n_totEvents_refTrigJetId += 1 

        # **************************************
        # Check trigger as a function of lead pT
        # **************************************
        if len(njetAcc) < 1:
            return False

        #sort_jets = sorted(njetAcc, key=lambda x: x.Pt(), reverse=True)
        #self.h_pt_leading_all.Fill(sort_jets[0].Pt()) 
        #if dst.PFScouting_JetHT == 1:
        #    self.h_pt_leading_passed.Fill(sort_jets[0].Pt()) 
        sort_jets = sorted(njetAcc, key=lambda x: x.pt, reverse=True)
        self.h_pt_leading_all.Fill(sort_jets[0].pt) 
        if dst.PFScouting_JetHT == 1:
            self.h_pt_leading_passed.Fill(sort_jets[0].pt) 
        
        return True
    

    def endJob(self):
        print(f"----------------------------------------------------------------")
        print("Summary after processing:")
        print(f"----------------------------------------------------------------")
        print(f"Total events after reference trigger: {self.n_totEvents_refTrig}")
        print(f"----------------------------------------------------------------")
        Module.endJob(self)

reference_paths = ["PFScouting_SingleMuon"]
signal_paths    = ["PFScouting_JetHT"]

preselection= "" #DST_PFScouting_SingleMuon == 1 && DST_PFScouting_JetHT == 1"

### --------------- ###
### Parse arguments ###
### --------------- ###
args = dict(arg.split('=') for arg in sys.argv[1:] if '=' in arg)
inputFile = args.get('inputFile')
outputFile = args.get('outputFile', 'histos_InclusiveTrigNanoAOD.root')

### ------- ###
### Running ###
### ------- ###
p = PostProcessor(
    ".",
    [inputFile],
    cut=preselection,
    branchsel=None,
    modules=[TrigDijetHTAnalysis()],
    noOut=True,
    histFileName=outputFile,
    histDirName="InclusiveTrigNanoAOD"
)

p.run()
