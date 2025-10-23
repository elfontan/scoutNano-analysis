#----------------------------------------------------------------------------
# Original code from Marina Kolosova: "TriggerEfficiency" for Collection Objs
# from helpers.utils import deltaPhi, deltaR, deltaEta
#----------------------------------------------------------------------------

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


class TrigDijetHTAnalysis(Module):
    def __init__(self):
        self.writeHistFile=True
        self.reference_paths=reference_paths
        print(reference_paths)
        self.signal_paths=signal_paths
        
    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)

        # Counters
        self.n_total_events = 0        # after refAccept
        self.n_VBF = 0
        self.n_Boosted = 0
        self.n_Resolved = 0
        self.n_Rest = 0

        self.h_passreftrig = ROOT.TH1F("h_passreftrig" , "; passed ref trigger", 2, 0. , 2.)
        self.h_minv_1_all = ROOT.TH1F("h_minv_1_all", "; M_{jj} VBF [GeV]; Efficency;", 100, 100., 1100.)
        self.h_minv_1_passed = ROOT.TH1F("h_minv_1_passed", "; M_{jj} VBF [GeV]; Efficency;", 100, 100., 1100.)
        self.h_minv_2_all = ROOT.TH1F("h_minv_2_all", "; M_{jj} Boosted ISR [GeV]; Efficiency;", 100, 100., 1100.)
        self.h_minv_2_passed = ROOT.TH1F("h_minv_2_passed", "; M_{jj} Boosted ISR [GeV]; Efficiency;", 100, 100., 1100.)
        self.h_minv_3_all = ROOT.TH1F("h_minv_3_all", "; M_{jj} Resolved ISR [GeV]; Efficency;", 100, 100., 1100.)
        self.h_minv_3_passed = ROOT.TH1F("h_minv_3_passed", "; M_{jj} Resolved ISR [GeV]; Efficency;", 100, 100., 1100.)
        self.h_minv_4_all = ROOT.TH1F("h_minv_4_all", "; M_{jj} pp [GeV]; Efficency;", 100, 100., 1100.)
        self.h_minv_4_passed = ROOT.TH1F("h_minv_4_passed", "; M_{jj} pp [GeV]; Efficency;", 100, 100., 1100.)
        
        self.addObject(self.h_passreftrig)
        self.addObject(self.h_minv_1_all)
        self.addObject(self.h_minv_1_passed)        
        self.addObject(self.h_minv_2_all)
        self.addObject(self.h_minv_2_passed)
        self.addObject(self.h_minv_3_all)
        self.addObject(self.h_minv_3_passed)
        self.addObject(self.h_minv_4_all)
        self.addObject(self.h_minv_4_passed)

        # for h in self.hList:
        #     self.addObject(h)

    def analyze(self, event):
        dst = Object(event, "DST")
        
        # --- Reference trigger ---
        refAccept = any(getattr(dst, path) for path in self.reference_paths)
        self.h_passreftrig.Fill(refAccept)
        if not refAccept:
            return False

        self.n_total_events += 1
                
        jets = Collection(event, "ScoutingPFJet")
        njets = getattr(event, "n" + jets._prefix)

        # --- Jet preselection ---
        njetAN = [j for j in jets if j.pt > 30 and abs(j.eta) < 5]

        if len(njetAN) < 2:
            return False

        sort_jets = sorted(njetAN, key=lambda x: x.pt, reverse=True)
        jetdef_ak4 = fj.JetDefinition(fj.antikt_algorithm, 0.4)
        jetdef_ak8 = fj.JetDefinition(fj.antikt_algorithm, 0.8)
        
        Mjj_1 = Mjj_2 = Mjj_3 = Mjj_4 = 0.0
        
        # -------------------------------
        # [1] VBF CATEGORY
        # -------------------------------
        if len(sort_jets) >= 2:
            dEta = abs(sort_jets[0].eta - sort_jets[1].eta)
        else:
            return False  # no valid dijet system

        if dEta > 4:
            Dijet_1 = sort_jets.copy()
            del Dijet_1[0:2]

            cluster_1 = fj.ClusterSequence(PseudoJ(Dijet_1), jetdef_ak4)
            if len(cluster_1.inclusive_jets()) > 1:
                wj_1 = Lorentz(fj.sorted_by_pt(cluster_1.inclusive_jets()))
                if abs(wj_1[0].Eta() - wj_1[1].Eta()) < 1.3:
                    Mjj_1 = (wj_1[0] + wj_1[1]).M()

            # Remove VBF events from further categories
            self.h_minv_1_all.Fill(Mjj_1)
            if dst.PFScouting_JetHT == 1:
                self.h_minv_1_passed.Fill(Mjj_1)

            self.n_VBF += 1
            return True
                
        # -------------------------------
        # [2] BOOSTED w/ ISR CATEGORY (AK8)
        # -------------------------------
        cluster_boosted = fj.ClusterSequence(PseudoJ(sort_jets), jetdef_ak8)
        ak8jets = fj.sorted_by_pt(cluster_boosted.inclusive_jets())
        
        if len(ak8jets) >= 1:
            boosted_lv = Lorentz(ak8jets)[0]
            if boosted_lv.Pt() > 170:
                for j in sort_jets:
                    j_lv = TLorentzVector()                
                    j_lv.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.m)
                    
                    dphi = abs(boosted_lv.DeltaPhi(j_lv))
                    dR = boosted_lv.DeltaR(j_lv)
                    pt_asym = abs(boosted_lv.Pt() - j_lv.Pt()) / (boosted_lv.Pt() + j_lv.Pt())
                    
                    if abs(dphi - math.pi) < 0.4 and dR > 1.0 and pt_asym < 0.3:
                        Mjj_2 = (boosted_lv + j_lv).M()
                        self.h_minv_2_all.Fill(Mjj_2)
                        if dst.PFScouting_JetHT == 1:
                            self.h_minv_2_passed.Fill(Mjj_2)

                        self.n_Boosted += 1
                        return True  # remove boosted events

        # -----------------------------------
        # [3] RESOLVED w/ ISR CATEGORY (AK4)
        # -----------------------------------
        cluster_resolved = fj.ClusterSequence(PseudoJ(sort_jets), jetdef_ak4)
        wj_ak4 = Lorentz(fj.sorted_by_pt(cluster_resolved.inclusive_jets()))

        if len(wj_ak4) > 2:
            if (abs(wj_ak4[0].Eta() - wj_ak4[1].Eta()) < 1.3) and wj_ak4[2].Pt()>70:
            #if (abs(wj_ak4[0].Eta() - wj_ak4[1].Eta()) < 1.3) and wj_ak4[0].Pt()>70 and wj_ak4[1].Pt()>70 and wj_ak4[2].Pt()>70:
                Mjj_3 = (wj_ak4[0] + wj_ak4[1]).M()
                self.h_minv_3_all.Fill(Mjj_3)
                if dst.PFScouting_JetHT == 1:
                    self.h_minv_3_passed.Fill(Mjj_3)

                self.n_Resolved += 1
                return True  # remove ISR resolved

        #if len(wj_ak4) < 3:
        #    return False

        # --- Identify resonance pair (closest DR)
        #best_pair = None
        #min_dR = 999
        #for i in range(len(wj_ak4)):
        #    for j in range(i+1, len(wj_ak4)):
        #        #dEta = abs(wj_ak4[i].Eta() - wj_ak4[j].Eta())
        #        #if dEta < min_dEta:
        #        #    min_dEta = dEta
        #        #    best_pair = (i, j)
        #        dR = wj_ak4[i].DeltaR(wj_ak4[j])
        #        if dR < min_dR:
        #            min_dR = dR
        #            best_pair = (i, j)

        #if not best_pair:
        #    return False

        #if min_dEta > 1.3:
        #    return False
        
        #if wj_ak4[best_pair[0]].Pt() < 30 or wj_ak4[best_pair[1]].Pt() < 30:
        #    return False

        #res_pair = wj_ak4[best_pair[0]] + wj_ak4[best_pair[1]]

        # --- Find ISR jet (back-to-back with resonance system) ---
        #ISR_jet = None
        #max_dphi = 0
        #for k in range(len(wj_ak4)):
        #    if k in best_pair:
        #        continue
        #    dphi = abs(res_pair.DeltaPhi(wj_ak4[k]))
        #    if dphi > max_dphi:
        #        max_dphi = dphi
        #        ISR_jet = wj_ak4[k]

        # --- Require back-to-back and ISR pT cut ---
        #if ISR_jet and ISR_jet.Pt() > 70:
        ##if ISR_jet and abs(max_dphi - math.pi) < 0.4 and ISR_jet.Pt() > 70:
        #    Mjj_3 = (res_pair).M()
        #    self.h_minv_3_all.Fill(Mjj_3)
        #    if dst.PFScouting_JetHT == 1:
        #        self.h_minv_3_passed.Fill(Mjj_3)
        #        return True # remove ISR resolved
        #Mjj_3 = (res_pair).M()
        #self.h_minv_3_all.Fill(Mjj_3)
        #if dst.PFScouting_JetHT == 1:
        #    self.h_minv_3_passed.Fill(Mjj_3)
        #    return True # remove ISR resolved

        # -------------------------------
        # [4] REST CATEGORY
        # -------------------------------
        if len(wj_ak4) >= 2:
            if (abs(wj_ak4[0].Eta() - wj_ak4[1].Eta()) < 1.3):
                Mjj_4 = (wj_ak4[0] + wj_ak4[1]).M()

        # Fill histos
        if Mjj_4 > 0:
            self.h_minv_4_all.Fill(Mjj_4)
            if dst.PFScouting_JetHT == 1:
                self.h_minv_4_passed.Fill(Mjj_4)

            self.n_Rest += 1
            return True

        # Nothing matched (rare), just return False to skip
        return False

    def endJob(self):
        print(f"------------------------------------------------------------")
        print("Summary after processing:")
        print(f"Total events after reference trigger: {self.n_total_events}")
        print(f"VBF:      {self.n_VBF}")
        print(f"Boosted:  {self.n_Boosted}")
        print(f"Resolved: {self.n_Resolved}")
        print(f"Rest:     {self.n_Rest}")
        print(f"Sum categories (should equal total): {self.n_VBF + self.n_Boosted + self.n_Resolved + self.n_Rest}")
        Module.endJob(self)

reference_paths = ["PFScouting_SingleMuon"]
signal_paths    = ["PFScouting_JetHT"]

# To make event loops faster ex. DST_PFScouting_SingleMuon == 1 && DST_PFScouting_JetHT == 1
preselection= ""#DST_PFScouting_SingleMuon == 1 && DST_PFScouting_JetHT == 1"

files=[
        "root://cms-xrd-global.cern.ch///store/data/Run2024I/ScoutingPFRun3/NANOAOD/PromptReco-v2/000/386/694/00000/005815b0-42c7-4968-9666-f11e22b18a13.root",
        ]

## Unused files. If you want to use all of them, you remove three lines here
file2=[##
        ]

### CHANGE NAME OF OUTPUT FILE PRODUCED AFTER SKIM AND SELECTION
p=PostProcessor(".",files,cut=preselection,branchsel=None,modules=[TrigDijetHTAnalysis()],noOut=True,histFileName="histos_DijetHTTrigNanoAOD.root",histDirName="DijethtTrigAnalyzerNanoAOD")

p.run()
