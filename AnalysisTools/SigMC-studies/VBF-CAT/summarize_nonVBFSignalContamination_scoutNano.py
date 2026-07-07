#!/usr/bin/env python3

"""
Run the central-widejet VBF-category selections on non-VBF ScoutNano signals
to quantify how often these signals are selected by the current VBF signal
strategies.

The script applies the same trigger baseline, JECs, jet ID, central-widejet
selection, and forward-pair requirements used in the current
studyCentralWideJetVBF_scoutNano.py flow, then reports per-mass event counts
and efficiencies for:
 - Central-to-VBFTag
 - VBFTag-to-Central
 - only Central-to-VBFTag
 - only VBFTag-to-Central
 - both strategies
 - either strategy
"""

import os
import re
import csv
import math
import glob
import argparse
import itertools
from types import SimpleNamespace

import ROOT
import correctionlib
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

ROOT.gROOT.SetBatch(True)

SAMPLE_PATTERNS = {
    "zprime_qq": [
        "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/ZPrimeToQQ-Central/ZprimeToQQ_Par-M-*.root",
    ],
    "glugluspin0_2b": [
        "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/GluGluSpin0-Central/GluGluSpin0To2B/GluGluSpin0To2B_Par-M-*.root",
    ],
    "rsgraviton_2q": [
        "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/RSGraviton-Central/RSGravTo2Q/RSGravTo2Q_Par-kMpl01-M-*.root",
    ],
    "rsgraviton_2glu": [
        "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/RSGraviton-Central/RSGravTo2Glu/RSGravTo2Glu_Par-kMpl01-M-*.root",
    ],
}

WORKING_POINT_PRESETS = {
    "fwdPt24_VBF-Deta4-Mjj500": {
        "forwardPairAlgo": "maxDEta",
        "forwardJetPtMin": 24.0,
        "forwardPairDEtaMin": 4.0,
        "forwardPairExtraMjjMin": 500.0,
    },
    "fwdPt21_VBF-Deta4p5-Mjj450": {
        "forwardPairAlgo": "maxDEta",
        "forwardJetPtMin": 21.0,
        "forwardPairDEtaMin": 4.5,
        "forwardPairExtraMjjMin": 450.0,
    },
    "fwdPt20_VBF-Deta5": {
        "forwardPairAlgo": "maxDEta",
        "forwardJetPtMin": 20.0,
        "forwardPairDEtaMin": 5.0,
        "forwardPairExtraMjjMin": 0.0,
    },
}

SINGLE_MU_L1_BITS = (
    "L1_SingleMu11_SQ14_BMTF",
    "L1_SingleMu10_SQ14_BMTF",
    "L1_SingleMu9_SQ14_BMTF",
    "L1_SingleMu8_SQ14_BMTF",
    "L1_SingleMu7_SQ14_BMTF",
    "L1_SingleMu6_SQ14_BMTF",
    "L1_SingleMu5_SQ14_BMTF",
)


parser = argparse.ArgumentParser(
    description="Summarize how often non-VBF ScoutNano signals pass the VBF-category selections."
)
parser.add_argument(
    "--inputs",
    nargs="+",
    default=None,
    help="Explicit input ROOT files or glob patterns. Overrides --sampleSet.",
)
parser.add_argument(
    "--sampleSet",
    nargs="+",
    choices=sorted(SAMPLE_PATTERNS.keys()) + ["all"],
    default=["all"],
    help="Named non-VBF ScoutNano sample sets to process when --inputs is not given.",
)
parser.add_argument(
    "--workingPointPreset",
    nargs="+",
    choices=sorted(WORKING_POINT_PRESETS.keys()) + ["all"],
    default=None,
    help="Named VBF working-point presets to run. Overrides the scan options when provided.",
)
parser.add_argument(
    "--outputCsv",
    type=str,
    default="summary_nonVBFSignalContamination_scoutNano.csv",
    help="CSV output path.",
)
parser.add_argument(
    "--outputRoot",
    type=str,
    default="summary_nonVBFSignalContamination_scoutNano.root",
    help="ROOT output path for forward-pair kinematic histograms.",
)
parser.add_argument(
    "--jecJson",
    type=str,
    default=None,
    help="Path to the HLT JEC JSON. If omitted, try CMSSW_BASE and known local fallbacks.",
)
parser.add_argument(
    "--maxFiles",
    type=int,
    default=None,
    help="Maximum number of files to process after expansion and sorting.",
)
parser.add_argument(
    "--maxEntries",
    type=int,
    default=None,
    help="Maximum number of entries to process per file.",
)
parser.add_argument(
    "--firstEntry",
    type=int,
    default=0,
    help="First entry to process.",
)
parser.add_argument(
    "--prefetch",
    action="store_true",
    help="Enable NanoAODTools prefetching.",
)
parser.add_argument(
    "--centralJetPtMin",
    type=float,
    default=30.0,
    help="Minimum pT [GeV] for central seed jets.",
)
parser.add_argument(
    "--centralEtaMax",
    type=float,
    default=2.5,
    help="Maximum |eta| for central seed jets.",
)
parser.add_argument(
    "--centralPairDEtaMin",
    type=float,
    default=0.0,
    help="Minimum |DeltaEta| required for the leading central pair.",
)
parser.add_argument(
    "--centralPairDEtaMax",
    type=float,
    default=1.5,
    help="Maximum |DeltaEta| allowed for the leading central pair.",
)
parser.add_argument(
    "--centralPairDPhiMin",
    type=float,
    default=1.5,
    help="Minimum |DeltaPhi| required for the leading central pair.",
)
parser.add_argument(
    "--leadingPtOverMMin",
    type=float,
    default=0.0,
    help="Minimum leading widejet pT / m(widejet pair).",
)
parser.add_argument(
    "--subleadingPtOverMMin",
    type=float,
    default=0.0,
    help="Minimum subleading widejet pT / m(widejet pair).",
)
parser.add_argument(
    "--widejetRadius",
    type=float,
    default=1.1,
    help="Maximum DeltaR used to absorb reco jets into the wide jets.",
)
parser.add_argument(
    "--forwardJetPtMin",
    type=float,
    default=20.0,
    help="Minimum pT [GeV] for forward VBF-tag jets.",
)
parser.add_argument(
    "--forwardEtaMin",
    type=float,
    default=2.5,
    help="Minimum |eta| for forward VBF-tag jets.",
)
parser.add_argument(
    "--forwardEtaMax",
    type=float,
    default=5.0,
    help="Maximum |eta| for forward VBF-tag jets.",
)
parser.add_argument(
    "--requireForwardEtaMinOnCentralFirst",
    action="store_true",
    help="On the Central-to-VBFTag strategy, require leftover VBF-tag jets to satisfy |eta| >= forwardEtaMin.",
)
parser.add_argument(
    "--forwardPairAlgo",
    type=str,
    choices=["maxDEta", "maxMjj"],
    default="maxDEta",
    help="How to choose the forward pair.",
)
parser.add_argument(
    "--forwardPairDEtaMin",
    type=float,
    default=5.0,
    help="Minimum |DeltaEta| required for the selected forward pair when --forwardPairAlgo=maxDEta.",
)
parser.add_argument(
    "--forwardPairMjjMin",
    type=float,
    default=750.0,
    help="Minimum mjj [GeV] required for the selected forward pair when --forwardPairAlgo=maxMjj.",
)
parser.add_argument(
    "--forwardPairExtraMjjMin",
    type=float,
    default=0.0,
    help="Optional extra mjj [GeV] requirement applied on top of the chosen forward-pair algorithm. Useful with --forwardPairAlgo=maxDEta. Set to 0 to disable.",
)
parser.add_argument(
    "--scanForwardJetPtMin",
    nargs="+",
    type=float,
    default=None,
    help="Optional list of forward-jet pT thresholds to scan in one pass.",
)
parser.add_argument(
    "--scanForwardPairDEtaMin",
    nargs="+",
    type=float,
    default=None,
    help="Optional list of forward-pair |DeltaEta| thresholds to scan in one pass for --forwardPairAlgo=maxDEta.",
)
parser.add_argument(
    "--scanForwardPairExtraMjjMin",
    nargs="+",
    type=float,
    default=None,
    help="Optional list of extra forward-pair mjj thresholds to scan in one pass.",
)
parser.add_argument(
    "--requireTriggerBaseline",
    dest="requireTriggerBaseline",
    action="store_true",
    help="Require the full trigger baseline.",
)
parser.add_argument(
    "--noRequireTriggerBaseline",
    dest="requireTriggerBaseline",
    action="store_false",
    help="Disable the trigger-baseline requirement.",
)
parser.set_defaults(requireTriggerBaseline=True)
args = parser.parse_args()


def build_working_points():
    if args.workingPointPreset:
        preset_names = (
            sorted(WORKING_POINT_PRESETS.keys())
            if "all" in args.workingPointPreset
            else args.workingPointPreset
        )
        working_points = []
        for idx, preset_name in enumerate(preset_names):
            preset = WORKING_POINT_PRESETS[preset_name]
            working_points.append(
                {
                    "index": idx,
                    "label": preset_name,
                    "forwardJetPtMin": float(preset["forwardJetPtMin"]),
                    "forwardPairDEtaMin": float(preset["forwardPairDEtaMin"]),
                    "forwardPairExtraMjjMin": float(preset["forwardPairExtraMjjMin"]),
                    "forwardPairAlgo": str(preset["forwardPairAlgo"]),
                    "forwardPairMjjMin": float(args.forwardPairMjjMin),
                }
            )
        return working_points

    forward_pts = args.scanForwardJetPtMin if args.scanForwardJetPtMin is not None else [args.forwardJetPtMin]
    forward_detas = args.scanForwardPairDEtaMin if args.scanForwardPairDEtaMin is not None else [args.forwardPairDEtaMin]
    extra_mjjs = (
        args.scanForwardPairExtraMjjMin
        if args.scanForwardPairExtraMjjMin is not None
        else [args.forwardPairExtraMjjMin]
    )

    working_points = []
    for idx, (forward_pt, forward_deta, extra_mjj) in enumerate(
        itertools.product(forward_pts, forward_detas, extra_mjjs)
    ):
        label = (
            f"{args.forwardPairAlgo}"
            f"_fwdPt{forward_pt:g}"
            f"_dEta{forward_deta:g}"
            f"_extraMjj{extra_mjj:g}"
        )
        working_points.append(
            {
                "index": idx,
                "label": label,
                "forwardJetPtMin": float(forward_pt),
                "forwardPairDEtaMin": float(forward_deta),
                "forwardPairExtraMjjMin": float(extra_mjj),
                "forwardPairAlgo": args.forwardPairAlgo,
                "forwardPairMjjMin": float(args.forwardPairMjjMin),
            }
        )
    return working_points


WORKING_POINTS = build_working_points()
MIN_FORWARD_JET_PT = min(wp["forwardJetPtMin"] for wp in WORKING_POINTS)


FORWARD_HIST_SPECS = [
    ("forward_pair_dEta", 80, 0.0, 10.0, "Forward pair |#Delta#eta|"),
    ("forward_pair_dPhi", 64, 0.0, 3.2, "Forward pair |#Delta#phi|"),
    ("forward_pair_dR", 80, 0.0, 10.0, "Forward pair #DeltaR"),
    ("forward_pair_mass", 120, 0.0, 4000.0, "Forward pair mass [GeV]"),
    ("forward_jet1_pt", 100, 0.0, 800.0, "Leading forward jet p_{T} [GeV]"),
    ("forward_jet2_pt", 100, 0.0, 600.0, "Subleading forward jet p_{T} [GeV]"),
    ("forward_jet1_abseta", 60, 0.0, 6.0, "Leading forward jet |#eta|"),
    ("forward_jet2_abseta", 60, 0.0, 6.0, "Subleading forward jet |#eta|"),
]


def sanitize_token(text):
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(text)).strip("_")


def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2.0 * math.pi
    while dphi <= -math.pi:
        dphi += 2.0 * math.pi
    return dphi


def delta_r(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = delta_phi(phi1, phi2)
    return math.sqrt(deta * deta + dphi * dphi)


def recojet_mass(jet):
    for attr in ["mass", "m"]:
        try:
            return float(getattr(jet, attr))
        except Exception:
            continue
    return 0.0


def corrected_recojet(jet, pt_corr):
    return SimpleNamespace(
        original=jet,
        pt=pt_corr,
        eta=jet.eta,
        phi=jet.phi,
        mass=recojet_mass(jet),
        neutralHadronEnergy=jet.neutralHadronEnergy,
        HFHadronEnergy=jet.HFHadronEnergy,
        photonEnergy=jet.photonEnergy,
        HFEMEnergy=jet.HFEMEnergy,
        muonEnergy=jet.muonEnergy,
        electronEnergy=jet.electronEnergy,
        chargedHadronEnergy=jet.chargedHadronEnergy,
        chargedHadronMultiplicity=jet.chargedHadronMultiplicity,
        HFHadronMultiplicity=jet.HFHadronMultiplicity,
        neutralHadronMultiplicity=jet.neutralHadronMultiplicity,
        HFEMMultiplicity=jet.HFEMMultiplicity,
        muonMultiplicity=jet.muonMultiplicity,
        electronMultiplicity=jet.electronMultiplicity,
        photonMultiplicity=jet.photonMultiplicity,
    )


def p4_from_recojet(jet):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, recojet_mass(jet))
    return vec


def get_trigger_bits(event):
    pass_dst = bool(getattr(event, "DST_PFScouting_JetHT", False))
    pass_htt280 = bool(getattr(event, "L1_HTT280er", False))
    pass_singlejet180 = bool(getattr(event, "L1_SingleJet180", False))
    pass_singlemu_veto = not any(
        bool(getattr(event, bit_name, False))
        for bit_name in SINGLE_MU_L1_BITS
    )
    pass_baseline = pass_dst and (pass_htt280 or pass_singlejet180) and pass_singlemu_veto
    return {
        "pass_baseline": pass_baseline,
    }


def jet_id(jet):
    aeta = abs(jet.eta)
    total_e = (
        jet.neutralHadronEnergy + jet.HFHadronEnergy +
        jet.photonEnergy + jet.HFEMEnergy +
        jet.muonEnergy + jet.electronEnergy +
        jet.chargedHadronEnergy
    )
    if total_e <= 0:
        return False

    nhf = (jet.neutralHadronEnergy + jet.HFHadronEnergy) / float(total_e)
    nemf = (jet.photonEnergy + jet.HFEMEnergy) / float(total_e)
    mu_frac = jet.muonEnergy / float(total_e)

    charged_mult = jet.chargedHadronMultiplicity + jet.HFHadronMultiplicity
    neutral_mult = jet.neutralHadronMultiplicity + jet.HFEMMultiplicity
    nconst = (
        jet.chargedHadronMultiplicity + jet.neutralHadronMultiplicity +
        jet.muonMultiplicity + jet.electronMultiplicity +
        jet.photonMultiplicity
    )

    if aeta < 2.6:
        if nhf >= 0.99 or nemf >= 0.90 or nconst <= 1:
            return False
        if charged_mult <= 0 or mu_frac >= 0.80:
            return False
        return True

    if 2.6 <= aeta < 2.7:
        if nemf >= 0.99 or mu_frac >= 0.80:
            return False
        return True

    if 2.7 <= aeta < 3.0:
        if nemf >= 0.99 or neutral_mult <= 1:
            return False
        return True

    if 3.0 <= aeta < 5.0:
        if nemf >= 0.10:
            return False
        return True

    return False


def build_widejets_from_seeds(jets, seed1, seed2, widejet_radius):
    widejet1 = ROOT.TLorentzVector(p4_from_recojet(seed1))
    widejet2 = ROOT.TLorentzVector(p4_from_recojet(seed2))
    assigned_to_widejet1_indices = [id(seed1)]
    assigned_to_widejet2_indices = [id(seed2)]
    leftover_jets = []

    for jet in jets:
        if id(jet) in {id(seed1), id(seed2)}:
            continue

        dR1 = delta_r(jet.eta, jet.phi, seed1.eta, seed1.phi)
        dR2 = delta_r(jet.eta, jet.phi, seed2.eta, seed2.phi)

        if dR1 < widejet_radius and dR1 < dR2:
            widejet1 += p4_from_recojet(jet)
            assigned_to_widejet1_indices.append(id(jet))
        elif dR2 < widejet_radius and dR2 < dR1:
            widejet2 += p4_from_recojet(jet)
            assigned_to_widejet2_indices.append(id(jet))
        else:
            leftover_jets.append(jet)

    return {
        "widejet1": widejet1,
        "widejet2": widejet2,
        "assigned_to_widejet1_ids": assigned_to_widejet1_indices,
        "assigned_to_widejet2_ids": assigned_to_widejet2_indices,
        "leftover_jets": leftover_jets,
    }


def choose_max_deta_pair(jets):
    if len(jets) < 2:
        return None

    best_pair = None
    best_deta = -1.0
    for i in range(len(jets)):
        for j in range(i + 1, len(jets)):
            j1, j2 = jets[i], jets[j]
            deta = abs(j1.eta - j2.eta)
            dphi = abs(delta_phi(j1.phi, j2.phi))
            dr = delta_r(j1.eta, j1.phi, j2.eta, j2.phi)
            mjj = (p4_from_recojet(j1) + p4_from_recojet(j2)).M()
            if deta > best_deta:
                best_deta = deta
                best_pair = (j1, j2, deta, dphi, dr, mjj)
    return best_pair


def choose_max_mjj_pair(jets):
    if len(jets) < 2:
        return None

    best_pair = None
    best_mjj = -1.0
    for i in range(len(jets)):
        for j in range(i + 1, len(jets)):
            j1, j2 = jets[i], jets[j]
            deta = abs(j1.eta - j2.eta)
            dphi = abs(delta_phi(j1.phi, j2.phi))
            dr = delta_r(j1.eta, j1.phi, j2.eta, j2.phi)
            mjj = (p4_from_recojet(j1) + p4_from_recojet(j2)).M()
            if mjj > best_mjj:
                best_mjj = mjj
                best_pair = (j1, j2, deta, dphi, dr, mjj)
    return best_pair


def choose_forward_pair_by_mode(jets, mode):
    if mode == "maxDEta":
        return choose_max_deta_pair(jets)
    if mode == "maxMjj":
        return choose_max_mjj_pair(jets)
    raise ValueError(f"Unsupported forward pair mode: {mode}")


def pass_vbf_requirement(pair, working_point=None):
    if pair is None:
        return False

    if working_point is None:
        working_point = {
            "forwardPairAlgo": args.forwardPairAlgo,
            "forwardPairDEtaMin": args.forwardPairDEtaMin,
            "forwardPairMjjMin": args.forwardPairMjjMin,
            "forwardPairExtraMjjMin": args.forwardPairExtraMjjMin,
        }

    _, _, deta, _, _, mjj = pair
    pass_pair_algo = False
    if working_point["forwardPairAlgo"] == "maxDEta":
        pass_pair_algo = deta >= working_point["forwardPairDEtaMin"]
    elif working_point["forwardPairAlgo"] == "maxMjj":
        pass_pair_algo = mjj >= working_point["forwardPairMjjMin"]

    if not pass_pair_algo:
        return False
    if working_point["forwardPairExtraMjjMin"] > 0.0 and mjj < working_point["forwardPairExtraMjjMin"]:
        return False
    return True


def book_forward_hist_set(tag):
    hists = {}
    for name, nbins, xmin, xmax, xtitle in FORWARD_HIST_SPECS:
        hist_name = f"h_{tag}_{name}"
        hists[name] = ROOT.TH1F(hist_name, f";{xtitle};Events", nbins, xmin, xmax)
        hists[name].Sumw2()
        hists[name].SetDirectory(0)
    return hists


def fill_forward_hist_set(hists, pair):
    if pair is None:
        return

    jet1, jet2, deta, dphi, dr, mjj = pair
    hists["forward_pair_dEta"].Fill(deta)
    hists["forward_pair_dPhi"].Fill(dphi)
    hists["forward_pair_dR"].Fill(dr)
    hists["forward_pair_mass"].Fill(mjj)

    jet1_pt, jet2_pt = sorted([jet1.pt, jet2.pt], reverse=True)
    jet1_abseta, jet2_abseta = sorted([abs(jet1.eta), abs(jet2.eta)], reverse=True)
    hists["forward_jet1_pt"].Fill(jet1_pt)
    hists["forward_jet2_pt"].Fill(jet2_pt)
    hists["forward_jet1_abseta"].Fill(jet1_abseta)
    hists["forward_jet2_abseta"].Fill(jet2_abseta)


def infer_mass(path):
    patterns = [
        r"_M(\d+)\.root$",
        r"Par-M-(\d+)\.root$",
        r"Par-[^-]+-M-(\d+)\.root$",
    ]
    base = os.path.basename(path)
    for pattern in patterns:
        match = re.search(pattern, base)
        if match:
            return int(match.group(1))
    return -1


def infer_sample_label(path):
    if "ZPrimeToQQ" in path or "ZprimeToQQ" in path:
        return "zprime_qq"
    if "GluGluSpin0To2B" in path:
        return "glugluspin0_2b"
    if "RSGravTo2Glu" in path or "RSGravitonToGluGlu" in path:
        return "rsgraviton_2glu"
    if "RSGravTo2Q" in path or "RSGravitonToQQ" in path:
        return "rsgraviton_2q"
    return "unknown"


def expand_input_files():
    files = []
    if args.inputs:
        for pattern in args.inputs:
            matches = sorted(glob.glob(pattern))
            if matches:
                files.extend(matches)
            elif pattern.endswith(".root"):
                files.append(pattern)
    else:
        sample_sets = list(SAMPLE_PATTERNS.keys()) if "all" in args.sampleSet else args.sampleSet
        for sample_set in sample_sets:
            for pattern in SAMPLE_PATTERNS[sample_set]:
                files.extend(sorted(glob.glob(pattern)))

    files = sorted(set(files))
    if args.maxFiles is not None:
        files = files[:args.maxFiles]
    if not files:
        raise RuntimeError("No input files found.")
    return files


def jec_json_path():
    candidates = []
    if args.jecJson:
        candidates.append(args.jecJson)

    cmssw_base = os.getenv("CMSSW_BASE", "")
    if cmssw_base:
        candidates.append(
            os.path.join(cmssw_base, "src", "2024_UtilsDataQuality", "jetHLT_jerc.json")
        )

    candidates.append(
        "/afs/cern.ch/work/e/elfontan/private/dijetAnalysis_ScoutingRun3/BKGModelling/CMSSW_14_1_0_pre4/src/2024_UtilsDataQuality/jetHLT_jerc.json"
    )

    for candidate in candidates:
        if candidate and os.path.exists(candidate):
            return candidate

    return candidates[0] if candidates else ""


class NonVBFSignalContaminationModule(Module):
    def __init__(self):
        self.writeHistFile = False
        self.working_points = WORKING_POINTS
        self.total = 0
        self.after_trigger = 0
        self.counters = {
            wp["label"]: {
                "n_central_pass_central": 0,
                "n_central_has_vbf_pair": 0,
                "n_central_selected": 0,
            }
            for wp in self.working_points
        }
        self.forward_hists = {
            wp["label"]: {
                "all_pairs": book_forward_hist_set(f"{sanitize_token(wp['label'])}_allPairs"),
                "selected_pairs": book_forward_hist_set(f"{sanitize_token(wp['label'])}_selectedPairs"),
                "rejected_pairs": book_forward_hist_set(f"{sanitize_token(wp['label'])}_rejectedPairs"),
            }
            for wp in self.working_points
        }

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)
        jec_json = jec_json_path()
        if not os.path.exists(jec_json):
            raise RuntimeError(f"HLT JEC JSON not found: {jec_json}")
        cset = correctionlib.CorrectionSet.from_file(jec_json)
        self.jec = cset.compound["HLT_Winter24_V1_MC_L1L2L3Res_AK4PFHLT"]

    def corrected_jets(self, event, jets_raw):
        rho = float(getattr(event, "ScoutingRho_fixedGridRhoFastjetAll", 0.0))
        run = int(getattr(event, "run", 1))
        corrected = []
        for jet in jets_raw:
            variables = {
                "JetPt": jet.pt,
                "JetEta": jet.eta,
                "JetPhi": jet.phi,
                "JetA": getattr(jet, "jetArea", 0.0),
                "Rho": rho,
                "run": run,
            }
            inputs = [variables[inp.name] for inp in self.jec.inputs]
            corr_factor = self.jec.evaluate(*inputs)
            corrected.append(corrected_recojet(jet, jet.pt * corr_factor))
        return corrected

    def analyze(self, event):
        self.total += 1
        trigger_bits = get_trigger_bits(event)
        if args.requireTriggerBaseline and not trigger_bits["pass_baseline"]:
            return True

        if trigger_bits["pass_baseline"]:
            self.after_trigger += 1

        jets_raw = Collection(event, "ScoutingPFJet")
        jets = self.corrected_jets(event, jets_raw)

        good_jets = [
            jet for jet in jets
            if jet_id(jet) and jet.pt > min(args.centralJetPtMin, MIN_FORWARD_JET_PT) and abs(jet.eta) < args.forwardEtaMax
        ]
        central_jets = [
            jet for jet in good_jets
            if jet.pt > args.centralJetPtMin and abs(jet.eta) < args.centralEtaMax
        ]
        central_jets = sorted(central_jets, key=lambda x: x.pt, reverse=True)
        pass_central_selection = False
        remaining_jets = None

        if len(central_jets) >= 2:
            jc1, jc2 = central_jets[0], central_jets[1]
            central_deta = abs(jc1.eta - jc2.eta)
            central_dphi = abs(delta_phi(jc1.phi, jc2.phi))

            widejet_info = build_widejets_from_seeds(good_jets, jc1, jc2, args.widejetRadius)
            widejet1 = widejet_info["widejet1"]
            widejet2 = widejet_info["widejet2"]
            widejet_pair_mass = (widejet1 + widejet2).M()
            leading_widejet_pt = max(widejet1.Pt(), widejet2.Pt())
            subleading_widejet_pt = min(widejet1.Pt(), widejet2.Pt())

            pass_central_selection = True
            if args.centralPairDEtaMin > 0.0 and central_deta < args.centralPairDEtaMin:
                pass_central_selection = False
            if args.centralPairDEtaMax > 0.0 and central_deta > args.centralPairDEtaMax:
                pass_central_selection = False
            if args.centralPairDPhiMin > 0.0 and central_dphi < args.centralPairDPhiMin:
                pass_central_selection = False
            if args.leadingPtOverMMin > 0.0:
                if widejet_pair_mass <= 0.0 or (leading_widejet_pt / widejet_pair_mass) < args.leadingPtOverMMin:
                    pass_central_selection = False
            if args.subleadingPtOverMMin > 0.0:
                if widejet_pair_mass <= 0.0 or (subleading_widejet_pt / widejet_pair_mass) < args.subleadingPtOverMMin:
                    pass_central_selection = False

            if pass_central_selection:
                remaining_jets = widejet_info["leftover_jets"]

        for working_point in self.working_points:
            counts = self.counters[working_point["label"]]
            central_selected = False

            if pass_central_selection:
                counts["n_central_pass_central"] += 1
                forward_remaining_jets = [
                    jet for jet in remaining_jets
                    if (
                        jet.pt > working_point["forwardJetPtMin"]
                        and abs(jet.eta) <= args.forwardEtaMax
                        and (
                            not args.requireForwardEtaMinOnCentralFirst
                            or abs(jet.eta) >= args.forwardEtaMin
                        )
                    )
                ]
                vbf_pair = choose_forward_pair_by_mode(
                    forward_remaining_jets,
                    working_point["forwardPairAlgo"],
                )
                if vbf_pair is not None:
                    counts["n_central_has_vbf_pair"] += 1
                    fill_forward_hist_set(self.forward_hists[working_point["label"]]["all_pairs"], vbf_pair)
                if vbf_pair is not None:
                    if pass_vbf_requirement(vbf_pair, working_point):
                        central_selected = True
                        fill_forward_hist_set(self.forward_hists[working_point["label"]]["selected_pairs"], vbf_pair)
                    else:
                        fill_forward_hist_set(self.forward_hists[working_point["label"]]["rejected_pairs"], vbf_pair)

            if central_selected:
                counts["n_central_selected"] += 1

        return True

    def summary(self):
        def frac(count, denom):
            return 100.0 * count / denom if denom > 0 else 0.0

        ref = self.after_trigger if args.requireTriggerBaseline else self.total
        rows = []
        for working_point in self.working_points:
            counts = self.counters[working_point["label"]]
            rows.append(
                {
                    "working_point": working_point["label"],
                    "forwardPairAlgo": working_point["forwardPairAlgo"],
                    "forwardJetPtMin": working_point["forwardJetPtMin"],
                    "forwardPairDEtaMin": working_point["forwardPairDEtaMin"],
                    "forwardPairMjjMin": working_point["forwardPairMjjMin"],
                    "forwardPairExtraMjjMin": working_point["forwardPairExtraMjjMin"],
                    "n_total": self.total,
                    "n_after_trigger": self.after_trigger,
                    "n_central_pass_central": counts["n_central_pass_central"],
                    "n_central_has_vbf_pair": counts["n_central_has_vbf_pair"],
                    "n_central_selected": counts["n_central_selected"],
                    "eff_central_pct": frac(counts["n_central_selected"], ref),
                }
            )
        return rows

    def cloned_forward_hists(self):
        cloned = {}
        for working_point in self.working_points:
            wp_label = working_point["label"]
            cloned[wp_label] = {}
            for hist_group, hist_map in self.forward_hists[wp_label].items():
                cloned[wp_label][hist_group] = {}
                for hist_name, hist in hist_map.items():
                    hist_clone = hist.Clone(f"{hist.GetName()}_clone")
                    hist_clone.SetDirectory(0)
                    cloned[wp_label][hist_group][hist_name] = hist_clone
        return cloned


def process_file(path):
    module = NonVBFSignalContaminationModule()
    processor = PostProcessor(
        ".",
        [path],
        modules=[module],
        noOut=True,
        maxEntries=args.maxEntries,
        firstEntry=args.firstEntry,
        prefetch=args.prefetch,
    )
    processor.run()

    rows = []
    for summary_row in module.summary():
        row = {
            "sample": infer_sample_label(path),
            "mass": infer_mass(path),
            "file": path,
        }
        row.update(summary_row)
        rows.append(row)
    return rows, module.cloned_forward_hists()


def write_root(path, payloads):
    def get_or_make_dir(parent_dir, name):
        existing = parent_dir.GetDirectory(name)
        if existing:
            return existing
        return parent_dir.mkdir(name)

    root_file = ROOT.TFile(path, "RECREATE")
    try:
        for payload in payloads:
            sample = payload["sample"]
            mass = payload["mass"]
            histograms = payload["histograms"]

            sample_dir = get_or_make_dir(root_file, sanitize_token(sample))
            mass_dir = get_or_make_dir(sample_dir, f"mass_{mass}")
            for working_point_label, hist_groups in histograms.items():
                wp_dir = get_or_make_dir(mass_dir, sanitize_token(working_point_label))
                for hist_group, hist_map in hist_groups.items():
                    group_dir = get_or_make_dir(wp_dir, hist_group)
                    group_dir.cd()
                    for hist_name, hist in hist_map.items():
                        hist.SetName(hist_name)
                        hist.Write()
                    wp_dir.cd()
            root_file.cd()
    finally:
        root_file.Close()


def write_csv(path, rows):
    fieldnames = [
        "sample",
        "mass",
        "file",
        "working_point",
        "forwardPairAlgo",
        "forwardJetPtMin",
        "forwardPairDEtaMin",
        "forwardPairMjjMin",
        "forwardPairExtraMjjMin",
        "n_total",
        "n_after_trigger",
        "n_central_pass_central",
        "n_central_has_vbf_pair",
        "n_central_selected",
        "eff_central_pct",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def print_rows(rows):
    print(
        f"{'sample':<16s} {'mass':>6s} {'fwdPt':>7s} {'dEta':>7s} {'xMjj':>7s} {'afterTrig':>10s} "
        f"{'Sel%':>8s}"
    )
    for row in rows:
        ref = row["n_after_trigger"] if args.requireTriggerBaseline else row["n_total"]
        print(
            f"{row['sample']:<16s} "
            f"{row['mass']:6d} "
            f"{row['forwardJetPtMin']:7.1f} "
            f"{row['forwardPairDEtaMin']:7.2f} "
            f"{row['forwardPairExtraMjjMin']:7.1f} "
            f"{ref:10d} "
            f"{row['eff_central_pct']:8.3f} "
        )


def main():
    files = expand_input_files()
    rows = []
    root_payloads = []
    for path in files:
        file_rows, file_hists = process_file(path)
        rows.extend(file_rows)
        root_payloads.append(
            {
                "sample": infer_sample_label(path),
                "mass": infer_mass(path),
                "histograms": file_hists,
            }
        )
    rows.sort(
        key=lambda row: (
            row["sample"],
            row["mass"],
            row["forwardJetPtMin"],
            row["forwardPairDEtaMin"],
            row["forwardPairExtraMjjMin"],
            row["file"],
        )
    )

    print_rows(rows)
    write_csv(args.outputCsv, rows)
    write_root(args.outputRoot, root_payloads)
    print(f"\n[INFO] Wrote CSV summary to: {args.outputCsv}")
    print(f"[INFO] Wrote ROOT histograms to: {args.outputRoot}")


if __name__ == "__main__":
    main()
