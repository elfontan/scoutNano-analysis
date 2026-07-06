#!/usr/bin/env python3

"""
Joint forward-cut optimization on one VBF and one non-VBF 500 GeV ScoutNano sample.

It scans the Central-to-VBFTag forward selection over:
 - forwardJetPtMin
 - forwardPairDEtaMin
 - forwardPairExtraMjjMin

and reports, for each working point:
 - VBF selected efficiency [%]
 - non-VBF selected efficiency [%]
 - simple ranking metrics to balance the two
"""

import os
import csv
import math
import argparse
import itertools
from types import SimpleNamespace

import ROOT
import correctionlib
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

ROOT.gROOT.SetBatch(True)


DEFAULT_VBF_INPUT = "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/VBFH-Central/VBFHTo2B/VBFH-Hto2B_Par-M-500.root"
DEFAULT_NONVBF_INPUT = "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/GluGluSpin0-Central/GluGluSpin0To2B_Par-M-500.root"

SINGLE_MU_L1_BITS = (
    "L1_SingleMu11_SQ14_BMTF",
    "L1_SingleMu10_SQ14_BMTF",
    "L1_SingleMu9_SQ14_BMTF",
    "L1_SingleMu8_SQ14_BMTF",
    "L1_SingleMu7_SQ14_BMTF",
    "L1_SingleMu6_SQ14_BMTF",
    "L1_SingleMu5_SQ14_BMTF",
)


parser = argparse.ArgumentParser(description="Jointly optimize forward VBF cuts on one VBF and one non-VBF 500 GeV ScoutNano sample")
parser.add_argument("--vbfInput", default=DEFAULT_VBF_INPUT, help="500 GeV VBF signal input ROOT file")
parser.add_argument("--nonVBFInput", default=DEFAULT_NONVBF_INPUT, help="500 GeV non-VBF signal input ROOT file")
parser.add_argument("--outputCsv", default="scan_jointVBF_vs_nonVBF_500_scoutNano.csv", help="Output CSV path")
parser.add_argument("--outputParetoCsv", default="scan_jointVBF_vs_nonVBF_500_pareto.csv", help="Output CSV path for Pareto-optimal working points")
parser.add_argument("--jecJson", default=None, help="Optional explicit JEC JSON path")
parser.add_argument("--maxEntries", type=int, default=None, help="Maximum entries per file")
parser.add_argument("--firstEntry", type=int, default=0, help="First entry per file")
parser.add_argument("--prefetch", action="store_true", help="Enable NanoAODTools prefetching")
parser.add_argument("--requireTriggerBaseline", dest="requireTriggerBaseline", action="store_true", help="Require the full trigger baseline")
parser.add_argument("--noRequireTriggerBaseline", dest="requireTriggerBaseline", action="store_false", help="Disable the trigger-baseline requirement")
parser.set_defaults(requireTriggerBaseline=True)

parser.add_argument("--centralJetPtMin", type=float, default=30.0)
parser.add_argument("--centralEtaMax", type=float, default=2.5)
parser.add_argument("--centralPairDEtaMin", type=float, default=0.0)
parser.add_argument("--centralPairDEtaMax", type=float, default=0.0)
parser.add_argument("--centralPairDPhiMin", type=float, default=0.0)
parser.add_argument("--widejetPairDEtaMin", type=float, default=0.0)
parser.add_argument("--widejetPairDEtaMax", type=float, default=1.3)
parser.add_argument("--widejetPairDPhiMin", type=float, default=0.0)
parser.add_argument("--leadingPtOverMMin", type=float, default=0.0)
parser.add_argument("--subleadingPtOverMMin", type=float, default=0.0)
parser.add_argument("--widejetRadius", type=float, default=1.1)
parser.add_argument("--forwardEtaMin", type=float, default=2.5)
parser.add_argument("--forwardEtaMax", type=float, default=5.0)
parser.add_argument("--requireForwardEtaMinOnCentralFirst", action="store_true")
parser.add_argument("--forwardPairAlgo", choices=["maxDEta"], default="maxDEta", help="Forward-pair choice mode. This joint optimizer is built for maxDEta.")

parser.add_argument("--scanForwardJetPtMinMin", type=float, default=20.0)
parser.add_argument("--scanForwardJetPtMinMax", type=float, default=30.0)
parser.add_argument("--scanForwardJetPtMinStep", type=float, default=1.0)
parser.add_argument("--scanForwardPairDEtaMinMin", type=float, default=4.0)
parser.add_argument("--scanForwardPairDEtaMinMax", type=float, default=6.0)
parser.add_argument("--scanForwardPairDEtaMinStep", type=float, default=0.25)
parser.add_argument("--scanForwardPairExtraMjjMinMin", type=float, default=0.0)
parser.add_argument("--scanForwardPairExtraMjjMinMax", type=float, default=500.0)
parser.add_argument("--scanForwardPairExtraMjjMinStep", type=float, default=50.0)
parser.add_argument("--nonVBFPenalty", type=float, default=1.0, help="Score = eff_vbf_pct - penalty * eff_nonvbf_pct")
args = parser.parse_args()


def frange(start, stop, step):
    values = []
    current = start
    while current <= stop + 1e-9:
        values.append(round(current, 6))
        current += step
    return values


def build_working_points():
    working_points = []
    for idx, (fwd_pt, deta, extra_mjj) in enumerate(
        itertools.product(
            frange(args.scanForwardJetPtMinMin, args.scanForwardJetPtMinMax, args.scanForwardJetPtMinStep),
            frange(args.scanForwardPairDEtaMinMin, args.scanForwardPairDEtaMinMax, args.scanForwardPairDEtaMinStep),
            frange(args.scanForwardPairExtraMjjMinMin, args.scanForwardPairExtraMjjMinMax, args.scanForwardPairExtraMjjMinStep),
        )
    ):
        working_points.append(
            {
                "index": idx,
                "label": f"maxDEta_fwdPt{fwd_pt:g}_dEta{deta:g}_extraMjj{extra_mjj:g}",
                "forwardJetPtMin": float(fwd_pt),
                "forwardPairDEtaMin": float(deta),
                "forwardPairExtraMjjMin": float(extra_mjj),
            }
        )
    return working_points


WORKING_POINTS = build_working_points()
MIN_FORWARD_JET_PT = min(wp["forwardJetPtMin"] for wp in WORKING_POINTS)


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
    pass_singlemu_veto = not any(bool(getattr(event, bit_name, False)) for bit_name in SINGLE_MU_L1_BITS)
    return {"pass_baseline": pass_dst and (pass_htt280 or pass_singlejet180) and pass_singlemu_veto}


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
        return not (nemf >= 0.99 or mu_frac >= 0.80)
    if 2.7 <= aeta < 3.0:
        return not (nemf >= 0.99 or neutral_mult <= 1)
    if 3.0 <= aeta < 5.0:
        return nemf < 0.10
    return False


def build_widejets_from_seeds(jets, seed1, seed2, widejet_radius):
    widejet1 = ROOT.TLorentzVector(p4_from_recojet(seed1))
    widejet2 = ROOT.TLorentzVector(p4_from_recojet(seed2))
    leftover_jets = []
    for jet in jets:
        if id(jet) in {id(seed1), id(seed2)}:
            continue
        dR1 = delta_r(jet.eta, jet.phi, seed1.eta, seed1.phi)
        dR2 = delta_r(jet.eta, jet.phi, seed2.eta, seed2.phi)
        if dR1 < widejet_radius and dR1 < dR2:
            widejet1 += p4_from_recojet(jet)
        elif dR2 < widejet_radius and dR2 < dR1:
            widejet2 += p4_from_recojet(jet)
        else:
            leftover_jets.append(jet)
    return {"widejet1": widejet1, "widejet2": widejet2, "leftover_jets": leftover_jets}


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


def pass_forward_selection(pair, working_point):
    if pair is None:
        return False
    _, _, deta, _, _, mjj = pair
    if deta < working_point["forwardPairDEtaMin"]:
        return False
    if working_point["forwardPairExtraMjjMin"] > 0.0 and mjj < working_point["forwardPairExtraMjjMin"]:
        return False
    return True


def jec_json_path():
    candidates = []
    if args.jecJson:
        candidates.append(args.jecJson)
    cmssw_base = os.getenv("CMSSW_BASE", "")
    if cmssw_base:
        candidates.append(os.path.join(cmssw_base, "src", "2024_UtilsDataQuality", "jetHLT_jerc.json"))
    candidates.append("/afs/cern.ch/work/e/elfontan/private/dijetAnalysis_ScoutingRun3/BKGModelling/CMSSW_14_1_0_pre4/src/2024_UtilsDataQuality/jetHLT_jerc.json")
    for candidate in candidates:
        if candidate and os.path.exists(candidate):
            return candidate
    return candidates[0] if candidates else ""


class JointScanModule(Module):
    def __init__(self, sample_label):
        self.writeHistFile = False
        self.sample_label = sample_label
        self.total = 0
        self.after_trigger = 0
        self.counters = {
            wp["label"]: {
                "n_step1": 0,
                "n_has_forward_pair": 0,
                "n_selected": 0,
            }
            for wp in WORKING_POINTS
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
        if len(central_jets) < 2:
            return True

        jc1, jc2 = central_jets[0], central_jets[1]
        central_deta = abs(jc1.eta - jc2.eta)
        central_dphi = abs(delta_phi(jc1.phi, jc2.phi))
        widejet_info = build_widejets_from_seeds(good_jets, jc1, jc2, args.widejetRadius)
        widejet1 = widejet_info["widejet1"]
        widejet2 = widejet_info["widejet2"]
        widejet_pair_deta = abs(widejet1.Eta() - widejet2.Eta())
        widejet_pair_dphi = abs(delta_phi(widejet1.Phi(), widejet2.Phi()))
        widejet_pair_mass = (widejet1 + widejet2).M()
        leading_widejet_pt = max(widejet1.Pt(), widejet2.Pt())
        subleading_widejet_pt = min(widejet1.Pt(), widejet2.Pt())

        pass_step1 = True
        if args.centralPairDEtaMin > 0.0 and central_deta < args.centralPairDEtaMin:
            pass_step1 = False
        if args.centralPairDEtaMax > 0.0 and central_deta > args.centralPairDEtaMax:
            pass_step1 = False
        if args.centralPairDPhiMin > 0.0 and central_dphi < args.centralPairDPhiMin:
            pass_step1 = False
        if args.leadingPtOverMMin > 0.0 and (widejet_pair_mass <= 0.0 or (leading_widejet_pt / widejet_pair_mass) < args.leadingPtOverMMin):
            pass_step1 = False
        if args.subleadingPtOverMMin > 0.0 and (widejet_pair_mass <= 0.0 or (subleading_widejet_pt / widejet_pair_mass) < args.subleadingPtOverMMin):
            pass_step1 = False
        if args.widejetPairDEtaMin > 0.0 and widejet_pair_deta < args.widejetPairDEtaMin:
            pass_step1 = False
        if args.widejetPairDEtaMax > 0.0 and widejet_pair_deta > args.widejetPairDEtaMax:
            pass_step1 = False
        if args.widejetPairDPhiMin > 0.0 and widejet_pair_dphi < args.widejetPairDPhiMin:
            pass_step1 = False
        if not pass_step1:
            return True

        remaining_jets = widejet_info["leftover_jets"]
        for wp in WORKING_POINTS:
            counters = self.counters[wp["label"]]
            counters["n_step1"] += 1
            forward_remaining_jets = [
                jet for jet in remaining_jets
                if (
                    jet.pt > wp["forwardJetPtMin"]
                    and abs(jet.eta) <= args.forwardEtaMax
                    and (
                        not args.requireForwardEtaMinOnCentralFirst
                        or abs(jet.eta) >= args.forwardEtaMin
                    )
                )
            ]
            pair = choose_max_deta_pair(forward_remaining_jets)
            if pair is not None:
                counters["n_has_forward_pair"] += 1
            if pair is not None and pass_forward_selection(pair, wp):
                counters["n_selected"] += 1
        return True

    def summary_rows(self):
        ref = self.after_trigger if args.requireTriggerBaseline else self.total
        rows = []
        for wp in WORKING_POINTS:
            counts = self.counters[wp["label"]]
            rows.append(
                {
                    "sample_kind": self.sample_label,
                    "working_point": wp["label"],
                    "forwardJetPtMin": wp["forwardJetPtMin"],
                    "forwardPairDEtaMin": wp["forwardPairDEtaMin"],
                    "forwardPairExtraMjjMin": wp["forwardPairExtraMjjMin"],
                    "n_total": self.total,
                    "n_after_trigger": self.after_trigger,
                    "n_step1": counts["n_step1"],
                    "n_has_forward_pair": counts["n_has_forward_pair"],
                    "n_selected": counts["n_selected"],
                    "eff_selected_pct": 100.0 * counts["n_selected"] / ref if ref > 0 else 0.0,
                }
            )
        return rows


def run_one_file(input_path, sample_label):
    module = JointScanModule(sample_label)
    processor = PostProcessor(
        ".",
        [input_path],
        modules=[module],
        noOut=True,
        maxEntries=args.maxEntries,
        firstEntry=args.firstEntry,
        prefetch=args.prefetch,
    )
    processor.run()
    return module.summary_rows()


def combine_rows(vbf_rows, nonvbf_rows):
    nonvbf_map = {row["working_point"]: row for row in nonvbf_rows}
    rows = []
    for vbf in vbf_rows:
        wp = vbf["working_point"]
        nonvbf = nonvbf_map[wp]
        eff_vbf = vbf["eff_selected_pct"]
        eff_nonvbf = nonvbf["eff_selected_pct"]
        score = eff_vbf - args.nonVBFPenalty * eff_nonvbf
        eff_ratio = eff_vbf / eff_nonvbf if eff_nonvbf > 0 else 999.0
        significance_like = eff_vbf / math.sqrt(max(eff_nonvbf, 1e-9))
        purity_like = eff_vbf / max(eff_vbf + eff_nonvbf, 1e-9)
        rows.append(
            {
                "working_point": wp,
                "forwardJetPtMin": vbf["forwardJetPtMin"],
                "forwardPairDEtaMin": vbf["forwardPairDEtaMin"],
                "forwardPairExtraMjjMin": vbf["forwardPairExtraMjjMin"],
                "vbf_n_total": vbf["n_total"],
                "vbf_n_after_trigger": vbf["n_after_trigger"],
                "vbf_n_step1": vbf["n_step1"],
                "vbf_n_has_forward_pair": vbf["n_has_forward_pair"],
                "vbf_n_selected": vbf["n_selected"],
                "vbf_eff_selected_pct": eff_vbf,
                "nonvbf_n_total": nonvbf["n_total"],
                "nonvbf_n_after_trigger": nonvbf["n_after_trigger"],
                "nonvbf_n_step1": nonvbf["n_step1"],
                "nonvbf_n_has_forward_pair": nonvbf["n_has_forward_pair"],
                "nonvbf_n_selected": nonvbf["n_selected"],
                "nonvbf_eff_selected_pct": eff_nonvbf,
                "eff_diff_pct": eff_vbf - eff_nonvbf,
                "eff_ratio": eff_ratio,
                "score": score,
                "significance_like": significance_like,
                "purity_like": purity_like,
                "score_definition": f"score = vbf_eff_selected_pct - {args.nonVBFPenalty:g} * nonvbf_eff_selected_pct",
                "significance_like_definition": "significance_like = vbf_eff_selected_pct / sqrt(nonvbf_eff_selected_pct)",
                "purity_like_definition": "purity_like = vbf_eff_selected_pct / (vbf_eff_selected_pct + nonvbf_eff_selected_pct)",
            }
        )
    return rows


def pareto_frontier(rows):
    frontier = []
    for row in rows:
        dominated = False
        for other in rows:
            if other is row:
                continue
            if (
                other["vbf_eff_selected_pct"] >= row["vbf_eff_selected_pct"]
                and other["nonvbf_eff_selected_pct"] <= row["nonvbf_eff_selected_pct"]
                and (
                    other["vbf_eff_selected_pct"] > row["vbf_eff_selected_pct"]
                    or other["nonvbf_eff_selected_pct"] < row["nonvbf_eff_selected_pct"]
                )
            ):
                dominated = True
                break
        if not dominated:
            frontier.append(row)
    return sorted(
        frontier,
        key=lambda row: (
            -row["vbf_eff_selected_pct"],
            row["nonvbf_eff_selected_pct"],
            -row["score"],
        ),
    )


def write_csv(path, rows):
    fieldnames = [
        "working_point",
        "forwardJetPtMin",
        "forwardPairDEtaMin",
        "forwardPairExtraMjjMin",
        "vbf_n_total",
        "vbf_n_after_trigger",
        "vbf_n_step1",
        "vbf_n_has_forward_pair",
        "vbf_n_selected",
        "vbf_eff_selected_pct",
        "nonvbf_n_total",
        "nonvbf_n_after_trigger",
        "nonvbf_n_step1",
        "nonvbf_n_has_forward_pair",
        "nonvbf_n_selected",
        "nonvbf_eff_selected_pct",
        "eff_diff_pct",
        "eff_ratio",
        "score",
        "significance_like",
        "purity_like",
        "score_definition",
        "significance_like_definition",
        "purity_like_definition",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def print_best(rows, key, reverse=True, title="Top working points", nmax=10):
    print(f"\n[INFO] {title}")
    if key == "score":
        print(f"[INFO] score = VBF% - {args.nonVBFPenalty:g} * NoVBF%")
        print("[INFO] Higher is better: it rewards VBF efficiency and penalizes non-VBF contamination.")
    elif key == "significance_like":
        print("[INFO] significance_like = VBF% / sqrt(NoVBF%). Higher values favor rejection more strongly.")
    elif key == "purity_like":
        print("[INFO] purity_like = VBF% / (VBF% + NoVBF%). Higher values mean a cleaner selected sample.")
    print(
        f"{'Rank':>4s} {'fwdPt':>7s} {'dEta':>7s} {'xMjj':>7s} "
        f"{'VBF%':>8s} {'NoVBF%':>8s} {'Diff%':>8s} {'Score':>8s}"
    )
    sorted_rows = sorted(rows, key=lambda row: row[key], reverse=reverse)[:nmax]
    for idx, row in enumerate(sorted_rows, start=1):
        print(
            f"{idx:4d} "
            f"{row['forwardJetPtMin']:7.1f} "
            f"{row['forwardPairDEtaMin']:7.2f} "
            f"{row['forwardPairExtraMjjMin']:7.1f} "
            f"{row['vbf_eff_selected_pct']:8.3f} "
            f"{row['nonvbf_eff_selected_pct']:8.3f} "
            f"{row['eff_diff_pct']:8.3f} "
            f"{row['score']:8.3f}"
        )


def print_pareto(rows, nmax=20):
    print(f"\n[INFO] Pareto-optimal working points (top {min(nmax, len(rows))} shown)")
    print(
        f"{'Rank':>4s} {'fwdPt':>7s} {'dEta':>7s} {'xMjj':>7s} "
        f"{'VBF%':>8s} {'NoVBF%':>8s} {'Diff%':>8s} {'Score':>8s}"
    )
    for idx, row in enumerate(rows[:nmax], start=1):
        print(
            f"{idx:4d} "
            f"{row['forwardJetPtMin']:7.1f} "
            f"{row['forwardPairDEtaMin']:7.2f} "
            f"{row['forwardPairExtraMjjMin']:7.1f} "
            f"{row['vbf_eff_selected_pct']:8.3f} "
            f"{row['nonvbf_eff_selected_pct']:8.3f} "
            f"{row['eff_diff_pct']:8.3f} "
            f"{row['score']:8.3f}"
        )


def main():
    vbf_rows = run_one_file(args.vbfInput, "vbf")
    nonvbf_rows = run_one_file(args.nonVBFInput, "nonvbf")
    rows = combine_rows(vbf_rows, nonvbf_rows)
    rows.sort(key=lambda row: (row["forwardJetPtMin"], row["forwardPairDEtaMin"], row["forwardPairExtraMjjMin"]))
    pareto_rows = pareto_frontier(rows)
    write_csv(args.outputCsv, rows)
    write_csv(args.outputParetoCsv, pareto_rows)

    print(f"[INFO] Wrote joint scan CSV to: {args.outputCsv}")
    print(f"[INFO] Wrote Pareto CSV to:     {args.outputParetoCsv}")
    print(f"[INFO] VBF input:     {args.vbfInput}")
    print(f"[INFO] non-VBF input: {args.nonVBFInput}")
    print(f"[INFO] Number of working points: {len(rows)}")
    print(f"[INFO] Pareto working points:   {len(pareto_rows)}")
    print_best(rows, "score", True, f"Best by score = eff_vbf - {args.nonVBFPenalty:g} * eff_nonvbf")
    print_best(rows, "eff_diff_pct", True, "Best by efficiency difference")
    print_best(rows, "significance_like", True, "Best by significance-like metric")
    print_best(rows, "purity_like", True, "Best by purity-like metric")
    print_best(rows, "nonvbf_eff_selected_pct", False, "Lowest non-VBF selected efficiency")
    print_pareto(pareto_rows)


if __name__ == "__main__":
    main()
