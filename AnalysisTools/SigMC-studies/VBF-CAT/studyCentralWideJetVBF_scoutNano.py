#!/usr/bin/env python3

"""
Reco-level ScoutNano study for a central dijet + wide-jet + VBF topology.

The script:
 - reads ScoutNano signal samples with NanoAODTools
 - selects the leading/subleading central reco jets as the AK4 dijet seeds
 - builds wide jets around those two seeds
 - applies the STEP 1 central selection on the wide-jet pair
 - removes all jets used by the wide-jet building from the reco-jet collection
 - searches the remaining jets for the best VBF-like pair
 - classifies STEP 2 events as isVBF or isPP
 - matches the selected central wide-jet pair to Higgs daughters from GenPart
 - matches the selected forward pair to LHE VBF quarks
 - stores event-level flags in an output tree and summary histograms
"""

import os
import math
import glob
import argparse
from array import array
from types import SimpleNamespace

import ROOT
import correctionlib
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

ROOT.gROOT.SetBatch(True)

BASE_DIR = "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/VBFH-Central/"
#BASE_DIR = "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/VBFHToXX"
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
    description="Central wide-jet + VBF reco study for ScoutNano signal samples"
)
parser.add_argument(
    "--input",
    type=str,
    default=None,
    help="Input ScoutNano ROOT file or glob pattern. Overrides --decay/--mass.",
)
parser.add_argument(
    "--decay",
    type=str,
    default="2B",
    help="Signal decay label used in the directory/file naming, e.g. 2B, 2C, 2Q, 2Glu.",
)
parser.add_argument(
    "--mass",
    type=int,
    default=300,
    help="Signal mass hypothesis used in the file naming.",
)
parser.add_argument(
    "--output",
    type=str,
    default=None,
    help="Output ROOT file. If omitted, a name is built from decay and mass.",
)
parser.add_argument(
    "--maxEntries",
    type=int,
    default=None,
    help="Maximum number of entries to process.",
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
    help="Minimum |DeltaEta| required for the leading central AK4 seed pair. Set to 0 to disable.",
)
parser.add_argument(
    "--centralPairDEtaMax",
    type=float,
    default=0.0,
    help="Maximum |DeltaEta| allowed for the leading central AK4 seed pair. Set to 0 to disable.",
)
parser.add_argument(
    "--centralPairDPhiMin",
    type=float,
    default=0.0,
    help="Minimum |DeltaPhi| required for the leading central AK4 seed pair. Set to 0 to disable.",
)
parser.add_argument(
    "--widejetPairDEtaMin",
    type=float,
    default=0.0,
    help="Minimum |DeltaEta| required for the selected wide-jet pair. Set to 0 to disable.",
)
parser.add_argument(
    "--widejetPairDEtaMax",
    type=float,
    default=1.3,
    help="Maximum |DeltaEta| allowed for the selected wide-jet pair. Set to 0 to disable.",
)
parser.add_argument(
    "--widejetPairDPhiMin",
    type=float,
    default=0.0,
    help="Minimum |DeltaPhi| required for the selected wide-jet pair. Set to 0 to disable.",
)
parser.add_argument(
    "--leadingPtOverMMin",
    type=float,
    default=0.0,
    help="Minimum leading widejet pT / m(widejet pair). Set to 0 to disable.",
)
parser.add_argument(
    "--subleadingPtOverMMin",
    type=float,
    default=0.0,
    help="Minimum subleading widejet pT / m(widejet pair). Set to 0 to disable.",
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
    help="Minimum |eta| for forward VBF-tag jets. Kept for optional bookkeeping.",
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
    help="On the Central-to-VBFTag strategy, require leftover VBF-tag jets to also satisfy |eta| >= forwardEtaMin.",
)
parser.add_argument(
    "--forwardPairAlgo",
    type=str,
    choices=["maxDEta", "maxMjj"],
    default="maxDEta",
    help="How to choose the forward leftover-jet pair: maximize |DeltaEta| or mjj.",
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
    "--requireMinimalVBF",
    action="store_true",
    help="Report the minimal VBF subset explicitly in the summary/output flags.",
)
parser.add_argument(
    "--matchDR",
    type=float,
    default=0.4,
    help="DeltaR threshold used for reco-to-truth matching decisions.",
)
parser.add_argument(
    "--widejetGenMatchDR",
    type=float,
    default=None,
    help="DeltaR threshold for matching Higgs daughters to the two selected wide jets. Defaults to --widejetRadius.",
)
parser.add_argument(
    "--requireTriggerBaseline",
    action="store_true",
    help="Require DST_PFScouting_JetHT and (L1_HTT280er or L1_SingleJet180), vetoing SingleMu SQ14_BMTF bits.",
)
args = parser.parse_args()


JET_MULTIPLICITY_PT_THRESHOLDS = (20, 30)
JET_MULTIPLICITY_ETA_SPLIT = 2.5


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


def p4_from_recojet_raw(jet):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, recojet_mass(jet))
    return vec


def p4_from_lhepart(lp):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(lp.pt, lp.eta, lp.phi, lp.mass)
    return vec


def count_recojets_by_region(jets, pt_threshold, eta_split=JET_MULTIPLICITY_ETA_SPLIT, eta_max=None):
    if eta_max is None:
        eta_max = args.forwardEtaMax

    selected_jets = [
        jet for jet in jets
        if jet_id(jet) and jet.pt > pt_threshold and abs(jet.eta) < eta_max
    ]
    central_jets = [jet for jet in selected_jets if abs(jet.eta) < eta_split]
    forward_jets = [jet for jet in selected_jets if eta_split <= abs(jet.eta) < eta_max]
    return {
        "total": len(selected_jets),
        "central_eta2p5": len(central_jets),
        "forward_eta2p5": len(forward_jets),
    }


def multiplicity_hist_title(branch_name):
    title = branch_name
    title = title.replace("n_recojets_", "")
    title = title.replace("central_eta2p5", "central |eta| < 2.5")
    title = title.replace("forward_eta2p5", "forward 2.5 <= |eta| < 5.0")
    title = title.replace("postwide", "after widejet")
    title = title.replace("raw", "raw")
    title = title.replace("jec", "JEC")
    title = title.replace("total", "inclusive")
    title = title.replace("_pt", ", pT > ")
    title = title.replace("_", " ")
    return f"{title};Jet multiplicity;Events"


def get_input_files():
    if args.input:
        files = sorted(glob.glob(args.input))
        if files:
            return files
        if os.path.exists(args.input):
            return [args.input]
        raise RuntimeError(f"No files found for explicit input pattern: {args.input}")

    decay_dir = f"VBFHTo{args.decay}"
    file_name = f"VBFH-Hto2B_Par-M-{args.mass}.root"
    #file_name = f"hadd_VBFTo{args.decay}_M{args.mass}.root"
    candidate = os.path.join(BASE_DIR, decay_dir, file_name)
    if not os.path.exists(candidate):
        raise RuntimeError(f"Default sample not found: {candidate}")
    return [candidate]


def default_output_name():
    if args.output:
        return args.output
    return (
        f"VBFHTo{args.decay}_M{args.mass}_centralWideJetVBF_scoutNano_"
        f"centPt{int(args.centralJetPtMin)}_fwdPt{int(args.forwardJetPtMin)}.root"
    )


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
        "pass_dst": pass_dst,
        "pass_htt280": pass_htt280,
        "pass_singlejet180": pass_singlejet180,
        "pass_singlemu_veto": pass_singlemu_veto,
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


def final_nonself_daughters(genparts, higgs_idx):
    daughters = []
    queue = [higgs_idx]
    visited = set()

    while queue:
        idx = queue.pop()
        if idx in visited:
            continue
        visited.add(idx)

        parent = genparts[idx]
        child_indices = [
            i for i, gp in enumerate(genparts)
            if getattr(gp, "genPartIdxMother", -1) == idx
        ]

        for child_idx in child_indices:
            child = genparts[child_idx]
            if child.pdgId == parent.pdgId:
                queue.append(child_idx)
            else:
                daughters.append(child)

    return daughters


def find_last_higgs(genparts):
    higgs_indices = [i for i, gp in enumerate(genparts) if abs(gp.pdgId) == 25]
    if not higgs_indices:
        return None, []

    higgs_indices.sort(key=lambda idx: genparts[idx].pt, reverse=True)
    for idx in higgs_indices:
        daughters = final_nonself_daughters(genparts, idx)
        if len(daughters) >= 2:
            return genparts[idx], daughters

    best_idx = higgs_indices[0]
    return genparts[best_idx], final_nonself_daughters(genparts, best_idx)


def is_final_state_lhepart(lp):
    return getattr(lp, "status", -999) == 1


def is_lhe_quark(lp):
    return abs(getattr(lp, "pdgId", 0)) in {1, 2, 3, 4, 5, 6}


def lhe_mothers(lp):
    return (
        getattr(lp, "firstMotherIdx", -1),
        getattr(lp, "lastMotherIdx", -1),
    )


def is_incoming_lhepart(lp):
    return getattr(lp, "status", 0) == -1


def find_higgs_lhe_daughters_from_mothers(lheparts):
    higgs_idx = None
    for idx, lp in enumerate(lheparts):
        if abs(getattr(lp, "pdgId", 0)) == 25:
            higgs_idx = idx
            if getattr(lp, "status", 0) == 2:
                break

    if higgs_idx is None:
        return None, []

    daughter_indices = []
    daughters = []
    for idx, lp in enumerate(lheparts):
        if not is_final_state_lhepart(lp):
            continue
        first_mother, last_mother = lhe_mothers(lp)
        if first_mother <= higgs_idx <= last_mother:
            daughter_indices.append(idx)
            daughters.append(lp)

    return higgs_idx, list(zip(daughter_indices, daughters))


def find_vbf_lhe_quarks_from_mothers(lheparts, excluded_indices=None):
    if excluded_indices is None:
        excluded_indices = set()

    vbf_quarks = []
    for idx, lp in enumerate(lheparts):
        if idx in excluded_indices:
            continue
        if not (is_final_state_lhepart(lp) and is_lhe_quark(lp)):
            continue
        first_mother, last_mother = lhe_mothers(lp)
        if first_mother < 0 and last_mother < 0:
            continue

        mother_indices = {m for m in {first_mother, last_mother} if m >= 0}
        mothers = [lheparts[m] for m in mother_indices if m < len(lheparts)]
        if any(is_incoming_lhepart(mother) for mother in mothers):
            vbf_quarks.append((idx, lp))

    return vbf_quarks


def choose_vbf_lhe_pair(lheparts, excluded_indices=None):
    if excluded_indices is None:
        excluded_indices = set()

    candidate_quarks = [
        (idx, lp)
        for idx, lp in enumerate(lheparts)
        if idx not in excluded_indices and is_final_state_lhepart(lp) and is_lhe_quark(lp)
    ]
    if len(candidate_quarks) < 2:
        return None

    best_pair = None
    best_deta = -1.0
    for i in range(len(candidate_quarks)):
        for j in range(i + 1, len(candidate_quarks)):
            idx1, lp1 = candidate_quarks[i]
            idx2, lp2 = candidate_quarks[j]
            deta = abs(lp1.eta - lp2.eta)
            dphi = abs(delta_phi(lp1.phi, lp2.phi))
            dr = delta_r(lp1.eta, lp1.phi, lp2.eta, lp2.phi)
            mjj = (p4_from_lhepart(lp1) + p4_from_lhepart(lp2)).M()
            if deta > best_deta:
                best_deta = deta
                best_pair = (idx1, lp1, idx2, lp2, deta, dphi, dr, mjj)
    return best_pair


def match_pair_by_delta_r(obj1, obj2, truth1, truth2):
    assignment_a = (
        delta_r(obj1.eta, obj1.phi, truth1.eta, truth1.phi),
        delta_r(obj2.eta, obj2.phi, truth2.eta, truth2.phi),
    )
    assignment_b = (
        delta_r(obj1.eta, obj1.phi, truth2.eta, truth2.phi),
        delta_r(obj2.eta, obj2.phi, truth1.eta, truth1.phi),
    )

    if max(assignment_a) <= max(assignment_b):
        return assignment_a
    return assignment_b


def build_widejets_from_seeds(jets, seed1, seed2, widejet_radius):
    seed1_p4 = p4_from_recojet(seed1)
    seed2_p4 = p4_from_recojet(seed2)
    widejet1 = ROOT.TLorentzVector(seed1_p4)
    widejet2 = ROOT.TLorentzVector(seed2_p4)

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


def choose_forward_pair(jets):
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


def choose_max_deta_pair(jets):
    return choose_forward_pair(jets)


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


def pass_vbf_requirement(pair):
    if pair is None:
        return False

    _, _, deta, _, _, mjj = pair
    pass_pair_algo = False
    if args.forwardPairAlgo == "maxDEta":
        pass_pair_algo = deta >= args.forwardPairDEtaMin
    elif args.forwardPairAlgo == "maxMjj":
        pass_pair_algo = mjj >= args.forwardPairMjjMin

    if not pass_pair_algo:
        return False
    if args.forwardPairExtraMjjMin > 0.0 and mjj < args.forwardPairExtraMjjMin:
        return False
    return True


def match_widejet_pair_to_higgs_daughters(widejet1, widejet2, daughters, max_dr):
    if len(daughters) < 2:
        return None

    daughters_sorted = sorted(daughters, key=lambda x: x.pt, reverse=True)[:2]
    widejet_obj1 = SimpleNamespace(eta=widejet1.Eta(), phi=widejet1.Phi())
    widejet_obj2 = SimpleNamespace(eta=widejet2.Eta(), phi=widejet2.Phi())
    dr1, dr2 = match_pair_by_delta_r(
        widejet_obj1,
        widejet_obj2,
        daughters_sorted[0],
        daughters_sorted[1],
    )
    is_matched = max(dr1, dr2) < max_dr
    return dr1, dr2, max(dr1, dr2), is_matched


class CentralWideJetVBFScoutNanoStudy(Module):
    def __init__(self, output_name):
        self.output_name = output_name
        self.writeHistFile = False
        self.total = 0
        self.after_trigger = 0
        self.n_with_central_pair = 0
        self.n_step1_pre_deta = 0
        self.n_step1 = 0
        self.n_with_central_pair_genmatch = 0
        self.n_with_vbf_pair = 0
        self.n_with_vbf_tag = 0
        self.n_with_pp_tag = 0
        self.n_with_vbf_lhe_match = 0

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        jec_json = os.path.join(
            os.getenv("CMSSW_BASE", ""),
            "src",
            "2024_UtilsDataQuality",
            "jetHLT_jerc.json"
        )
        if not os.path.exists(jec_json):
            raise RuntimeError(f"HLT JEC JSON not found: {jec_json}")

        cset = correctionlib.CorrectionSet.from_file(jec_json)
        self.jec = cset.compound["HLT_Winter24_V1_MC_L1L2L3Res_AK4PFHLT"]

        self.outfile = ROOT.TFile(self.output_name, "RECREATE")
        self.tree = ROOT.TTree("Events", "Central wide-jet + VBF reco study tree from ScoutNano")
        self.mult_dir = self.outfile.mkdir("Multiplicity")

        self.run = array("I", [0])
        self.lumi = array("I", [0])
        self.event = array("L", [0])

        self.pass_dst_jetht = array("i", [0])
        self.pass_l1_htt280 = array("i", [0])
        self.pass_l1_singlejet180 = array("i", [0])
        self.pass_trigger_baseline = array("i", [0])

        self.n_recojets_total = array("i", [0])
        self.n_recojets_central = array("i", [0])
        self.n_recojets_remaining = array("i", [0])
        self.n_recojets_forward_remaining = array("i", [0])
        self.n_lhe_vbf_quarks = array("i", [0])
        self.recojet_mult_branches = {}
        self.recojet_mult_hists = {}
        for level in ["raw", "jec", "postwide"]:
            for region in ["total", "central_eta2p5", "forward_eta2p5"]:
                for pt_threshold in JET_MULTIPLICITY_PT_THRESHOLDS:
                    branch_name = f"n_recojets_{level}_{region}_pt{pt_threshold}"
                    self.recojet_mult_branches[branch_name] = array("i", [0])
                    self.recojet_mult_hists[branch_name] = ROOT.TH1F(
                        f"h_{branch_name}",
                        multiplicity_hist_title(branch_name),
                        31,
                        -0.5,
                        30.5,
                    )

        self.has_gen_higgs = array("i", [0])
        self.has_two_higgs_daughters = array("i", [0])
        self.has_lhe_vbf_pair = array("i", [0])
        self.used_lhe_vbf_topology_fallback = array("i", [0])

        self.has_central_pair = array("i", [0])
        self.pass_step1_kinematics = array("i", [0])
        self.pass_central_pair_deta = array("i", [0])
        self.pass_preclean_selection = array("i", [0])
        self.pass_central_pair_selection = array("i", [0])
        self.central_pair_genmatched = array("i", [0])
        self.has_widejet_pair = array("i", [0])
        self.has_vbf_pair = array("i", [0])
        self.pass_minimal_vbf = array("i", [0])
        self.pass_final_selection = array("i", [0])
        self.isVBF = array("i", [0])
        self.isPP = array("i", [0])
        self.vbf_pair_lhematched = array("i", [0])

        self.inv_n_recojets_remaining = array("i", [0])
        self.inv_n_recojets_central = array("i", [0])
        self.inv_has_vbf_pair = array("i", [0])
        self.inv_pass_minimal_vbf = array("i", [0])
        self.inv_has_central_pair = array("i", [0])
        self.inv_pass_central_pair_selection = array("i", [0])
        self.inv_has_widejet_pair = array("i", [0])
        self.inv_pass_final_selection = array("i", [0])
        self.inv_isVBF = array("i", [0])
        self.inv_central_pair_genmatched = array("i", [0])
        self.inv_vbf_pair_lhematched = array("i", [0])

        self.central_jet1_pt = array("f", [0.0])
        self.central_jet1_eta = array("f", [0.0])
        self.central_jet1_phi = array("f", [0.0])
        self.central_jet1_mass = array("f", [0.0])
        self.central_jet1_pt_raw = array("f", [0.0])
        self.central_jet2_pt = array("f", [0.0])
        self.central_jet2_eta = array("f", [0.0])
        self.central_jet2_phi = array("f", [0.0])
        self.central_jet2_mass = array("f", [0.0])
        self.central_jet2_pt_raw = array("f", [0.0])
        self.central_pair_dEta = array("f", [0.0])
        self.central_pair_dPhi = array("f", [0.0])
        self.central_pair_dR = array("f", [0.0])
        self.central_pair_mass = array("f", [0.0])
        self.central_pair_mass_raw = array("f", [0.0])
        self.central_pair_match_dR1 = array("f", [0.0])
        self.central_pair_match_dR2 = array("f", [0.0])
        self.central_pair_match_maxdr = array("f", [0.0])

        self.widejet1_pt = array("f", [0.0])
        self.widejet1_eta = array("f", [0.0])
        self.widejet1_phi = array("f", [0.0])
        self.widejet1_mass = array("f", [0.0])
        self.widejet2_pt = array("f", [0.0])
        self.widejet2_eta = array("f", [0.0])
        self.widejet2_phi = array("f", [0.0])
        self.widejet2_mass = array("f", [0.0])
        self.widejet_pair_dEta = array("f", [0.0])
        self.widejet_pair_dPhi = array("f", [0.0])
        self.widejet_pair_dR = array("f", [0.0])
        self.widejet_pair_mass = array("f", [0.0])
        self.widejet_leading_pt_over_m = array("f", [0.0])
        self.widejet_subleading_pt_over_m = array("f", [0.0])
        self.n_widejet1_constituents = array("i", [0])
        self.n_widejet2_constituents = array("i", [0])

        self.forward_jet1_pt = array("f", [0.0])
        self.forward_jet1_eta = array("f", [0.0])
        self.forward_jet1_phi = array("f", [0.0])
        self.forward_jet1_mass = array("f", [0.0])
        self.forward_jet2_pt = array("f", [0.0])
        self.forward_jet2_eta = array("f", [0.0])
        self.forward_jet2_phi = array("f", [0.0])
        self.forward_jet2_mass = array("f", [0.0])
        self.forward_pair_dEta = array("f", [0.0])
        self.forward_pair_dPhi = array("f", [0.0])
        self.forward_pair_dR = array("f", [0.0])
        self.forward_pair_mass = array("f", [0.0])
        self.vbf_pair_match_dR1 = array("f", [0.0])
        self.vbf_pair_match_dR2 = array("f", [0.0])
        self.vbf_pair_match_maxdr = array("f", [0.0])

        self.inv_central_pair_dEta = array("f", [0.0])
        self.inv_central_pair_dPhi = array("f", [0.0])
        self.inv_central_pair_dR = array("f", [0.0])
        self.inv_central_pair_mass = array("f", [0.0])
        self.inv_central_pair_mass_raw = array("f", [0.0])
        self.inv_central_pair_match_dR1 = array("f", [0.0])
        self.inv_central_pair_match_dR2 = array("f", [0.0])
        self.inv_central_pair_match_maxdr = array("f", [0.0])
        self.inv_widejet_pair_dEta = array("f", [0.0])
        self.inv_widejet_pair_dPhi = array("f", [0.0])
        self.inv_widejet_pair_dR = array("f", [0.0])
        self.inv_widejet_pair_mass = array("f", [0.0])
        self.inv_forward_pair_dEta = array("f", [0.0])
        self.inv_forward_pair_dPhi = array("f", [0.0])
        self.inv_forward_pair_dR = array("f", [0.0])
        self.inv_forward_pair_mass = array("f", [0.0])
        self.inv_vbf_pair_match_dR1 = array("f", [0.0])
        self.inv_vbf_pair_match_dR2 = array("f", [0.0])
        self.inv_vbf_pair_match_maxdr = array("f", [0.0])

        self.lhe_vbf1_pdgId = array("i", [0])
        self.lhe_vbf1_pt = array("f", [0.0])
        self.lhe_vbf1_eta = array("f", [0.0])
        self.lhe_vbf1_phi = array("f", [0.0])
        self.lhe_vbf1_mass = array("f", [0.0])
        self.lhe_vbf2_pdgId = array("i", [0])
        self.lhe_vbf2_pt = array("f", [0.0])
        self.lhe_vbf2_eta = array("f", [0.0])
        self.lhe_vbf2_phi = array("f", [0.0])
        self.lhe_vbf2_mass = array("f", [0.0])
        self.lhe_vbf_pair_dEta = array("f", [0.0])
        self.lhe_vbf_pair_dPhi = array("f", [0.0])
        self.lhe_vbf_pair_dR = array("f", [0.0])
        self.lhe_vbf_pair_mass = array("f", [0.0])

        self.tree.Branch("run", self.run, "run/i")
        self.tree.Branch("lumi", self.lumi, "lumi/i")
        self.tree.Branch("event", self.event, "event/l")
        self.tree.Branch("pass_dst_jetht", self.pass_dst_jetht, "pass_dst_jetht/I")
        self.tree.Branch("pass_l1_htt280", self.pass_l1_htt280, "pass_l1_htt280/I")
        self.tree.Branch("pass_l1_singlejet180", self.pass_l1_singlejet180, "pass_l1_singlejet180/I")
        self.tree.Branch("pass_trigger_baseline", self.pass_trigger_baseline, "pass_trigger_baseline/I")

        self.tree.Branch("n_recojets_total", self.n_recojets_total, "n_recojets_total/I")
        self.tree.Branch("n_recojets_central", self.n_recojets_central, "n_recojets_central/I")
        self.tree.Branch("n_recojets_remaining", self.n_recojets_remaining, "n_recojets_remaining/I")
        self.tree.Branch("n_recojets_forward_remaining", self.n_recojets_forward_remaining, "n_recojets_forward_remaining/I")
        self.tree.Branch("n_lhe_vbf_quarks", self.n_lhe_vbf_quarks, "n_lhe_vbf_quarks/I")
        for branch_name, branch_array in self.recojet_mult_branches.items():
            self.tree.Branch(branch_name, branch_array, f"{branch_name}/I")

        self.tree.Branch("has_gen_higgs", self.has_gen_higgs, "has_gen_higgs/I")
        self.tree.Branch("has_two_higgs_daughters", self.has_two_higgs_daughters, "has_two_higgs_daughters/I")
        self.tree.Branch("has_lhe_vbf_pair", self.has_lhe_vbf_pair, "has_lhe_vbf_pair/I")
        self.tree.Branch(
            "used_lhe_vbf_topology_fallback",
            self.used_lhe_vbf_topology_fallback,
            "used_lhe_vbf_topology_fallback/I",
        )

        self.tree.Branch("has_central_pair", self.has_central_pair, "has_central_pair/I")
        self.tree.Branch("pass_step1_kinematics", self.pass_step1_kinematics, "pass_step1_kinematics/I")
        self.tree.Branch("pass_central_pair_deta", self.pass_central_pair_deta, "pass_central_pair_deta/I")
        self.tree.Branch("pass_preclean_selection", self.pass_preclean_selection, "pass_preclean_selection/I")
        self.tree.Branch("pass_central_pair_selection", self.pass_central_pair_selection, "pass_central_pair_selection/I")
        self.tree.Branch("central_pair_genmatched", self.central_pair_genmatched, "central_pair_genmatched/I")
        self.tree.Branch("has_widejet_pair", self.has_widejet_pair, "has_widejet_pair/I")
        self.tree.Branch("has_vbf_pair", self.has_vbf_pair, "has_vbf_pair/I")
        self.tree.Branch("pass_minimal_vbf", self.pass_minimal_vbf, "pass_minimal_vbf/I")
        self.tree.Branch("pass_final_selection", self.pass_final_selection, "pass_final_selection/I")
        self.tree.Branch("isVBF", self.isVBF, "isVBF/I")
        self.tree.Branch("isPP", self.isPP, "isPP/I")
        self.tree.Branch("vbf_pair_lhematched", self.vbf_pair_lhematched, "vbf_pair_lhematched/I")
        self.tree.Branch("inv_n_recojets_remaining", self.inv_n_recojets_remaining, "inv_n_recojets_remaining/I")
        self.tree.Branch("inv_n_recojets_central", self.inv_n_recojets_central, "inv_n_recojets_central/I")
        self.tree.Branch("inv_has_vbf_pair", self.inv_has_vbf_pair, "inv_has_vbf_pair/I")
        self.tree.Branch("inv_pass_minimal_vbf", self.inv_pass_minimal_vbf, "inv_pass_minimal_vbf/I")
        self.tree.Branch("inv_has_central_pair", self.inv_has_central_pair, "inv_has_central_pair/I")
        self.tree.Branch(
            "inv_pass_central_pair_selection",
            self.inv_pass_central_pair_selection,
            "inv_pass_central_pair_selection/I",
        )
        self.tree.Branch("inv_has_widejet_pair", self.inv_has_widejet_pair, "inv_has_widejet_pair/I")
        self.tree.Branch("inv_pass_final_selection", self.inv_pass_final_selection, "inv_pass_final_selection/I")
        self.tree.Branch("inv_isVBF", self.inv_isVBF, "inv_isVBF/I")
        self.tree.Branch("inv_central_pair_genmatched", self.inv_central_pair_genmatched, "inv_central_pair_genmatched/I")
        self.tree.Branch("inv_vbf_pair_lhematched", self.inv_vbf_pair_lhematched, "inv_vbf_pair_lhematched/I")

        self.tree.Branch("central_jet1_pt", self.central_jet1_pt, "central_jet1_pt/F")
        self.tree.Branch("central_jet1_eta", self.central_jet1_eta, "central_jet1_eta/F")
        self.tree.Branch("central_jet1_phi", self.central_jet1_phi, "central_jet1_phi/F")
        self.tree.Branch("central_jet1_mass", self.central_jet1_mass, "central_jet1_mass/F")
        self.tree.Branch("central_jet1_pt_raw", self.central_jet1_pt_raw, "central_jet1_pt_raw/F")
        self.tree.Branch("central_jet2_pt", self.central_jet2_pt, "central_jet2_pt/F")
        self.tree.Branch("central_jet2_eta", self.central_jet2_eta, "central_jet2_eta/F")
        self.tree.Branch("central_jet2_phi", self.central_jet2_phi, "central_jet2_phi/F")
        self.tree.Branch("central_jet2_mass", self.central_jet2_mass, "central_jet2_mass/F")
        self.tree.Branch("central_jet2_pt_raw", self.central_jet2_pt_raw, "central_jet2_pt_raw/F")
        self.tree.Branch("central_pair_dEta", self.central_pair_dEta, "central_pair_dEta/F")
        self.tree.Branch("central_pair_dPhi", self.central_pair_dPhi, "central_pair_dPhi/F")
        self.tree.Branch("central_pair_dR", self.central_pair_dR, "central_pair_dR/F")
        self.tree.Branch("central_pair_mass", self.central_pair_mass, "central_pair_mass/F")
        self.tree.Branch("central_pair_mass_raw", self.central_pair_mass_raw, "central_pair_mass_raw/F")
        self.tree.Branch("central_pair_match_dR1", self.central_pair_match_dR1, "central_pair_match_dR1/F")
        self.tree.Branch("central_pair_match_dR2", self.central_pair_match_dR2, "central_pair_match_dR2/F")
        self.tree.Branch("central_pair_match_maxdr", self.central_pair_match_maxdr, "central_pair_match_maxdr/F")

        self.tree.Branch("widejet1_pt", self.widejet1_pt, "widejet1_pt/F")
        self.tree.Branch("widejet1_eta", self.widejet1_eta, "widejet1_eta/F")
        self.tree.Branch("widejet1_phi", self.widejet1_phi, "widejet1_phi/F")
        self.tree.Branch("widejet1_mass", self.widejet1_mass, "widejet1_mass/F")
        self.tree.Branch("widejet2_pt", self.widejet2_pt, "widejet2_pt/F")
        self.tree.Branch("widejet2_eta", self.widejet2_eta, "widejet2_eta/F")
        self.tree.Branch("widejet2_phi", self.widejet2_phi, "widejet2_phi/F")
        self.tree.Branch("widejet2_mass", self.widejet2_mass, "widejet2_mass/F")
        self.tree.Branch("widejet_pair_dEta", self.widejet_pair_dEta, "widejet_pair_dEta/F")
        self.tree.Branch("widejet_pair_dPhi", self.widejet_pair_dPhi, "widejet_pair_dPhi/F")
        self.tree.Branch("widejet_pair_dR", self.widejet_pair_dR, "widejet_pair_dR/F")
        self.tree.Branch("widejet_pair_mass", self.widejet_pair_mass, "widejet_pair_mass/F")
        self.tree.Branch("widejet_leading_pt_over_m", self.widejet_leading_pt_over_m, "widejet_leading_pt_over_m/F")
        self.tree.Branch("widejet_subleading_pt_over_m", self.widejet_subleading_pt_over_m, "widejet_subleading_pt_over_m/F")
        self.tree.Branch("n_widejet1_constituents", self.n_widejet1_constituents, "n_widejet1_constituents/I")
        self.tree.Branch("n_widejet2_constituents", self.n_widejet2_constituents, "n_widejet2_constituents/I")

        self.tree.Branch("forward_jet1_pt", self.forward_jet1_pt, "forward_jet1_pt/F")
        self.tree.Branch("forward_jet1_eta", self.forward_jet1_eta, "forward_jet1_eta/F")
        self.tree.Branch("forward_jet1_phi", self.forward_jet1_phi, "forward_jet1_phi/F")
        self.tree.Branch("forward_jet1_mass", self.forward_jet1_mass, "forward_jet1_mass/F")
        self.tree.Branch("forward_jet2_pt", self.forward_jet2_pt, "forward_jet2_pt/F")
        self.tree.Branch("forward_jet2_eta", self.forward_jet2_eta, "forward_jet2_eta/F")
        self.tree.Branch("forward_jet2_phi", self.forward_jet2_phi, "forward_jet2_phi/F")
        self.tree.Branch("forward_jet2_mass", self.forward_jet2_mass, "forward_jet2_mass/F")
        self.tree.Branch("forward_pair_dEta", self.forward_pair_dEta, "forward_pair_dEta/F")
        self.tree.Branch("forward_pair_dPhi", self.forward_pair_dPhi, "forward_pair_dPhi/F")
        self.tree.Branch("forward_pair_dR", self.forward_pair_dR, "forward_pair_dR/F")
        self.tree.Branch("forward_pair_mass", self.forward_pair_mass, "forward_pair_mass/F")
        self.tree.Branch("vbf_pair_match_dR1", self.vbf_pair_match_dR1, "vbf_pair_match_dR1/F")
        self.tree.Branch("vbf_pair_match_dR2", self.vbf_pair_match_dR2, "vbf_pair_match_dR2/F")
        self.tree.Branch("vbf_pair_match_maxdr", self.vbf_pair_match_maxdr, "vbf_pair_match_maxdr/F")
        self.tree.Branch("inv_central_pair_dEta", self.inv_central_pair_dEta, "inv_central_pair_dEta/F")
        self.tree.Branch("inv_central_pair_dPhi", self.inv_central_pair_dPhi, "inv_central_pair_dPhi/F")
        self.tree.Branch("inv_central_pair_dR", self.inv_central_pair_dR, "inv_central_pair_dR/F")
        self.tree.Branch("inv_central_pair_mass", self.inv_central_pair_mass, "inv_central_pair_mass/F")
        self.tree.Branch("inv_central_pair_mass_raw", self.inv_central_pair_mass_raw, "inv_central_pair_mass_raw/F")
        self.tree.Branch("inv_central_pair_match_dR1", self.inv_central_pair_match_dR1, "inv_central_pair_match_dR1/F")
        self.tree.Branch("inv_central_pair_match_dR2", self.inv_central_pair_match_dR2, "inv_central_pair_match_dR2/F")
        self.tree.Branch("inv_central_pair_match_maxdr", self.inv_central_pair_match_maxdr, "inv_central_pair_match_maxdr/F")
        self.tree.Branch("inv_widejet_pair_dEta", self.inv_widejet_pair_dEta, "inv_widejet_pair_dEta/F")
        self.tree.Branch("inv_widejet_pair_dPhi", self.inv_widejet_pair_dPhi, "inv_widejet_pair_dPhi/F")
        self.tree.Branch("inv_widejet_pair_dR", self.inv_widejet_pair_dR, "inv_widejet_pair_dR/F")
        self.tree.Branch("inv_widejet_pair_mass", self.inv_widejet_pair_mass, "inv_widejet_pair_mass/F")
        self.tree.Branch("inv_forward_pair_dEta", self.inv_forward_pair_dEta, "inv_forward_pair_dEta/F")
        self.tree.Branch("inv_forward_pair_dPhi", self.inv_forward_pair_dPhi, "inv_forward_pair_dPhi/F")
        self.tree.Branch("inv_forward_pair_dR", self.inv_forward_pair_dR, "inv_forward_pair_dR/F")
        self.tree.Branch("inv_forward_pair_mass", self.inv_forward_pair_mass, "inv_forward_pair_mass/F")
        self.tree.Branch("inv_vbf_pair_match_dR1", self.inv_vbf_pair_match_dR1, "inv_vbf_pair_match_dR1/F")
        self.tree.Branch("inv_vbf_pair_match_dR2", self.inv_vbf_pair_match_dR2, "inv_vbf_pair_match_dR2/F")
        self.tree.Branch("inv_vbf_pair_match_maxdr", self.inv_vbf_pair_match_maxdr, "inv_vbf_pair_match_maxdr/F")

        self.tree.Branch("lhe_vbf1_pdgId", self.lhe_vbf1_pdgId, "lhe_vbf1_pdgId/I")
        self.tree.Branch("lhe_vbf1_pt", self.lhe_vbf1_pt, "lhe_vbf1_pt/F")
        self.tree.Branch("lhe_vbf1_eta", self.lhe_vbf1_eta, "lhe_vbf1_eta/F")
        self.tree.Branch("lhe_vbf1_phi", self.lhe_vbf1_phi, "lhe_vbf1_phi/F")
        self.tree.Branch("lhe_vbf1_mass", self.lhe_vbf1_mass, "lhe_vbf1_mass/F")
        self.tree.Branch("lhe_vbf2_pdgId", self.lhe_vbf2_pdgId, "lhe_vbf2_pdgId/I")
        self.tree.Branch("lhe_vbf2_pt", self.lhe_vbf2_pt, "lhe_vbf2_pt/F")
        self.tree.Branch("lhe_vbf2_eta", self.lhe_vbf2_eta, "lhe_vbf2_eta/F")
        self.tree.Branch("lhe_vbf2_phi", self.lhe_vbf2_phi, "lhe_vbf2_phi/F")
        self.tree.Branch("lhe_vbf2_mass", self.lhe_vbf2_mass, "lhe_vbf2_mass/F")
        self.tree.Branch("lhe_vbf_pair_dEta", self.lhe_vbf_pair_dEta, "lhe_vbf_pair_dEta/F")
        self.tree.Branch("lhe_vbf_pair_dPhi", self.lhe_vbf_pair_dPhi, "lhe_vbf_pair_dPhi/F")
        self.tree.Branch("lhe_vbf_pair_dR", self.lhe_vbf_pair_dR, "lhe_vbf_pair_dR/F")
        self.tree.Branch("lhe_vbf_pair_mass", self.lhe_vbf_pair_mass, "lhe_vbf_pair_mass/F")

        self.cutflow = ROOT.TH1F("cutflow", "Selection cutflow;Selection step;Events", 6, 0.5, 6.5)
        self.cutflow.GetXaxis().SetBinLabel(1, "All")
        self.cutflow.GetXaxis().SetBinLabel(2, "Trigger")
        self.cutflow.GetXaxis().SetBinLabel(3, "STEP1 pre-#Delta#eta")
        self.cutflow.GetXaxis().SetBinLabel(4, "STEP1")
        self.cutflow.GetXaxis().SetBinLabel(5, "VBF cand")
        self.cutflow.GetXaxis().SetBinLabel(6, "isVBF")

        self.h_central_mass_selected = ROOT.TH1F(
            "h_central_mass_selected",
            "Central wide-jet pair mass after STEP 1;m_{JJ}^{wide} [GeV];Events",
            240,
            0.0,
            2400.0,
        )
        self.h_central_deta_selected = ROOT.TH1F(
            "h_central_deta_selected",
            "Central wide-jet pair |#Delta#eta| after STEP 1;|#Delta#eta(J_{1}^{wide},J_{2}^{wide})|;Events",
            60,
            0.0,
            6.0,
        )
        self.h_forward_mass_isvbf = ROOT.TH1F(
            "h_forward_mass_isvbf",
            "Forward-pair mass for isVBF;m_{jj}^{VBF} [GeV];Events",
            300,
            0.0,
            6000.0,
        )
        self.h_forward_deta_isvbf = ROOT.TH1F(
            "h_forward_deta_isvbf",
            "Forward-pair |#Delta#eta| for isVBF;|#Delta#eta(j_{1}^{VBF},j_{2}^{VBF})|;Events",
            80,
            0.0,
            8.0,
        )
        self.h_forward_mass_ispp = ROOT.TH1F(
            "h_forward_mass_ispp",
            "Forward-pair mass for isPP;m_{jj}^{VBF cand} [GeV];Events",
            300,
            0.0,
            6000.0,
        )
        self.h_forward_deta_ispp = ROOT.TH1F(
            "h_forward_deta_ispp",
            "Forward-pair |#Delta#eta| for isPP;|#Delta#eta(j_{1}^{VBF cand},j_{2}^{VBF cand})|;Events",
            80,
            0.0,
            8.0,
        )
        self.h_central_mass_matched = ROOT.TH1F(
            "h_central_mass_matched",
            "Central wide-jet pair mass, matched;m_{JJ}^{wide} [GeV];Events",
            240,
            0.0,
            2400.0,
        )
        self.h_central_mass_unmatched = ROOT.TH1F(
            "h_central_mass_unmatched",
            "Central wide-jet pair mass, not matched;m_{JJ}^{wide} [GeV];Events",
            240,
            0.0,
            2400.0,
        )

    def reset(self):
        int_branches = [
            self.pass_dst_jetht,
            self.pass_l1_htt280,
            self.pass_l1_singlejet180,
            self.pass_trigger_baseline,
            self.n_recojets_total,
            self.n_recojets_central,
            self.n_recojets_remaining,
            self.n_recojets_forward_remaining,
            self.n_lhe_vbf_quarks,
            *self.recojet_mult_branches.values(),
            self.has_gen_higgs,
            self.has_two_higgs_daughters,
            self.has_lhe_vbf_pair,
            self.used_lhe_vbf_topology_fallback,
            self.has_central_pair,
            self.pass_step1_kinematics,
            self.pass_central_pair_deta,
            self.pass_preclean_selection,
            self.pass_central_pair_selection,
            self.central_pair_genmatched,
            self.has_widejet_pair,
            self.has_vbf_pair,
            self.pass_minimal_vbf,
            self.pass_final_selection,
            self.isVBF,
            self.isPP,
            self.vbf_pair_lhematched,
            self.inv_n_recojets_remaining,
            self.inv_n_recojets_central,
            self.inv_has_vbf_pair,
            self.inv_pass_minimal_vbf,
            self.inv_has_central_pair,
            self.inv_pass_central_pair_selection,
            self.inv_has_widejet_pair,
            self.inv_pass_final_selection,
            self.inv_isVBF,
            self.inv_central_pair_genmatched,
            self.inv_vbf_pair_lhematched,
            self.n_widejet1_constituents,
            self.n_widejet2_constituents,
            self.lhe_vbf1_pdgId,
            self.lhe_vbf2_pdgId,
        ]
        for branch in int_branches:
            branch[0] = 0

        float_branches = [
            self.central_jet1_pt,
            self.central_jet1_eta,
            self.central_jet1_phi,
            self.central_jet1_mass,
            self.central_jet1_pt_raw,
            self.central_jet2_pt,
            self.central_jet2_eta,
            self.central_jet2_phi,
            self.central_jet2_mass,
            self.central_jet2_pt_raw,
            self.central_pair_dEta,
            self.central_pair_dPhi,
            self.central_pair_dR,
            self.central_pair_mass,
            self.central_pair_mass_raw,
            self.central_pair_match_dR1,
            self.central_pair_match_dR2,
            self.central_pair_match_maxdr,
            self.widejet1_pt,
            self.widejet1_eta,
            self.widejet1_phi,
            self.widejet1_mass,
            self.widejet2_pt,
            self.widejet2_eta,
            self.widejet2_phi,
            self.widejet2_mass,
            self.widejet_pair_dEta,
            self.widejet_pair_dPhi,
            self.widejet_pair_dR,
            self.widejet_pair_mass,
            self.widejet_leading_pt_over_m,
            self.widejet_subleading_pt_over_m,
            self.forward_jet1_pt,
            self.forward_jet1_eta,
            self.forward_jet1_phi,
            self.forward_jet1_mass,
            self.forward_jet2_pt,
            self.forward_jet2_eta,
            self.forward_jet2_phi,
            self.forward_jet2_mass,
            self.forward_pair_dEta,
            self.forward_pair_dPhi,
            self.forward_pair_dR,
            self.forward_pair_mass,
            self.vbf_pair_match_dR1,
            self.vbf_pair_match_dR2,
            self.vbf_pair_match_maxdr,
            self.inv_central_pair_dEta,
            self.inv_central_pair_dPhi,
            self.inv_central_pair_dR,
            self.inv_central_pair_mass,
            self.inv_central_pair_mass_raw,
            self.inv_central_pair_match_dR1,
            self.inv_central_pair_match_dR2,
            self.inv_central_pair_match_maxdr,
            self.inv_widejet_pair_dEta,
            self.inv_widejet_pair_dPhi,
            self.inv_widejet_pair_dR,
            self.inv_widejet_pair_mass,
            self.inv_forward_pair_dEta,
            self.inv_forward_pair_dPhi,
            self.inv_forward_pair_dR,
            self.inv_forward_pair_mass,
            self.inv_vbf_pair_match_dR1,
            self.inv_vbf_pair_match_dR2,
            self.inv_vbf_pair_match_maxdr,
            self.lhe_vbf1_pt,
            self.lhe_vbf1_eta,
            self.lhe_vbf1_phi,
            self.lhe_vbf1_mass,
            self.lhe_vbf2_pt,
            self.lhe_vbf2_eta,
            self.lhe_vbf2_phi,
            self.lhe_vbf2_mass,
            self.lhe_vbf_pair_dEta,
            self.lhe_vbf_pair_dPhi,
            self.lhe_vbf_pair_dR,
            self.lhe_vbf_pair_mass,
        ]
        for branch in float_branches:
            branch[0] = -999.0

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

    def fill_multiplicity_hists(self):
        for branch_name, branch_array in self.recojet_mult_branches.items():
            self.recojet_mult_hists[branch_name].Fill(branch_array[0])

    def analyze(self, event):
        self.total += 1

        self.run[0] = getattr(event, "run", 0)
        self.lumi[0] = getattr(event, "luminosityBlock", 0)
        self.event[0] = getattr(event, "event", 0)
        self.reset()

        trigger_bits = get_trigger_bits(event)
        self.pass_dst_jetht[0] = int(trigger_bits["pass_dst"])
        self.pass_l1_htt280[0] = int(trigger_bits["pass_htt280"])
        self.pass_l1_singlejet180[0] = int(trigger_bits["pass_singlejet180"])
        self.pass_trigger_baseline[0] = int(trigger_bits["pass_baseline"])

        if args.requireTriggerBaseline and not trigger_bits["pass_baseline"]:
            self.fill_multiplicity_hists()
            self.tree.Fill()
            return True

        if trigger_bits["pass_baseline"]:
            self.after_trigger += 1

        genparts = Collection(event, "GenPart")
        lheparts = Collection(event, "LHEPart") if hasattr(event, "nLHEPart") else []
        jets_raw = Collection(event, "ScoutingPFJet")
        jets = self.corrected_jets(event, jets_raw)

        for pt_threshold in JET_MULTIPLICITY_PT_THRESHOLDS:
            raw_counts = count_recojets_by_region(jets_raw, pt_threshold)
            jec_counts = count_recojets_by_region(jets, pt_threshold)
            for region, count in raw_counts.items():
                self.recojet_mult_branches[f"n_recojets_raw_{region}_pt{pt_threshold}"][0] = count
            for region, count in jec_counts.items():
                self.recojet_mult_branches[f"n_recojets_jec_{region}_pt{pt_threshold}"][0] = count

        higgs, daughters = find_last_higgs(genparts)
        if higgs is not None:
            self.has_gen_higgs[0] = 1
        if len(daughters) >= 2:
            self.has_two_higgs_daughters[0] = 1

        good_jets = [
            jet for jet in jets
            if jet_id(jet) and jet.pt > min(args.centralJetPtMin, args.forwardJetPtMin) and abs(jet.eta) < args.forwardEtaMax
        ]
        self.n_recojets_total[0] = len(good_jets)

        central_jets = [
            jet for jet in good_jets
            if jet.pt > args.centralJetPtMin and abs(jet.eta) < args.centralEtaMax
        ]
        central_jets = sorted(central_jets, key=lambda x: x.pt, reverse=True)
        self.n_recojets_central[0] = len(central_jets)

        lhe_higgs_idx, lhe_higgs_daughters = find_higgs_lhe_daughters_from_mothers(lheparts)
        excluded_lhe_indices = {idx for idx, _ in lhe_higgs_daughters}
        mother_based_vbf_quarks = find_vbf_lhe_quarks_from_mothers(lheparts, excluded_lhe_indices)
        self.n_lhe_vbf_quarks[0] = len(mother_based_vbf_quarks)

        lhe_vbf_pair = None
        if len(mother_based_vbf_quarks) >= 2:
            best_deta = -1.0
            for i in range(len(mother_based_vbf_quarks)):
                for j in range(i + 1, len(mother_based_vbf_quarks)):
                    idx1, lp1 = mother_based_vbf_quarks[i]
                    idx2, lp2 = mother_based_vbf_quarks[j]
                    deta = abs(lp1.eta - lp2.eta)
                    dphi = abs(delta_phi(lp1.phi, lp2.phi))
                    dr = delta_r(lp1.eta, lp1.phi, lp2.eta, lp2.phi)
                    mjj = (p4_from_lhepart(lp1) + p4_from_lhepart(lp2)).M()
                    if deta > best_deta:
                        best_deta = deta
                        lhe_vbf_pair = (idx1, lp1, idx2, lp2, deta, dphi, dr, mjj)
        else:
            lhe_vbf_pair = choose_vbf_lhe_pair(lheparts, excluded_lhe_indices)
            if lhe_vbf_pair is not None:
                self.used_lhe_vbf_topology_fallback[0] = 1

        if lhe_vbf_pair is not None:
            self.has_lhe_vbf_pair[0] = 1
            idx1, lp1, idx2, lp2, deta, dphi, dr, mjj = lhe_vbf_pair
            self.lhe_vbf1_pdgId[0] = getattr(lp1, "pdgId", 0)
            self.lhe_vbf1_pt[0] = lp1.pt
            self.lhe_vbf1_eta[0] = lp1.eta
            self.lhe_vbf1_phi[0] = lp1.phi
            self.lhe_vbf1_mass[0] = lp1.mass
            self.lhe_vbf2_pdgId[0] = getattr(lp2, "pdgId", 0)
            self.lhe_vbf2_pt[0] = lp2.pt
            self.lhe_vbf2_eta[0] = lp2.eta
            self.lhe_vbf2_phi[0] = lp2.phi
            self.lhe_vbf2_mass[0] = lp2.mass
            self.lhe_vbf_pair_dEta[0] = deta
            self.lhe_vbf_pair_dPhi[0] = dphi
            self.lhe_vbf_pair_dR[0] = dr
            self.lhe_vbf_pair_mass[0] = mjj

        self.cutflow.Fill(1)
        if trigger_bits["pass_baseline"]:
            self.cutflow.Fill(2)

        if len(central_jets) >= 2:
            self.has_central_pair[0] = 1
            self.n_with_central_pair += 1

            jc1, jc2 = central_jets[0], central_jets[1]
            central_pair_p4 = p4_from_recojet(jc1) + p4_from_recojet(jc2)
            central_pair_raw_p4 = (
                p4_from_recojet_raw(jc1.original) + p4_from_recojet_raw(jc2.original)
            )
            central_deta = abs(jc1.eta - jc2.eta)
            central_dphi = abs(delta_phi(jc1.phi, jc2.phi))
            central_dr = delta_r(jc1.eta, jc1.phi, jc2.eta, jc2.phi)

            self.central_jet1_pt[0] = jc1.pt
            self.central_jet1_eta[0] = jc1.eta
            self.central_jet1_phi[0] = jc1.phi
            self.central_jet1_mass[0] = recojet_mass(jc1)
            self.central_jet1_pt_raw[0] = getattr(jc1.original, "pt", -999.0)
            self.central_jet2_pt[0] = jc2.pt
            self.central_jet2_eta[0] = jc2.eta
            self.central_jet2_phi[0] = jc2.phi
            self.central_jet2_mass[0] = recojet_mass(jc2)
            self.central_jet2_pt_raw[0] = getattr(jc2.original, "pt", -999.0)
            self.central_pair_dEta[0] = central_deta
            self.central_pair_dPhi[0] = central_dphi
            self.central_pair_dR[0] = central_dr
            self.central_pair_mass[0] = central_pair_p4.M()
            self.central_pair_mass_raw[0] = central_pair_raw_p4.M()

            if args.centralPairDEtaMin <= 0.0 or central_deta >= args.centralPairDEtaMin:
                self.pass_central_pair_deta[0] = 1

            widejet_info = build_widejets_from_seeds(good_jets, jc1, jc2, args.widejetRadius)
            widejet1 = widejet_info["widejet1"]
            widejet2 = widejet_info["widejet2"]
            self.has_widejet_pair[0] = 1
            self.widejet1_pt[0] = widejet1.Pt()
            self.widejet1_eta[0] = widejet1.Eta()
            self.widejet1_phi[0] = widejet1.Phi()
            self.widejet1_mass[0] = widejet1.M()
            self.widejet2_pt[0] = widejet2.Pt()
            self.widejet2_eta[0] = widejet2.Eta()
            self.widejet2_phi[0] = widejet2.Phi()
            self.widejet2_mass[0] = widejet2.M()
            self.widejet_pair_dEta[0] = abs(widejet1.Eta() - widejet2.Eta())
            self.widejet_pair_dPhi[0] = abs(delta_phi(widejet1.Phi(), widejet2.Phi()))
            self.widejet_pair_dR[0] = delta_r(widejet1.Eta(), widejet1.Phi(), widejet2.Eta(), widejet2.Phi())
            self.widejet_pair_mass[0] = (widejet1 + widejet2).M()
            widejet_pair_mass = self.widejet_pair_mass[0]
            leading_widejet_pt = max(widejet1.Pt(), widejet2.Pt())
            subleading_widejet_pt = min(widejet1.Pt(), widejet2.Pt())
            if widejet_pair_mass > 0.0:
                self.widejet_leading_pt_over_m[0] = leading_widejet_pt / widejet_pair_mass
                self.widejet_subleading_pt_over_m[0] = subleading_widejet_pt / widejet_pair_mass
            else:
                self.widejet_leading_pt_over_m[0] = -999.0
                self.widejet_subleading_pt_over_m[0] = -999.0
            self.n_widejet1_constituents[0] = len(widejet_info["assigned_to_widejet1_ids"])
            self.n_widejet2_constituents[0] = len(widejet_info["assigned_to_widejet2_ids"])

            widejet_match_dr = args.widejetGenMatchDR if args.widejetGenMatchDR is not None else args.widejetRadius
            widejet_match_result = match_widejet_pair_to_higgs_daughters(
                widejet1,
                widejet2,
                daughters,
                widejet_match_dr,
            )
            if widejet_match_result is not None:
                dr1, dr2, maxdr, is_matched = widejet_match_result
                self.central_pair_match_dR1[0] = dr1
                self.central_pair_match_dR2[0] = dr2
                self.central_pair_match_maxdr[0] = maxdr
                if is_matched:
                    self.central_pair_genmatched[0] = 1
                    self.n_with_central_pair_genmatch += 1

            pass_step1_pre_deta = True
            if args.centralPairDEtaMin > 0.0 and central_deta < args.centralPairDEtaMin:
                pass_step1_pre_deta = False
            if args.centralPairDEtaMax > 0.0 and central_deta > args.centralPairDEtaMax:
                pass_step1_pre_deta = False
            if args.centralPairDPhiMin > 0.0 and central_dphi < args.centralPairDPhiMin:
                pass_step1_pre_deta = False
            if args.leadingPtOverMMin > 0.0 and self.widejet_leading_pt_over_m[0] < args.leadingPtOverMMin:
                pass_step1_pre_deta = False
            if args.subleadingPtOverMMin > 0.0 and self.widejet_subleading_pt_over_m[0] < args.subleadingPtOverMMin:
                pass_step1_pre_deta = False

            if pass_step1_pre_deta:
                self.pass_step1_kinematics[0] = 1
                self.n_step1_pre_deta += 1
                self.cutflow.Fill(3)

            remaining_jets = widejet_info["leftover_jets"]
            for pt_threshold in JET_MULTIPLICITY_PT_THRESHOLDS:
                postwide_counts = count_recojets_by_region(remaining_jets, pt_threshold)
                for region, count in postwide_counts.items():
                    self.recojet_mult_branches[f"n_recojets_postwide_{region}_pt{pt_threshold}"][0] = count
            forward_remaining_jets = [
                jet for jet in remaining_jets
                if (
                    jet.pt > args.forwardJetPtMin
                    and abs(jet.eta) <= args.forwardEtaMax
                    and (
                        not args.requireForwardEtaMinOnCentralFirst
                        or abs(jet.eta) >= args.forwardEtaMin
                    )
                )
            ]
            preclean_vbf_pair = choose_forward_pair_by_mode(forward_remaining_jets, args.forwardPairAlgo)
            if preclean_vbf_pair is not None and pass_vbf_requirement(preclean_vbf_pair):
                self.pass_preclean_selection[0] = 1

            pass_step1 = pass_step1_pre_deta
            if args.widejetPairDEtaMin > 0.0 and self.widejet_pair_dEta[0] < args.widejetPairDEtaMin:
                pass_step1 = False
            if args.widejetPairDEtaMax > 0.0 and self.widejet_pair_dEta[0] > args.widejetPairDEtaMax:
                pass_step1 = False
            if args.widejetPairDPhiMin > 0.0 and self.widejet_pair_dPhi[0] < args.widejetPairDPhiMin:
                pass_step1 = False

            if not pass_step1:
                self.fill_multiplicity_hists()
                self.tree.Fill()
                return True

            self.pass_central_pair_selection[0] = 1
            self.n_step1 += 1
            self.cutflow.Fill(4)
            self.h_central_mass_selected.Fill(self.widejet_pair_mass[0])
            self.h_central_deta_selected.Fill(self.widejet_pair_dEta[0])
            if self.central_pair_genmatched[0] == 1:
                self.h_central_mass_matched.Fill(self.widejet_pair_mass[0])
            else:
                self.h_central_mass_unmatched.Fill(self.widejet_pair_mass[0])

            self.n_recojets_remaining[0] = len(remaining_jets)
            self.n_recojets_forward_remaining[0] = len(forward_remaining_jets)

            vbf_pair = preclean_vbf_pair
            if vbf_pair is not None:
                self.has_vbf_pair[0] = 1
                self.n_with_vbf_pair += 1

                jf1, jf2, deta, dphi, dr, mjj = vbf_pair
                self.forward_jet1_pt[0] = jf1.pt
                self.forward_jet1_eta[0] = jf1.eta
                self.forward_jet1_phi[0] = jf1.phi
                self.forward_jet1_mass[0] = recojet_mass(jf1)
                self.forward_jet2_pt[0] = jf2.pt
                self.forward_jet2_eta[0] = jf2.eta
                self.forward_jet2_phi[0] = jf2.phi
                self.forward_jet2_mass[0] = recojet_mass(jf2)
                self.forward_pair_dEta[0] = deta
                self.forward_pair_dPhi[0] = dphi
                self.forward_pair_dR[0] = dr
                self.forward_pair_mass[0] = mjj
                self.cutflow.Fill(5)

                if pass_vbf_requirement(vbf_pair):
                    self.pass_minimal_vbf[0] = 1
                    self.isVBF[0] = 1
                    self.n_with_vbf_tag += 1
                    self.cutflow.Fill(6)
                    self.h_forward_mass_isvbf.Fill(mjj)
                    self.h_forward_deta_isvbf.Fill(deta)
                else:
                    self.isPP[0] = 1
                    self.n_with_pp_tag += 1
                    self.h_forward_mass_ispp.Fill(mjj)
                    self.h_forward_deta_ispp.Fill(deta)

                if lhe_vbf_pair is not None:
                    _, lp1, _, lp2, _, _, _, _ = lhe_vbf_pair
                    dr1, dr2 = match_pair_by_delta_r(jf1, jf2, lp1, lp2)
                    self.vbf_pair_match_dR1[0] = dr1
                    self.vbf_pair_match_dR2[0] = dr2
                    self.vbf_pair_match_maxdr[0] = max(dr1, dr2)
                    if max(dr1, dr2) < args.matchDR:
                        self.vbf_pair_lhematched[0] = 1
                        self.n_with_vbf_lhe_match += 1
            else:
                self.isPP[0] = 1
                self.n_with_pp_tag += 1

        if self.has_central_pair[0] == 1 and self.pass_central_pair_selection[0] == 1:
            self.pass_final_selection[0] = 1

        self.fill_multiplicity_hists()
        self.tree.Fill()
        return True

    def endJob(self):
        n_after_trigger = self.tree.GetEntries("pass_trigger_baseline == 1")
        n_central_pair = self.tree.GetEntries("pass_trigger_baseline == 1 && has_central_pair == 1")
        n_step1_pre_deta = self.tree.GetEntries("pass_trigger_baseline == 1 && pass_step1_kinematics == 1")
        n_central_pair_pass = self.tree.GetEntries("pass_trigger_baseline == 1 && pass_central_pair_selection == 1")
        n_central_pair_genmatched = self.tree.GetEntries(
            "pass_trigger_baseline == 1 && pass_central_pair_selection == 1 && central_pair_genmatched == 1"
        )
        n_central_pair_not_matched = self.tree.GetEntries(
            "pass_trigger_baseline == 1 && pass_central_pair_selection == 1 && central_pair_genmatched == 0"
        )
        n_vbf_pair = self.tree.GetEntries("pass_trigger_baseline == 1 && pass_central_pair_selection == 1 && has_vbf_pair == 1")
        n_is_vbf = self.tree.GetEntries("pass_trigger_baseline == 1 && isVBF == 1")
        n_is_pp = self.tree.GetEntries("pass_trigger_baseline == 1 && isPP == 1")
        n_truth_vbf = self.tree.GetEntries("pass_trigger_baseline == 1 && has_lhe_vbf_pair == 1")
        n_vbf_pair_matched = self.tree.GetEntries(
            "pass_trigger_baseline == 1 && has_vbf_pair == 1 && vbf_pair_lhematched == 1"
        )
        n_vbf_lhematched = self.tree.GetEntries(
            "pass_trigger_baseline == 1 && isVBF == 1 && vbf_pair_lhematched == 1"
        )
        central_purity = 100.0 * n_central_pair_genmatched / n_central_pair_pass if n_central_pair_pass > 0 else 0.0
        vbf_pair_match_purity = 100.0 * n_vbf_pair_matched / n_vbf_pair if n_vbf_pair > 0 else 0.0
        vbf_tag_match_purity = 100.0 * n_vbf_lhematched / n_is_vbf if n_is_vbf > 0 else 0.0
        vbf_tag_match_eff = 100.0 * n_vbf_lhematched / n_truth_vbf if n_truth_vbf > 0 else 0.0

        print("\n=========== ScoutNano Central WideJet + VBF Study ==========")
        print(f"Total events:                     {self.total}")
        print(f"Passing trigger:                  {n_after_trigger}")
        print(f"Central jet pT min:               {args.centralJetPtMin}")
        print(f"Central jet |eta| max:            {args.centralEtaMax}")
        if args.centralPairDEtaMin > 0.0:
            print(f"Central pair |DeltaEta| min:      {args.centralPairDEtaMin}")
        else:
            print("Central pair |DeltaEta| min:      disabled")
        if args.centralPairDEtaMax > 0.0:
            print(f"Central pair |DeltaEta| max:      {args.centralPairDEtaMax}")
        else:
            print("Central pair |DeltaEta| max:      disabled")
        if args.centralPairDPhiMin > 0.0:
            print(f"Central pair |DeltaPhi| min:      {args.centralPairDPhiMin}")
        else:
            print("Central pair |DeltaPhi| min:      disabled")
        if args.leadingPtOverMMin > 0.0:
            print(f"Leading widejet pT/m min:         {args.leadingPtOverMMin}")
        else:
            print("Leading widejet pT/m min:         disabled")
        if args.subleadingPtOverMMin > 0.0:
            print(f"Subleading widejet pT/m min:      {args.subleadingPtOverMMin}")
        else:
            print("Subleading widejet pT/m min:      disabled")
        print(f"Wide-jet radius:                  {args.widejetRadius}")
        if args.widejetPairDEtaMin > 0.0:
            print(f"Wide-jet pair |DeltaEta| min:     {args.widejetPairDEtaMin}")
        else:
            print("Wide-jet pair |DeltaEta| min:     disabled")
        if args.widejetPairDEtaMax > 0.0:
            print(f"Wide-jet pair |DeltaEta| max:     {args.widejetPairDEtaMax}")
        else:
            print("Wide-jet pair |DeltaEta| max:     disabled")
        if args.widejetPairDPhiMin > 0.0:
            print(f"Wide-jet pair |DeltaPhi| min:     {args.widejetPairDPhiMin}")
        else:
            print("Wide-jet pair |DeltaPhi| min:     disabled")
        print(f"VBF leftover-jet pT min:          {args.forwardJetPtMin}")
        print(f"Central-to-VBFTag require |eta| min: {args.requireForwardEtaMinOnCentralFirst}")
        if args.requireForwardEtaMinOnCentralFirst:
            print(f"Central-to-VBFTag jet |eta| min:  {args.forwardEtaMin}")
        print(f"VBF leftover-jet |eta| max:       {args.forwardEtaMax}")
        print(f"Forward pair algorithm:           {args.forwardPairAlgo}")
        if args.forwardPairAlgo == "maxDEta":
            print(f"VBF pair |DeltaEta| min:          {args.forwardPairDEtaMin}")
        elif args.forwardPairAlgo == "maxMjj":
            print(f"VBF pair mjj min [GeV]:           {args.forwardPairMjjMin}")
        if args.forwardPairExtraMjjMin > 0.0:
            print(f"Extra forward-pair mjj min [GeV]: {args.forwardPairExtraMjjMin}")
        else:
            print("Extra forward-pair mjj min [GeV]: disabled")
        print(f"Require minimal VBF selection:    {args.requireMinimalVBF}")
        print(f"Passing central seed pair:        {n_central_pair}")
        print(f"Passing STEP 1 before wide #Delta#eta: {n_step1_pre_deta}")
        print(f"Passing STEP 1 after wide #Delta#eta:  {n_central_pair_pass}")
        print(f"Events with a VBF pair candidate found: {n_vbf_pair}")
        print(f"Events passing VBF selection:     {n_is_vbf}")
        print(f"Number of events passing STEP 1:  {n_central_pair_pass}")
        print(f"STEP 2 split: isVBF / isPP        {n_is_vbf} / {n_is_pp}")
        print(f"STEP 1 matched central pairs:     {n_central_pair_genmatched}")
        print(f"STEP 1 unmatched central pairs:   {n_central_pair_not_matched}")
        print(f"Leftover pairs matched to LHE:    {n_vbf_pair_matched}")
        print(f"Leftover-pair match purity [%]:   {vbf_pair_match_purity:.2f}")
        print(f"VBF pairs matched to LHE quarks:  {n_vbf_lhematched}")
        print(f"VBF-tag match efficiency [%]:     {vbf_tag_match_eff:.2f}")
        print(f"VBF-tag match purity [%]:         {vbf_tag_match_purity:.2f}")
        print(f"Central-pair match purity [%]:    {central_purity:.2f}")
        print("============================================================\n")

        self.outfile.cd()
        self.tree.Write()
        self.cutflow.Write()
        self.mult_dir.cd()
        for hist in self.recojet_mult_hists.values():
            hist.Write()
        self.outfile.cd()
        self.h_central_mass_selected.Write()
        self.h_central_deta_selected.Write()
        self.h_forward_mass_isvbf.Write()
        self.h_forward_deta_isvbf.Write()
        self.h_forward_mass_ispp.Write()
        self.h_forward_deta_ispp.Write()
        self.h_central_mass_matched.Write()
        self.h_central_mass_unmatched.Write()
        self.outfile.Close()
        Module.endJob(self)


if __name__ == "__main__":
    input_files = get_input_files()
    output_name = default_output_name()

    print(f"[INFO] Using {len(input_files)} input file(s)")
    print(f"[INFO] First file: {input_files[0]}")
    print(f"[INFO] Output file: {output_name}")

    processor = PostProcessor(
        ".",
        input_files,
        modules=[CentralWideJetVBFScoutNanoStudy(output_name)],
        noOut=True,
        maxEntries=args.maxEntries,
        firstEntry=args.firstEntry,
        prefetch=args.prefetch,
    )
    processor.run()
