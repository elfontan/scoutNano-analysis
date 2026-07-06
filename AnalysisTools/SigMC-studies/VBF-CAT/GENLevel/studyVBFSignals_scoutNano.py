#!/usr/bin/env python3

"""
Minimal generator-level Higgs study for ScoutNano signal samples.

The script:
 - reads ScoutNano files with NanoAODTools
 - finds a generated Higgs boson (pdgId = 25)
 - stores the Higgs pt/eta/phi/mass
 - identifies its final non-self decay products
 - stores a two-body decay summary using the two daughters sorted by pt

Input can be provided either explicitly with --input or via --decay and --mass.
"""

import os
import math
import glob
import argparse
from array import array

import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

ROOT.gROOT.SetBatch(True)

BASE_DIR = "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/VBFH-Central/"
#BASE_DIR = "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/VBFHToXX"


parser = argparse.ArgumentParser(description="Generator-level ScoutNano Higgs study")
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
    "--forwardEtaMin",
    type=float,
    default=3.0,
    help="Minimum |eta| for forward gen jets.",
)
parser.add_argument(
    "--forwardEtaMax",
    type=float,
    default=5.0,
    help="Maximum |eta| for forward gen jets.",
)
parser.add_argument(
    "--forwardPairDEtaMin",
    type=float,
    default=5.0,
    help="Minimum |DeltaEta| for the chosen VBF gen-jet pair when required.",
)
parser.add_argument(
    "--vbfTagMode",
    type=str,
    default="twoForward_maxDEta_dEtaGtMin",
    choices=[
        "twoForward_maxDEta",
        "twoForward_maxDEta_dEtaGtMin",
        "oneForward_maxDEta",
        "oneForward_maxDEta_dEtaGtMin",
        "oneForward_maxMjj",
        "maxDEta",
        "maxMjj",
    ],
    help=(
        "Strategy used to choose the VBF gen-jet pair. "
        "'twoForward_maxDEta_dEtaGtMin' reproduces the baseline."
    ),
)
parser.add_argument(
    "--genJetLheMatchDRMax",
    type=float,
    default=0.4,
    help="Maximum DeltaR used to declare a selected GenJet matched to an LHE VBF quark.",
)
args = parser.parse_args()


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


def p4_from_genpart(gp):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(gp.pt, gp.eta, gp.phi, gp.mass)
    return vec


def p4_from_genjet(gj):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(gj.pt, gj.eta, gj.phi, gj.mass)
    return vec


def p4_from_lhepart(lp):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(lp.pt, lp.eta, lp.phi, lp.mass)
    return vec


def is_final_state_lhepart(lp):
    return getattr(lp, "status", -999) == 1


def is_lhe_quark(lp):
    return abs(getattr(lp, "pdgId", 0)) in {1, 2, 3, 4, 5, 6}


def is_lhe_gluon(lp):
    return abs(getattr(lp, "pdgId", 0)) == 21


def lhe_mothers(lp):
    return (
        getattr(lp, "firstMotherIdx", -1),
        getattr(lp, "lastMotherIdx", -1),
    )


def lhe_has_mother_info(lp):
    first_mother, last_mother = lhe_mothers(lp)
    return first_mother >= 0 or last_mother >= 0


def find_lhe_higgs_index(lheparts):
    for idx, lp in enumerate(lheparts):
        if abs(getattr(lp, "pdgId", 0)) == 25 and getattr(lp, "status", 0) == 2:
            return idx
    for idx, lp in enumerate(lheparts):
        if abs(getattr(lp, "pdgId", 0)) == 25:
            return idx
    return None


def find_higgs_lhe_daughters_from_mothers(lheparts):
    higgs_idx = find_lhe_higgs_index(lheparts)
    if higgs_idx is None:
        return None, None

    daughters = []
    daughter_indices = []
    for idx, lp in enumerate(lheparts):
        if not is_final_state_lhepart(lp):
            continue
        first_mother, last_mother = lhe_mothers(lp)
        if first_mother <= higgs_idx <= last_mother:
            daughter_indices.append(idx)
            daughters.append(lp)

    if len(daughters) != 2:
        return higgs_idx, None
    return higgs_idx, daughters


def is_incoming_lhepart(lp):
    return getattr(lp, "status", 0) == -1


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
        f"VBFHTo{args.decay}_M{args.mass}_genHiggs_scoutNano_"
        f"VBFTagMode-{args.vbfTagMode}.root"
    )


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


def match_lhe_particles_to_genparts(lheparts, genparts_to_match, max_dr=0.4):
    matched = []
    used_indices = set()

    for gp in genparts_to_match:
        best_idx = None
        best_dr = max_dr
        for idx, lp in enumerate(lheparts):
            if idx in used_indices:
                continue
            if abs(getattr(lp, "pdgId", 0)) != abs(gp.pdgId):
                continue
            dr = delta_r(lp.eta, lp.phi, gp.eta, gp.phi)
            if dr < best_dr:
                best_dr = dr
                best_idx = idx

        if best_idx is None:
            return None

        used_indices.add(best_idx)
        matched.append((best_idx, lheparts[best_idx]))

    return matched


def choose_vbf_lhe_pair(lheparts, excluded_indices=None):
    if excluded_indices is None:
        excluded_indices = set()

    candidate_quarks = [
        (idx, lp)
        for idx, lp in enumerate(lheparts)
        if idx not in excluded_indices and is_final_state_lhepart(lp) and is_lhe_quark(lp)
    ]
    if len(candidate_quarks) < 2:
        return None, candidate_quarks

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

    return best_pair, candidate_quarks


def match_genjet_pair_to_lhe_pair(genjet1, genjet2, lhe1_eta, lhe1_phi, lhe2_eta, lhe2_phi):
    assignment_a = (
        delta_r(genjet1.eta, genjet1.phi, lhe1_eta, lhe1_phi),
        delta_r(genjet2.eta, genjet2.phi, lhe2_eta, lhe2_phi),
    )
    assignment_b = (
        delta_r(genjet1.eta, genjet1.phi, lhe2_eta, lhe2_phi),
        delta_r(genjet2.eta, genjet2.phi, lhe1_eta, lhe1_phi),
    )

    score_a = max(assignment_a)
    score_b = max(assignment_b)
    if score_a <= score_b:
        return assignment_a
    return assignment_b


def best_forward_genjet_pair(genjets):
    forward = [gj for gj in genjets if args.forwardEtaMin <= abs(gj.eta) <= args.forwardEtaMax]
    if len(forward) < 2:
        return None, forward

    best_pair = None
    best_deta = -1.0
    for i in range(len(forward)):
        for j in range(i + 1, len(forward)):
            deta = abs(forward[i].eta - forward[j].eta)
            if deta > best_deta:
                best_deta = deta
                best_pair = (forward[i], forward[j])

    return best_pair, forward


def choose_vbf_genjet_pair(genjets, forward_genjets, mode):
    if mode in {"twoForward_maxDEta", "twoForward_maxDEta_dEtaGtMin"}:
        candidate_jets = forward_genjets
        require_one_forward = False
    elif mode in {"oneForward_maxDEta", "oneForward_maxDEta_dEtaGtMin", "oneForward_maxMjj"}:
        candidate_jets = genjets
        require_one_forward = True
    elif mode in {"maxDEta", "maxMjj"}:
        candidate_jets = genjets
        require_one_forward = False
    else:
        raise ValueError(f"Unknown vbfTagMode: {mode}")

    if len(candidate_jets) < 2:
        return None

    forward_ids = {id(jet) for jet in forward_genjets}
    candidates = []
    for i in range(len(candidate_jets)):
        for j in range(i + 1, len(candidate_jets)):
            gj1, gj2 = candidate_jets[i], candidate_jets[j]
            n_forward = int(id(gj1) in forward_ids) + int(id(gj2) in forward_ids)
            if require_one_forward and n_forward < 1:
                continue

            deta = abs(gj1.eta - gj2.eta)
            dphi = abs(delta_phi(gj1.phi, gj2.phi))
            dr = delta_r(gj1.eta, gj1.phi, gj2.eta, gj2.phi)
            mjj = (p4_from_genjet(gj1) + p4_from_genjet(gj2)).M()
            candidates.append((gj1, gj2, deta, dphi, dr, mjj))

    if not candidates:
        return None

    if mode in {
        "twoForward_maxDEta",
        "twoForward_maxDEta_dEtaGtMin",
        "oneForward_maxDEta",
        "oneForward_maxDEta_dEtaGtMin",
        "maxDEta",
    }:
        return max(candidates, key=lambda item: item[2])
    if mode in {"oneForward_maxMjj", "maxMjj"}:
        return max(candidates, key=lambda item: item[5])
    raise ValueError(f"Unknown vbfTagMode: {mode}")


def vbf_genjet_pair_passes_selection(pair_info):
    if pair_info is None:
        return False

    _, _, deta, _, _, _ = pair_info
    if args.vbfTagMode in {"twoForward_maxDEta_dEtaGtMin", "oneForward_maxDEta_dEtaGtMin"}:
        return deta > args.forwardPairDEtaMin
    return True


def describe_vbf_tag_mode():
    if args.vbfTagMode == "twoForward_maxDEta":
        return (
            f">=2 forward gen jets with {args.forwardEtaMin} <= |eta| <= {args.forwardEtaMax}; "
            "choose the pair with largest |DeltaEta|"
        )
    if args.vbfTagMode == "twoForward_maxDEta_dEtaGtMin":
        return (
            f">=2 forward gen jets with {args.forwardEtaMin} <= |eta| <= {args.forwardEtaMax}; "
            f"choose the pair with largest |DeltaEta| and require |DeltaEta| > {args.forwardPairDEtaMin}"
        )
    if args.vbfTagMode == "oneForward_maxDEta":
        return (
            f">=1 forward gen jet with {args.forwardEtaMin} <= |eta| <= {args.forwardEtaMax}; "
            "choose the pair with largest |DeltaEta| and leave the second jet unrestricted"
        )
    if args.vbfTagMode == "oneForward_maxDEta_dEtaGtMin":
        return (
            f">=1 forward gen jet with {args.forwardEtaMin} <= |eta| <= {args.forwardEtaMax}; "
            f"choose the pair with largest |DeltaEta|, second jet unrestricted, require |DeltaEta| > {args.forwardPairDEtaMin}"
        )
    if args.vbfTagMode == "oneForward_maxMjj":
        return (
            f">=1 forward gen jet with {args.forwardEtaMin} <= |eta| <= {args.forwardEtaMax}; "
            "choose the pair with largest m_jj and leave the second jet unrestricted"
        )
    if args.vbfTagMode == "maxDEta":
        return "Choose the pair with largest |DeltaEta| among all gen jets"
    if args.vbfTagMode == "maxMjj":
        return "Choose the pair with largest m_jj among all gen jets"
    return args.vbfTagMode


class GenHiggsScoutNanoStudy(Module):
    def __init__(self, output_name):
        self.output_name = output_name
        self.writeHistFile = False
        self.total = 0
        self.n_higgs_found = 0
        self.n_two_body = 0
        self.n_events_with_lhe = 0
        self.n_events_with_mother_based_vbf_quarks = 0
        self.n_events_with_topology_fallback = 0
        self.n_events_with_selected_genjet_pair = 0
        self.n_events_with_truth_lhe_vbf_pair = 0
        self.n_events_with_genjet_lhe_comparison = 0
        self.n_events_with_genjet_lhe_match = 0
        self.n_events_with_tagged_genjet_lhe_comparison = 0
        self.n_events_with_tagged_genjet_lhe_match = 0

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        self.outfile = ROOT.TFile(self.output_name, "RECREATE")
        self.tree = ROOT.TTree("Events", "Generator-level Higgs tree from ScoutNano")

        self.run = array("I", [0])
        self.lumi = array("I", [0])
        self.event = array("L", [0])
        self.n_decay_products = array("i", [0])
        self.n_genjets_total = array("i", [0])
        self.n_genjets_central = array("i", [0])
        self.n_genjets_forward = array("i", [0])
        self.has_forward_genjet_pair = array("i", [0])
        self.n_lheparts_total = array("i", [0])
        self.n_lheparts_final_state = array("i", [0])
        self.n_lhe_final_state_partons = array("i", [0])
        self.n_lhe_final_state_quarks = array("i", [0])
        self.n_lhe_final_state_gluons = array("i", [0])
        self.n_lhe_vbf_quarks = array("i", [0])
        self.has_lhe_vbf_pair = array("i", [0])
        self.used_lhe_vbf_topology_fallback = array("i", [0])
        self.genjet_lhe_match_n = array("i", [0])
        self.has_genjet_lhe_match_pair = array("i", [0])
        self.genjet_lhe_match_dr1 = array("f", [0.0])
        self.genjet_lhe_match_dr2 = array("f", [0.0])
        self.genjet_lhe_match_maxdr = array("f", [0.0])
        self.genjet_lhe_match_sumdr = array("f", [0.0])

        self.higgs_pt = array("f", [0.0])
        self.higgs_eta = array("f", [0.0])
        self.higgs_phi = array("f", [0.0])
        self.higgs_mass = array("f", [0.0])

        self.decay1_pdgId = array("i", [0])
        self.decay2_pdgId = array("i", [0])
        self.decay_dR = array("f", [0.0])
        self.decay_dEta = array("f", [0.0])
        self.decay_dPhi = array("f", [0.0])
        self.decay_mass = array("f", [0.0])

        self.lead_decay_pt = array("f", [0.0])
        self.lead_decay_eta = array("f", [0.0])
        self.lead_decay_phi = array("f", [0.0])
        self.sublead_decay_pt = array("f", [0.0])
        self.sublead_decay_eta = array("f", [0.0])
        self.sublead_decay_phi = array("f", [0.0])

        self.forward_genjet1_pt = array("f", [0.0])
        self.forward_genjet1_eta = array("f", [0.0])
        self.forward_genjet1_phi = array("f", [0.0])
        self.forward_genjet1_mass = array("f", [0.0])
        self.forward_genjet2_pt = array("f", [0.0])
        self.forward_genjet2_eta = array("f", [0.0])
        self.forward_genjet2_phi = array("f", [0.0])
        self.forward_genjet2_mass = array("f", [0.0])
        self.forward_genjet_pair_dEta = array("f", [0.0])
        self.forward_genjet_pair_dPhi = array("f", [0.0])
        self.forward_genjet_pair_dR = array("f", [0.0])
        self.forward_genjet_pair_mass = array("f", [0.0])

        self.lhe_vbf1_pdgId = array("i", [0])
        self.lhe_vbf2_pdgId = array("i", [0])
        self.lhe_vbf1_pt = array("f", [0.0])
        self.lhe_vbf1_eta = array("f", [0.0])
        self.lhe_vbf1_phi = array("f", [0.0])
        self.lhe_vbf1_mass = array("f", [0.0])
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
        self.tree.Branch("n_decay_products", self.n_decay_products, "n_decay_products/I")
        self.tree.Branch("n_genjets_total", self.n_genjets_total, "n_genjets_total/I")
        self.tree.Branch("n_genjets_central", self.n_genjets_central, "n_genjets_central/I")
        self.tree.Branch("n_genjets_forward", self.n_genjets_forward, "n_genjets_forward/I")
        self.tree.Branch("has_forward_genjet_pair", self.has_forward_genjet_pair, "has_forward_genjet_pair/I")
        self.tree.Branch("n_lheparts_total", self.n_lheparts_total, "n_lheparts_total/I")
        self.tree.Branch("n_lheparts_final_state", self.n_lheparts_final_state, "n_lheparts_final_state/I")
        self.tree.Branch("n_lhe_final_state_partons", self.n_lhe_final_state_partons, "n_lhe_final_state_partons/I")
        self.tree.Branch("n_lhe_final_state_quarks", self.n_lhe_final_state_quarks, "n_lhe_final_state_quarks/I")
        self.tree.Branch("n_lhe_final_state_gluons", self.n_lhe_final_state_gluons, "n_lhe_final_state_gluons/I")
        self.tree.Branch("n_lhe_vbf_quarks", self.n_lhe_vbf_quarks, "n_lhe_vbf_quarks/I")
        self.tree.Branch("has_lhe_vbf_pair", self.has_lhe_vbf_pair, "has_lhe_vbf_pair/I")
        self.tree.Branch(
            "used_lhe_vbf_topology_fallback",
            self.used_lhe_vbf_topology_fallback,
            "used_lhe_vbf_topology_fallback/I",
        )
        self.tree.Branch("genjet_lhe_match_n", self.genjet_lhe_match_n, "genjet_lhe_match_n/I")
        self.tree.Branch(
            "has_genjet_lhe_match_pair",
            self.has_genjet_lhe_match_pair,
            "has_genjet_lhe_match_pair/I",
        )
        self.tree.Branch("genjet_lhe_match_dr1", self.genjet_lhe_match_dr1, "genjet_lhe_match_dr1/F")
        self.tree.Branch("genjet_lhe_match_dr2", self.genjet_lhe_match_dr2, "genjet_lhe_match_dr2/F")
        self.tree.Branch("genjet_lhe_match_maxdr", self.genjet_lhe_match_maxdr, "genjet_lhe_match_maxdr/F")
        self.tree.Branch("genjet_lhe_match_sumdr", self.genjet_lhe_match_sumdr, "genjet_lhe_match_sumdr/F")

        self.tree.Branch("higgs_pt", self.higgs_pt, "higgs_pt/F")
        self.tree.Branch("higgs_eta", self.higgs_eta, "higgs_eta/F")
        self.tree.Branch("higgs_phi", self.higgs_phi, "higgs_phi/F")
        self.tree.Branch("higgs_mass", self.higgs_mass, "higgs_mass/F")

        self.tree.Branch("decay1_pdgId", self.decay1_pdgId, "decay1_pdgId/I")
        self.tree.Branch("decay2_pdgId", self.decay2_pdgId, "decay2_pdgId/I")
        self.tree.Branch("decay_dR", self.decay_dR, "decay_dR/F")
        self.tree.Branch("decay_dEta", self.decay_dEta, "decay_dEta/F")
        self.tree.Branch("decay_dPhi", self.decay_dPhi, "decay_dPhi/F")
        self.tree.Branch("decay_mass", self.decay_mass, "decay_mass/F")

        self.tree.Branch("lead_decay_pt", self.lead_decay_pt, "lead_decay_pt/F")
        self.tree.Branch("lead_decay_eta", self.lead_decay_eta, "lead_decay_eta/F")
        self.tree.Branch("lead_decay_phi", self.lead_decay_phi, "lead_decay_phi/F")
        self.tree.Branch("sublead_decay_pt", self.sublead_decay_pt, "sublead_decay_pt/F")
        self.tree.Branch("sublead_decay_eta", self.sublead_decay_eta, "sublead_decay_eta/F")
        self.tree.Branch("sublead_decay_phi", self.sublead_decay_phi, "sublead_decay_phi/F")

        self.tree.Branch("forward_genjet1_pt", self.forward_genjet1_pt, "forward_genjet1_pt/F")
        self.tree.Branch("forward_genjet1_eta", self.forward_genjet1_eta, "forward_genjet1_eta/F")
        self.tree.Branch("forward_genjet1_phi", self.forward_genjet1_phi, "forward_genjet1_phi/F")
        self.tree.Branch("forward_genjet1_mass", self.forward_genjet1_mass, "forward_genjet1_mass/F")
        self.tree.Branch("forward_genjet2_pt", self.forward_genjet2_pt, "forward_genjet2_pt/F")
        self.tree.Branch("forward_genjet2_eta", self.forward_genjet2_eta, "forward_genjet2_eta/F")
        self.tree.Branch("forward_genjet2_phi", self.forward_genjet2_phi, "forward_genjet2_phi/F")
        self.tree.Branch("forward_genjet2_mass", self.forward_genjet2_mass, "forward_genjet2_mass/F")
        self.tree.Branch("forward_genjet_pair_dEta", self.forward_genjet_pair_dEta, "forward_genjet_pair_dEta/F")
        self.tree.Branch("forward_genjet_pair_dPhi", self.forward_genjet_pair_dPhi, "forward_genjet_pair_dPhi/F")
        self.tree.Branch("forward_genjet_pair_dR", self.forward_genjet_pair_dR, "forward_genjet_pair_dR/F")
        self.tree.Branch("forward_genjet_pair_mass", self.forward_genjet_pair_mass, "forward_genjet_pair_mass/F")

        self.tree.Branch("lhe_vbf1_pdgId", self.lhe_vbf1_pdgId, "lhe_vbf1_pdgId/I")
        self.tree.Branch("lhe_vbf2_pdgId", self.lhe_vbf2_pdgId, "lhe_vbf2_pdgId/I")
        self.tree.Branch("lhe_vbf1_pt", self.lhe_vbf1_pt, "lhe_vbf1_pt/F")
        self.tree.Branch("lhe_vbf1_eta", self.lhe_vbf1_eta, "lhe_vbf1_eta/F")
        self.tree.Branch("lhe_vbf1_phi", self.lhe_vbf1_phi, "lhe_vbf1_phi/F")
        self.tree.Branch("lhe_vbf1_mass", self.lhe_vbf1_mass, "lhe_vbf1_mass/F")
        self.tree.Branch("lhe_vbf2_pt", self.lhe_vbf2_pt, "lhe_vbf2_pt/F")
        self.tree.Branch("lhe_vbf2_eta", self.lhe_vbf2_eta, "lhe_vbf2_eta/F")
        self.tree.Branch("lhe_vbf2_phi", self.lhe_vbf2_phi, "lhe_vbf2_phi/F")
        self.tree.Branch("lhe_vbf2_mass", self.lhe_vbf2_mass, "lhe_vbf2_mass/F")
        self.tree.Branch("lhe_vbf_pair_dEta", self.lhe_vbf_pair_dEta, "lhe_vbf_pair_dEta/F")
        self.tree.Branch("lhe_vbf_pair_dPhi", self.lhe_vbf_pair_dPhi, "lhe_vbf_pair_dPhi/F")
        self.tree.Branch("lhe_vbf_pair_dR", self.lhe_vbf_pair_dR, "lhe_vbf_pair_dR/F")
        self.tree.Branch("lhe_vbf_pair_mass", self.lhe_vbf_pair_mass, "lhe_vbf_pair_mass/F")

    def reset(self):
        self.n_decay_products[0] = 0
        self.n_genjets_total[0] = 0
        self.n_genjets_central[0] = 0
        self.n_genjets_forward[0] = 0
        self.has_forward_genjet_pair[0] = 0
        self.n_lheparts_total[0] = 0
        self.n_lheparts_final_state[0] = 0
        self.n_lhe_final_state_partons[0] = 0
        self.n_lhe_final_state_quarks[0] = 0
        self.n_lhe_final_state_gluons[0] = 0
        self.n_lhe_vbf_quarks[0] = 0
        self.has_lhe_vbf_pair[0] = 0
        self.used_lhe_vbf_topology_fallback[0] = 0
        self.genjet_lhe_match_n[0] = 0
        self.has_genjet_lhe_match_pair[0] = 0
        self.genjet_lhe_match_dr1[0] = -999.0
        self.genjet_lhe_match_dr2[0] = -999.0
        self.genjet_lhe_match_maxdr[0] = -999.0
        self.genjet_lhe_match_sumdr[0] = -999.0

        self.higgs_pt[0] = -999.0
        self.higgs_eta[0] = -999.0
        self.higgs_phi[0] = -999.0
        self.higgs_mass[0] = -999.0

        self.decay1_pdgId[0] = 0
        self.decay2_pdgId[0] = 0
        self.decay_dR[0] = -999.0
        self.decay_dEta[0] = -999.0
        self.decay_dPhi[0] = -999.0
        self.decay_mass[0] = -999.0

        self.lead_decay_pt[0] = -999.0
        self.lead_decay_eta[0] = -999.0
        self.lead_decay_phi[0] = -999.0
        self.sublead_decay_pt[0] = -999.0
        self.sublead_decay_eta[0] = -999.0
        self.sublead_decay_phi[0] = -999.0

        self.forward_genjet1_pt[0] = -999.0
        self.forward_genjet1_eta[0] = -999.0
        self.forward_genjet1_phi[0] = -999.0
        self.forward_genjet1_mass[0] = -999.0
        self.forward_genjet2_pt[0] = -999.0
        self.forward_genjet2_eta[0] = -999.0
        self.forward_genjet2_phi[0] = -999.0
        self.forward_genjet2_mass[0] = -999.0
        self.forward_genjet_pair_dEta[0] = -999.0
        self.forward_genjet_pair_dPhi[0] = -999.0
        self.forward_genjet_pair_dR[0] = -999.0
        self.forward_genjet_pair_mass[0] = -999.0

        self.lhe_vbf1_pdgId[0] = 0
        self.lhe_vbf2_pdgId[0] = 0
        self.lhe_vbf1_pt[0] = -999.0
        self.lhe_vbf1_eta[0] = -999.0
        self.lhe_vbf1_phi[0] = -999.0
        self.lhe_vbf1_mass[0] = -999.0
        self.lhe_vbf2_pt[0] = -999.0
        self.lhe_vbf2_eta[0] = -999.0
        self.lhe_vbf2_phi[0] = -999.0
        self.lhe_vbf2_mass[0] = -999.0
        self.lhe_vbf_pair_dEta[0] = -999.0
        self.lhe_vbf_pair_dPhi[0] = -999.0
        self.lhe_vbf_pair_dR[0] = -999.0
        self.lhe_vbf_pair_mass[0] = -999.0

    def analyze(self, event):
        self.total += 1

        self.run[0] = getattr(event, "run", 0)
        self.lumi[0] = getattr(event, "luminosityBlock", 0)
        self.event[0] = getattr(event, "event", 0)
        self.reset()

        genparts = Collection(event, "GenPart")
        genjets = Collection(event, "GenJet")
        higgs, daughters = find_last_higgs(genparts)
        lheparts = Collection(event, "LHEPart") if hasattr(event, "nLHEPart") else []

        self.n_genjets_total[0] = len(genjets)
        self.n_genjets_central[0] = len([gj for gj in genjets if abs(gj.eta) < args.forwardEtaMin])
        _, forward_genjets = best_forward_genjet_pair(genjets)
        pair = choose_vbf_genjet_pair(genjets, forward_genjets, args.vbfTagMode)
        self.n_genjets_forward[0] = len(forward_genjets)
        if pair is not None:
            self.n_events_with_selected_genjet_pair += 1
        self.n_lheparts_total[0] = len(lheparts)
        final_state_lheparts = [lp for lp in lheparts if is_final_state_lhepart(lp)]
        final_state_lhe_quarks = [lp for lp in final_state_lheparts if is_lhe_quark(lp)]
        final_state_lhe_gluons = [lp for lp in final_state_lheparts if is_lhe_gluon(lp)]
        self.n_lheparts_final_state[0] = len(final_state_lheparts)
        self.n_lhe_final_state_partons[0] = len(final_state_lhe_quarks) + len(final_state_lhe_gluons)
        self.n_lhe_final_state_quarks[0] = len(final_state_lhe_quarks)
        self.n_lhe_final_state_gluons[0] = len(final_state_lhe_gluons)
        if lheparts:
            self.n_events_with_lhe += 1

        if higgs is not None:
            self.n_higgs_found += 1
            self.n_decay_products[0] = len(daughters)
            self.higgs_pt[0] = higgs.pt
            self.higgs_eta[0] = higgs.eta
            self.higgs_phi[0] = higgs.phi
            self.higgs_mass[0] = higgs.mass

        if higgs is not None and len(daughters) >= 2:
            daughters = sorted(daughters, key=lambda x: x.pt, reverse=True)
            lead = daughters[0]
            sublead = daughters[1]

            self.n_two_body += 1
            self.decay1_pdgId[0] = lead.pdgId
            self.decay2_pdgId[0] = sublead.pdgId
            self.decay_dEta[0] = abs(lead.eta - sublead.eta)
            self.decay_dPhi[0] = abs(delta_phi(lead.phi, sublead.phi))
            self.decay_dR[0] = delta_r(lead.eta, lead.phi, sublead.eta, sublead.phi)
            self.decay_mass[0] = (p4_from_genpart(lead) + p4_from_genpart(sublead)).M()

            self.lead_decay_pt[0] = lead.pt
            self.lead_decay_eta[0] = lead.eta
            self.lead_decay_phi[0] = lead.phi
            self.sublead_decay_pt[0] = sublead.pt
            self.sublead_decay_eta[0] = sublead.eta
            self.sublead_decay_phi[0] = sublead.phi

            lhe_vbf_quarks = find_vbf_lhe_quarks_from_mothers(lheparts, excluded_indices=set())
            if lhe_vbf_quarks:
                self.n_lhe_vbf_quarks[0] = len(lhe_vbf_quarks)
                self.n_events_with_mother_based_vbf_quarks += 1
            else:
                self.used_lhe_vbf_topology_fallback[0] = 1
                self.n_events_with_topology_fallback += 1
                lhe_vbf_pair, lhe_vbf_quarks = choose_vbf_lhe_pair(lheparts, excluded_indices=set())
                self.n_lhe_vbf_quarks[0] = len(lhe_vbf_quarks)
                if lhe_vbf_pair is not None:
                    _, vbf1, _, vbf2, _, _, _, _ = lhe_vbf_pair
                else:
                    vbf1, vbf2 = None, None

            if len(lhe_vbf_quarks) == 2:
                vbf1, vbf2 = lhe_vbf_quarks[0][1], lhe_vbf_quarks[1][1]
            elif lhe_vbf_quarks and "vbf1" in locals() and vbf1 is not None:
                pass
            else:
                vbf1, vbf2 = None, None

            if vbf1 is not None and vbf2 is not None:
                ordered_vbf = sorted([vbf1, vbf2], key=lambda x: x.pt, reverse=True)
                vbf1, vbf2 = ordered_vbf[0], ordered_vbf[1]
                self.has_lhe_vbf_pair[0] = 1
                self.lhe_vbf1_pdgId[0] = vbf1.pdgId
                self.lhe_vbf2_pdgId[0] = vbf2.pdgId
                self.lhe_vbf1_pt[0] = vbf1.pt
                self.lhe_vbf1_eta[0] = vbf1.eta
                self.lhe_vbf1_phi[0] = vbf1.phi
                self.lhe_vbf1_mass[0] = vbf1.mass
                self.lhe_vbf2_pt[0] = vbf2.pt
                self.lhe_vbf2_eta[0] = vbf2.eta
                self.lhe_vbf2_phi[0] = vbf2.phi
                self.lhe_vbf2_mass[0] = vbf2.mass
                self.lhe_vbf_pair_dEta[0] = abs(vbf1.eta - vbf2.eta)
                self.lhe_vbf_pair_dPhi[0] = abs(delta_phi(vbf1.phi, vbf2.phi))
                self.lhe_vbf_pair_dR[0] = delta_r(vbf1.eta, vbf1.phi, vbf2.eta, vbf2.phi)
                self.lhe_vbf_pair_mass[0] = (p4_from_lhepart(vbf1) + p4_from_lhepart(vbf2)).M()
                self.n_events_with_truth_lhe_vbf_pair += 1
        else:
            lhe_vbf_quarks = find_vbf_lhe_quarks_from_mothers(lheparts, excluded_indices=set())
            if lhe_vbf_quarks:
                self.n_lhe_vbf_quarks[0] = len(lhe_vbf_quarks)
                self.n_events_with_mother_based_vbf_quarks += 1
                if len(lhe_vbf_quarks) == 2:
                    vbf1, vbf2 = lhe_vbf_quarks[0][1], lhe_vbf_quarks[1][1]
                else:
                    vbf1, vbf2 = None, None
            else:
                self.used_lhe_vbf_topology_fallback[0] = 1
                self.n_events_with_topology_fallback += 1
                lhe_vbf_pair, lhe_vbf_quarks = choose_vbf_lhe_pair(lheparts, excluded_indices=set())
                self.n_lhe_vbf_quarks[0] = len(lhe_vbf_quarks)
                if lhe_vbf_pair is not None:
                    _, vbf1, _, vbf2, _, _, _, _ = lhe_vbf_pair
                else:
                    vbf1, vbf2 = None, None

            if vbf1 is not None and vbf2 is not None:
                ordered_vbf = sorted([vbf1, vbf2], key=lambda x: x.pt, reverse=True)
                vbf1, vbf2 = ordered_vbf[0], ordered_vbf[1]
                self.has_lhe_vbf_pair[0] = 1
                self.lhe_vbf1_pdgId[0] = vbf1.pdgId
                self.lhe_vbf2_pdgId[0] = vbf2.pdgId
                self.lhe_vbf1_pt[0] = vbf1.pt
                self.lhe_vbf1_eta[0] = vbf1.eta
                self.lhe_vbf1_phi[0] = vbf1.phi
                self.lhe_vbf1_mass[0] = vbf1.mass
                self.lhe_vbf2_pt[0] = vbf2.pt
                self.lhe_vbf2_eta[0] = vbf2.eta
                self.lhe_vbf2_phi[0] = vbf2.phi
                self.lhe_vbf2_mass[0] = vbf2.mass
                self.lhe_vbf_pair_dEta[0] = abs(vbf1.eta - vbf2.eta)
                self.lhe_vbf_pair_dPhi[0] = abs(delta_phi(vbf1.phi, vbf2.phi))
                self.lhe_vbf_pair_dR[0] = delta_r(vbf1.eta, vbf1.phi, vbf2.eta, vbf2.phi)
                self.lhe_vbf_pair_mass[0] = (p4_from_lhepart(vbf1) + p4_from_lhepart(vbf2)).M()
                self.n_events_with_truth_lhe_vbf_pair += 1

        if pair is not None:
            gj1, gj2, deta, dphi, dr, mjj = pair
            self.has_forward_genjet_pair[0] = 1 if vbf_genjet_pair_passes_selection(pair) else 0
            self.forward_genjet1_pt[0] = gj1.pt
            self.forward_genjet1_eta[0] = gj1.eta
            self.forward_genjet1_phi[0] = gj1.phi
            self.forward_genjet1_mass[0] = gj1.mass
            self.forward_genjet2_pt[0] = gj2.pt
            self.forward_genjet2_eta[0] = gj2.eta
            self.forward_genjet2_phi[0] = gj2.phi
            self.forward_genjet2_mass[0] = gj2.mass
            self.forward_genjet_pair_dEta[0] = deta
            self.forward_genjet_pair_dPhi[0] = dphi
            self.forward_genjet_pair_dR[0] = dr
            self.forward_genjet_pair_mass[0] = mjj

            if self.has_lhe_vbf_pair[0]:
                self.n_events_with_genjet_lhe_comparison += 1
                if self.has_forward_genjet_pair[0]:
                    self.n_events_with_tagged_genjet_lhe_comparison += 1

                dr1, dr2 = match_genjet_pair_to_lhe_pair(
                    gj1,
                    gj2,
                    self.lhe_vbf1_eta[0],
                    self.lhe_vbf1_phi[0],
                    self.lhe_vbf2_eta[0],
                    self.lhe_vbf2_phi[0],
                )
                self.genjet_lhe_match_dr1[0] = dr1
                self.genjet_lhe_match_dr2[0] = dr2
                self.genjet_lhe_match_maxdr[0] = max(dr1, dr2)
                self.genjet_lhe_match_sumdr[0] = dr1 + dr2
                self.genjet_lhe_match_n[0] = int(dr1 < args.genJetLheMatchDRMax) + int(dr2 < args.genJetLheMatchDRMax)

                if self.genjet_lhe_match_n[0] == 2:
                    self.has_genjet_lhe_match_pair[0] = 1
                    self.n_events_with_genjet_lhe_match += 1
                    if self.has_forward_genjet_pair[0]:
                        self.n_events_with_tagged_genjet_lhe_match += 1

        self.tree.Fill()
        return True

    def endJob(self):
        print("\n================ ScoutNano Gen Higgs Study ===============")
        print(f"Total events:          {self.total}")
        print(f"Events with gen Higgs: {self.n_higgs_found}")
        print(f"Events with >=2 daughters: {self.n_two_body}")
        print(f"Events with LHEPart info: {self.n_events_with_lhe}")
        print(f"Events with mother-based LHE VBF quarks: {self.n_events_with_mother_based_vbf_quarks}")
        print(f"Events using topology fallback for LHE VBF pair: {self.n_events_with_topology_fallback}")
        print(f"Events with selected GenJet pair: {self.n_events_with_selected_genjet_pair}")
        print(f"Events with truth LHE VBF pair: {self.n_events_with_truth_lhe_vbf_pair}")
        print(f"Events with GenJet/LHE comparison: {self.n_events_with_genjet_lhe_comparison}")
        if self.n_events_with_genjet_lhe_comparison > 0:
            eff = 100.0 * self.n_events_with_genjet_lhe_match / self.n_events_with_genjet_lhe_comparison
            print(f"Pair-match efficiency (all comparable events, DR<{args.genJetLheMatchDRMax:g}): {self.n_events_with_genjet_lhe_match} ({eff:.2f}%)")
        if self.n_events_with_tagged_genjet_lhe_comparison > 0:
            eff_tag = 100.0 * self.n_events_with_tagged_genjet_lhe_match / self.n_events_with_tagged_genjet_lhe_comparison
            print(f"Pair-match efficiency (tagged events only, DR<{args.genJetLheMatchDRMax:g}): {self.n_events_with_tagged_genjet_lhe_match} ({eff_tag:.2f}%)")
        print(f"VBF tag definition:    {describe_vbf_tag_mode()}")
        print("==========================================================\n")

        self.outfile.cd()
        self.tree.Write()
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
        modules=[GenHiggsScoutNanoStudy(output_name)],
        noOut=True,
        maxEntries=args.maxEntries,
        firstEntry=args.firstEntry,
        prefetch=args.prefetch,
    )
    processor.run()
