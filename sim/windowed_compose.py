#!/usr/bin/env python3
"""Host-window composition legality. Thesis #29.

The catalogue is a set of named observation-channel compositions.
Two evidence schedules score the same files. The ledger schedule
carries immunometabolic evidence identifiers. The host schedule
adds the two declared infection windows. A file is accepted when
the structural rule list is empty and every identifier in the
active schedule is hosted.

Theta is a frozen seven-key object. The ranker and the composition
validator do not write it. The script records SHA-256 of the
canonical JSON before either schedule and after the refused
promotion.

Seed 20260929 governs the permutation audit of the rank control.
Composition decisions do not draw random numbers.

Kinetic literals and window endpoints are declared inputs. They
are not fitted. No results file from Thesis #22 or Thesis #23 is
read.
"""

from __future__ import annotations

import copy
import hashlib
import inspect
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from jsonschema import Draft7Validator

ROOT = Path(__file__).resolve().parents[1]
SIM = Path(__file__).resolve().parent
PROF = ROOT / "profiles"
FIG = SIM / "figures"
SCHEMA_PATH = ROOT / "schema" / "windowed_composition.schema.json"
SNAPSHOT_PATH = SIM / "crossref_snapshot.json"

SEED = 20260929
DISEASE = "named-channel-class-unspecified"
DISCLAIMER = (
    "Research object only. Not a medical device, not clinical decision "
    "support, not a diagnostic, not a dose, and not a cure. No document DOI."
)

LEDGER_IDS = ["E_lac", "E_ck", "E_sup"]
WINDOW_IDS = ["W_I", "W_M"]
SCHEDULES = {
    "ledger": list(LEDGER_IDS),
    "host": list(LEDGER_IDS) + list(WINDOW_IDS),
}

# Declared kinetic object. The canonical string is the digest input.
# Roles are the cartoon names from the parent ledger. This script does
# not evaluate that cartoon.
THETA = {
    "r": 0.3,
    "K": 1.2,
    "kappa": 0.5,
    "sigma": 0.12,
    "delta": 0.25,
    "pi": 0.7,
    "lam": 0.55,
}
CANONICAL_THETA = (
    '{"K":1.2,"delta":0.25,"kappa":0.5,"lam":0.55,'
    '"pi":0.7,"r":0.3,"sigma":0.12}'
)

WINDOWS = {
    "W_I": {"label": "infection", "interval": [1.5, 2.5]},
    "W_M": {"label": "marrow_stress", "interval": [10.0, 12.0]},
}

# Names that would smuggle a window, a ledger constant, or an evidence
# identifier into a map coefficient. Hosting those identifiers is a
# different act from naming a coefficient after them.
WINDOW_DERIVED = {
    "k_inf",
    "k_host",
    "sigma_host",
    "phi",
    "L_lo",
    "L_hi",
    "checkpoint",
    "W_I",
    "W_M",
    "E_lac",
    "E_ck",
    "E_sup",
}

SOURCES = {"T08", "T15", "declared_map"}
IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def canonical_theta(theta: dict) -> str:
    return json.dumps(theta, sort_keys=True, separators=(",", ":"))


def digest(theta: dict) -> str:
    return hashlib.sha256(canonical_theta(theta).encode("utf-8")).hexdigest()


def load_snapshot() -> dict:
    snap = json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))
    return snap["works"]


def channel(
    cid: str,
    name: str,
    symbol: str,
    coeff: str,
    observes: str,
    factor: str,
    source: str,
    status: str = "not_in_theta",
) -> dict:
    return {
        "channel_id": cid,
        "name": name,
        "symbol": symbol,
        "map_coefficient": coeff,
        "observes": observes,
        "parameter_status": status,
        "factor_id": factor,
        "source_prior": source,
        "refusal": (
            "This coefficient names an observation map. "
            "It is not a therapeutic parameter."
        ),
    }


CH_A = channel(
    "ch_amyl_map",
    "amylase-shaped map",
    "y_A",
    "alpha",
    "A declared saturating map of a substrate coordinate. Not an inhibition percentage.",
    "F_map",
    "T08",
)
CH_S = channel(
    "ch_spectral",
    "spectral contrast",
    "y_S",
    "beta",
    "A declared spectral contrast. Not a copied reporter loading.",
    "F_multi",
    "T15",
)
CH_C = channel(
    "ch_compartment",
    "compartment indicator",
    "y_C",
    "gamma",
    "A declared compartment indicator. Not a vesicle count.",
    "F_multi",
    "T15",
)


def factor(fid: str, source: str, channel_ids: list[str], disease: str = DISEASE) -> dict:
    return {
        "factor_id": fid,
        "source_prior": source,
        "disease_id": disease,
        "channel_ids": list(channel_ids),
    }


def mechanisms(one_theta: str = "refused") -> list[dict]:
    rows = [
        {
            "id": "M-DISTINCT",
            "statement": (
                "The channel symbols stay distinct under disjoint union. "
                "A host window is an evidence identifier, not a new coefficient."
            ),
            "evidence_class": "in_silico",
            "falsifier": (
                "A preparation series whose loadings differ and which is "
                "nevertheless fitted by one coefficient."
            ),
            "status": "admitted",
        },
        {
            "id": "M-ONE-THETA",
            "statement": "One therapeutic parameter accounts for every channel and for both windows.",
            "evidence_class": "refused_merge",
            "falsifier": (
                "A file whose operator is disjoint union, whose theta is empty, "
                "and whose window identifiers are hosted as names."
            ),
            "status": one_theta,
            "refusal_rule": "one symbol used as every map coefficient, or a window written into theta",
        },
    ]
    return rows


def hypotheses() -> list[dict]:
    return [
        {
            "id": "H-COEFF-OUT",
            "statement": "Map coefficients alpha, beta, and gamma stay out of theta.",
            "falsifier": "A legal file that lists alpha, beta, or gamma inside theta.",
            "parameter_status": "forbidden_to_enter_theta",
        },
        {
            "id": "H-WINDOW-NAME",
            "statement": (
                "Hosting W_I or W_M is membership of an identifier. "
                "It does not define k_inf or k_host."
            ),
            "falsifier": "A file accepted while theta contains k_inf or k_host.",
            "parameter_status": "forbidden_to_enter_theta",
        },
    ]


def non_parameters() -> list[dict]:
    rows = [
        ("merged_therapeutic_theta", "One symbol used as every map coefficient."),
        ("composition_quotient", "A join that identifies distinct coefficients."),
        ("product_of_map_coefficients", "A product written as a parameter."),
        ("arithmetic_mean_of_channels", "A mean written as a parameter."),
        ("channel_promoted_into_theta", "One channel promoted by the act of joining."),
        ("cross_disease_join", "A union across two disease-class identifiers."),
        ("self_composition", "A repeated channel identifier."),
        ("assay_scalar_as_treatment", "An assay-shaped coefficient copied in as a treatment."),
        ("host_window_endpoint", "An endpoint of W_I or W_M copied into theta."),
        ("window_rate_k_inf", "k_inf = 1/2, from the infection-window midpoint, copied into theta."),
        ("window_rate_k_host", "k_host = 1/9, from the gap between window midpoints, copied into theta."),
        ("ledger_constant_phi", "The checkpoint scale phi copied into theta."),
    ]
    return [{"id": i, "statement": s} for i, s in rows]


def citations() -> list[dict]:
    return [
        {
            "id": "c_may",
            "vancouver": (
                "May RM. Uses and abuses of mathematics in biology. "
                "Science. 2004;303(5659):790-793."
            ),
            "doi": "10.1126/science.1094442",
        },
        {
            "id": "c_saltelli",
            "vancouver": (
                "Saltelli A, Bammer G, Bruno I, Charters E, Di Fiore M, Didier E, et al. "
                "Five ways to ensure that models serve society: a manifesto. "
                "Nature. 2020;582(7813):482-484."
            ),
            "doi": "10.1038/d41586-020-01812-9",
        },
        {
            "id": "c_fair",
            "vancouver": (
                "Wilkinson MD, Dumontier M, Aalbersberg IJ, Appleton G, Axton M, Baak A, et al. "
                "The FAIR Guiding Principles for scientific data management and stewardship. "
                "Sci Data. 2016;3:160018."
            ),
            "doi": "10.1038/sdata.2016.18",
        },
        {
            "id": "c_fusion",
            "vancouver": (
                "Bareinboim E, Pearl J. Causal inference and the data-fusion problem. "
                "Proc Natl Acad Sci U S A. 2016;113(27):7345-7352."
            ),
            "doi": "10.1073/pnas.1510507113",
        },
    ]


def shell(
    profile_id: str,
    channels: list[dict],
    factors: list[dict],
    hosted: list[str],
    operator: str = "disjoint_union",
    quotient=None,
    combiner=None,
    theta: list[str] | None = None,
    one_theta: str = "refused",
    disease: str = DISEASE,
) -> dict:
    return {
        "schema_version": "1.3.0-windowed-compose",
        "profile_id": profile_id,
        "extends": {
            "object": "DiseaseProfile",
            "schema_version": "1.0.0",
            "source": "https://github.com/cloudynirvana/thesis-03-disease-profile",
        },
        "disease_id": disease,
        "disease_label": "Named observation-channel class, unspecified. Not a patient.",
        "created_at": "2026-09-21",
        "asker_role": "methods_researcher",
        "answers": [
            {
                "id": "A1",
                "question": "Who asks?",
                "label": "A methods researcher comparing two evidence schedules.",
            },
            {
                "id": "A2",
                "question": "Which class?",
                "label": DISEASE,
            },
            {
                "id": "A3",
                "question": "Which regime?",
                "label": "In silico. Declared records. No cohort.",
            },
            {
                "id": "A4",
                "question": "Where is the portfolio stuck?",
                "label": (
                    "The channels and the host windows both exist. "
                    "The joint legality of the composition under those windows does not."
                ),
            },
        ],
        "observables": [
            {
                "id": "obs_channels",
                "statement": "Named maps remain separate symbols under a legal join.",
                "citation_ids": ["c_may", "c_fusion"],
            },
            {
                "id": "obs_windows",
                "statement": "A host window is an evidence identifier with a declared interval.",
                "citation_ids": ["c_saltelli", "c_fair"],
            },
        ],
        "observation_channels": channels,
        "composition": {
            "operator": operator,
            "quotient": quotient,
            "combiner": combiner,
            "factors": factors,
            "hosted_evidence": list(hosted),
        },
        "candidate_mechanisms": mechanisms(one_theta),
        "non_parameters": non_parameters(),
        "theta": [] if theta is None else list(theta),
        "admitted_hypotheses": hypotheses(),
        "citations": citations(),
        "disclaimer": DISCLAIMER,
    }


def clone(ch: dict, **over) -> dict:
    row = copy.deepcopy(ch)
    row.update(over)
    return row


def build_catalogue() -> list[dict]:
    """Sixteen files. Illegal operator files host every evidence identifier.

    That choice is deliberate. A later schedule cannot be what refuses
    them. Files that are structurally legal host only the evidence named
    in the row, so the schedule clause is the thing that can move.
    """
    a, s, c = CH_A, CH_S, CH_C
    both = LEDGER_IDS + WINDOW_IDS
    ledger = list(LEDGER_IDS)
    partial = LEDGER_IDS + ["W_I"]
    windows = list(WINDOW_IDS)

    a_map = clone(a, factor_id="F_map")
    s_multi = clone(s, factor_id="F_multi")
    c_multi = clone(c, factor_id="F_multi")
    s_spec = clone(s, factor_id="F_spec")
    c_comp = clone(c, factor_id="F_comp")

    rows = []

    def add(doc: dict, filename: str) -> None:
        doc["_filename"] = filename
        rows.append(doc)

    add(
        shell(
            "bundle",
            [a_map, s_multi, c_multi],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            ledger,
        ),
        "ledger_only_bundle.profile.yaml",
    )
    add(
        shell(
            "singletons",
            [clone(a, factor_id="F_map"), s_spec, c_comp],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_spec", "T15", ["ch_spectral"]),
                factor("F_comp", "T15", ["ch_compartment"]),
            ],
            ledger,
        ),
        "ledger_only_singletons.profile.yaml",
    )
    add(
        shell(
            "binary",
            [s_spec, c_comp],
            [
                factor("F_spec", "T15", ["ch_spectral"]),
                factor("F_comp", "T15", ["ch_compartment"]),
            ],
            ledger,
        ),
        "ledger_only_binary.profile.yaml",
    )
    add(
        shell(
            "windowed",
            [a_map, s_multi, c_multi],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            both,
        ),
        "both_schedules_windowed.profile.yaml",
    )
    add(
        shell(
            "partial_window",
            [a_map, s_multi, c_multi],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            partial,
        ),
        "ledger_only_partial_window.profile.yaml",
    )
    add(
        shell(
            "windows_only",
            [a_map, s_multi, c_multi],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            windows,
        ),
        "refuse_windows_only.profile.yaml",
    )
    add(
        shell(
            "empty_hosts",
            [a_map, s_multi, c_multi],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            [],
        ),
        "refuse_empty_hosts.profile.yaml",
    )

    identified = [
        clone(a, symbol="theta", map_coefficient="theta", parameter_status="in_theta"),
        clone(s, symbol="theta", map_coefficient="theta", parameter_status="in_theta"),
        clone(c, symbol="theta", map_coefficient="theta", parameter_status="in_theta"),
    ]
    add(
        shell(
            "identify",
            identified,
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            both,
            operator="identify",
            quotient="theta",
            combiner="identify",
            theta=["theta"],
        ),
        "refuse_identify.profile.yaml",
    )
    add(
        shell(
            "product",
            [a_map, s_multi, c_multi],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            both,
            operator="product",
            combiner="product",
            theta=["pi"],
        ),
        "refuse_product.profile.yaml",
    )
    add(
        shell(
            "mean",
            [a_map, s_multi, c_multi],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            both,
            operator="arithmetic_mean",
            combiner="arithmetic_mean",
            theta=["theta_mean"],
        ),
        "refuse_mean.profile.yaml",
    )
    promoted = [
        clone(a, parameter_status="in_theta"),
        s_multi,
        c_multi,
    ]
    add(
        shell(
            "promote_alpha",
            promoted,
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            both,
            theta=["alpha"],
        ),
        "refuse_promote_alpha.profile.yaml",
    )
    s_left = clone(s, factor_id="F_left")
    s_right = clone(s, factor_id="F_right")
    c_right = clone(c, factor_id="F_right")
    add(
        shell(
            "repeat",
            [s_left, s_right, c_right],
            [
                factor("F_left", "T15", ["ch_spectral"]),
                factor("F_right", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            both,
        ),
        "refuse_repeat.profile.yaml",
    )
    add(
        shell(
            "cross_class",
            [a_map, s_multi, c_multi],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"], "other-channel-class"),
            ],
            both,
        ),
        "refuse_cross_class.profile.yaml",
    )
    rate = clone(a, map_coefficient="k_inf")
    add(
        shell(
            "window_coeff",
            [rate, s_multi, c_multi],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            both,
        ),
        "refuse_window_coeff.profile.yaml",
    )
    add(
        shell(
            "window_theta",
            [a_map, s_multi, c_multi],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            both,
            theta=["k_inf", "k_host"],
        ),
        "refuse_window_theta.profile.yaml",
    )
    add(
        shell(
            "one_theta",
            [a_map, s_multi, c_multi],
            [
                factor("F_map", "T08", ["ch_amyl_map"]),
                factor("F_multi", "T15", ["ch_spectral", "ch_compartment"]),
            ],
            both,
            one_theta="admitted",
        ),
        "refuse_one_theta.profile.yaml",
    )
    return rows


def structural_codes(doc: dict, snapshot: dict) -> list[str]:
    """Schedule-blind rule list. Order is fixed so the tables are stable."""
    codes: list[str] = []
    comp = doc["composition"]
    channels = doc["observation_channels"]
    if comp.get("operator") != "disjoint_union":
        codes.append("OPERATOR")
    if comp.get("quotient") is not None:
        codes.append("QUOTIENT")
    if comp.get("combiner") is not None:
        codes.append("COMBINER")
    if doc.get("theta") != []:
        codes.append("THETA_NONEMPTY")
    if any(ch.get("parameter_status") != "not_in_theta" for ch in channels):
        codes.append("CHANNEL_IN_THETA")
    ids = [ch.get("channel_id") for ch in channels]
    if len(ids) != len(set(ids)):
        codes.append("REPEATED_CHANNEL")
    flat: list[str] = []
    for fac in comp.get("factors", []):
        flat.extend(fac.get("channel_ids", []))
    if len(flat) != len(set(flat)) or set(flat) != set(ids) or len(flat) != len(ids):
        codes.append("PARTITION")
    if len(channels) < 2:
        codes.append("TOO_FEW_CHANNELS")
    for fac in comp.get("factors", []):
        if fac.get("source_prior") not in SOURCES:
            codes.append("SOURCE")
            break
    if any(fac.get("disease_id") != doc.get("disease_id") for fac in comp.get("factors", [])):
        codes.append("CROSS_CLASS")
    disease = str(doc.get("disease_id", "")).lower()
    if any(tok in disease for tok in ("patient", "mrn", "subject")):
        codes.append("DISEASE_TOKEN")
    for ch in channels:
        if not IDENT.match(str(ch.get("symbol", ""))) or not IDENT.match(
            str(ch.get("map_coefficient", ""))
        ):
            codes.append("COEFF_EXPRESSION")
            break
    derived = []
    for ch in channels:
        if ch.get("map_coefficient") in WINDOW_DERIVED or ch.get("symbol") in WINDOW_DERIVED:
            derived.append(ch.get("map_coefficient"))
    if derived:
        codes.append("WINDOW_COEFF")
    for mech in doc.get("candidate_mechanisms", []):
        if mech.get("id") == "M-ONE-THETA" and mech.get("status") != "refused":
            codes.append("ONE_THETA_ADMITTED")
        if mech.get("status") == "admitted" and not str(mech.get("falsifier", "")).strip():
            codes.append("EMPTY_FALSIFIER")
    for hyp in doc.get("admitted_hypotheses", []):
        if hyp.get("parameter_status") != "forbidden_to_enter_theta":
            codes.append("HYPOTHESIS_IN_THETA")
            break
    if doc.get("disclaimer") != DISCLAIMER:
        codes.append("DISCLAIMER")
    for cit in doc.get("citations", []):
        doi = cit.get("doi")
        if doi is not None and doi not in snapshot:
            codes.append("DOI_NOT_IN_SNAPSHOT")
            break
    return codes


def schedule_codes(doc: dict, evidence: list[str]) -> list[dict]:
    hosted = set(doc["composition"].get("hosted_evidence", []))
    missing = [item for item in evidence if item not in hosted]
    if not missing:
        return []
    return [{"code": "UNHOSTED", "missing": missing}]


def code_labels(structural: list[str], scheduled: list[dict]) -> list[str]:
    labels = list(structural)
    for row in scheduled:
        if row["code"] == "UNHOSTED":
            labels.append("UNHOSTED:" + ",".join(row["missing"]))
        else:
            labels.append(row["code"])
    return labels


def schema_report(doc: dict, validator: Draft7Validator) -> dict:
    errors = sorted(validator.iter_errors(doc), key=lambda e: list(e.path))
    paths = ["/".join(str(p) for p in err.path) or "(root)" for err in errors]
    return {
        "valid": not errors,
        "n_errors": len(errors),
        "paths": paths[:8],
    }


def public_doc(doc: dict) -> dict:
    out = copy.deepcopy(doc)
    out.pop("_filename", None)
    return out


def associativity(docs: dict[str, dict]) -> dict:
    bundle = docs["bundle"]
    single = docs["singletons"]
    b_ids = [ch["channel_id"] for ch in bundle["observation_channels"]]
    s_ids = [ch["channel_id"] for ch in single["observation_channels"]]
    # (spectral ⊔ compartment) ⊔ amylase, as a set, against the bundle.
    left = {"ch_spectral", "ch_compartment"}
    reassoc = left | {"ch_amyl_map"}
    return {
        "bundle_id_list": b_ids,
        "singleton_id_list": s_ids,
        "lists_equal": b_ids == s_ids,
        "reassociated_set_equals_bundle": reassoc == set(b_ids),
    }


def guard_theta(proposal: dict, theta: dict) -> dict:
    """Refuse a proposal that is not exactly the seven legal pairs.

    The ranker is not called. Theta is returned unchanged.
    """
    legal = set(theta)
    extra = sorted(set(proposal) - legal)
    missing = sorted(legal - set(proposal))
    changed = sorted(k for k in legal & set(proposal) if proposal[k] != theta[k])
    if extra or missing or changed:
        return {
            "status": "REFUSED",
            "code": "ILLEGAL_PROMOTION",
            "extra_keys": extra,
            "missing_keys": missing,
            "changed_keys": changed,
            "theta_returned": canonical_theta(theta),
        }
    return {
        "status": "NAMES_UNCHANGED",
        "code": "NAMES_UNCHANGED",
        "extra_keys": [],
        "missing_keys": [],
        "changed_keys": [],
        "theta_returned": canonical_theta(theta),
    }


def promotion_proposal() -> dict:
    """The map Thesis #23 refuses, rebuilt from the declared literals.

    Midpoints are t_I = 2 and t_M = 11. k_inf = 1/t_I. k_host = 1/(t_M - t_I).
    sigma is overwritten with the host supply budget. Ledger constants are added.
    """
    t_i = 0.5 * (WINDOWS["W_I"]["interval"][0] + WINDOWS["W_I"]["interval"][1])
    t_m = 0.5 * (WINDOWS["W_M"]["interval"][0] + WINDOWS["W_M"]["interval"][1])
    proposal = dict(THETA)
    proposal["sigma"] = 0.20
    proposal["phi"] = 0.55
    proposal["L_lo"] = 0
    proposal["L_hi"] = 0.45
    proposal["sigma_host"] = 0.20
    proposal["checkpoint"] = "present"
    proposal["k_inf"] = 1.0 / t_i
    proposal["k_host"] = 1.0 / (t_m - t_i)
    return proposal


# --- rank control: declared structure records, recomputed here ------------

STRUCTURES = [
    {
        "id": "U0",
        "name": "open_baseline",
        "immuno": [],
        "windows": [],
        "inside_band": True,
        "within_budget": True,
    },
    {
        "id": "U1",
        "name": "lactate_masked",
        "immuno": ["E_lac", "E_sup"],
        "windows": [],
        "inside_band": True,
        "within_budget": True,
    },
    {
        "id": "U2",
        "name": "checkpoint_scaled",
        "immuno": ["E_lac", "E_ck", "E_sup"],
        "windows": [],
        "inside_band": True,
        "within_budget": True,
    },
    {
        "id": "U3",
        "name": "windowed_checkpoint",
        "immuno": ["E_lac", "E_ck", "E_sup"],
        "windows": ["W_I", "W_M"],
        "inside_band": True,
        "within_budget": True,
    },
    {
        "id": "U4",
        "name": "band_bypass",
        "immuno": ["E_ck"],
        "windows": ["W_I", "W_M"],
        "inside_band": False,
        "within_budget": True,
    },
]


def enrich(row: dict, evidence: list[str]) -> dict:
    immuno_ev = [i for i in evidence if i in LEDGER_IDS]
    window_ev = [i for i in evidence if i in WINDOW_IDS]
    admissible = bool(row["inside_band"] and row["within_budget"])
    # Checkpoint flag is true for every record in this catalogue:
    # none claims a checkpoint slot while the label is absent.
    out = dict(row)
    out["admissible"] = admissible
    out["unhosted_windows"] = [i for i in window_ev if i not in row["windows"]]
    out["unhosted_immuno"] = [i for i in immuno_ev if i not in row["immuno"]]
    out["slots"] = len(row["immuno"]) + len(row["windows"])
    return out


def key_primary(row: dict) -> tuple:
    return (
        0 if row["admissible"] else 1,
        len(row["unhosted_windows"]),
        len(row["unhosted_immuno"]),
        row["slots"],
        row["id"],
    )


def key_immuno_terms(row: dict) -> tuple:
    return (
        0 if row["admissible"] else 1,
        len(row["unhosted_immuno"]),
        row["slots"],
        row["id"],
    )


def key_tie_break(row: dict) -> tuple:
    # Slot count is read before unhosted windows.
    return (
        0 if row["admissible"] else 1,
        row["slots"],
        len(row["unhosted_windows"]),
        len(row["unhosted_immuno"]),
        row["id"],
    )


def key_window_only(row: dict) -> tuple:
    return (
        len(row["unhosted_windows"]),
        row["slots"],
        row["id"],
    )


def rank_records(records: list[dict], evidence: list[str], key) -> list[dict]:
    rows = [enrich(r, evidence) for r in records]
    rows.sort(key=key)
    ranked = []
    for i, row in enumerate(rows, start=1):
        ranked.append(
            {
                "rank": i,
                "id": row["id"],
                "name": row["name"],
                "admissible": row["admissible"],
                "unhosted_windows": len(row["unhosted_windows"]),
                "unhosted_immuno": len(row["unhosted_immuno"]),
                "slots": row["slots"],
            }
        )
    return ranked


def permutation_audit(evidence: list[str], n: int, rng: np.random.Generator) -> dict:
    base = [r["id"] for r in rank_records(STRUCTURES, evidence, key_primary)]
    discordant = 0
    order = np.arange(len(STRUCTURES))
    for _ in range(n):
        rng.shuffle(order)
        shuffled = [STRUCTURES[i] for i in order]
        got = [r["id"] for r in rank_records(shuffled, evidence, key_primary)]
        if got != base:
            discordant += 1
    return {"n": n, "discordant": discordant, "order": base}


def draw_matrix(files: list[dict], path: Path) -> None:
    labels = [f["profile_id"] for f in files]
    cols = ["Ledger rule", "Host rule", "Schema"]
    grid = []
    for f in files:
        grid.append(
            [
                1 if f["accept"]["ledger"] else 0,
                1 if f["accept"]["host"] else 0,
                1 if f["schema"]["valid"] else 0,
            ]
        )
    arr = np.array(grid, dtype=float)
    fig, ax = plt.subplots(figsize=(7.4, 8.6))
    cmap = matplotlib.colors.ListedColormap(["#9b2c2c", "#1b7f4e"])
    ax.imshow(arr, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(3), cols)
    ax.set_yticks(range(len(labels)), labels)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            ax.text(
                j,
                i,
                "accept" if arr[i, j] == 1 else "refuse",
                ha="center",
                va="center",
                color="white",
                fontsize=8,
            )
    ax.set_title("Composition accept and refuse sets")
    ax.tick_params(length=0)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def draw_movement(movement: dict, path: Path) -> None:
    columns = [
        ("Accepted on both", movement["stay_accept"]),
        ("Leave when windows attach", movement["leave_on_host"]),
        ("Enter when windows attach", movement["enter_on_host"]),
        ("Refused on both", movement["stay_refuse"]),
    ]
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 14)
    ax.axis("off")
    ax.set_title("What the host schedule does to the rule-list accept set")
    for i, (title, members) in enumerate(columns):
        ax.text(i + 0.5, 13.2, title, ha="center", va="top", fontsize=9, fontweight="bold")
        ax.text(
            i + 0.5,
            12.4,
            f"n = {len(members)}",
            ha="center",
            va="top",
            fontsize=8,
            color="#333",
        )
        body = members if members else ["(empty)"]
        for k, name in enumerate(body):
            ax.text(i + 0.5, 11.4 - 0.62 * k, name, ha="center", va="top", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def draw_ranks(ledger_order: list[dict], host_order: list[dict], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    rank_l = {r["id"]: r["rank"] for r in ledger_order}
    rank_h = {r["id"]: r["rank"] for r in host_order}
    xs = [0, 1]
    colors = {
        "U0": "#4c4c4c",
        "U1": "#1f4e79",
        "U2": "#b85c38",
        "U3": "#1b7f4e",
        "U4": "#6b3fa0",
    }
    for sid in ["U0", "U1", "U4", "U2", "U3"]:
        ys = [rank_l[sid], rank_h[sid]]
        ax.plot(xs, ys, marker="o", lw=1.8, color=colors[sid], label=sid, zorder=3)
        ax.text(-0.08, ys[0], sid, ha="right", va="center", fontsize=8, color=colors[sid])
        ax.text(1.08, ys[1], sid, ha="left", va="center", fontsize=8, color=colors[sid])
    ax.annotate(
        "U2 and U3 exchange rank 1",
        xy=(0.5, 1.5),
        xytext=(0.50, 3.15),
        ha="center",
        fontsize=8,
        arrowprops={"arrowstyle": "-", "color": "#444"},
        color="#222",
    )
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(5.5, 0.6)
    ax.set_xticks(xs, ["Ledger schedule", "Host schedule"])
    ax.set_ylabel("Rank (1 is preferred)")
    ax.set_title("Rank control: the windows still move first place")
    ax.legend(frameon=False, ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.28))
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    assert "theta" not in inspect.signature(rank_records).parameters
    assert canonical_theta(THETA) == CANONICAL_THETA
    theta = dict(THETA)
    d0 = digest(theta)

    snapshot = load_snapshot()
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    validator = Draft7Validator(schema)

    catalogue = build_catalogue()
    by_id = {doc["profile_id"]: doc for doc in catalogue}

    files_out = []
    for doc in catalogue:
        structural = structural_codes(doc, snapshot)
        per_schedule = {}
        accept = {}
        for name, evidence in SCHEDULES.items():
            scheduled = schedule_codes(doc, evidence)
            labels = code_labels(structural, scheduled)
            per_schedule[name] = {
                "codes": labels,
                "accept": not labels,
                "missing": next(
                    (row["missing"] for row in scheduled if row["code"] == "UNHOSTED"),
                    [],
                ),
            }
            accept[name] = per_schedule[name]["accept"]
        report = schema_report(public_doc(doc), validator)
        # Schema does not receive the schedule. Score it again under the
        # other name and require the same boolean.
        report_again = schema_report(public_doc(doc), validator)
        if report["valid"] != report_again["valid"]:
            raise SystemExit("schema result depended on a repeated call")
        files_out.append(
            {
                "profile_id": doc["profile_id"],
                "filename": doc["_filename"],
                "operator": doc["composition"]["operator"],
                "theta": list(doc["theta"]),
                "hosted_evidence": list(doc["composition"]["hosted_evidence"]),
                "channel_ids": [ch["channel_id"] for ch in doc["observation_channels"]],
                "map_coefficients": [
                    ch["map_coefficient"] for ch in doc["observation_channels"]
                ],
                "structural_codes": structural,
                "schedules": per_schedule,
                "accept": accept,
                "schema": report,
            }
        )

    d1 = digest(theta)
    accept_ledger = [f["profile_id"] for f in files_out if f["accept"]["ledger"]]
    accept_host = [f["profile_id"] for f in files_out if f["accept"]["host"]]
    structural_accept = [f["profile_id"] for f in files_out if not f["structural_codes"]]
    schema_accept = [f["profile_id"] for f in files_out if f["schema"]["valid"]]
    movement = {
        "stay_accept": [i for i in accept_ledger if i in accept_host],
        "leave_on_host": [i for i in accept_ledger if i not in accept_host],
        "enter_on_host": [i for i in accept_host if i not in accept_ledger],
        "stay_refuse": [
            f["profile_id"]
            for f in files_out
            if not f["accept"]["ledger"] and not f["accept"]["host"]
        ],
    }

    # In-memory DOI that matches the schema pattern and misses the snapshot.
    injected = public_doc(by_id["windowed"])
    injected["citations"] = copy.deepcopy(injected["citations"])
    injected["citations"][0]["doi"] = "10.1234/not-a-real-record"
    injected_schema = schema_report(injected, validator)
    injected_codes = structural_codes(injected, snapshot)

    proposal = promotion_proposal()
    refused = guard_theta(proposal, theta)
    d2 = digest(theta)
    unchanged = guard_theta(dict(theta), theta)
    d3 = digest(theta)

    rng = np.random.default_rng(SEED)
    audit_ledger = permutation_audit(SCHEDULES["ledger"], 200, rng)
    audit_host = permutation_audit(SCHEDULES["host"], 200, rng)
    rank_ledger = rank_records(STRUCTURES, SCHEDULES["ledger"], key_primary)
    rank_host = rank_records(STRUCTURES, SCHEDULES["host"], key_primary)
    rank_immuno_on_host = rank_records(STRUCTURES, SCHEDULES["host"], key_immuno_terms)
    rank_tie = rank_records(STRUCTURES, SCHEDULES["host"], key_tie_break)
    rank_window_only = rank_records(STRUCTURES, SCHEDULES["host"], key_window_only)

    assoc = associativity(by_id)

    digests = {
        "start": d0,
        "after_catalogue": d1,
        "after_refused_promotion": d2,
        "after_names_unchanged": d3,
    }
    if len(set(digests.values())) != 1:
        raise SystemExit(f"theta digest moved: {digests}")
    if d0 != hashlib.sha256(CANONICAL_THETA.encode("utf-8")).hexdigest():
        raise SystemExit("digest is not the hash of the declared canonical string")
    if accept_ledger == accept_host:
        raise SystemExit("accept sets did not change")
    if movement["enter_on_host"]:
        raise SystemExit("a refused file entered under the host schedule")
    if schema_accept != [f["profile_id"] for f in files_out if f["schema"]["valid"]]:
        raise SystemExit("schema list drifted")
    if not injected_schema["valid"]:
        raise SystemExit("schema rejected the unsnapshotted DOI; the pattern should allow it")
    if "DOI_NOT_IN_SNAPSHOT" not in injected_codes:
        raise SystemExit("rule list accepted an unsnapshotted DOI")
    if refused["status"] != "REFUSED" or unchanged["status"] != "NAMES_UNCHANGED":
        raise SystemExit("theta guard failed")
    if audit_ledger["discordant"] or audit_host["discordant"]:
        raise SystemExit("permutation audit discordant")
    if not assoc["lists_equal"] or not assoc["reassociated_set_equals_bundle"]:
        raise SystemExit("associativity failed")
    for f in files_out:
        if f["accept"]["ledger"] or f["accept"]["host"]:
            if f["theta"] != []:
                raise SystemExit(f"accepted file with theta: {f['profile_id']}")
            if any(c in WINDOW_DERIVED for c in f["map_coefficients"]):
                raise SystemExit(f"accepted window coefficient: {f['profile_id']}")

    PROF.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    for doc in catalogue:
        target = PROF / doc["_filename"]
        target.write_text(
            yaml.safe_dump(public_doc(doc), sort_keys=False, allow_unicode=True),
            encoding="utf-8",
        )

    draw_matrix(files_out, FIG / "fig_accept_refuse.png")
    draw_movement(movement, FIG / "fig_set_movement.png")
    draw_ranks(rank_ledger, rank_host, FIG / "fig_rank_control.png")

    results = {
        "thesis": 29,
        "seed": SEED,
        "seed_governs": "permutation audit of the rank control only",
        "inputs_not_results": {
            "theta_canonical": CANONICAL_THETA,
            "windows": WINDOWS,
            "ledger_ids": LEDGER_IDS,
            "note": (
                "Kinetic literals and window endpoints are declared inputs. "
                "No results file from Thesis #22 or Thesis #23 was read. "
                "Fisher diagonals from Thesis #22 were not recomputed."
            ),
        },
        "theta_sha256": d0,
        "digest_checkpoints": digests,
        "digest_unchanged": True,
        "n_files": len(files_out),
        "accept": {"ledger": accept_ledger, "host": accept_host},
        "accept_counts": {"ledger": len(accept_ledger), "host": len(accept_host)},
        "structural_accept": structural_accept,
        "structural_accept_count": len(structural_accept),
        "schema_accept": schema_accept,
        "schema_accept_count": len(schema_accept),
        "schema_schedule_blind": True,
        "movement": movement,
        "files": files_out,
        "associativity": assoc,
        "unsnapshotted_doi": {
            "doi": "10.1234/not-a-real-record",
            "written_to_disk": False,
            "schema_valid": injected_schema["valid"],
            "rule_codes": injected_codes,
        },
        "promotion": {
            "status": refused["status"],
            "code": refused["code"],
            "extra_keys": refused["extra_keys"],
            "changed_keys": refused["changed_keys"],
            "k_inf": proposal["k_inf"],
            "k_host": proposal["k_host"],
            "k_host_exact": "1/9",
            "sigma_proposed": proposal["sigma"],
            "phi": proposal["phi"],
            "stored_as": "proposal_not_a_result",
        },
        "names_unchanged": unchanged,
        "rank_control": {
            "primary_ledger": rank_ledger,
            "primary_host": rank_host,
            "immuno_terms_on_host": rank_immuno_on_host,
            "tie_break_on_host": rank_tie,
            "window_only_on_host": rank_window_only,
            "permutation": {"ledger": audit_ledger, "host": audit_host},
            "note": (
                "Recomputed from the declared structure records. "
                "The ranker does not take theta as an argument."
            ),
        },
    }
    (SIM / "results.json").write_text(
        json.dumps(results, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    print(f"theta_sha256 {d0}")
    print(f"accept ledger ({len(accept_ledger)}): {accept_ledger}")
    print(f"accept host   ({len(accept_host)}): {accept_host}")
    print(f"schema accept ({len(schema_accept)}): {schema_accept}")
    print(f"leave {movement['leave_on_host']}")
    print(f"enter {movement['enter_on_host']}")
    print("primary ledger", [r["id"] for r in rank_ledger])
    print("primary host  ", [r["id"] for r in rank_host])
    print("immuno terms  ", [r["id"] for r in rank_immuno_on_host])
    print("tie break     ", [r["id"] for r in rank_tie])
    print("window only   ", [r["id"] for r in rank_window_only])
    print("promotion", refused["status"], refused["extra_keys"], refused["changed_keys"])
    print(f"wrote {SIM / 'results.json'}")


if __name__ == "__main__":
    main()
