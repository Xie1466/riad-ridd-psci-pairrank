"""Frozen Literature-prior and PSCI-PairRank reproduction pipeline."""

from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

from common.model import ArchitecturePrior, FoldPCA, RidgeLogistic, Standardizer, correlation_metrics, sigmoid

SYSTEMS = ("D-PA", "GI-DAE", "GP-PGM", "RHLA-RHLB")
FEATURES = ("D", "C", "F1", "F4", "F5", "L")
MODELS = ("M0_GLOBAL_PSCI_PAIRWISE", "M1_LITERATURE_PRIOR", "M2_PRIOR_PLUS_PSCI", "M3_SYSTEM_CONDITIONED_PSCI")


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields=None):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def topology_code(construct_id):
    match = re.search(r"T0([1-4])$", construct_id)
    if not match:
        raise ValueError(f"unrecognized topology ID: {construct_id}")
    return f"T0{match.group(1)}"


def internal_architecture(construct_id):
    code = topology_code(construct_id)
    terminus = "C" if code in {"T01", "T02"} else "N"
    recruitment = ("RIAD", "RIDD") if code in {"T01", "T03"} else ("RIDD", "RIAD")
    return {"construct_id": construct_id, "attachment_A": terminus, "attachment_B": terminus,
            "N/C_topology": f"{terminus}/{terminus}", "recruitment_A": recruitment[0],
            "recruitment_B": recruitment[1], "scaffold_platform": "RIAD_RIDD",
            "linker_length_A": "15", "linker_length_B": "15"}


def load_literature(data_root):
    root = Path(data_root) / "literature_prior"
    constructs = {r["construct_id"]: r for r in read_csv(root / "construct_architecture_91.csv")}
    all_pairs, valid = read_csv(root / "validated_pairs_381.csv"), []
    for source in all_pairs:
        if source["tie"].lower() == "true" or source["construct_i"] not in constructs or source["construct_j"] not in constructs:
            continue
        row = dict(source); direction = row["rank_direction"]
        if direction in {"i>j", "I>J"}: direction = "I_GT_J"
        if direction in {"i<j", "I<J"}: direction = "J_GT_I"
        if direction not in {"I_GT_J", "J_GT_I"}: continue
        row["rank_direction"] = direction; valid.append(row)
    if len(all_pairs) != 381 or len(valid) != 377:
        raise RuntimeError(f"literature freeze mismatch: {len(all_pairs)} total/{len(valid)} eligible")
    return constructs, valid


def load_internal(data_root):
    root = Path(data_root) / "internal16"
    features = {r["topology_id"]: r for r in read_csv(root / "INTERNAL16_PSCI_FINAL_VALUES.csv")}
    phenotypes = {r["topology_id"]: r for r in read_csv(root / "phenotype_rank_freeze.csv")}
    contexts = {}
    for r in read_csv(root / "system_context.csv"):
        contexts[r["system"]] = np.asarray([float(r[x]) for x in (
            "mean_log_enzyme_length", "absolute_log_length_ratio", "amino_acid_composition_cosine",
            "mean_hydrophobic_fraction", "Q_context_known_fraction", "Q_context_mean_log_oligomer_size")])
    constructs = {}
    for cid, row in features.items():
        values = np.asarray([float(row[f]) for f in FEATURES])
        constructs[cid] = {"construct_id": cid, "system": row["system"], "features": values,
                           "phenotype": float(phenotypes[cid]["rank_score_higher_is_better"]),
                           "within_system_rank": int(phenotypes[cid]["within_system_rank"]),
                           "architecture": internal_architecture(cid)}
    pairs = read_csv(root / "pair_definitions_24.csv")
    if len(constructs) != 16 or len(pairs) != 24 or set(contexts) != set(SYSTEMS):
        raise RuntimeError("Internal16 shape gate failed")
    return constructs, pairs, contexts


def exclusions(data_root):
    rows = read_csv(Path(data_root) / "literature_prior/overlap_exclusion_audit.csv")
    result = {s: set() for s in SYSTEMS}
    for r in rows:
        if r["exclude_from_fold_prior_training"].lower() == "true":
            result[r["held_out_internal_system"]].add(r["literature_system"])
    return result


def pair_matrix(rows):
    return np.asarray([[float(r[f"delta_{f}"]) for f in FEATURES] for r in rows], float)


def fit_scaled(train_x, y, test_x):
    scaler = Standardizer.fit_scale_only(train_x)
    model = RidgeLogistic.fit(scaler.transform(train_x), y, penalty=1.0)
    return model.predict_proba(scaler.transform(test_x))


def fit_priors(data_root, output_dir):
    constructs, pairs, _ = load_internal(data_root)
    lit_constructs, lit_pairs = load_literature(data_root)
    excluded = exclusions(data_root)
    out, artifacts = [], {}
    for heldout in SYSTEMS:
        allowed = [r for r in lit_pairs if r["system_id"] not in excluded[heldout]]
        prior = ArchitecturePrior.fit(lit_constructs, allowed, penalty=1.0)
        artifacts[heldout] = prior.to_dict()
        for r in [x for x in pairs if x["system"] == heldout]:
            value = prior.logit(constructs[r["construct_i"]]["architecture"], constructs[r["construct_j"]]["architecture"])
            out.append({"fold": heldout, "pair_id": r["pair_id"], "construct_i": r["construct_i"],
                        "construct_j": r["construct_j"], "literature_prior_logit": value,
                        "probability_i_gt_j": sigmoid(value), "eligible_literature_pair_N": len(allowed)})
    output_dir = Path(output_dir); output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "literature_prior_internal16.csv", out)
    (output_dir / "literature_prior_folds.json").write_text(json.dumps(artifacts, indent=2, sort_keys=True) + "\n")
    return out


def run_loso(data_root, output_dir, prior_artifact):
    constructs, pairs, contexts = load_internal(data_root)
    artifacts = json.loads(Path(prior_artifact).read_text())
    predictions = []
    for heldout in SYSTEMS:
        train, test = [r for r in pairs if r["system"] != heldout], [r for r in pairs if r["system"] == heldout]
        y = np.asarray([1.0 if int(r["target_direction"]) == 1 else 0.0 for r in train])
        train_dx, test_dx = pair_matrix(train), pair_matrix(test)
        prior = ArchitecturePrior.from_dict(artifacts[heldout])
        train_prior = np.asarray([prior.logit(constructs[r["construct_i"]]["architecture"], constructs[r["construct_j"]]["architecture"]) for r in train])
        test_prior = np.asarray([prior.logit(constructs[r["construct_i"]]["architecture"], constructs[r["construct_j"]]["architecture"]) for r in test])
        p0 = fit_scaled(train_dx, y, test_dx)
        p1 = sigmoid(test_prior)
        p2 = fit_scaled(np.column_stack([train_prior, train_dx]), y, np.column_stack([test_prior, test_dx]))
        train_systems = [s for s in SYSTEMS if s != heldout]
        pca = FoldPCA.fit(np.vstack([contexts[s] for s in train_systems]), 1)
        z_lookup = {s: pca.transform(contexts[s][None, :])[0] for s in train_systems}
        z_train = np.vstack([z_lookup[r["system"]] for r in train])
        z_test = np.repeat(pca.transform(contexts[heldout][None, :]), len(test), axis=0)
        x3_train = np.column_stack([train_prior, train_dx, train_dx * z_train])
        x3_test = np.column_stack([test_prior, test_dx, test_dx * z_test])
        p3 = fit_scaled(x3_train, y, x3_test)
        for model, probs in zip(MODELS, (p0, p1, p2, p3)):
            for r, probability, prior_logit in zip(test, probs, test_prior):
                pred = 1 if probability > 0.5 else -1
                predictions.append({"model": model, "fold": heldout, "system": heldout, "pair_id": r["pair_id"],
                                    "construct_i": r["construct_i"], "construct_j": r["construct_j"],
                                    "true_direction": r["target_direction"], "predicted_direction": pred,
                                    "probability_i_gt_j": float(probability), "probability_j_gt_i": float(1-probability),
                                    "pairwise_margin": float(abs(probability-.5)*2), "correct": pred == int(r["target_direction"]),
                                    "literature_prior_logit": float(prior_logit), "unordered_pair_counted_once": True})
    output_dir = Path(output_dir); output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "pairwise_probabilities.csv", predictions)
    return predictions


def rankings_and_metrics(data_root, prediction_path, output_dir):
    constructs, _, _ = load_internal(data_root); predictions = read_csv(prediction_path)
    rankings, system_metrics, aggregate = [], [], []
    for model in MODELS:
        for system in SYSTEMS:
            rows = [r for r in predictions if r["model"] == model and r["system"] == system]
            scores = defaultdict(list)
            for r in rows:
                scores[r["construct_i"]].append(float(r["probability_i_gt_j"])); scores[r["construct_j"]].append(float(r["probability_j_gt_i"]))
            scores = {k: float(np.mean(v)) for k, v in scores.items()}
            order = sorted(scores, key=lambda cid: (-scores[cid], cid))
            true_order = sorted(order, key=lambda cid: -constructs[cid]["phenotype"])
            for rank, cid in enumerate(order, 1):
                rankings.append({"model": model, "system": system, "construct": cid, "predicted_win_probability": scores[cid],
                                 "predicted_rank": rank, "observed_rank": constructs[cid]["within_system_rank"],
                                 "observed_rank_score": constructs[cid]["phenotype"]})
            rho, tau = correlation_metrics([constructs[c]["phenotype"] for c in order], [scores[c] for c in order])
            top2, true_top2 = set(order[:2]), set(true_order[:2])
            system_metrics.append({"model": model, "system": system, "pair_N": 6,
                                   "pair_accuracy": np.mean([r["correct"].lower()=="true" for r in rows]),
                                   "rho": rho, "tau": tau, "Top1": order[0] == true_order[0],
                                   "true_best_in_predicted_top2": true_order[0] in top2, "Top2_exact_set": top2 == true_top2})
        all_pairs = [r for r in predictions if r["model"] == model]
        all_rank = [r for r in rankings if r["model"] == model]
        rho, tau = correlation_metrics([float(r["observed_rank_score"]) for r in all_rank], [float(r["predicted_win_probability"]) for r in all_rank])
        sm = [r for r in system_metrics if r["model"] == model]
        aggregate.append({"model": model, "pair_N": 24, "system_N": 4,
                          "pair_accuracy": np.mean([r["correct"].lower()=="true" for r in all_pairs]), "rho": rho, "tau": tau,
                          "Top1_exact": np.mean([r["Top1"] for r in sm]),
                          "true_best_in_predicted_top2": np.mean([r["true_best_in_predicted_top2"] for r in sm]),
                          "Top2_exact_set": np.mean([r["Top2_exact_set"] for r in sm])})
    output_dir = Path(output_dir); output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "construct_rankings.csv", rankings); write_csv(output_dir / "system_metrics.csv", system_metrics); write_csv(output_dir / "metrics_summary.csv", aggregate)
    return aggregate
