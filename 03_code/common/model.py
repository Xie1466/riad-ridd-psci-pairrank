#!/usr/bin/env python3
"""Small, dependency-light primitives for the frozen PSCI pairwise V2 analysis.

Only NumPy/SciPy are required.  The estimators deliberately have no intercept:
for every fitted score f(dx), f(-dx) = -f(dx), hence sigmoid(f(-dx)) is exactly
1 - sigmoid(f(dx)).  Reverse orientations are never added as biological rows.
"""

from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy.stats import kendalltau, rankdata, spearmanr


EPS = 1e-12
AA = "ACDEFGHIKLMNPQRSTVWY"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[dict], fields: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_float(value: object) -> float:
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NONE", "NOT_IDENTIFIABLE"}:
        return math.nan
    return float(text)


def sigmoid(x: np.ndarray | float) -> np.ndarray | float:
    arr = np.asarray(x, dtype=float)
    out = np.empty_like(arr)
    pos = arr >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-arr[pos]))
    exp_x = np.exp(arr[~pos])
    out[~pos] = exp_x / (1.0 + exp_x)
    return float(out) if out.ndim == 0 else out


@dataclass
class Standardizer:
    mean: np.ndarray
    scale: np.ndarray

    @classmethod
    def fit(cls, x: np.ndarray) -> "Standardizer":
        mean = np.nanmean(x, axis=0)
        scale = np.nanstd(x, axis=0)
        mean = np.where(np.isfinite(mean), mean, 0.0)
        scale = np.where(np.isfinite(scale) & (scale > EPS), scale, 1.0)
        return cls(mean, scale)

    @classmethod
    def fit_scale_only(cls, x: np.ndarray) -> "Standardizer":
        """Scale antisymmetric delta designs without creating an effective intercept."""
        scale = np.nanstd(x, axis=0)
        scale = np.where(np.isfinite(scale) & (scale > EPS), scale, 1.0)
        return cls(np.zeros(x.shape[1]), scale)

    def transform(self, x: np.ndarray) -> np.ndarray:
        filled = np.where(np.isfinite(x), x, self.mean)
        return (filled - self.mean) / self.scale

    def to_dict(self) -> dict:
        return {"mean": self.mean.tolist(), "scale": self.scale.tolist()}

    @classmethod
    def from_dict(cls, value: dict) -> "Standardizer":
        return cls(np.asarray(value["mean"], float), np.asarray(value["scale"], float))


@dataclass
class FoldPCA:
    standardizer: Standardizer
    components: np.ndarray

    @classmethod
    def fit(cls, x: np.ndarray, rank: int) -> "FoldPCA":
        scaler = Standardizer.fit(x)
        z = scaler.transform(x)
        _, _, vt = np.linalg.svd(z, full_matrices=False)
        available = min(rank, vt.shape[0])
        components = vt[:available]
        if available < rank:
            components = np.vstack([components, np.zeros((rank - available, z.shape[1]))])
        # Deterministic component sign.
        for i in range(len(components)):
            nz = np.flatnonzero(np.abs(components[i]) > EPS)
            if len(nz) and components[i, nz[0]] < 0:
                components[i] *= -1
        return cls(scaler, components)

    def transform(self, x: np.ndarray) -> np.ndarray:
        return self.standardizer.transform(x) @ self.components.T

    def to_dict(self) -> dict:
        return {"standardizer": self.standardizer.to_dict(), "components": self.components.tolist()}

    @classmethod
    def from_dict(cls, value: dict) -> "FoldPCA":
        return cls(Standardizer.from_dict(value["standardizer"]), np.asarray(value["components"], float))


@dataclass
class RidgeLogistic:
    coefficient: np.ndarray
    penalty: float = 1.0

    @classmethod
    def fit(cls, x: np.ndarray, y: np.ndarray, penalty: float = 1.0) -> "RidgeLogistic":
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        if x.ndim != 2 or len(y) != len(x) or not len(x):
            raise ValueError("non-empty two-dimensional X and matching y are required")
        if not set(np.unique(y)).issubset({0.0, 1.0}) or len(np.unique(y)) != 2:
            raise ValueError("both binary classes are required")

        def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
            score = x @ beta
            # stable logistic negative log likelihood
            loss = np.logaddexp(0.0, score).sum() - np.dot(y, score)
            loss += 0.5 * penalty * np.dot(beta, beta)
            prob = sigmoid(score)
            grad = x.T @ (prob - y) + penalty * beta
            return float(loss), grad

        beta = np.zeros(x.shape[1])
        for _ in range(200):
            loss, grad = objective(beta)
            prob = sigmoid(x @ beta)
            weight = np.clip(prob * (1.0 - prob), 1e-9, None)
            hessian = x.T @ (x * weight[:, None]) + penalty * np.eye(x.shape[1])
            step = np.linalg.solve(hessian, grad)
            if np.linalg.norm(step) < 1e-9:
                break
            factor = 1.0
            while factor > 1e-8 and objective(beta - factor * step)[0] > loss:
                factor *= 0.5
            beta -= factor * step
            if factor <= 1e-8:
                break
        if not np.isfinite(beta).all():
            raise RuntimeError("logistic IRLS produced non-finite coefficients")
        return cls(beta, penalty)

    def decision_function(self, x: np.ndarray) -> np.ndarray:
        return np.asarray(x, float) @ self.coefficient

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        return sigmoid(self.decision_function(x))

    def to_dict(self) -> dict:
        return {"coefficient": self.coefficient.tolist(), "penalty": self.penalty, "fit_intercept": False}

    @classmethod
    def from_dict(cls, value: dict) -> "RidgeLogistic":
        if value.get("fit_intercept", False):
            raise ValueError("antisymmetric artifact must not contain an intercept")
        return cls(np.asarray(value["coefficient"], float), float(value["penalty"]))


@dataclass
class RidgeRanker:
    coefficient: np.ndarray
    penalty: float = 1.0

    @classmethod
    def fit(cls, x: np.ndarray, delta: np.ndarray, penalty: float = 1.0) -> "RidgeRanker":
        x = np.asarray(x, float)
        delta = np.asarray(delta, float)
        gram = x.T @ x + penalty * np.eye(x.shape[1])
        return cls(np.linalg.solve(gram, x.T @ delta), penalty)

    def decision_function(self, x: np.ndarray) -> np.ndarray:
        return np.asarray(x, float) @ self.coefficient

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        return sigmoid(self.decision_function(x))


def sequence_context(sequence_a: str, sequence_b: str, q_values: tuple[float, float]) -> np.ndarray:
    """Six predefined, label-free native sequence/Q-context summaries."""
    a = "".join(c for c in sequence_a.upper() if c in AA)
    b = "".join(c for c in sequence_b.upper() if c in AA)
    if not a or not b:
        raise ValueError("two non-empty amino-acid sequences are required")
    va = np.asarray([a.count(c) for c in AA], float) / len(a)
    vb = np.asarray([b.count(c) for c in AA], float) / len(b)
    cosine = float(np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb) + EPS))
    hydrophobic = set("AVILMFWY")
    hydro_mean = 0.5 * (sum(c in hydrophobic for c in a) / len(a) + sum(c in hydrophobic for c in b) / len(b))
    return np.asarray([
        0.5 * (math.log(len(a)) + math.log(len(b))),
        abs(math.log(len(a) / len(b))),
        cosine,
        hydro_mean,
        float(q_values[0]),
        float(q_values[1]),
    ])


ARCH_CATEGORICAL = (
    "attachment_A", "attachment_B", "N/C_topology", "recruitment_A",
    "recruitment_B", "scaffold_platform",
)
ARCH_NUMERIC = ("linker_length_A", "linker_length_B")


@dataclass
class ArchitectureEncoder:
    levels: dict[str, list[str]]

    @classmethod
    def fit(cls, constructs: Iterable[dict[str, str]]) -> "ArchitectureEncoder":
        rows = list(constructs)
        levels = {
            col: sorted({(r.get(col) or "<MISSING>").strip() for r in rows})
            for col in ARCH_CATEGORICAL
        }
        return cls(levels)

    @property
    def feature_names(self) -> list[str]:
        names: list[str] = []
        for col in ARCH_CATEGORICAL:
            names.extend(f"{col}={level}" for level in self.levels[col])
        names.extend(ARCH_NUMERIC)
        return names

    def transform_one(self, row: dict[str, str]) -> np.ndarray:
        values: list[float] = []
        for col in ARCH_CATEGORICAL:
            value = (row.get(col) or "<MISSING>").strip()
            values.extend(float(value == level) for level in self.levels[col])
        for col in ARCH_NUMERIC:
            text = str(row.get(col, "")).strip()
            try:
                raw = clean_float(text)
            except ValueError:
                numbers = [float(value) for value in re.findall(r"\d+(?:\.\d+)?", text)]
                raw = sum(numbers) if "+" in text and numbers else (numbers[0] if numbers else math.nan)
            values.append(math.log1p(raw) if math.isfinite(raw) and raw >= 0 else 0.0)
        return np.asarray(values)

    def novelty(self, row: dict[str, str]) -> int:
        return sum((row.get(col) or "<MISSING>").strip() not in self.levels[col] for col in ARCH_CATEGORICAL)

    def to_dict(self) -> dict:
        return {"levels": self.levels, "feature_names": self.feature_names}

    @classmethod
    def from_dict(cls, value: dict) -> "ArchitectureEncoder":
        return cls({key: list(items) for key, items in value["levels"].items()})


@dataclass
class ArchitecturePrior:
    encoder: ArchitectureEncoder
    scaler: Standardizer
    model: RidgeLogistic

    @classmethod
    def fit(
        cls, constructs: dict[str, dict[str, str]], pairs: Sequence[dict[str, str]], penalty: float = 1.0
    ) -> "ArchitecturePrior":
        used_ids = {r[k] for r in pairs for k in ("construct_i", "construct_j")}
        encoder = ArchitectureEncoder.fit(constructs[cid] for cid in sorted(used_ids))
        x = np.vstack([
            encoder.transform_one(constructs[r["construct_i"]]) - encoder.transform_one(constructs[r["construct_j"]])
            for r in pairs
        ])
        y = np.asarray([1.0 if r["rank_direction"] == "I_GT_J" else 0.0 for r in pairs])
        scaler = Standardizer.fit_scale_only(x)
        model = RidgeLogistic.fit(scaler.transform(x), y, penalty)
        return cls(encoder, scaler, model)

    def logit(self, left: dict[str, str], right: dict[str, str]) -> float:
        delta = self.encoder.transform_one(left) - self.encoder.transform_one(right)
        return float(self.model.decision_function(self.scaler.transform(delta[None, :]))[0])

    def probability(self, left: dict[str, str], right: dict[str, str]) -> float:
        return float(sigmoid(self.logit(left, right)))

    def to_dict(self) -> dict:
        return {"encoder": self.encoder.to_dict(), "scaler": self.scaler.to_dict(), "model": self.model.to_dict()}

    @classmethod
    def from_dict(cls, value: dict) -> "ArchitecturePrior":
        return cls(
            ArchitectureEncoder.from_dict(value["encoder"]),
            Standardizer.from_dict(value["scaler"]),
            RidgeLogistic.from_dict(value["model"]),
        )


def correlation_metrics(observed: Sequence[float], predicted: Sequence[float]) -> tuple[float, float]:
    if len(observed) < 2 or len(set(predicted)) < 2:
        return 0.0, 0.0
    rho = spearmanr(observed, predicted).statistic
    tau = kendalltau(observed, predicted).statistic
    return float(0.0 if not np.isfinite(rho) else rho), float(0.0 if not np.isfinite(tau) else tau)


def json_dump(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
