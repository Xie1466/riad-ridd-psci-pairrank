#!/usr/bin/env python3
"""Reanalyse the 16 existing trajectories using the phenotype-blind PSCI freeze.

No MD is run.  Existing protein-only, PBC-corrected trajectories and existing
uniform D/L/0.45-nm-contact transforms are reused read-only.  New calculations
are sampled at the frozen 100-ps stride and are restartable per topology.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(os.environ.get("PSCI_MD_SOURCE_ROOT", "05_data/md_reduced_or_manifest"))
GMX = Path(os.environ.get("GMX_BIN", "gmx"))
MANIFEST = ROOT / "metadata/analysis_trajectory_manifest.csv"
COMPONENTS = ROOT / "metadata/component_group_map.csv"
ACTIVE = ROOT / "metadata/ACTIVE_SITE_RESIDUE_REGISTRY_INTERNAL16_V2.csv"
FREEZE = ROOT / "metadata/PSCI_PHENOTYPE_BLIND_FEATURE_FREEZE.json"
OUT = ROOT / "results/PSCI_FINAL"
LOGROOT = ROOT / "logs/psci_internal16_reanalysis"
VERSION = "PSCI_INTERNAL16_REANALYSIS_V1"
WINDOWS = (("primary_0_75", 0.0, 75.0), ("sensitivity_0_50", 0.0, 50.0),
           ("sensitivity_0_100", 0.0, 100.0))


def args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--only", action="append")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--aggregate-only", action="store_true")
    return p.parse_args()


def utc():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields=None):
    if not rows:
        raise RuntimeError(f"refusing empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0])
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    tmp.replace(path)


def sha(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def xvg(path: Path):
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.strip() and not line.lstrip().startswith(("#", "@")):
            rows.append([float(x) for x in line.split()])
    out = np.asarray(rows, float)
    if out.ndim != 2 or not len(out):
        raise RuntimeError(f"invalid XVG {path}")
    return out


def run(log: Path, command: list[str], stdin=""):
    log.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    proc = subprocess.run(command, input=stdin or None, text=True, capture_output=True)
    elapsed = time.monotonic() - started
    log.write_text(f"started_utc={utc()}\ncommand={shlex.join(command)}\nstdin={stdin!r}\n"
                   f"version={VERSION}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}\n"
                   f"exit_status={proc.returncode}\nwall_seconds={elapsed:.3f}\n", encoding="utf-8")
    if proc.returncode:
        raise RuntimeError(f"command failed; see {log}")


def gro_atoms(path: Path):
    lines = path.read_text().splitlines(); n = int(lines[1])
    if len(lines) < n + 3:
        raise RuntimeError(f"truncated GRO {path}")
    out = {}
    for atom_id, line in enumerate(lines[2:2+n], 1):
        out[atom_id] = {"resnr": int(line[:5]), "resname": line[5:10].strip(),
                        "atomname": line[10:15].strip(),
                        "xyz": np.asarray([float(line[20:28]), float(line[28:36]), float(line[36:44])])}
    return out


def write_ndx(path: Path, groups: dict[str, list[int]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for name, values in groups.items():
            if not values:
                raise RuntimeError(f"empty group {name}")
            handle.write(f"[ {name} ]\n")
            for i in range(0, len(values), 15):
                handle.write(" ".join(f"{x:6d}" for x in values[i:i+15]) + "\n")


def active_position(value: str):
    m = re.fullmatch(r"[A-Z](\d+)", value)
    return int(m.group(1)) if m else None


def build_groups(tid: str, reference: Path, component_rows, active_rows):
    atoms = gro_atoms(reference)
    rows = [x for x in component_rows if x["topology_id"] == tid and x["chain"] in ("1", "2")]
    if len(rows) != 6:
        raise RuntimeError(f"{tid}: expected 6 component rows, got {len(rows)}")
    groups = {"Protein": sorted(atoms)}
    meta = {}
    for chain in (1, 2):
        cr = [x for x in rows if int(x["chain"]) == chain]
        enzyme = next(x for x in cr if x["component"] == f"enzyme_{chain}")
        linker = next(x for x in cr if x["component"] == f"linker_{chain}")
        tag = next(x for x in cr if x["component"] in ("RIAD", "RIDD"))
        ids = lambda r: list(range(int(r["atom_start_global"]), int(r["atom_end_global"]) + 1))
        groups[f"enzyme_{chain}"] = ids(enzyme)
        groups[f"enzyme_{chain}_heavy"] = [i for i in groups[f"enzyme_{chain}"] if not atoms[i]["atomname"].startswith("H")]
        groups[f"enzyme_{chain}_ca"] = [i for i in groups[f"enzyme_{chain}"] if atoms[i]["atomname"] == "CA"]
        groups[f"chain_{chain}_heavy"] = [i for r in cr for i in ids(r) if not atoms[i]["atomname"].startswith("H")]
        groups[f"linker_{chain}_ca"] = [i for i in ids(linker) if atoms[i]["atomname"] == "CA"]
        if len(groups[f"linker_{chain}_ca"]) != int(linker["residue_count"]):
            raise RuntimeError(f"{tid}: linker CA count mismatch")
        fused_c = int(linker["residue_start_local"]) == int(enzyme["residue_end_local"]) + 1
        fused_n = int(linker["residue_end_local"]) + 1 == int(enzyme["residue_start_local"])
        if fused_c == fused_n:
            raise RuntimeError(f"{tid}: cannot determine fused terminus")
        erange = list(range(int(enzyme["residue_start_local"]), int(enzyme["residue_end_local"]) + 1))
        terminal_res = set(erange[-5:] if fused_c else erange[:5])
        groups[f"terminal_{chain}"] = [i for i in ids(enzyme) if atoms[i]["resnr"] in terminal_res]
        groups[f"terminal_{chain}_heavy"] = [i for i in groups[f"terminal_{chain}"] if not atoms[i]["atomname"].startswith("H")]
        groups[f"terminal_{chain}_ca"] = [i for i in groups[f"terminal_{chain}"] if atoms[i]["atomname"] == "CA"]
        site_positions = {active_position(x["local_residue"]) for x in active_rows
                          if x["topology_id"] == tid and int(x["chain"]) == chain and active_position(x["local_residue"]) is not None}
        if site_positions:
            groups[f"active_{chain}"] = [i for i in ids(enzyme) if atoms[i]["resnr"] in site_positions]
            groups[f"active_{chain}_ca"] = [i for i in groups[f"active_{chain}"] if atoms[i]["atomname"] == "CA"]
            if len(groups[f"active_{chain}_ca"]) != len(site_positions):
                raise RuntimeError(f"{tid}: active CA mismatch chain {chain}: {len(groups[f'active_{chain}_ca'])}/{len(site_positions)}")
        coords = np.asarray([atoms[i]["xyz"] for i in groups[f"enzyme_{chain}_ca"]])
        rg = float(np.sqrt(np.mean(np.sum((coords - coords.mean(axis=0)) ** 2, axis=1))))
        meta[chain] = {"enzyme": next(x["enzyme"] for x in active_rows if x["topology_id"] == tid and int(x["chain"]) == chain),
                       "rg_nm": rg, "active_n": len(site_positions), "terminal": "C" if fused_c else "N",
                       "terminal_residues": sorted(terminal_res), "linker_n": int(linker["residue_count"]),
                       "tag": tag["component"]}
    groups["terminal_both_ca"] = groups["terminal_1_ca"] + groups["terminal_2_ca"]
    return groups, meta


def ensure_analysis(row, groups, meta):
    tid, system = row["topology_id"], row["canonical_system_id"]
    work = OUT / "timeseries" / tid; logs = LOGROOT / tid
    work.mkdir(parents=True, exist_ok=True); logs.mkdir(parents=True, exist_ok=True)
    xtc = SOURCE / row["analysis_xtc"]; ref = SOURCE / row["analysis_reference"]
    ndx = work / "psci_groups.ndx"
    write_ndx(ndx, groups)
    # Reuse validated 10-ps source transforms, but uniformly select 100-ps frames.
    spatial = SOURCE / f"results/{system}/{tid}/spatial_metrics"
    base = read_csv(spatial / "spatial_metrics_timeseries.csv")
    base100 = [x for x in base if abs((float(x["time_ns"]) * 10) - round(float(x["time_ns"]) * 10)) < 1e-6]
    if len(base100) != 1001:
        raise RuntimeError(f"{tid}: expected 1001 spatial points, got {len(base100)}")
    # Heavy-atom contact pairs at 0.45 nm already exist from the same validated transform.
    c45 = xvg(spatial / "gromacs/enzyme_atom_contacts_0p45.xvg")[::10]
    if len(c45) != 1001:
        raise RuntimeError(f"{tid}: bad reused 0.45 contact grid")
    for cutoff in (0.40, 0.50):
        path = work / f"enzyme_atom_contacts_{cutoff:.2f}.xvg"
        if path.exists() and len(xvg(path)) != 1001:
            path.unlink()
        if not path.exists():
            run(logs / f"contacts_{cutoff:.2f}.log",
                [str(GMX), "mindist", "-f", str(xtc), "-s", str(ref), "-n", str(ndx),
                 "-on", str(path), "-od", str(work / f"enzyme_min_distance_{cutoff:.2f}.xvg"),
                 "-d", f"{cutoff:.2f}", "-dt", "0.1", "-tu", "ns"],
                "enzyme_1_heavy\nenzyme_2_heavy\n")
    c40, c50 = xvg(work / "enzyme_atom_contacts_0.40.xvg"), xvg(work / "enzyme_atom_contacts_0.50.xvg")
    if len(c40) != 1001 or len(c50) != 1001:
        raise RuntimeError(f"{tid}: contact grid failure")
    have_active = "active_1" in groups and "active_2" in groups
    geom = work / "geometry.xvg"; vec = work / "geometry_vectors.xvg"
    if have_active and (not geom.exists() or not vec.exists()):
        selections = "; ".join([
            'cog of group "enzyme_1_ca" plus cog of group "enzyme_2_ca"',
            'cog of group "enzyme_1_ca" plus cog of group "active_1_ca"',
            'cog of group "enzyme_2_ca" plus cog of group "active_2_ca"',
            'cog of group "terminal_1_ca" plus cog of group "active_1_ca"',
            'cog of group "terminal_2_ca" plus cog of group "active_2_ca"'])
        run(logs / "geometry.log", [str(GMX), "distance", "-f", str(xtc), "-s", str(ref), "-n", str(ndx),
            "-oall", str(geom), "-oxyz", str(vec), "-dt", "0.1", "-tu", "ns", "-select", selections])
    sasa = work / "sasa.xvg"
    if not sasa.exists():
        selections = ['group "terminal_1"', 'group "terminal_2"']
        if have_active: selections = ['group "active_1"', 'group "active_2"'] + selections
        run(logs / "sasa.log", [str(GMX), "sasa", "-f", str(xtc), "-s", str(ref), "-n", str(ndx),
            "-o", str(sasa), "-dt", "0.1", "-tu", "ns", "-probe", "0.14",
            "-surface", 'group "Protein"', "-output", "; ".join(selections)])
    for chain in (1, 2):
        path = work / f"terminal_{chain}_nonself_mindist.xvg"
        if not path.exists():
            other = 2 if chain == 1 else 1
            run(logs / f"terminal_{chain}_nonself.log", [str(GMX), "mindist", "-f", str(xtc), "-s", str(ref),
                "-n", str(ndx), "-od", str(path), "-dt", "0.1", "-tu", "ns"],
                f"terminal_{chain}_heavy\nchain_{other}_heavy\n")
    for label, _, end in WINDOWS:
        path = work / f"terminal_rmsf_{label}.xvg"
        if not path.exists():
            run(logs / f"terminal_rmsf_{label}.log", [str(GMX), "rmsf", "-f", str(xtc), "-s", str(ref),
                "-n", str(ndx), "-o", str(path), "-b", "0", "-e", str(int(end * 1000)), "-res"],
                "terminal_both_ca\n")
    time_ns = np.asarray([float(x["time_ns"]) for x in base100])
    values = {
        "time_ns": time_ns,
        "D_com_nm": np.asarray([float(x["enzyme_com_distance_nm"]) for x in base100]),
        "D_min_nm": np.asarray([float(x["enzyme_min_distance_nm"]) for x in base100]),
        "C_atom_contacts_0p40": c40[:, 1], "C_atom_contacts_0p45": c45[:, 1], "C_atom_contacts_0p50": c50[:, 1],
        "L_extension_ratio": np.asarray([(float(x["linker_1_extension_ratio"]) + float(x["linker_2_extension_ratio"])) / 2 for x in base100]),
        "Q_riad_ridd_contacts_0p45": np.asarray([float(x["riad_ridd_contacts_0p45"]) for x in base100]),
        "Q_riad_ridd_com_nm": np.asarray([float(x["riad_ridd_com_distance_nm"]) for x in base100]),
    }
    characteristic_rg = (meta[1]["rg_nm"] + meta[2]["rg_nm"]) / 2
    values["D_normalized"] = values["D_com_nm"] / characteristic_rg
    denom = math.sqrt(len(groups["enzyme_1_heavy"]) * len(groups["enzyme_2_heavy"]))
    values["C_density_0p45"] = values["C_atom_contacts_0p45"] / denom
    sd = xvg(sasa)
    # gmx sasa writes time, total surface, then one column per output selection.
    expected_cols = 6 if have_active else 4
    if sd.shape != (1001, expected_cols):
        raise RuntimeError(f"{tid}: unexpected SASA shape {sd.shape}")
    offset = 2 if have_active else 0
    f1 = np.column_stack([sd[:, 2 + offset] / max(sd[0, 2 + offset], 1e-12),
                          sd[:, 3 + offset] / max(sd[0, 3 + offset], 1e-12)])
    values["F1_terminal_accessibility"] = f1.mean(axis=1)
    mind1, mind2 = xvg(work / "terminal_1_nonself_mindist.xvg"), xvg(work / "terminal_2_nonself_mindist.xvg")
    values["F4_terminal_nonself_contact"] = ((mind1[:, 1] <= .45) | (mind2[:, 1] <= .45)).astype(float)
    if have_active:
        active_reference_valid = bool(sd[0, 2] > 1e-12 and sd[0, 3] > 1e-12)
        if active_reference_valid:
            active_ratio = np.column_stack([sd[:, 2] / sd[0, 2], sd[:, 3] / sd[0, 3]])
            values["A_normalized_accessibility"] = active_ratio.mean(axis=1)
        gd, gv = xvg(geom), xvg(vec)
        if gd.shape != (1001, 6) or gv.shape != (1001, 16):
            raise RuntimeError(f"{tid}: geometry shape mismatch {gd.shape}/{gv.shape}")
        vectors = [gv[:, 1 + 3*i:4 + 3*i] for i in range(5)]
        unit = lambda a: a / np.maximum(np.linalg.norm(a, axis=1)[:, None], 1e-12)
        o1 = np.sum(unit(vectors[1]) * unit(vectors[0]), axis=1)
        o2 = np.sum(unit(vectors[2]) * -unit(vectors[0]), axis=1)
        values["O1_cosine"] = o1; values["O2_cosine"] = o2
        values["O_mutually_facing"] = ((o1 >= .5) & (o2 >= .5)).astype(float)
        values["F2_terminal_active_normalized"] = np.column_stack([gd[:, 4] / meta[1]["rg_nm"],
                                                                    gd[:, 5] / meta[2]["rg_nm"]]).mean(axis=1)
    fields = list(values)
    rows = [{k: (f"{values[k][i]:.10g}" if k != "time_ns" else f"{values[k][i]:.1f}") for k in fields}
            for i in range(len(time_ns))]
    write_csv(work / "PSCI_timeseries_100ps.csv", rows, fields)
    qc = {"topology_id": tid, "status": "PASS", "source_xtc": str(xtc), "source_xtc_sha256": row.get("analysis_sha256", "not_recorded"),
          "source_reference": str(ref), "frames": 1001, "start_ns": 0, "end_ns": 100, "stride_ps": 100,
          "active_pair_status": "EVIDENCE_ELIGIBLE" if have_active else "SCIENTIFICALLY_NOT_IDENTIFIABLE",
          "active_frame0_sasa_status": "PASS" if have_active and active_reference_valid else
                                       "NUMERIC_MISSING_ZERO_REFERENCE" if have_active else "NOT_APPLICABLE",
          "characteristic_rg_nm": characteristic_rg, "contact_density_denominator": denom,
          "groups": {k: len(v) for k, v in groups.items()}, "component_meta": meta,
          "timeseries_sha256": sha(work / "PSCI_timeseries_100ps.csv"), "completed_utc": utc(), "version": VERSION}
    (work / "qc.json").write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return qc


def summary_stats(v):
    v = np.asarray(v, float)
    return {"mean": float(np.mean(v)), "sd": float(np.std(v, ddof=1)) if len(v) > 1 else 0.0, "median": float(np.median(v)),
            "q25": float(np.quantile(v, .25)), "q75": float(np.quantile(v, .75)),
            "min": float(np.min(v)), "max": float(np.max(v))}


def aggregate(manifest):
    raw, master, confidence, complete = [], [], [], []
    metrics = ["D_normalized", "D_com_nm", "D_min_nm", "C_atom_contacts_0p40", "C_atom_contacts_0p45",
               "C_atom_contacts_0p50", "C_density_0p45", "A_normalized_accessibility", "O1_cosine", "O2_cosine",
               "O_mutually_facing", "F1_terminal_accessibility", "F2_terminal_active_normalized",
               "F4_terminal_nonself_contact", "L_extension_ratio", "Q_riad_ridd_contacts_0p45", "Q_riad_ridd_com_nm"]
    for src in manifest:
        tid, system = src["topology_id"], src["canonical_system_id"]
        tsp = OUT / "timeseries" / tid / "PSCI_timeseries_100ps.csv"
        qcp = OUT / "timeseries" / tid / "qc.json"
        if not tsp.exists() or not qcp.exists():
            raise RuntimeError(f"missing completed topology {tid}")
        ts = read_csv(tsp); qc = json.loads(qcp.read_text())
        for label, start, end in WINDOWS:
            use = [x for x in ts if float(x["time_ns"]) >= start and float(x["time_ns"]) <= end]
            for metric in metrics:
                vals = [float(x[metric]) for x in use if x.get(metric, "") not in ("", "NA")]
                if not vals:
                    raw.append({"topology_id": tid, "system_id": system, "window": label, "metric": metric,
                                "mean": "NA", "sd": "NA", "median": "NA", "q25": "NA", "q75": "NA",
                                "min": "NA", "max": "NA", "n_frames": 0,
                                "missingness": ("SCIENTIFICALLY_NOT_IDENTIFIABLE" if system == "GI-DAE" and metric.startswith(("A_", "O", "F2"))
                                                else "NUMERIC_MISSING_ZERO_REFERENCE" if metric == "A_normalized_accessibility"
                                                else "NUMERIC_MISSING"),
                                "source": str(tsp.relative_to(ROOT)), "confidence": "not_identifiable"})
                else:
                    stat = summary_stats(vals)
                    conf = "low_proxy" if (system in ("GP-PGM", "RHLA-RHLB") and metric.startswith(("A_", "O", "F2"))) else "high_simulation_descriptor"
                    raw.append({"topology_id": tid, "system_id": system, "window": label, "metric": metric,
                                **stat, "n_frames": len(vals), "missingness": "NONE", "source": str(tsp.relative_to(ROOT)), "confidence": conf})
            rmsf = xvg(OUT / "timeseries" / tid / f"terminal_rmsf_{label}.xvg")[:, 1]
            f5 = 1 / (1 + float(np.mean(rmsf)))
            raw.append({"topology_id": tid, "system_id": system, "window": label, "metric": "F5_terminal_constraint",
                        **summary_stats([f5]), "n_frames": len(rmsf), "missingness": "NONE",
                        "source": str((OUT / "timeseries" / tid / f"terminal_rmsf_{label}.xvg").relative_to(ROOT)),
                        "confidence": "high_simulation_descriptor"})
        primary = {(x["metric"]): x for x in raw if x["topology_id"] == tid and x["window"] == "primary_0_75"}
        get = lambda k: primary[k]["mean"]
        row = {"dataset_version": VERSION, "topology_id": tid, "system_id": system, "window": "0-75 ns",
               "D": get("D_normalized"), "C": get("C_density_0p45"),
               "A": get("A_normalized_accessibility"), "O": get("O_mutually_facing"),
               "F_scalar": "NA", "F1": get("F1_terminal_accessibility"),
               "F2": get("F2_terminal_active_normalized"), "F3": "NA", "F4": get("F4_terminal_nonself_contact"),
               "F5": get("F5_terminal_constraint"), "L": get("L_extension_ratio"),
               "Q_context": "SYSTEM_CONTEXT_ONLY", "Q_topology": "NA",
               "representation": "STRUCTURED_PSCI_[D,C,A,O,F_VECTOR,L,Q_CONTEXT,Q_TOPOLOGY]",
               "source_timeseries": str(tsp.relative_to(ROOT)), "qc": qc["status"]}
        master.append(row)
        dim_status = {"D": "NUMERIC", "C": "NUMERIC", "A": "NUMERIC" if row["A"] != "NA" else "SCIENTIFICALLY_NOT_IDENTIFIABLE",
                      "O": "NUMERIC" if row["O"] != "NA" else "SCIENTIFICALLY_NOT_IDENTIFIABLE",
                      "F": "VECTOR_PARTIAL_F3_NOT_IDENTIFIABLE", "L": "NUMERIC", "Q": "SYSTEM_CONTEXT_ONLY"}
        for dim, status in dim_status.items():
            confidence.append({"topology_id": tid, "system_id": system, "dimension": dim, "status": status,
                               "confidence": "LOW_PROXY" if dim in ("A", "O") and system in ("GP-PGM", "RHLA-RHLB") else
                                             "NOT_IDENTIFIABLE" if "NOT_IDENTIFIABLE" in status else "HIGH_DESCRIPTOR",
                               "source": str(qcp.relative_to(ROOT))})
        numeric = sum(status == "NUMERIC" for status in dim_status.values())
        complete.append({"topology_id": tid, "system_id": system, "numeric_scalar_dimensions": numeric,
                         "full_7D_complete": "false", "reduced_D_C_F1_F4_F5_L_complete": "true",
                         "A_available": str(row["A"] != "NA").lower(), "O_available": str(row["O"] != "NA").lower(),
                         "A_O_both_available": str(row["A"] != "NA" and row["O"] != "NA").lower(), "F_scalar": "F_SCALAR_NOT_IDENTIFIABLE",
                         "Q_status": "SYSTEM_CONTEXT_ONLY", "missingness": ";".join(f"{d}:{s}" for d, s in dim_status.items() if s != "NUMERIC")})
    write_csv(OUT / "INTERNAL16_PSCI_RAW_DESCRIPTORS_FINAL.csv", raw)
    write_csv(OUT / "INTERNAL16_PSCI_MASTER_FINAL.csv", master)
    write_csv(OUT / "INTERNAL16_PSCI_CONFIDENCE_FINAL.csv", confidence)
    write_csv(OUT / "INTERNAL16_PSCI_COMPLETENESS_FINAL.csv", complete)
    return master


def main():
    a = args()
    freeze = json.loads(FREEZE.read_text())
    if freeze["phenotype_opened_by_this_pipeline"] is not False or freeze["primary_window_ns"] != [0, 75]:
        raise RuntimeError("phenotype-blind freeze gate failed")
    manifest = read_csv(MANIFEST)
    if len(manifest) != 16:
        raise RuntimeError("trajectory manifest is not Internal16")
    if a.only:
        manifest_run = [x for x in manifest if x["topology_id"] in set(a.only)]
    else:
        manifest_run = manifest
    components, active = read_csv(COMPONENTS), read_csv(ACTIVE)
    if a.dry_run:
        print(json.dumps({"topologies": [x["topology_id"] for x in manifest_run], "new_md": False,
                          "primary_window_ns": [0, 75], "stride_ps": 100, "expected_persistent_storage_gb": "<1"}))
        return 0
    version = subprocess.run([str(GMX), "--version"], capture_output=True, text=True)
    if version.returncode:
        raise RuntimeError("GROMACS unavailable")
    if not a.aggregate_only:
        for row in manifest_run:
            tid = row["topology_id"]; tic = time.monotonic()
            ref = SOURCE / row["analysis_reference"]
            groups, meta = build_groups(tid, ref, components, active)
            qc = ensure_analysis(row, groups, meta)
            print(json.dumps({"topology_id": tid, "status": qc["status"], "seconds": round(time.monotonic()-tic, 2)}, sort_keys=True), flush=True)
    if not a.only and all((OUT / "timeseries" / x["topology_id"] / "qc.json").exists() for x in manifest):
        master = aggregate(manifest)
        print(json.dumps({"status": "PASS_INTERNAL16_REANALYSIS", "topologies": len(master),
                          "full_7D": 0, "reduced_common": 16}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(2)
