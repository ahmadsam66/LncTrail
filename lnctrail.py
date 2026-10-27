#!/usr/bin/env python3
"""LncTrail: reproducible lncRNA-to-transcript RNA-RNA candidate discovery.

Primary engine: IntaRNA (external executable). Fallback: transparent seed-and-extend
sequence heuristic. The fallback does NOT calculate thermodynamic free energy and its
score is not a probability of biological interaction.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import datetime as dt
import gzip
import hashlib
import html
import io
import json
import math
import os
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

VERSION = "0.2.3"
DNA_RNA_CHARS = set("ACGUTNRYSWKMBDHV")
CANONICAL = {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")}
GU_WOBBLE = {("G", "U"), ("U", "G")}
COMP = str.maketrans({"A": "U", "U": "A", "T": "A", "C": "G", "G": "C",
                      "N": "N", "R": "Y", "Y": "R", "S": "S", "W": "W",
                      "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D"})

TSV_COLUMNS = [
    "lncrna_id", "lncrna_description", "target_id", "target_description",
    "lncrna_start", "lncrna_end", "target_start", "target_end",
    "lncrna_length", "target_length", "lncrna_gc_percent", "target_gc_percent",
    "engine", "evidence_type", "prediction_status", "target_rank", "site_rank",
    "energy_kcal_mol", "hybridization_energy_kcal_mol", "energy_normalized",
    "accessibility_target_probability", "accessibility_lncrna_probability", "interaction_probability",
    "seed_start_lncrna", "seed_end_lncrna", "seed_start_target", "seed_end_target", "seed_energy_kcal_mol",
    "seed_length", "canonical_pairs", "gu_wobble_pairs", "mismatches", "binding_score",
    "lncrna_interacting_sequence", "target_interacting_sequence", "interaction_structure", "interaction_pattern",
    "raw_engine_fields",
    "target_gene", "target_isoform_count", "priority_score", "integrated_priority_score", "priority_tier",
    "seed_transcriptome_occurrences", "seed_specificity_score",
    "lncrna_low_complexity_fraction", "target_low_complexity_fraction",
    "canonical_fraction", "alignment_length", "rime_max_score", "rime_windows_above_0_5", "global_priority_rank"
]

@dataclass(frozen=True)
class SequenceRecord:
    id: str
    description: str
    sequence: str
    original_had_t: bool = False


def open_text(path: Path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fasta(path: Path, label: str) -> tuple[list[SequenceRecord], list[str]]:
    if not path.exists():
        raise FileNotFoundError(f"{label} FASTA not found: {path}")
    records: list[SequenceRecord] = []
    warnings: list[str] = []
    header: str | None = None
    chunks: list[str] = []

    def add_record(h: str, seq_chunks: list[str]) -> None:
        if not h:
            raise ValueError(f"{label}: empty FASTA header in {path}")
        seq_raw = "".join(seq_chunks).upper().replace(" ", "").replace("\t", "")
        seq_raw = seq_raw.replace("-", "").replace(".", "")
        if not seq_raw:
            raise ValueError(f"{label}: no sequence for FASTA header >{h}")
        invalid = sorted(set(seq_raw) - DNA_RNA_CHARS)
        if invalid:
            raise ValueError(f"{label}: unsupported character(s) {invalid} in record >{h}. Accepted nucleotide/IUPAC codes: {''.join(sorted(DNA_RNA_CHARS))}")
        had_t = "T" in seq_raw
        seq = seq_raw.replace("T", "U")
        seq_id = h.split()[0]
        if any(rec.id == seq_id for rec in records):
            raise ValueError(f"{label}: duplicated FASTA ID '{seq_id}'. IDs (first word after >) must be unique.")
        if len(seq) < 5:
            warnings.append(f"{label} record '{seq_id}' is very short ({len(seq)} nt); predictions may be uninformative.")
        unknown_fraction = sum(base not in "ACGU" for base in seq) / len(seq)
        if unknown_fraction > 0.10:
            warnings.append(f"{label} record '{seq_id}' contains {unknown_fraction:.1%} ambiguous nucleotides.")
        if len(set(seq) & set("ACGU")) < 3:
            warnings.append(f"{label} record '{seq_id}' has very low nucleotide diversity.")
        records.append(SequenceRecord(seq_id, h, seq, had_t))

    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    add_record(header, chunks)
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"{label}: sequence text before first FASTA header at line {line_no}")
                chunks.append(line)
    if header is not None:
        add_record(header, chunks)
    if not records:
        raise ValueError(f"{label}: no FASTA records found in {path}")
    return records, warnings


def gc_percent(seq: str) -> float:
    known = [b for b in seq if b in "ACGU"]
    if not known:
        return 0.0
    return 100.0 * sum(b in "GC" for b in known) / len(known)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def first_available_engine(requested: str) -> str:
    inta = shutil.which("IntaRNA") or shutil.which("intarna")
    if requested == "auto":
        return "intarna" if inta else "builtin"
    if requested == "intarna" and not inta:
        raise RuntimeError("--engine intarna was requested, but IntaRNA is not on PATH. Install it with: conda install -c conda-forge -c bioconda intarna")
    return requested


def parse_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        v = str(value).strip()
        if not v or v.lower() in {"nan", "none", "na", "n/a", "-"}:
            return None
        return float(v)
    except (ValueError, TypeError):
        return None


def parse_int(value: Any) -> int | None:
    v = parse_float(value)
    return int(v) if v is not None else None


def write_normalized_fasta(records: list[SequenceRecord], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        for rec in records:
            handle.write(f">{rec.id}\n")
            for i in range(0, len(rec.sequence), 80):
                handle.write(rec.sequence[i:i+80] + "\n")


def run_intarna(lncs: list[SequenceRecord], targets: list[SequenceRecord], args: argparse.Namespace,
                temp_dir: Path, batch_index: int | None = None) -> tuple[list[dict[str, Any]], str, list[str]]:
    """Run IntaRNA for one target batch and parse its tabular output.

    Batching is handled by the caller. This function intentionally executes one
    bounded batch at a time so completed batches can be persisted immediately.
    """
    executable = shutil.which("IntaRNA") or shutil.which("intarna")
    if not executable:
        raise RuntimeError("IntaRNA executable not found on PATH")
    tag = f".batch_{batch_index:05d}" if batch_index is not None else ""
    qfile = temp_dir / "lncrna.normalized.fasta"
    tfile = temp_dir / f"targets.normalized{tag}.fasta"
    write_normalized_fasta(lncs, qfile)
    write_normalized_fasta(targets, tfile)
    cols = ["id1", "id2", "start1", "end1", "start2", "end2", "subseq1", "subseq2",
            "hybridDP", "hybridDB", "E", "ED1", "ED2", "E_hybrid", "E_norm",
            "seedStart1", "seedEnd1", "seedStart2", "seedEnd2", "seedE", "Pu1", "Pu2", "P_E"]
    cmd = [executable, "-t", str(tfile), "-q", str(qfile), "--outMode=C",
           "--outCsvCols=" + ",".join(cols), "--outBestSeedOnly", "--threads=" + str(args.threads),
           "-n", str(args.max_sites_per_pair), "--seedBP=" + str(args.seed_length)]
    if args.allow_overlapping_sites:
        cmd += ["--outOverlap=B"]
    else:
        cmd += ["--outOverlap=N"]
    if args.temperature is not None:
        cmd += ["--temperature=" + str(args.temperature)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    combined = result.stdout
    raw_prefix = temp_dir / (f"intarna{tag}" if tag else "intarna")
    if result.returncode != 0:
        (Path(str(raw_prefix) + ".stdout.txt")).write_text(result.stdout, encoding="utf-8")
        (Path(str(raw_prefix) + ".stderr.txt")).write_text(result.stderr, encoding="utf-8")
        err = result.stderr.strip()[-6000:]
        raise RuntimeError("IntaRNA execution failed. Command: " + " ".join(cmd) + "\n" + err)
    if args.keep_intermediate:
        Path(str(raw_prefix) + ".stdout.txt").write_text(result.stdout, encoding="utf-8")
        Path(str(raw_prefix) + ".stderr.txt").write_text(result.stderr, encoding="utf-8")
    records_by_id = {r.id: r for r in lncs}
    target_by_id = {r.id: r for r in targets}
    lines = combined.splitlines()
    header_i = None
    delimiter = ";"
    for i, line in enumerate(lines):
        stripped = line.lstrip("# ")
        if "id1" in stripped.split(";") or "id1" in stripped.split("\t"):
            header_i = i
            delimiter = ";" if ";" in line else "\t"
            break
    if header_i is None:
        if not combined.strip():
            if not args.keep_intermediate:
                qfile.unlink(missing_ok=True)
                tfile.unlink(missing_ok=True)
            return [], " ".join(cmd), []
        raise RuntimeError("Could not identify IntaRNA CSV header. Check the saved intermediate output. First output:\n" + combined[:2500])
    csv_text = "\n".join(lines[header_i:])
    reader = csv.DictReader(io.StringIO(csv_text), delimiter=delimiter)
    output: list[dict[str, Any]] = []
    for row in reader:
        row = {(k or "").strip().lstrip("# "): (v or "").strip() for k, v in row.items()}
        target_id = row.get("id1", "")
        lnc_id = row.get("id2", "")
        if target_id not in target_by_id or lnc_id not in records_by_id:
            continue
        energy = parse_float(row.get("E"))
        if energy is not None and args.min_energy is not None and energy > args.min_energy:
            continue
        lnc, target = records_by_id[lnc_id], target_by_id[target_id]
        q_start, q_end = parse_int(row.get("start2")), parse_int(row.get("end2"))
        t_start, t_end = parse_int(row.get("start1")), parse_int(row.get("end1"))
        if None in (q_start, q_end, t_start, t_end):
            continue
        q_sub = row.get("subseq2", "")
        t_sub = row.get("subseq1", "")
        seed_start_q, seed_end_q = parse_int(row.get("seedStart2")), parse_int(row.get("seedEnd2"))
        seed_start_t, seed_end_t = parse_int(row.get("seedStart1")), parse_int(row.get("seedEnd1"))
        seed_len = abs(seed_end_q - seed_start_q) + 1 if seed_start_q is not None and seed_end_q is not None else None
        output.append({
            "lncrna_id": lnc.id, "lncrna_description": lnc.description,
            "target_id": target.id, "target_description": target.description,
            "lncrna_start": q_start, "lncrna_end": q_end, "target_start": t_start, "target_end": t_end,
            "lncrna_length": len(lnc.sequence), "target_length": len(target.sequence),
            "lncrna_gc_percent": round(gc_percent(lnc.sequence), 3), "target_gc_percent": round(gc_percent(target.sequence), 3),
            "engine": "IntaRNA", "evidence_type": "thermodynamic_prediction_accessibility_and_seed",
            "target_rank": "", "site_rank": "", "energy_kcal_mol": energy,
            "hybridization_energy_kcal_mol": parse_float(row.get("E_hybrid")), "energy_normalized": parse_float(row.get("E_norm")),
            "accessibility_target_probability": parse_float(row.get("Pu1")), "accessibility_lncrna_probability": parse_float(row.get("Pu2")),
            "interaction_probability": parse_float(row.get("P_E")),
            "seed_start_lncrna": seed_start_q, "seed_end_lncrna": seed_end_q,
            "seed_start_target": seed_start_t, "seed_end_target": seed_end_t,
            "seed_energy_kcal_mol": parse_float(row.get("seedE")), "seed_length": seed_len,
            "canonical_pairs": "", "gu_wobble_pairs": "", "mismatches": "", "binding_score": "",
            "lncrna_interacting_sequence": q_sub, "target_interacting_sequence": t_sub,
            "interaction_structure": row.get("hybridDP", ""), "interaction_pattern": row.get("hybridDB", ""),
            "raw_engine_fields": json.dumps(row, ensure_ascii=False, sort_keys=True)
        })
    if not args.keep_intermediate:
        qfile.unlink(missing_ok=True)
        tfile.unlink(missing_ok=True)
    return output, " ".join(cmd), []


def pair_type(q: str, t: str) -> str:
    pair = (q, t)
    if pair in CANONICAL:
        return "canonical"
    if pair in GU_WOBBLE:
        return "wobble"
    return "mismatch"


def pair_score(q: str, t: str) -> int:
    typ = pair_type(q, t)
    return 2 if typ == "canonical" else (1 if typ == "wobble" else -2)


def gc_counts_for_alignment(q: str, t: str, q_start: int, q_end: int, t_start: int) -> tuple[int, int, int, str]:
    canonical = wobble = mismatch = 0
    pattern_chars: list[str] = []
    # Query is paired antiparallel to target: target left end contacts query's right end.
    q_indices = range(q_end - 1, q_start - 1, -1)
    t_indices = range(t_start, t_start + (q_end - q_start))
    for qi, ti in zip(q_indices, t_indices):
        if ti < 0 or ti >= len(t):
            break
        typ = pair_type(q[qi], t[ti])
        if typ == "canonical":
            canonical += 1; pattern_chars.append("|")
        elif typ == "wobble":
            wobble += 1; pattern_chars.append(":")
        else:
            mismatch += 1; pattern_chars.append(".")
    return canonical, wobble, mismatch, "".join(pattern_chars)


def builtin_pair_task(task: tuple[Any, ...]) -> list[dict[str, Any]]:
    (lnc_id, lnc_desc, q, q_gc, target_id, target_desc, t, t_gc,
     seed_len, max_span, min_pairs, min_score, max_sites, max_occurrences) = task
    if seed_len > len(q) or seed_len > len(t):
        return []
    q_rc = q.translate(COMP)[::-1]
    seed_index: dict[str, list[int]] = defaultdict(list)
    for r in range(0, len(q_rc) - seed_len + 1):
        kmer = q_rc[r:r+seed_len]
        if all(c in "ACGU" for c in kmer):
            if len(seed_index[kmer]) < max_occurrences:
                seed_index[kmer].append(r)
    candidates: dict[tuple[int, int, int, int], dict[str, Any]] = {}
    if not seed_index:
        return []
    for t_start_seed in range(0, len(t) - seed_len + 1):
        target_seed = t[t_start_seed:t_start_seed + seed_len]
        if not all(c in "ACGU" for c in target_seed):
            continue
        q_rc_positions = seed_index.get(target_seed)
        if not q_rc_positions:
            continue
        for rc_start in q_rc_positions:
            q_start_seed = len(q) - (rc_start + seed_len)
            q_end_seed = len(q) - rc_start
            # Extend target-left / query-right and target-right / query-left independently.
            left_len = 0; left_score = 0; running = 0
            for d in range(1, max_span + 1):
                ti = t_start_seed - d
                qi = q_end_seed + d - 1
                if ti < 0 or qi >= len(q): break
                running += pair_score(q[qi], t[ti])
                if running > left_score:
                    left_score, left_len = running, d
            right_len = 0; right_score = 0; running = 0
            for d in range(1, max_span + 1):
                ti = t_start_seed + seed_len + d - 1
                qi = q_start_seed - d
                if ti >= len(t) or qi < 0: break
                running += pair_score(q[qi], t[ti])
                if running > right_score:
                    right_score, right_len = running, d
            q_start = max(0, q_start_seed - right_len)
            q_end = min(len(q), q_end_seed + left_len)
            t_start = max(0, t_start_seed - left_len)
            t_end = min(len(t), t_start_seed + seed_len + right_len)
            # Ensure coordinates represent a contiguous antiparallel alignment of equal lengths.
            span = min(q_end - q_start, t_end - t_start)
            q_start = q_end - span
            t_end = t_start + span
            canonical, gu, mismatches, pattern = gc_counts_for_alignment(q, t, q_start, q_end, t_start)
            actual_score = 2 * canonical + gu - 2 * mismatches
            pair_count = canonical + gu + mismatches
            if pair_count < min_pairs or actual_score < min_score:
                continue
            key = (q_start, q_end, t_start, t_end)
            item = {
                "lncrna_id": lnc_id, "lncrna_description": lnc_desc,
                "target_id": target_id, "target_description": target_desc,
                "lncrna_start": q_start + 1, "lncrna_end": q_end,
                "target_start": t_start + 1, "target_end": t_end,
                "lncrna_length": len(q), "target_length": len(t),
                "lncrna_gc_percent": round(q_gc, 3), "target_gc_percent": round(t_gc, 3),
                "engine": "LncTrail-builtin", "evidence_type": "sequence_heuristic_exact_seed_plus_contiguous_extension",
                "target_rank": "", "site_rank": "", "energy_kcal_mol": "",
                "hybridization_energy_kcal_mol": "", "energy_normalized": "",
                "accessibility_target_probability": "", "accessibility_lncrna_probability": "", "interaction_probability": "",
                "seed_start_lncrna": q_start_seed + 1, "seed_end_lncrna": q_end_seed,
                "seed_start_target": t_start_seed + 1, "seed_end_target": t_start_seed + seed_len,
                "seed_energy_kcal_mol": "", "seed_length": seed_len,
                "canonical_pairs": canonical, "gu_wobble_pairs": gu, "mismatches": mismatches,
                "binding_score": actual_score,
                "lncrna_interacting_sequence": q[q_start:q_end],
                "target_interacting_sequence": t[t_start:t_end],
                "interaction_structure": "", "interaction_pattern": pattern,
                "raw_engine_fields": json.dumps({"seed_start_rc_0based": rc_start, "exact_seed_len": seed_len,
                                                 "extension_limit_per_side": max_span,
                                                 "score_definition": "2*canonical_pairs + GU_wobbles - 2*mismatches"}, sort_keys=True)
            }
            old = candidates.get(key)
            if old is None or actual_score > float(old["binding_score"]):
                candidates[key] = item
    ordered = sorted(candidates.values(), key=lambda r: (-int(r["binding_score"]), -int(r["canonical_pairs"]), r["target_start"], r["lncrna_start"]))
    return ordered[:max_sites]


def run_builtin(lncs: list[SequenceRecord], targets: list[SequenceRecord], args: argparse.Namespace) -> tuple[list[dict[str, Any]], str, list[str]]:
    tasks = []
    for lnc in lncs:
        for target in targets:
            tasks.append((lnc.id, lnc.description, lnc.sequence, gc_percent(lnc.sequence),
                          target.id, target.description, target.sequence, gc_percent(target.sequence),
                          args.seed_length, args.extension, args.min_pairs, args.min_score,
                          args.max_sites_per_pair, args.max_seed_occurrences))
    rows: list[dict[str, Any]] = []
    # ProcessPoolExecutor sidesteps the Python GIL for the pure-Python fallback scan.
    if args.threads > 1 and len(tasks) > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as pool:
            for hits in pool.map(builtin_pair_task, tasks, chunksize=max(1, min(32, len(tasks) // (args.threads * 4) or 1))):
                rows.extend(hits)
    else:
        for task in tasks:
            rows.extend(builtin_pair_task(task))
    notes = ["Built-in mode is a sequence-only seed-and-extend heuristic; scores are not thermodynamic energies or probabilities."]
    if not rows:
        notes.append("No candidates passed the builtin heuristic thresholds. Consider changing --seed-length/--min-pairs/--min-score or using IntaRNA.")
    return rows, "builtin seed-index + exact-seed extension; process workers=" + str(args.threads), notes


def run_builtin_batch(lncs: list[SequenceRecord], targets: list[SequenceRecord],
                      args: argparse.Namespace,
                      pool: concurrent.futures.ProcessPoolExecutor | None = None) -> list[dict[str, Any]]:
    """Run the built-in engine for one target batch, optionally using a persistent worker pool."""
    tasks = []
    for lnc in lncs:
        for target in targets:
            tasks.append((lnc.id, lnc.description, lnc.sequence, gc_percent(lnc.sequence),
                          target.id, target.description, target.sequence, gc_percent(target.sequence),
                          args.seed_length, args.extension, args.min_pairs, args.min_score,
                          args.max_sites_per_pair, args.max_seed_occurrences))
    rows: list[dict[str, Any]] = []
    if pool is not None and tasks:
        chunksize = max(1, min(32, len(tasks) // (args.threads * 4) or 1))
        for hits in pool.map(builtin_pair_task, tasks, chunksize=chunksize):
            rows.extend(hits)
    else:
        for task in tasks:
            rows.extend(builtin_pair_task(task))
    return rows


def infer_target_gene(rec: SequenceRecord) -> str:
    """Infer an HGNC-like symbol from common RefSeq FASTA descriptions.

    This is intentionally conservative: use explicit [gene=...] metadata first,
    otherwise the first short parenthesized token, and finally the accession.
    A separate annotation table is preferable when exact gene mapping is essential.
    """
    import re
    desc = rec.description
    m = re.search(r"\[gene=([^\]]+)\]", desc, flags=re.IGNORECASE)
    if m and m.group(1).strip():
        return m.group(1).strip()
    # Common NCBI RefSeq style: "... transmembrane protein 173 (TMEM173), mRNA"
    for m in re.finditer(r"\(([A-Za-z][A-Za-z0-9.-]{0,24})\)", desc):
        token = m.group(1)
        if token.lower() not in {"homo", "sapiens", "predicted", "xenopus", "putative", "isoform"}:
            return token
    return rec.id.split(".")[0]


def low_complexity_fraction(seq: str, window: int = 10) -> float:
    """Fast low-complexity proxy: fraction of windows with <=2 bases or >=80% one base."""
    seq = ''.join(c for c in seq.upper().replace('T', 'U') if c in 'ACGU')
    if not seq:
        return 1.0
    if len(seq) < window:
        from collections import Counter
        counts = Counter(seq)
        return 1.0 if len(counts) <= 2 or max(counts.values()) / len(seq) >= 0.8 else 0.0
    counts: dict[str, int] = defaultdict(int)
    for c in seq[:window]:
        counts[c] += 1
    low = 1 if len(counts) <= 2 or max(counts.values()) / window >= 0.8 else 0
    total = 1
    for i in range(window, len(seq)):
        old, new = seq[i-window], seq[i]
        counts[old] -= 1
        if counts[old] == 0:
            del counts[old]
        counts[new] += 1
        total += 1
        if len(counts) <= 2 or max(counts.values()) / window >= 0.8:
            low += 1
    return low / total


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def annotate_and_score(rows: list[dict[str, Any]], lncs: list[SequenceRecord],
                       targets: list[SequenceRecord], engine: str) -> None:
    """Add interpretable sequence/seed specificity features and a heuristic priority index."""
    if not rows:
        return
    lnc_by_id = {r.id: r for r in lncs}
    target_by_id = {r.id: r for r in targets}
    genes = {r.id: infer_target_gene(r) for r in targets}
    target_lc = {r.id: low_complexity_fraction(r.sequence) for r in targets}
    lnc_lc = {r.id: low_complexity_fraction(r.sequence) for r in lncs}
    isoforms_per_gene: dict[str, int] = defaultdict(int)
    for target in targets:
        isoforms_per_gene[genes[target.id]] += 1

    # Find seed words observed in predictions, then count how many distinct target genes
    # contain each word anywhere in the target transcriptome. This is a motif-promiscuity
    # penalty, not a statistical p-value.
    motifs: set[str] = set()
    seed_lengths: set[int] = set()
    for row in rows:
        tid = str(row.get("target_id", ""))
        t = target_by_id.get(tid)
        ss, se = parse_int(row.get("seed_start_target")), parse_int(row.get("seed_end_target"))
        if t and ss and se and 1 <= ss <= se <= len(t.sequence):
            motif = t.sequence[ss-1:se]
            if motif and all(x in "ACGU" for x in motif):
                motifs.add(motif); seed_lengths.add(len(motif))
    motif_genes: dict[str, set[str]] = {m: set() for m in motifs}
    if motifs:
        by_len: dict[int, set[str]] = defaultdict(set)
        for motif in motifs:
            by_len[len(motif)].add(motif)
        for target in targets:
            gene = genes[target.id]
            seq = target.sequence
            for k, kset in by_len.items():
                if len(seq) < k:
                    continue
                for i in range(len(seq)-k+1):
                    mer = seq[i:i+k]
                    if mer in kset:
                        motif_genes[mer].add(gene)

    for row in rows:
        tid, qid = str(row.get("target_id", "")), str(row.get("lncrna_id", ""))
        target = target_by_id.get(tid)
        lnc = lnc_by_id.get(qid)
        row["target_gene"] = genes.get(tid, tid)
        row["target_isoform_count"] = isoforms_per_gene.get(row["target_gene"], 1)
        canonical = parse_float(row.get("canonical_pairs")) or 0.0
        gu = parse_float(row.get("gu_wobble_pairs")) or 0.0
        mismatches = parse_float(row.get("mismatches")) or 0.0
        n_pairs = int(canonical + gu + mismatches)
        n_pairs = max(n_pairs, (parse_int(row.get("lncrna_end")) or 0) - (parse_int(row.get("lncrna_start")) or 0) + 1)
        row["alignment_length"] = n_pairs
        row["canonical_fraction"] = round(canonical / n_pairs, 4) if n_pairs else ""
        seed_len = parse_int(row.get("seed_length")) or 0
        ss, se = parse_int(row.get("seed_start_target")), parse_int(row.get("seed_end_target"))
        motif = ""
        if target and ss and se and 1 <= ss <= se <= len(target.sequence):
            motif = target.sequence[ss-1:se]
        gene_freq = len(motif_genes.get(motif, set())) if motif else 0
        row["seed_transcriptome_occurrences"] = gene_freq
        specificity = 1.0 / (1.0 + math.log10(max(1, gene_freq)))
        row["seed_specificity_score"] = round(100 * specificity, 3)
        q_lc = lnc_lc.get(qid, 0.0)
        t_lc = target_lc.get(tid, 0.0)
        row["lncrna_low_complexity_fraction"] = round(q_lc, 4)
        row["target_low_complexity_fraction"] = round(t_lc, 4)
        paired = max(1, n_pairs)
        pair_quality = clamp01((canonical + 0.5 * gu) / paired)
        seed_quality = clamp01(seed_len / 10.0)
        complexity_quality = clamp01(1.0 - max(q_lc, t_lc))
        energy = parse_float(row.get("energy_kcal_mol"))
        if energy is not None:
            energy_quality = clamp01((-energy - 4.0) / 16.0)
            score = 100 * (0.35 * energy_quality + 0.20 * seed_quality +
                           0.20 * pair_quality + 0.15 * specificity + 0.10 * complexity_quality)
        else:
            # Transparent sequence-only ranking index. Not a probability or ΔG.
            score = 100 * (0.30 * pair_quality + 0.25 * seed_quality +
                           0.25 * specificity + 0.20 * complexity_quality)
        row["priority_score"] = round(score, 3)
        row["integrated_priority_score"] = round(score, 3)
        row["rime_max_score"] = ""
        row["rime_windows_above_0_5"] = ""
        row["priority_tier"] = "A" if score >= 70 else ("B" if score >= 50 else "C")


def rank_rows(rows: list[dict[str, Any]], engine: str) -> list[dict[str, Any]]:
    """Sort sites by priority index, falling back to energy/heuristic score for legacy rows."""
    def metric(r: dict[str, Any]) -> float:
        score = parse_float(r.get("integrated_priority_score"))
        if score is None:
            score = parse_float(r.get("priority_score"))
        if score is not None:
            return score
        energy = parse_float(r.get("energy_kcal_mol"))
        if energy is not None:
            return -energy  # more negative energy is better
        binding = parse_float(r.get("binding_score"))
        return binding if binding is not None else -1e9

    rows.sort(key=lambda r: (str(r.get("lncrna_id", "")), -metric(r),
                             str(r.get("target_gene", r.get("target_id", ""))),
                             parse_int(r.get("target_start")) or 0))
    targets_by_lnc: dict[str, dict[str, float]] = defaultdict(dict)
    for r in rows:
        q, tid = str(r.get("lncrna_id", "")), str(r.get("target_id", ""))
        score = metric(r)
        if tid not in targets_by_lnc[q] or score > targets_by_lnc[q][tid]:
            targets_by_lnc[q][tid] = score
    ranks: dict[tuple[str, str], int] = {}
    for q, mapping in targets_by_lnc.items():
        for i, (tid, _) in enumerate(sorted(mapping.items(), key=lambda kv: (-kv[1], kv[0])), 1):
            ranks[(q, tid)] = i
    site_counts: dict[tuple[str, str], int] = defaultdict(int)
    for r in rows:
        q, tid = str(r.get("lncrna_id", "")), str(r.get("target_id", ""))
        r["target_rank"] = ranks.get((q, tid), "")
        pairkey = (q, tid); site_counts[pairkey] += 1
        r["site_rank"] = site_counts[pairkey]
    return rows


def build_gene_summary(rows: list[dict[str, Any]], targets: list[SequenceRecord]) -> list[dict[str, Any]]:
    """Collapse transcript isoforms to one ranked record per lncRNA × inferred target gene."""
    if not rows:
        return []
    targets_by_gene: dict[str, set[str]] = defaultdict(set)
    for t in targets:
        targets_by_gene[infer_target_gene(t)].add(t.id)
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row.get("lncrna_id", "")), str(row.get("target_gene") or row.get("target_id", "")))].append(row)
    summary: list[dict[str, Any]] = []
    for (qid, gene), group in grouped.items():
        best = max(group, key=lambda r: parse_float(r.get("integrated_priority_score")) or 0.0)
        transcripts_hit = sorted({str(r.get("target_id", "")) for r in group})
        target_coords = sorted({(str(r.get("target_id", "")), parse_int(r.get("target_start")), parse_int(r.get("target_end"))) for r in group})
        item = {
            "lncrna_id": qid, "target_gene": gene,
            "best_transcript_id": best.get("target_id", ""),
            "best_transcript_description": best.get("target_description", ""),
            "priority_score": best.get("priority_score", ""),
            "integrated_priority_score": best.get("integrated_priority_score", best.get("priority_score", "")),
            "priority_tier": best.get("priority_tier", ""),
            "engine": best.get("engine", ""),
            "lncrna_start": best.get("lncrna_start", ""), "lncrna_end": best.get("lncrna_end", ""),
            "target_start": best.get("target_start", ""), "target_end": best.get("target_end", ""),
            "energy_kcal_mol": best.get("energy_kcal_mol", ""),
            "binding_score": best.get("binding_score", ""),
            "seed_length": best.get("seed_length", ""),
            "seed_specificity_score": best.get("seed_specificity_score", ""),
            "seed_transcriptome_occurrences": best.get("seed_transcriptome_occurrences", ""),
            "canonical_fraction": best.get("canonical_fraction", ""),
            "rime_max_score": best.get("rime_max_score", ""),
            "rime_windows_above_0_5": best.get("rime_windows_above_0_5", ""),
            "candidate_site_count": len(group),
            "transcripts_with_hits": len(transcripts_hit),
            "target_isoforms_in_input": len(targets_by_gene.get(gene, set())),
            "distinct_transcript_sites": len(target_coords),
            "lncrna_interacting_sequence": best.get("lncrna_interacting_sequence", ""),
            "target_interacting_sequence": best.get("target_interacting_sequence", ""),
            "interaction_pattern": best.get("interaction_pattern", ""),
            "transcript_ids_with_hits": ";".join(transcripts_hit)
        }
        summary.append(item)
    summary.sort(key=lambda r: (str(r["lncrna_id"]), -(parse_float(r.get("integrated_priority_score")) or 0.0), str(r["target_gene"])))
    for qid in sorted({str(r["lncrna_id"]) for r in summary}):
        i = 0
        for item in summary:
            if str(item["lncrna_id"]) == qid:
                i += 1; item["gene_rank"] = i
    return summary


def build_target_scan(rows: list[dict[str, Any]], lncs: list[SequenceRecord],
                     targets: list[SequenceRecord], engine: str,
                     target_universe: list[SequenceRecord] | None = None) -> list[dict[str, Any]]:
    """Return exactly one audit row for every lncRNA × input-target pair.

    If a pair has one or more sites passing the configured predictor filters, retain
    its best-ranked site. Otherwise retain the pair with explicit status and blank
    site-level scores. A missing candidate is not interpreted as a negative interaction.
    """
    def metric(row: dict[str, Any]) -> float:
        for key in ("integrated_priority_score", "priority_score"):
            val = parse_float(row.get(key))
            if val is not None:
                return val
        energy = parse_float(row.get("energy_kcal_mol"))
        if energy is not None:
            return -energy
        binding = parse_float(row.get("binding_score"))
        return binding if binding is not None else -1e12

    best_by_pair: dict[tuple[str, str], dict[str, Any]] = {}
    for row in rows:
        key = (str(row.get("lncrna_id", "")), str(row.get("target_id", "")))
        previous = best_by_pair.get(key)
        if previous is None or metric(row) > metric(previous):
            best_by_pair[key] = row

    universe = target_universe if target_universe is not None else targets
    genes = {t.id: infer_target_gene(t) for t in universe}
    isoforms: dict[str, int] = defaultdict(int)
    for target in universe:
        isoforms[genes[target.id]] += 1
    engine_label = "IntaRNA" if engine == "intarna" else "LncTrail-builtin"
    out: list[dict[str, Any]] = []
    for lnc in lncs:
        for target in targets:
            candidate = best_by_pair.get((lnc.id, target.id))
            if candidate is not None:
                item = dict(candidate)
                item["prediction_status"] = "candidate_found"
                item["target_gene"] = genes[target.id]
                item["target_isoform_count"] = isoforms[genes[target.id]]
                out.append(item)
                continue
            item = {column: "" for column in TSV_COLUMNS}
            item.update({
                "lncrna_id": lnc.id,
                "lncrna_description": lnc.description,
                "target_id": target.id,
                "target_description": target.description,
                "lncrna_length": len(lnc.sequence),
                "target_length": len(target.sequence),
                "lncrna_gc_percent": round(gc_percent(lnc.sequence), 3),
                "target_gc_percent": round(gc_percent(target.sequence), 3),
                "engine": engine_label,
                "evidence_type": "no_candidate_passed_configured_filters",
                "prediction_status": "no_candidate_passed_filters",
                "target_gene": genes[target.id],
                "target_isoform_count": isoforms[genes[target.id]],
                "rime_max_score": "",
                "rime_windows_above_0_5": "",
            })
            out.append(item)
    return out


def write_dict_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str] | None = None) -> None:
    if columns is None:
        cols: list[str] = []
        for row in rows:
            for k in row:
                if k not in cols:
                    cols.append(k)
        columns = cols or ["no_results"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") if row.get(k, "") is not None else "" for k in columns})


def run_rime_shortlist(lncs: list[SequenceRecord], targets: list[SequenceRecord], rows: list[dict[str, Any]],
                       args: argparse.Namespace, outdir: Path) -> tuple[list[dict[str, Any]], list[str]]:
    """Optionally run the published RIME deep-learning model on a prefiltered shortlist.

    RIME installation uses two named conda environments (download_embeddings and rime),
    a checkpoint, BEDTools, and a large Nucleotide Transformer embedding model. Only the
    top target-gene representatives are passed to prevent an accidental all-transcript run.
    """
    if not args.rime_repo:
        return rows, []
    repo = args.rime_repo.resolve()
    wrapper = repo / "rime_infer.py"
    if not wrapper.exists():
        raise RuntimeError(f"--rime-repo must point to a RIME repository containing rime_infer.py: {repo}")
    if not args.rime_bedtools_path or not Path(args.rime_bedtools_path).exists():
        raise RuntimeError("RIME requires a valid --rime-bedtools-path (the bedtools executable).")
    if not rows:
        return rows, ["RIME skipped because no candidate rows passed primary filters."]
    summaries = build_gene_summary(rows, targets)
    selected_genes: list[str] = []
    selected_ids: list[str] = []
    targets_by_id = {t.id: t for t in targets}
    for rec in summaries:
        gene = str(rec["target_gene"])
        tid = str(rec["best_transcript_id"])
        if gene not in selected_genes and tid in targets_by_id:
            selected_genes.append(gene); selected_ids.append(tid)
        if len(selected_ids) >= args.rime_top_targets:
            break
    if not selected_ids:
        return rows, ["RIME skipped: no target transcripts were selected for the shortlist."]
    rime_root = outdir / "rime_shortlist"
    input_dir = rime_root / "input"
    output_dir = rime_root / "output"
    input_dir.mkdir(parents=True, exist_ok=True); output_dir.mkdir(parents=True, exist_ok=True)
    write_normalized_fasta(lncs, input_dir / "lncrna.fa")
    shortlist = [targets_by_id[x] for x in selected_ids]
    write_normalized_fasta(shortlist, input_dir / "targets.fa")
    cmd = [sys.executable, str(wrapper), "--input_dir", str(input_dir), "--query", "lncrna.fa",
           "--target", "targets.fa", "--output_dir", str(output_dir),
           "--bedtools_path", str(Path(args.rime_bedtools_path).resolve()), "--model", "RIMEfull"]
    log = subprocess.run(cmd, capture_output=True, text=True, check=False, cwd=str(repo))
    (rime_root / "rime_stdout.log").write_text(log.stdout, encoding="utf-8")
    (rime_root / "rime_stderr.log").write_text(log.stderr, encoding="utf-8")
    if log.returncode != 0:
        raise RuntimeError("RIME failed. See " + str(rime_root / "rime_stderr.log") + " and ensure the RIME conda environments/checkpoints are installed.")
    bedpe = output_dir / "output_table.bedpe"
    if not bedpe.exists():
        matches = list(output_dir.rglob("output_table.bedpe"))
        if matches: bedpe = matches[0]
    if not bedpe.exists():
        raise RuntimeError(f"RIME completed but output_table.bedpe was not found under {output_dir}")
    with bedpe.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        windows = list(reader)
    pair_windows: dict[tuple[str, str], list[float]] = defaultdict(list)
    for w in windows:
        # RIME BEDPE uses sequence IDs for window_id1 (query RNA) and window_id2 (target RNA).
        qid, tid = str(w.get("window_id1", "")), str(w.get("window_id2", ""))
        val = parse_float(w.get("RIME_score"))
        if qid and tid and val is not None:
            pair_windows[(qid, tid)].append(val)
    write_dict_tsv(rime_root / "rime_window_scores.tsv", windows)
    for row in rows:
        pair_key = (str(row.get("lncrna_id", "")), str(row.get("target_id", "")))
        vals = pair_windows.get(pair_key, [])
        if vals:
            rmax = max(vals)
            row["rime_max_score"] = round(rmax, 6)
            row["rime_windows_above_0_5"] = sum(v > 0.5 for v in vals)
            # RIME is a complementary deep model; this is a ranking heuristic, not a calibrated probability.
            row["integrated_priority_score"] = round(0.70 * (parse_float(row.get("priority_score")) or 0.0) + 30.0 * rmax, 3)
            val = row["integrated_priority_score"]
            row["priority_tier"] = "A" if val >= 70 else ("B" if val >= 50 else "C")
    return rows, [f"RIME deep-learning model evaluated {len(selected_ids)} target transcript representatives; its window scores were added as complementary evidence."]


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") if row.get(k, "") is not None else "" for k in TSV_COLUMNS})


def initialize_tsv(path: Path) -> None:
    """Create a valid empty TSV immediately so it is visible during analysis."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        handle.flush()
        os.fsync(handle.fileno())


def append_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Append a chunk to an already initialized TSV and flush it for live viewing."""
    if not rows:
        return
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        for row in rows:
            writer.writerow({k: row.get(k, "") if row.get(k, "") is not None else "" for k in TSV_COLUMNS})
        handle.flush()
        os.fsync(handle.fileno())


def atomic_write_json(path: Path, data: dict[str, Any]) -> None:
    """Write metadata atomically so readers never see half-written JSON."""
    temp_path = path.with_name(path.name + ".tmp")
    temp_path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    temp_path.replace(path)


def write_progress_line(path: Path, message: str) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def make_html(path: Path, scan_name: str, sites_name: str, rows: list[dict[str, Any]], metadata: dict[str, Any], html_rows: int) -> None:
    top = sorted(rows, key=lambda r: -(parse_float(r.get("integrated_priority_score", r.get("priority_score"))) or 0.0))[:html_rows]
    def fmt(v: Any) -> str:
        if v is None or v == "": return "—"
        if isinstance(v, float): return f"{v:.4g}"
        return str(v)
    heads = ["lncRNA", "Target gene", "Target transcript", "Priority index", "lncRNA site", "Target site", "Engine", "ΔG (kcal/mol)", "Seed nt", "RIME max score", "Interaction pattern"]
    body_rows = []
    for r in top:
        cells = [r.get("lncrna_id", ""), r.get("target_gene", ""), r.get("target_id", ""),
                 fmt(r.get("integrated_priority_score", r.get("priority_score"))),
                 f"{fmt(r.get('lncrna_start'))}–{fmt(r.get('lncrna_end'))}",
                 f"{fmt(r.get('target_start'))}–{fmt(r.get('target_end'))}",
                 r.get("engine", ""), fmt(r.get("energy_kcal_mol")), fmt(r.get("seed_length")),
                 fmt(r.get("rime_max_score")), r.get("interaction_pattern") or r.get("interaction_structure") or ""]
        body_rows.append("<tr>" + "".join(f"<td>{html.escape(fmt(c))}</td>" for c in cells) + "</tr>")
    warnings_html = "".join(f"<li>{html.escape(str(w))}</li>" for w in metadata.get("warnings", [])) or "<li>None recorded.</li>"
    engine = html.escape(str(metadata.get("engine_used", "unknown")))
    date = html.escape(str(metadata.get("started_at", "")))
    duration = metadata.get("runtime_seconds", 0)
    n_lnc = metadata.get("lncrna_count", 0); n_targets = metadata.get("target_count", 0)
    n_pairs = metadata.get("pair_count", 0); n_rows = len(rows)
    method_explanation = (
        "IntaRNA rows report the program's interaction energy and available seed/accessibility values. "
        "The built-in mode reports a sequence-only heuristic. The priority index combines pairing, seed length, seed promiscuity and low-complexity penalties; it is not a probability or ΔG. "
        "RIME (when enabled) contributes a published deep-learning window score. Scores are for within-run prioritization, not experimentally validated target probabilities."
    )
    doc = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>LncTrail results</title><style>
:root{{--ink:#172033;--muted:#667085;--line:#dce2ea;--paper:#f5f7fa;--accent:#176b87;--warning:#fff4dc}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}}
header{{background:#fff;border-bottom:1px solid var(--line);padding:28px max(22px,calc((100vw - 1400px)/2))}}h1{{margin:0;font-size:28px;letter-spacing:-.4px}}.sub{{color:var(--muted);margin-top:6px}}main{{max-width:1400px;margin:24px auto;padding:0 20px}}.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:18px 0}}.stat,.panel{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:16px}}.stat b{{display:block;font-size:26px}}.stat span{{color:var(--muted);font-size:13px}}.note{{background:var(--warning);border:1px solid #efd79c;border-radius:10px;padding:14px 16px;margin:16px 0}}.panel h2{{font-size:18px;margin:0 0 12px}}.tabletools{{display:flex;gap:12px;align-items:center;justify-content:space-between;flex-wrap:wrap;margin:14px 0}}input{{padding:10px 12px;border:1px solid var(--line);border-radius:8px;min-width:260px;max-width:100%;font:inherit}}a{{color:var(--accent)}}.tablewrap{{overflow:auto;border:1px solid var(--line);border-radius:10px;background:#fff}}table{{border-collapse:collapse;min-width:100%;font-size:13px}}th,td{{text-align:left;padding:10px 12px;border-bottom:1px solid #edf0f4;white-space:nowrap;max-width:360px}}th{{background:#edf4f7;position:sticky;top:0}}tr:hover td{{background:#f8fbfc}}footer{{color:var(--muted);font-size:12px;margin:24px 0 38px}}code{{font-family:ui-monospace,SFMono-Regular,Consolas,monospace}}@media(max-width:600px){{h1{{font-size:23px}}}}
</style></head><body><header><h1>LncTrail</h1><div class="sub">lncRNA–transcript interaction candidate discovery report</div></header>
<main><div class="stats"><div class="stat"><b>{n_lnc:,}</b><span>lncRNA sequences</span></div><div class="stat"><b>{n_targets:,}</b><span>target sequences</span></div><div class="stat"><b>{n_pairs:,}</b><span>query–target pairs screened</span></div><div class="stat"><b>{n_rows:,}</b><span>candidate interaction sites</span></div><div class="stat"><b>{metadata.get("target_pairs_with_candidates", 0):,}</b><span>query–target pairs with candidate(s)</span></div><div class="stat"><b>{duration:.2f}s</b><span>runtime</span></div><div class="stat"><b>{engine}</b><span>prediction engine</span></div></div>
<div class="note"><b>Interpretation and limitations.</b> {html.escape(method_explanation)} Target rankings are within each lncRNA and engine, not estimates of clinical or biological certainty. FASTA-only pairing analysis will not detect all cis-regulatory, protein-mediated, DNA-triplex, or miRNA-sponge mechanisms.</div>
<section class="panel"><h2>Run details</h2><p>Started: {date}<br>Parameters: <code>{html.escape(json.dumps(metadata.get('parameters', {}), ensure_ascii=False, sort_keys=True))}</code></p><p><a href="{html.escape(scan_name)}">Complete target-scan TSV</a> (one row per lncRNA × target transcript, including pairs with no candidate passing filters) ·
<a href="{html.escape(sites_name)}">All retained candidate-site rows</a> ·
<a href="lnctrail_ranked_targets.tsv">All ranked target genes</a> ·
<a href="lnctrail_top_targets.tsv">Top target genes</a> ·
<a href="lnctrail_top_sites.tsv">Top-ranked sites</a></p></section>
<section class="panel" style="margin-top:16px"><h2>Input quality notes</h2><ul>{warnings_html}</ul></section>
<section class="panel" style="margin-top:16px"><h2>Candidate interaction sites</h2><div class="tabletools"><input id="search" type="search" placeholder="Filter by lncRNA, target, sequence, engine…"><span class="sub">Showing up to {min(html_rows, n_rows):,} of {n_rows:,} rows</span></div>
<div class="tablewrap"><table id="results"><thead><tr>{''.join('<th>'+html.escape(h)+'</th>' for h in heads)}</tr></thead><tbody>{''.join(body_rows) if body_rows else '<tr><td colspan="11">No candidates passed the selected filters.</td></tr>'}</tbody></table></div></section>
<footer>Generated by LncTrail v{VERSION}. Start and end coordinates are 1-based inclusive relative to the submitted FASTA sequence. Verify genome/transcript accession versions and validate priority candidates experimentally. Run started {date}.</footer></main>
<script>const s=document.getElementById('search');s.addEventListener('input',()=>{{const q=s.value.toLowerCase();document.querySelectorAll('#results tbody tr').forEach(r=>{{r.style.display=r.innerText.toLowerCase().includes(q)?'':'none'}})}});</script></body></html>'''
    path.write_text(doc, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="lnctrail",
        description="LncTrail: lncRNA–transcript RNA-RNA interaction candidate discovery from FASTA inputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--lncrna", required=True, type=Path, help="FASTA of lncRNA query sequences (.gz supported)")
    parser.add_argument("--targets", required=True, type=Path, help="FASTA of target RNA transcripts (prefer spliced transcript FASTA; .gz supported)")
    parser.add_argument("--outdir", type=Path, default=Path("lnctrail_results"), help="Output directory")
    parser.add_argument("--engine", choices=["auto", "intarna", "builtin"], default="auto", help="auto uses IntaRNA if available, otherwise the transparent built-in heuristic")
    parser.add_argument("--threads", type=int, default=max(1, min(8, os.cpu_count() or 1)), help="CPU workers; IntaRNA uses native threading, builtin uses worker processes")
    parser.add_argument("--seed-length", type=int, default=None, help="Minimum seed length; auto default is 12 nt in builtin mode and 7 nt with IntaRNA")
    parser.add_argument("--max-sites-per-pair", type=int, default=None, help="Maximum sites per lncRNA–transcript pair; default is 1 for a focused first-pass ranking")
    parser.add_argument("--max-pairs", type=int, default=250000, help="Safety limit on all-vs-all lncRNA x target pairs")
    parser.add_argument("--batch-size", type=int, default=250, help="Target transcripts processed per batch; outputs are flushed after each batch (smaller is more responsive but adds overhead)")
    parser.add_argument("--force-large-search", action="store_true", help="Proceed if the number of all-vs-all pairs exceeds --max-pairs")
    parser.add_argument("--min-energy", type=float, default=None, help="IntaRNA post-filter: retain interactions with energy <= this kcal/mol value (e.g. -10). Ignored in builtin mode.")
    parser.add_argument("--min-score", type=int, default=26, help="Builtin-only minimum heuristic score: 2*canonical pairs + GU wobble pairs - 2*mismatches")
    parser.add_argument("--min-pairs", type=int, default=12, help="Builtin-only minimum contiguous paired positions including mismatches")
    parser.add_argument("--extension", type=int, default=18, help="Builtin-only maximum contiguous extension length on each side of an exact seed")
    parser.add_argument("--max-seed-occurrences", type=int, default=200, help="Builtin-only cap on occurrences indexed per query seed to control low-complexity explosion")
    parser.add_argument("--top-targets", type=int, default=200, help="Maximum rows in lnctrail_top_targets.tsv")
    parser.add_argument("--top-sites", type=int, default=1000, help="Maximum rows in lnctrail_top_sites.tsv")
    parser.add_argument("--target-biotype", choices=["all", "protein-coding"], default="all", help="Optionally restrict RefSeq target inputs to NM_/XM_ mRNA accessions")
    parser.add_argument("--rime-repo", type=Path, default=None, help="Optional local RIME repository for deep-learning rescoring of a shortlist")
    parser.add_argument("--rime-top-targets", type=int, default=30, help="Number of top distinct target genes passed to RIME when --rime-repo is set")
    parser.add_argument("--rime-bedtools-path", type=str, default=None, help="Path to the bedtools executable required by the RIME wrapper")
    parser.add_argument("--allow-overlapping-sites", action="store_true", help="Permit overlapping predictions in IntaRNA")
    parser.add_argument("--temperature", type=float, default=37.0, help="IntaRNA temperature in Celsius")
    parser.add_argument("--html-rows", type=int, default=1000, help="Maximum candidate rows embedded in HTML; TSV always contains all rows")
    parser.add_argument("--keep-intermediate", action="store_true", help="Keep per-batch normalized FASTA and raw IntaRNA outputs in outdir/intermediate")
    parser.add_argument("--version", action="version", version=f"LncTrail {VERSION}")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if args.seed_length is not None and not 4 <= args.seed_length <= 30:
        parser.error("--seed-length must be between 4 and 30")
    if args.max_pairs < 1 or args.min_pairs < 1 or args.extension < 0:
        parser.error("pair/score thresholds must be positive; --extension must be >= 0")
    if args.batch_size < 1:
        parser.error("--batch-size must be >= 1")
    if args.top_targets < 1 or args.top_sites < 0 or args.rime_top_targets < 1:
        parser.error("--top-targets and --rime-top-targets must be >= 1; --top-sites must be >= 0")
    if args.max_sites_per_pair is not None and args.max_sites_per_pair < 1:
        parser.error("--max-sites-per-pair must be >= 1")
    if args.html_rows < 0:
        parser.error("--html-rows must be >= 0")

    args.outdir.mkdir(parents=True, exist_ok=True)
    intermediate = args.outdir / "intermediate"
    intermediate.mkdir(exist_ok=True)
    start = time.time()
    started_at = dt.datetime.now().astimezone().isoformat(timespec="seconds")
    tsv_path = args.outdir / "lnctrail_results.tsv"
    candidate_sites_path = args.outdir / "lnctrail_candidate_sites.tsv"
    html_path = args.outdir / "lnctrail_report.html"
    meta_path = args.outdir / "run_metadata.json"
    ranked_targets_path = args.outdir / "lnctrail_ranked_targets.tsv"
    top_targets_path = args.outdir / "lnctrail_top_targets.tsv"
    top_sites_path = args.outdir / "lnctrail_top_sites.tsv"
    run_log_path = args.outdir / "run.log"
    progress_meta: dict[str, Any] | None = None
    warnings: list[str] = []
    try:
        lncs, w1 = read_fasta(args.lncrna, "lncRNA")
        targets, w2 = read_fasta(args.targets, "target")
        warnings.extend(w1); warnings.extend(w2)
        if args.target_biotype == "protein-coding":
            before = len(targets)
            targets = [t for t in targets if t.id.startswith(("NM_", "XM_")) and "non-coding RNA" not in t.description.lower()]
            warnings.append(f"Target biotype filter protein-coding retained {len(targets):,}/{before:,} FASTA records (RefSeq NM_/XM_ heuristic). Verify transcript annotations for your assembly.")
            if not targets:
                raise ValueError("--target-biotype protein-coding removed all target records; inspect IDs or use --target-biotype all.")
        pair_count = len(lncs) * len(targets)
        if pair_count > args.max_pairs and not args.force_large_search:
            raise ValueError(f"All-vs-all input contains {pair_count:,} query-target pairs, above --max-pairs={args.max_pairs:,}. Use a smaller target set / indexed search, increase --max-pairs, or explicitly pass --force-large-search.")
        if any(r.original_had_t for r in targets):
            warnings.append("One or more target records contained T and were converted to U. Confirm that the input represents RNA transcripts, not unspliced genomic DNA.")
        warnings.append("Target FASTA should preferably contain spliced, transcript-oriented RNA sequences (including relevant UTRs). Genomic DNA or intron-containing loci are not an appropriate substitute for transcript FASTA in RNA-RNA pairing prediction.")
        warnings.append("FASTA-only predictions do not establish expression co-occurrence, subcellular co-localization, direct in-cell binding, causal regulation, or biological validation.")
        engine = first_available_engine(args.engine)
        if args.seed_length is None:
            args.seed_length = 12 if engine == "builtin" else 7
        if args.max_sites_per_pair is None:
            args.max_sites_per_pair = 1
        if engine == "builtin" and args.min_energy is not None:
            warnings.append("--min-energy was supplied but is not applicable to the builtin heuristic; use --min-score instead.")

        # These files are deliberately created before the expensive prediction stage.
        # The main TSV grows one completed target batch at a time.
        initialize_tsv(tsv_path)
        initialize_tsv(candidate_sites_path)
        progress_meta = {
            "tool": "LncTrail", "version": VERSION, "status": "running", "started_at": started_at,
            "runtime_seconds": 0.0, "engine_requested": args.engine, "engine_used": engine,
            "lncrna_file": str(args.lncrna.resolve()), "target_file": str(args.targets.resolve()),
            "lncrna_sha256": sha256_file(args.lncrna), "target_sha256": sha256_file(args.targets),
            "lncrna_count": len(lncs), "target_count": len(targets), "pair_count": pair_count,
            "processed_targets": 0, "processed_pairs": 0, "progress_fraction": 0.0,
            "candidate_sites_so_far": 0, "target_pairs_with_candidates_so_far": 0,
            "threads": args.threads, "batch_size_targets": args.batch_size,
            "parameters": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
            "warnings": warnings,
            "note": "status=running means TSV rows are provisional and contain completed batches only. After success, the main TSV is rewritten with all input pairs and final global ranking fields.",
        }
        atomic_write_json(meta_path, progress_meta)
        run_log_path.write_text(
            f"LncTrail {VERSION}\nStarted: {started_at}\nStatus: running\nEngine: {engine}\n"
            f"Input pairs: {pair_count}\nTarget batch size: {args.batch_size}\n",
            encoding="utf-8")

        print(f"LncTrail {VERSION} | {len(lncs)} lncRNA(s) x {len(targets)} target(s) = {pair_count:,} pairs", flush=True)
        print(f"Engine: {engine} | Workers: {args.threads} | Target batch size: {args.batch_size}", flush=True)
        print(f"Live target scan: {tsv_path}", flush=True)
        print("Progress is flushed after every completed target batch; final ranking tables and HTML are written at completion.", flush=True)

        rows: list[dict[str, Any]] = []
        commands: list[str] = []
        processed_targets = 0
        progress_sites = 0
        progress_pairs_with_hits: set[tuple[str, str]] = set()
        batch_count = (len(targets) + args.batch_size - 1) // args.batch_size
        pool: concurrent.futures.ProcessPoolExecutor | None = None
        if engine == "builtin" and args.threads > 1 and pair_count > 1:
            pool = concurrent.futures.ProcessPoolExecutor(max_workers=args.threads)
        try:
            for batch_no, offset in enumerate(range(0, len(targets), args.batch_size), 1):
                batch_targets = targets[offset: offset + args.batch_size]
                if engine == "intarna":
                    batch_rows, command, engine_notes = run_intarna(lncs, batch_targets, args, intermediate, batch_index=batch_no)
                    commands.append(command)
                    warnings.extend(engine_notes)
                else:
                    batch_rows = run_builtin_batch(lncs, batch_targets, args, pool=pool)
                # Do not run expensive transcriptome-wide ranking for every batch.
                # The live TSV provides raw predictor evidence; final global scores/ranks are computed once at the end.
                batch_scan = build_target_scan(batch_rows, lncs, batch_targets, engine, target_universe=targets)
                append_tsv(tsv_path, batch_scan)
                append_tsv(candidate_sites_path, batch_rows)
                rows.extend(batch_rows)
                processed_targets += len(batch_targets)
                progress_sites += len(batch_rows)
                progress_pairs_with_hits.update((str(r.get("lncrna_id", "")), str(r.get("target_id", ""))) for r in batch_rows)
                processed_pairs = processed_targets * len(lncs)
                fraction = processed_targets / len(targets) if targets else 1.0
                elapsed = time.time() - start
                msg = (f"Progress: batch {batch_no}/{batch_count}; targets {processed_targets:,}/{len(targets):,} "
                       f"({fraction:.1%}); pairs {processed_pairs:,}/{pair_count:,}; candidate sites {progress_sites:,}; "
                       f"elapsed {elapsed:.1f}s")
                print(msg, flush=True)
                write_progress_line(run_log_path, msg)
                progress_meta.update({
                    "processed_targets": processed_targets, "processed_pairs": processed_pairs,
                    "progress_fraction": round(fraction, 6), "candidate_sites_so_far": progress_sites,
                    "target_pairs_with_candidates_so_far": len(progress_pairs_with_hits),
                    "runtime_seconds": round(elapsed, 3), "last_completed_batch": batch_no,
                })
                atomic_write_json(meta_path, progress_meta)
        finally:
            if pool is not None:
                pool.shutdown(wait=True, cancel_futures=False)

        if engine == "builtin":
            commands.append("builtin seed-index + exact-seed extension; process workers=" + str(args.threads))
            warnings.append("Built-in mode is a sequence-only seed-and-extend heuristic; scores are not thermodynamic energies or probabilities.")
            if not rows:
                warnings.append("No candidates passed the builtin heuristic thresholds. Consider changing --seed-length/--min-pairs/--min-score or using IntaRNA.")
        if engine == "intarna" and not rows:
            warnings.append("IntaRNA completed without any predictions passing the configured filters.")
        # The expensive pair scan is finished; tell the user before global ranking/report finalization.
        progress_meta.update({
            "status": "finalizing",
            "runtime_seconds": round(time.time() - start, 3),
            "processed_targets": len(targets),
            "processed_pairs": pair_count,
            "progress_fraction": 1.0,
            "note": "All target batches are complete. Global scoring, ranking, and report generation are in progress.",
        })
        atomic_write_json(meta_path, progress_meta)
        write_progress_line(run_log_path, "All target batches complete; finalizing global scores, rankings and reports.")
        print("All target batches complete; finalizing global scores, rankings and reports...", flush=True)
        # Recalculate and rank globally once all batches are collected; this finalizes the live TSV.
        annotate_and_score(rows, lncs, targets, engine)
        if args.rime_repo:
            rows, rime_notes = run_rime_shortlist(lncs, targets, rows, args, args.outdir)
            warnings.extend(rime_notes)
        rows = rank_rows(rows, engine)
        global_rows = sorted(rows, key=lambda r: -(parse_float(r.get("integrated_priority_score", r.get("priority_score"))) or 0.0))
        for rank_i, row in enumerate(global_rows, 1):
            row["global_priority_rank"] = rank_i
        gene_rows = build_gene_summary(rows, targets)
        target_scan_rows = build_target_scan(rows, lncs, targets, engine, target_universe=targets)
        pair_hits = {(str(r.get("lncrna_id", "")), str(r.get("target_id", ""))) for r in rows}
        runtime = time.time() - start
        final_command = "\n--- batch ---\n".join(commands)
        metadata: dict[str, Any] = {
            "tool": "LncTrail", "version": VERSION, "status": "completed", "started_at": started_at,
            "completed_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
            "runtime_seconds": round(runtime, 3), "engine_requested": args.engine, "engine_used": engine,
            "lncrna_file": str(args.lncrna.resolve()), "target_file": str(args.targets.resolve()),
            "lncrna_sha256": progress_meta["lncrna_sha256"], "target_sha256": progress_meta["target_sha256"],
            "lncrna_count": len(lncs), "target_count": len(targets), "pair_count": pair_count,
            "candidate_sites": len(rows), "ranked_target_gene_pairs": len(gene_rows),
            "target_scan_rows": len(target_scan_rows), "target_pairs_with_candidates": len(pair_hits),
            "target_pairs_without_candidates": len(target_scan_rows) - len(pair_hits),
            "unique_target_gene_symbols_in_input": len({infer_target_gene(t) for t in targets}),
            "processed_targets": len(targets), "processed_pairs": pair_count, "progress_fraction": 1.0,
            "threads": args.threads, "batch_size_targets": args.batch_size,
            "engine_command": final_command,
            "engine_commands_by_batch": commands,
            "ranking_note": "priority_score is a heuristic ranking index, not a calibrated probability; RIME is optional complementary ML evidence.",
            "parameters": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
            "warnings": warnings,
            "scientific_scope": "Sequence-based candidate prediction for RNA-RNA interaction between the supplied lncRNA and target transcript sequences; not a comprehensive predictor of all lncRNA mechanisms."
        }
        # The live main TSV is replaced with the finalized full input-pair table (same filename).
        write_tsv(tsv_path, target_scan_rows)
        write_tsv(candidate_sites_path, rows)
        gene_columns = list(gene_rows[0].keys()) if gene_rows else ["lncrna_id", "target_gene", "gene_rank", "integrated_priority_score"]
        write_dict_tsv(ranked_targets_path, gene_rows, gene_columns)
        write_dict_tsv(top_targets_path, gene_rows[:args.top_targets], gene_columns)
        write_tsv(top_sites_path, global_rows[:args.top_sites])
        make_html(html_path, tsv_path.name, candidate_sites_path.name, rows, metadata, args.html_rows)
        atomic_write_json(meta_path, metadata)
        write_progress_line(run_log_path, f"Completed successfully in {runtime:.2f}s. Candidate sites: {len(rows):,}; input pairs: {len(target_scan_rows):,}.")
        if not args.keep_intermediate:
            for pth in intermediate.iterdir():
                if pth.is_file() and (pth.name.startswith(("lncrna.normalized", "targets.normalized")) or pth.name.startswith("intarna.batch_") or pth.name.startswith("intarna.stdout") or pth.name.startswith("intarna.stderr")):
                    pth.unlink()
        print(f"Done in {runtime:.2f}s. Candidate sites: {len(rows):,}; target-scan rows: {len(target_scan_rows):,}; pairs with candidates: {len(pair_hits):,}")
        print(f"Full target scan: {tsv_path}")
        print(f"Candidate site details: {candidate_sites_path}")
        print(f"HTML: {html_path}")
        print(f"Meta/progress: {meta_path}")
        print(f"Ranked target genes: {ranked_targets_path}")
        print(f"Top target genes: {top_targets_path}")
        print(f"Top sites: {top_sites_path}")
        if args.rime_repo:
            print(f"RIME outputs: {args.outdir / 'rime_shortlist'}")
        return 0
    except Exception as exc:
        if progress_meta is not None:
            progress_meta["status"] = "failed"
            progress_meta["failed_at"] = dt.datetime.now().astimezone().isoformat(timespec="seconds")
            progress_meta["error"] = str(exc)
            progress_meta["runtime_seconds"] = round(time.time() - start, 3)
            try:
                atomic_write_json(meta_path, progress_meta)
                write_progress_line(run_log_path, f"FAILED after {progress_meta['runtime_seconds']:.2f}s: {exc}")
            except Exception:
                pass
        print(f"LncTrail error: {exc}", file=sys.stderr)
        return 2

if __name__ == "__main__":
    raise SystemExit(main())
