#!/usr/bin/env python3
"""Plot reference-anchored gene and assay-site conservation on a phylogeny.

Search genome assemblies with BLASTn, evaluate all qualifying HSP chains,
group alternative chains into putative loci, and display the best three-site
assay match per genome. Bootstrap support is preserved by unrooted taxon split.
Minimum-variance rooting is the default; --rooting midpoint or input overrides it.
Site detection is a sequence comparison, not a prediction of assay activity;
Four-base PAMs are compared with TTTN in the figure and genome-site tables.
PAM differences are descriptive and do not affect locus selection.
The gRNA/crRNA FASTA record defaults to PAM+guide (4-base 5-prime PAM).
Use --guide-input guide-only if the record already excludes its PAM.
Dependencies: Biopython, Matplotlib, BLASTn,
and MAFFT. Use --help for inputs and thresholds.
"""

import argparse
import csv
from copy import copy
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.Phylo.Newick import Clade

VERSION = "1.6.2"
REF = "__GENE_REFERENCE__"
DNA = set("ACGTRYSWKMBDHVN")
IUPAC = dict(
    zip(
        "ACGTRYSWKMBDHVN",
        [
            "A",
            "C",
            "G",
            "T",
            "AG",
            "CT",
            "GC",
            "AT",
            "GT",
            "AC",
            "CGT",
            "AGT",
            "ACT",
            "ACG",
            "ACGT",
        ],
    )
)
COMPONENTS = [
    ("fwd", "fwdPrimer", "Forward RPA", "fwdRPA_aligned.fasta"),
    ("rev", "revPrimer", "Reverse RPA", "revRPA_aligned.fasta"),
    ("gRNA", "gRNA", "crRNA", "gRNA_aligned.fasta"),
]
REPORT_COMPONENTS = [
    COMPONENTS[0],
    COMPONENTS[1],
    ("PAM", "PAM", "PAM", "PAM_aligned.fasta"),
    COMPONENTS[2],
]


# Input, output and sequence helpers


def log(msg):
    print(msg, flush=True)


def rc(s):
    return str(Seq(s).reverse_complement())


def read_fasta(path, aligned=False):
    """Validate all records; IDs are the first whitespace-delimited header token."""
    records = {}
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt") as handle:
        for r in SeqIO.parse(handle, "fasta"):
            if r.id in records:
                raise ValueError(f"Duplicate FASTA identifier {r.id!r} in {path}")
            s = str(r.seq).upper()
            if not s:
                raise ValueError(f"Empty sequence: {r.id} in {path}")
            bad = set(s) - DNA - ({"-"} if aligned else set())
            if bad:
                raise ValueError(
                    f"Invalid DNA characters {sorted(bad)}: {r.id} in {path}; use DNA (T), not RNA (U)."
                )
            records[r.id] = s
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    if aligned and len({len(s) for s in records.values()}) != 1:
        raise ValueError(f"Inconsistent alignment lengths in {path}")
    return records


def fasta(path, records):
    with open(path, "w") as out:
        for name, seq in records.items():
            if not seq:
                raise ValueError(f"Refusing to write empty sequence: {name}")
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")


def table(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run(cmd, stdout_path, stderr_path):
    with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
        r = subprocess.run([str(x) for x in cmd], stdout=out, stderr=err)
    if r.returncode:
        raise RuntimeError(
            f"Command failed ({r.returncode}): {cmd[0]}; see {stderr_path}"
        )


def resolve_threads(value):
    if value is not None:
        if value < 1:
            raise ValueError("--threads must be a positive integer.")
        return value
    if not sys.stdin.isatty():
        raise ValueError(
            "No interactive terminal available; supply --threads N to specify the CPU budget."
        )
    available = (
        len(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else (os.cpu_count() or 1)
    )
    while True:
        try:
            answer = input(
                f"How many CPU cores/threads can this run use? (available: {available}): "
            ).strip()
        except EOFError as exc:
            raise ValueError("No thread count entered; supply --threads N.") from exc
        if answer.isdecimal() and int(answer) > 0:
            return int(answer)
        log("Please enter a positive whole number of CPU threads.")


def mismatch_positions(ref, seq):
    # A degenerate oligo accepts its concrete members. Ambiguous genome bases
    # are unresolved, never silently considered matches (checked by caller).
    return [i + 1 for i, (a, b) in enumerate(zip(ref, seq)) if b not in IUPAC[a]]


def locate_oligo(gene, oligo, max_mismatches=0):
    candidates = []
    for strand in ["+", "-"]:
        for i in range(len(gene) - len(oligo) + 1):
            site = gene[i : i + len(oligo)]
            if strand == "-":
                site = rc(site)
            if set(site) - set("ACGT"):
                continue
            mm = len(mismatch_positions(oligo, site))
            if mm <= max_mismatches:
                candidates.append((mm, i, i + len(oligo), strand))
    if not candidates:
        raise ValueError(
            "Oligo does not map to the reference gene at the allowed reference mismatch threshold."
        )
    candidates.sort()
    best = [x for x in candidates if x[0] == candidates[0][0]]
    if len(best) != 1:
        raise ValueError(
            f"Oligo maps ambiguously to {len(best)} equally good reference locations. Supply a unique gene/oligo reference."
        )
    mm, start, end, strand = best[0]
    return dict(start=start, end=end, strand=strand, reference_mismatches=mm)


# Tree validation, support preservation and root placement


def separate_pam(reference, oligos, anchors, args):
    """Separate the PAM and guide in the supplied sequence's 5-prime direction.

    Anchor PAM+guide together before splitting to preserve the original locus
    and strand. Guide-only inputs use the adjacent reference bases for the PAM.
    This reports the observed PAM; it does not impose a Cas-specific motif.
    """
    full = oligos["gRNA"]
    original = anchors["gRNA"].copy()
    n = args.pam_length
    if args.guide_input == "pam-guide":
        if len(full) <= n:
            raise ValueError("The PAM+guide input must be longer than --pam-length.")
        pam, guide = (
            (full[:n], full[n:]) if args.pam_side == "5" else (full[-n:], full[:-n])
        )
        pam_at_low = (args.pam_side == "5") == (original["strand"] == "+")
        if pam_at_low:
            pam_start, pam_end = original["start"], original["start"] + n
            guide_start, guide_end = pam_end, original["end"]
        else:
            pam_start, pam_end = original["end"] - n, original["end"]
            guide_start, guide_end = original["start"], pam_start
    else:
        guide = full
        guide_start, guide_end = original["start"], original["end"]
        pam_at_low = (args.pam_side == "5") == (original["strand"] == "+")
        pam_start, pam_end = (
            (guide_start - n, guide_start) if pam_at_low else (guide_end, guide_end + n)
        )
        if pam_start < 0 or pam_end > len(reference):
            raise ValueError(
                "The reference does not contain the complete adjacent PAM region."
            )
        pam = reference[pam_start:pam_end]
        if original["strand"] == "-":
            pam = rc(pam)
    oligos.update(gRNA=guide, PAM=pam)
    for key, start, end, sequence in [
        ("gRNA", guide_start, guide_end, guide),
        ("PAM", pam_start, pam_end, pam),
    ]:
        ref_site = reference[start:end]
        if original["strand"] == "-":
            ref_site = rc(ref_site)
        anchors[key] = dict(
            start=start,
            end=end,
            strand=original["strand"],
            reference_mismatches=len(mismatch_positions(sequence, ref_site)),
        )


def comparison_sequence(key, sequence):
    # Compare four-base PAMs to the assay motif in guide orientation.
    # Other PAM lengths retain their supplied or extracted comparison sequence.
    return "TTTN" if key == "PAM" and len(sequence) == 4 else sequence


def component_threshold(args, key, oligos):
    # PAM substitutions are descriptive: accept all concrete, full-length PAMs.
    if key == "PAM":
        return len(oligos["PAM"])
    value = getattr(args, key.lower() + "_max_mismatches")
    return args.max_mismatches if value is None else value


def validate_newick_lengths(text):
    """Check actual length tokens, skipping quoted labels and bracket comments."""
    i = 0
    while i < len(text):
        if text[i] == "'":
            i += 1
            while i < len(text):
                if text[i] == "'":
                    i += 1
                    if i < len(text) and text[i] == "'":
                        i += 1
                        continue
                    break
                i += 1
            else:
                raise ValueError("Unterminated quoted Newick label.")
        elif text[i] == "[":
            depth = 1
            i += 1
            while i < len(text) and depth:
                depth += (text[i] == "[") - (text[i] == "]")
                i += 1
            if depth:
                raise ValueError("Unterminated Newick comment.")
        elif text[i] == ":":
            i += 1
            while i < len(text) and text[i].isspace():
                i += 1
            start = i
            while i < len(text) and text[i] not in ",();[]" and not text[i].isspace():
                i += 1
            token = text[start:i]
            if not re.fullmatch(
                r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", token
            ):
                raise ValueError(f"Invalid Newick branch length: {token!r}")
        else:
            i += 1


def split_supports(tree):
    """Snapshot branch support before rooting, keyed by an unrooted split."""
    all_tips = frozenset(t.name for t in tree.get_terminals())
    anchor = min(all_tips)
    supports = {}
    for node in tree.get_nonterminals():
        if node is tree.root:
            continue
        value = bootstrap_support(node)
        if value is None:
            continue
        inside = frozenset(t.name for t in node.get_terminals())
        key = inside if anchor not in inside else all_tips - inside
        if key in supports and not math.isclose(value, supports[key], abs_tol=1e-9):
            raise ValueError(
                f"Conflicting support labels for the same taxon split: {supports[key]} and {value}."
            )
        supports[key] = value
    scale = 100 if supports and max(supports.values()) <= 1 else 1
    if any(not math.isfinite(s) or s < 0 or s * scale > 100 for s in supports.values()):
        raise ValueError(
            "Bootstrap support must be percentages (0–100) or proportions (0–1), not replicate counts."
        )
    return {key: value * scale for key, value in supports.items()}


def restore_split_supports(tree, supports):
    """Clear stale node labels and assign support to the same splits after rooting."""
    all_tips = frozenset(t.name for t in tree.get_terminals())
    anchor = min(all_tips)
    seen = set()
    for node in tree.find_clades():
        node.confidence = None
        if node.is_terminal():
            continue
        node.name = None
        node.comment = None
        if node is tree.root:
            continue
        inside = frozenset(t.name for t in node.get_terminals())
        key = inside if anchor not in inside else all_tips - inside
        # The two branches incident to a bifurcating root represent one split.
        if key in supports and key not in seen:
            node.confidence = supports[key]
            seen.add(key)
    if seen != set(supports):
        raise ValueError("Rooting changed a supported taxon split.")


def root_minimum_variance(tree):
    """Find the global minimum of root-to-tip distance variance over every edge.

    Along an edge, each tip distance has slope +1 or -1, so variance is a
    quadratic. Subtree counts and distance moments give each bounded optimum
    in constant time after two traversals (Mai et al., 2017,
    doi:10.1371/journal.pone.0182238). No extra rooting software is needed.
    """
    nodes = list(tree.find_clades(order="preorder"))
    depths = tree.depths()
    unit = max(depths[t] for t in tree.get_terminals())
    depths = {node: depth / unit for node, depth in depths.items()}
    count, sums = {}, {}
    for node in reversed(nodes):
        count[node] = sum(count[c] for c in node.clades) if node.clades else 1
        sums[node] = sum(sums[c] for c in node.clades) if node.clades else depths[node]
    n = count[tree.root]
    total = {tree.root: sums[tree.root]}
    squares = {tree.root: math.fsum(depths[t] ** 2 for t in tree.get_terminals())}
    best = None
    graph = {node: [] for node in nodes}
    for parent in nodes:
        for child in parent.clades:
            length = child.branch_length / unit
            k = count[child]
            subtree_sum = sums[child] - k * depths[parent]
            slope = (n - 2 * k) / n
            covariance = (total[parent] - 2 * subtree_sum) / n - total[
                parent
            ] / n * slope
            curvature = 4 * k * (n - k) / (n * n)
            x = min(length, max(0.0, -covariance / curvature)) if curvature else 0.0
            variance = squares[parent] / n - (total[parent] / n) ** 2
            variance += 2 * covariance * x + curvature * x * x
            if best is None or variance < best[0]:
                best = (variance, parent, child, x * unit)
            total[child] = total[parent] + (n - 2 * k) * length
            squares[child] = (
                squares[parent]
                + 2 * length * (total[parent] - 2 * subtree_sum)
                + n * length * length
            )
            graph[parent].append((child, child.branch_length))
            graph[child].append((parent, child.branch_length))

    _, parent, child, offset = best

    def orient(node, previous, incoming):
        result = copy(node)
        result.branch_length = incoming
        result.clades = [
            orient(other, node, length)
            for other, length in graph[node]
            if other is not previous
        ]
        # Suppress the old artificial root (or other unary internal nodes).
        if len(result.clades) == 1:
            result.clades[0].branch_length += incoming
            return result.clades[0]
        return result

    length = next(length for other, length in graph[parent] if other is child)
    tree.root = Clade(
        branch_length=0.0,
        clades=[orient(parent, child, offset), orient(child, parent, length - offset)],
    )
    tree.rooted = True
    distances = [tree.distance(t) for t in tree.get_terminals()]
    mean = math.fsum(distances) / len(distances)
    return math.fsum((d - mean) ** 2 for d in distances) / len(distances)


def inspect_tree(path, rooting="minvar"):
    validate_newick_lengths(Path(path).read_text())
    tree = Phylo.read(path, "newick")
    tips = [x.name for x in tree.get_terminals()]
    if not tips or any(not x for x in tips) or len(tips) != len(set(tips)):
        raise ValueError("Tree tips must have nonempty, unique identifiers.")
    for node in tree.find_clades():
        if node.branch_length is not None and (
            not math.isfinite(node.branch_length) or node.branch_length < 0
        ):
            raise ValueError("Tree branch lengths must be finite and nonnegative.")
    # Validate all tree information before any BLAST or MAFFT work.
    edges = [node for node in tree.find_clades() if node is not tree.root]
    if (
        len(tips) < 3
        or any(node.branch_length is None for node in edges)
        or not any(node.branch_length > 0 for node in edges)
    ):
        raise ValueError(
            "The figure requires at least three tips and finite nonnegative lengths on every edge (with at least one positive length)."
        )
    supports = split_supports(tree)
    # A root stem is not an edge of the unrooted tree.
    tree.root.branch_length = 0.0
    if rooting == "minvar":
        root_minimum_variance(tree)
    elif rooting == "midpoint":
        tree.root_at_midpoint()
    elif rooting != "input":
        raise ValueError(f"Unknown rooting method: {rooting}")
    restore_split_supports(tree, supports)
    tree.rooting_method = rooting
    return tree


def bootstrap_support(clade):
    """Read UFBoot from an IQ-TREE SH-aLRT/UFBoot pair, or a single support."""
    for raw in (clade.name, clade.confidence, clade.comment):
        if raw is None:
            continue
        if isinstance(raw, (int, float)):
            return float(raw)
        parts = str(raw).strip().split("/")
        if parts and all(
            re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", part)
            for part in parts
        ):
            return float(parts[-1])
    return None


# Conservative reconciliation of tree tips and assembly filenames


def identifier_parts(identifier):
    """Return an accession-aware key, keeping any supplied assembly version."""
    name = re.sub(r"\.(?:fna|fasta|fa)(?:\.gz)?$", "", identifier.strip(), flags=re.I)
    name = re.sub(r"\s+", "_", name)
    match = re.fullmatch(
        r"(?P<species>.+)_(?P<accession>GC[AF]_\d+)(?:\.(?P<version>\d+))?", name
    )
    if match:
        return (match["species"], match["accession"], match["version"])
    return (name, "", None)


def reconcile_identifiers(tipnames, genome_names):
    """Resolve only unique one-to-one matches; never guess between duplicates."""
    unmatched_tips, unmatched_genomes = set(tipnames), set(genome_names)
    matched = {}
    methods = {}
    ambiguous = []

    def stage(method, tip_key, genome_key, compatible=lambda t, g: True):
        for tip in tipnames:
            if tip not in unmatched_tips:
                continue
            key = tip_key(tip)
            candidates = [
                g
                for g in unmatched_genomes
                if genome_key(g) == key and compatible(tip, g)
            ]
            competing = [
                t
                for t in unmatched_tips
                if tip_key(t) == key
                and any(
                    genome_key(g) == key and compatible(t, g) for g in unmatched_genomes
                )
            ]
            if len(candidates) == 1 and len(competing) == 1:
                genome = candidates[0]
                matched[tip] = genome
                methods[tip] = method
                unmatched_tips.remove(tip)
                unmatched_genomes.remove(genome)
            elif candidates:
                ambiguous.append((method, tip, tuple(sorted(candidates))))

    stage("exact", lambda s: s, lambda s: s)
    stage("normalized label / FASTA extension", identifier_parts, identifier_parts)
    # If one side omits an accession version, match the species and accession
    # only when they identify exactly one tree tip and one genome. Two explicit
    # but different version numbers never match.
    stage(
        "unique accession with version missing on one side",
        lambda s: identifier_parts(s)[:2],
        lambda s: identifier_parts(s)[:2],
        lambda t, g: identifier_parts(t)[2] is None or identifier_parts(g)[2] is None,
    )
    unresolved = [(m, t, c) for m, t, c in ambiguous if t in unmatched_tips]
    return (
        matched,
        methods,
        sorted(unmatched_tips),
        sorted(unmatched_genomes),
        unresolved,
    )


# BLAST HSP chaining and grouping into putative physical loci


def trim_hsp_overlap(hit, previous, max_overlap):
    """Trim a small overlapping HSP prefix using its actual aligned sequences.

    Every retained alignment column contributes once to identity and coverage.
    Bit score is prorated for the retained part and used only for ranking.
    Legacy ten-column BLAST input still works, without overlap trimming.
    """
    query_overlap = previous["qe"] - hit["qs"] + 1
    subject_overlap = previous["b"] - hit["a"] + 1
    if max(query_overlap, subject_overlap) <= 0:
        return hit
    if max(query_overlap, subject_overlap) > max_overlap or "qseq" not in hit:
        return None
    if hit["qe"] <= previous["qe"] or hit["b"] <= previous["b"]:
        return None
    query, subject = hit["qseq"], hit["sseq"]
    qs, start, cut = hit["qs"], hit["a"], 0
    while cut < len(query) and (qs <= previous["qe"] or start <= previous["b"]):
        qs += query[cut] != "-"
        start += subject[cut] != "-"
        cut += 1
    if cut == len(query):
        return None
    result = hit.copy()
    result.update(
        qs=qs,
        a=start,
        qseq=query[cut:],
        sseq=subject[cut:],
        alen=len(query) - cut,
        score=hit["score"] * (len(query) - cut) / len(query),
        nident=sum(q == s and q != "-" for q, s in zip(query[cut:], subject[cut:])),
    )
    return result


def select_hits(lines, contigs, reference_length, args):
    """Chain collinear HSPs only within a bounded locus.
    Dynamic programming retains score/coverage Pareto states per start/end HSP.
    Every member independently passes identity; chain coverage never double-counts.
    """
    groups = {}
    for line in lines:
        if not line.strip():
            continue
        fields = line.rstrip().split("\t")
        if len(fields) not in (10, 12):
            raise ValueError("Unexpected BLAST output columns.")
        sid, ident, alen, qs, qe, ss, se, bits, evalue, nident = fields[:10]
        ident, bits, evalue = float(ident), float(bits), float(evalue)
        alen, qs, qe, ss, se, nident = map(int, (alen, qs, qe, ss, se, nident))
        if sid not in contigs:
            raise ValueError(f"Unknown BLAST contig identifier: {sid}")
        if not (1 <= min(qs, qe) <= max(qs, qe) <= reference_length):
            raise ValueError("BLAST query coordinates are out of range.")
        low, high = sorted((ss, se))
        if not (1 <= low <= high <= len(contigs[sid])):
            raise ValueError(f"Invalid BLAST coordinates on {sid}: {ss}, {se}")
        if ident < args.min_identity:
            continue
        strand = "+" if (qe - qs) * (se - ss) >= 0 else "-"
        a, b = (low, high) if strand == "+" else (-high, -low)
        h = dict(
            qs=min(qs, qe),
            qe=max(qs, qe),
            a=a,
            b=b,
            score=bits,
            alen=alen,
            nident=nident,
            evalue=evalue,
        )
        if len(fields) == 12:
            qseq, sseq = fields[10:]
            if qs > qe:
                qseq, sseq = rc(qseq), rc(sseq)
            h.update(qseq=qseq.upper(), sseq=sseq.upper())
            if (
                len(qseq) != alen
                or len(sseq) != alen
                or len(qseq.replace("-", "")) != h["qe"] - h["qs"] + 1
                or len(sseq.replace("-", "")) != b - a + 1
            ):
                raise ValueError(
                    "BLAST aligned sequences disagree with their coordinates."
                )
        groups.setdefault((sid, strand), []).append(h)
    candidates = []
    max_gap = getattr(args, "max_hsp_gap", 1000)
    max_qgap = getattr(args, "max_query_gap", 200)
    for (sid, strand), hsps in groups.items():
        hsps.sort(key=lambda h: (h["qs"], h["qe"], h["a"]))
        states = []
        for j, h in enumerate(hsps):
            options = [
                dict(
                    first=j,
                    path=[j],
                    parts=[h],
                    score=h["score"],
                    covered=h["qe"] - h["qs"] + 1,
                )
            ]
            for i, prev in enumerate(hsps[:j]):
                extension = trim_hsp_overlap(
                    h, prev, getattr(args, "max_hsp_overlap", 10)
                )
                if extension is None:
                    continue
                if not (
                    0 <= extension["qs"] - prev["qe"] - 1 <= max_qgap
                    and 0 <= extension["a"] - prev["b"] - 1 <= max_gap
                ):
                    continue
                for state in states[i]:
                    if (
                        h["b"] - hsps[state["first"]]["a"] + 1
                        > args.max_length_factor * reference_length
                    ):
                        continue
                    options.append(
                        dict(
                            first=state["first"],
                            path=state["path"] + [j],
                            parts=state["parts"] + [extension],
                            score=state["score"] + extension["score"],
                            covered=state["covered"]
                            + extension["qe"]
                            - extension["qs"]
                            + 1,
                        )
                    )
            # Different starts cannot dominate each other because their spans differ.
            frontier = []
            for option in sorted(options, key=lambda v: (-v["score"], -v["covered"])):
                if any(
                    x["first"] == option["first"]
                    and x["score"] >= option["score"]
                    and x["covered"] >= option["covered"]
                    for x in frontier
                ):
                    continue
                frontier.append(option)
            states.append(frontier)
            for state in frontier:
                coverage = 100 * state["covered"] / reference_length
                first = hsps[state["first"]]
                low, high = (
                    (first["a"], h["b"]) if strand == "+" else (-h["b"], -first["a"])
                )
                span = high - low + 1
                if (
                    coverage < args.min_coverage
                    or not args.min_length_factor * reference_length
                    <= span
                    <= args.max_length_factor * reference_length
                ):
                    continue
                chain = state["parts"]
                identity = (
                    100
                    * sum(v["nident"] for v in chain)
                    / sum(v["alen"] for v in chain)
                )
                if identity < args.min_identity:
                    continue
                seq = contigs[sid][low - 1 : high]
                if strand == "-":
                    seq = rc(seq)
                if not seq or set(seq) <= {"N"}:
                    continue
                candidates.append(
                    dict(
                        identity=identity,
                        coverage=coverage,
                        contig=sid,
                        start=low,
                        end=high,
                        strand=strand,
                        sequence=seq,
                        bitscore=state["score"],
                        evalue=min(v["evalue"] for v in chain),
                        query_start=first["qs"],
                        query_end=h["qe"],
                        extracted_length=span,
                        hsp_count=len(chain),
                        _hsp_ids=frozenset(state["path"]),
                    )
                )
    candidates.sort(
        key=lambda h: (
            -h["bitscore"],
            -h["coverage"],
            -h["identity"],
            h["contig"],
            h["start"],
        )
    )
    # Alternative HSP chains can resolve to precisely the same locus. Keep
    # the highest-scoring chain per coordinate span, but retain other loci.
    distinct = []
    seen = set()
    for hit in candidates:
        locus = (hit["contig"], hit["strand"], hit["start"], hit["end"])
        if locus not in seen:
            seen.add(locus)
            distinct.append(hit)
    # Group chain alternatives without discarding any before assay evaluation.
    # Shared HSPs identify one locus; substantially overlapping spans also group
    # if their inferred reference starts agree. A one-base overlap is insufficient.
    # These are putative loci, not validated biological copy counts.
    parents = list(range(len(distinct)))

    def find(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i

    for i, hit in enumerate(distinct):
        for j, other in enumerate(distinct[:i]):
            if (hit["contig"], hit["strand"]) != (other["contig"], other["strand"]):
                continue
            overlap = (
                min(hit["end"], other["end"]) - max(hit["start"], other["start"]) + 1
            )
            origin = (
                lambda h: (h["start"] if h["strand"] == "+" else -h["end"])
                - h["query_start"]
            )
            same_span = overlap >= 0.8 * min(
                hit["extracted_length"], other["extracted_length"]
            ) and abs(origin(hit) - origin(other)) <= max(10, 0.1 * reference_length)
            if hit["_hsp_ids"] & other["_hsp_ids"] or same_span:
                parents[find(i)] = find(j)
    loci = {}
    for i, hit in enumerate(distinct):
        hit["locus_id"] = loci.setdefault(find(i), len(loci) + 1)
    return distinct


def search_genome(index, path, reference_path, reference_length, work, logs, args):
    """Search one assembly with isolated temporary files; return its output row."""
    contigs = read_fasta(path)
    idmap = {f"contig_{i}": name for i, name in enumerate(contigs)}
    safe = {sid: contigs[name] for sid, name in idmap.items()}
    subject = work / f"{index:04d}_genome.fasta"
    hits_path = work / f"{index:04d}_blast.tsv"
    try:
        fasta(subject, safe)
        cmd = [
            args.blastn,
            "-task",
            "blastn",
            "-query",
            reference_path,
            "-subject",
            subject,
            "-dust",
            "no",
            "-soft_masking",
            "false",
            "-word_size",
            "7",
            "-evalue",
            str(args.evalue),
            "-max_target_seqs",
            str(max(100, len(safe))),
            "-outfmt",
            "6 sseqid pident length qstart qend sstart send bitscore evalue nident qseq sseq",
        ]
        run(cmd, hits_path, logs / f"{index:04d}_blast.log")
        with hits_path.open() as handle:
            candidates = select_hits(handle, safe, reference_length, args)
        row = dict(
            genome_id=path.name,
            gene_present=int(bool(candidates)),
            identity="",
            coverage="",
            contig="",
            start="",
            end="",
            strand="",
            bitscore="",
            evalue="",
            query_start="",
            query_end="",
            extracted_length="",
            hsp_count="",
            valid_candidate_count=len({h["locus_id"] for h in candidates}),
            candidate_chain_count=len(candidates),
            selected_candidate_id="",
            assay_compatible_candidate_count=0,
            assay_compatible=0,
            site_geometry="binding_site_unresolved",
            amplicon_length="",
        )
        for rank, hit in enumerate(candidates, 1):
            hit["candidate_id"] = f"HIT_{index:04d}_{rank:05d}"
            hit["genome_id"] = path.name
            hit["contig"] = idmap[hit["contig"]]
        return row, candidates
    finally:
        subject.unlink(missing_ok=True)


# Reference-anchored alignments, site calls and assay ranking


def align_candidate_pairs(reference, hits, work, logs, args):
    """Align each candidate independently, then merge only by reference columns.

    Padding insertion columns permits FASTA export without realigning candidates
    against one another. It preserves each pair's indel and mismatch calls.
    """

    def align(item):
        cid, hit = item
        source, target = work / f"{cid}_pair.fasta", work / f"{cid}_pair_aligned.fasta"
        fasta(source, {REF: reference, cid: hit["sequence"]})
        run(
            [args.mafft, "--auto", "--thread", "1", "--inputorder", source],
            target,
            logs / f"{cid}_mafft.log",
        )
        pair = read_fasta(target, aligned=True)
        if (
            set(pair) != {REF, cid}
            or pair[REF].replace("-", "") != reference
            or pair[cid].replace("-", "") != hit["sequence"]
        ):
            raise ValueError(f"MAFFT pairwise integrity check failed: {cid}")
        bases = []
        insertions = [""] * (len(reference) + 1)
        position = 0
        for ref_base, base in zip(pair[REF], pair[cid]):
            if ref_base == "-":
                insertions[position] += base
            else:
                bases.append(base)
                position += 1
        return cid, bases, insertions

    workers = min(len(hits), args.threads, args.genome_jobs or args.threads)
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        pairs = list(pool.map(align, hits.items()))
    widths = [
        max((len(insertions[i]) for _, _, insertions in pairs), default=0)
        for i in range(len(reference) + 1)
    ]

    def merge(bases, insertions):
        parts = []
        for i, base in enumerate(bases):
            parts.extend((insertions[i].ljust(widths[i], "-"), base))
        parts.append(insertions[-1].ljust(widths[-1], "-"))
        return "".join(parts)

    aligned = {REF: merge(reference, [""] * (len(reference) + 1))}
    aligned.update({cid: merge(bases, insertions) for cid, bases, insertions in pairs})
    return aligned


def site_result(aligned_ref, aligned_gene, oligo, anchor, threshold):
    """Project only reference columns, while rejecting real internal indels."""
    cols = [i for i, c in enumerate(aligned_ref) if c != "-"]
    selected = cols[anchor["start"] : anchor["end"]]
    result = dict(
        present=0,
        status="not_detected",
        mismatches="",
        mismatch_positions="",
        sequence="",
        candidate_sequence="",
        candidate_mismatches="",
        insertions=0,
        deletions=0,
        start="",
        end="",
        strand=anchor["strand"],
        threshold=threshold,
    )
    if len(selected) != len(oligo):
        raise ValueError("Reference anchor length is inconsistent.")
    first, last = selected[0], selected[-1]
    seq = "".join(aligned_gene[i] for i in selected)
    result["insertions"] = sum(
        aligned_ref[i] == "-" and aligned_gene[i] != "-" for i in range(first, last + 1)
    )
    result["deletions"] = seq.count("-")
    non_gap = [i for i, b in enumerate(aligned_gene) if b != "-"]
    if not non_gap or first < non_gap[0] or last > non_gap[-1]:
        result["status"] = "outside_extracted_region"
        return result
    if result["insertions"] or result["deletions"]:
        result["status"] = "indel"
        return result
    result["start"] = sum(c != "-" for c in aligned_gene[:first]) + 1
    result["end"] = sum(c != "-" for c in aligned_gene[: last + 1])
    if anchor["strand"] == "-":
        seq = rc(seq)
    result["candidate_sequence"] = seq
    if set(seq) - set("ACGT"):
        result.update(present="", status="unresolved_bases")
        return result
    mm = mismatch_positions(oligo, seq)
    result["candidate_mismatches"] = len(mm)
    result["mismatch_positions"] = ";".join(map(str, mm))
    if len(mm) > threshold:
        result["status"] = "above_mismatch_threshold"
        return result
    result.update(present=1, status="detected", mismatches=len(mm), sequence=seq)
    return result


def evaluate_candidate(hit, aligned_ref, aligned_gene, oligos, anchors, args):
    """Evaluate all three sites on one physical locus, never across hits."""
    sites = {}
    for key, alias, _, _ in REPORT_COMPONENTS:
        if key not in anchors:
            continue
        threshold = component_threshold(args, key, oligos)
        site = site_result(
            aligned_ref, aligned_gene, comparison_sequence(key, oligos[alias]), anchors[key], threshold
        )
        site["genomic_start"] = ""
        site["genomic_end"] = ""
        site["genomic_strand"] = ""
        if site["start"] != "":
            if hit["strand"] == "+":
                site["genomic_start"] = hit["start"] + site["start"] - 1
                site["genomic_end"] = hit["start"] + site["end"] - 1
            else:
                site["genomic_start"] = hit["end"] - site["end"] + 1
                site["genomic_end"] = hit["end"] - site["start"] + 1
            site["genomic_strand"] = "+" if hit["strand"] == site["strand"] else "-"
        sites[key] = site
    f, r, g = sites["fwd"], sites["rev"], sites["gRNA"]
    amplicon_length = ""
    geometry = "binding_site_unresolved"
    if all(s["start"] != "" for s in (f, r, g)):
        if f["strand"] == r["strand"]:
            geometry = "primers_same_strand"
        elif f["strand"] == "+" and f["end"] < r["start"]:
            amplicon_length = r["end"] - f["start"] + 1
            geometry = (
                "valid"
                if f["end"] < g["start"] and g["end"] < r["start"]
                else "guide_outside_primers"
            )
        elif r["strand"] == "+" and r["end"] < f["start"]:
            amplicon_length = f["end"] - r["start"] + 1
            geometry = (
                "valid"
                if r["end"] < g["start"] and g["end"] < f["start"]
                else "guide_outside_primers"
            )
        else:
            geometry = "primers_not_inward_facing"
        if (
            geometry == "valid"
            and not args.min_amplicon_len <= amplicon_length <= args.max_amplicon_len
        ):
            geometry = "amplicon_out_of_range"
    detected = sum(sites[key]["present"] == 1 for key, _, _, _ in COMPONENTS)
    compatible = geometry == "valid" and detected == 3
    assessed = sum(
        isinstance(sites[key]["candidate_mismatches"], int)
        for key, _, _, _ in COMPONENTS
    )
    mismatch_penalty = sum(
        (
            s["candidate_mismatches"]
            if isinstance(s["candidate_mismatches"], int)
            else len(oligos[alias]) + 1
        )
        for key, alias, _, _ in COMPONENTS
        for s in [sites[key]]
    )
    rank = (
        int(compatible),
        int(geometry == "valid"),
        detected,
        assessed,
        -mismatch_penalty,
        hit["bitscore"],
        hit["coverage"],
        hit["identity"],
    )
    return (
        sites,
        dict(
            assay_compatible=int(compatible),
            site_geometry=geometry,
            amplicon_length=amplicon_length,
            sites_detected=detected,
            total_mismatches=mismatch_penalty if assessed == 3 else "",
        ),
        rank,
    )


# Figure rendering


def plot_figure(tree, genes, sites, oligos, output, args, tip_to_genome=None):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.font_manager import FontProperties

    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42}
    )
    tips = tree.get_terminals()
    n = len(tips)
    if tip_to_genome is None:
        tip_to_genome = {t.name: t.name for t in tips if t.name in genes}
    ys = {t: float(n - 1 - i) for i, t in enumerate(tips)}

    def set_y(node):
        if node not in ys:
            child_rows = [set_y(child) for child in node.clades]
            ys[node] = (child_rows[0] + child_rows[-1]) / 2
        return ys[node]

    set_y(tree.root)
    depths = tree.depths()
    mx = max(depths[t] for t in tips)

    def label(name):
        species, accession, _ = identifier_parts(name)
        words = species.split("_")
        if len(words) >= 2 and accession:
            # Genus-independent display, preserving strain/species punctuation.
            species = words[0][0] + ". " + " ".join(words[1:])
            full_accession = re.search(r"GC[AF]_\d+(?:\.\d+)?", name).group()
            return species + "  " + full_accession
        return species.replace("_", " ")

    labels = {t: label(t.name) for t in tips}
    left = 0.35
    tree_width = 3.0
    label_x = left + tree_width + 0.14
    label_width = max(3.7, max(len(s) for s in labels.values()) * 0.056)
    gene_x = label_x + label_width + 0.25
    step = 0.092
    panel_x = gene_x + 0.9
    panels = []
    for key, alias, title, filename in REPORT_COMPONENTS:
        L = len(oligos[alias])
        panels.append((key, alias, title, panel_x, L))
        panel_x += max(1.5, L * step + 0.50) if key == "PAM" else L * step + 0.72
    width = panel_x + 0.18
    row_h = args.row_height
    bottom = 0.75
    top = 0.8
    height = max(4.5, (n - 1) * row_h + bottom + top)
    fig = plt.figure(figsize=(width, height), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(-bottom / row_h, n - 1 + top / row_h)
    ax.axis("off")
    # x uses inches; y uses shared leaf rows for branches, labels and all panels.
    ink = "#111111"
    muted = "#737373"
    teal = "#087F8C"
    red = "#C23B45"
    # Subtle alternating rows connect labels, gene symbols and assay panels.
    for i, t in enumerate(tips):
        if i % 2 == 0:
            ax.add_patch(
                Rectangle(
                    (label_x - 0.06, ys[t] - 0.5),
                    width - label_x - 0.18,
                    1.0,
                    facecolor="#EEF0F2",
                    edgecolor="none",
                    zorder=-1,
                )
            )
    dx = {c: left + depths[c] / mx * tree_width for c in tree.find_clades()}
    supported_branches = []

    def draw_node(node):
        if node.clades:
            ax.plot(
                [dx[node]] * 2,
                [ys[node.clades[0]], ys[node.clades[-1]]],
                color=ink,
                lw=args.branch_width,
            )
            for child in node.clades:
                ax.plot(
                    [dx[node], dx[child]],
                    [ys[child]] * 2,
                    color=ink,
                    lw=args.branch_width,
                )
                if child.clades and bootstrap_support(child) is not None:
                    supported_branches.append((node, child))
                draw_node(child)

    draw_node(tree.root)
    if args.show_support:
        # Put support above the branch it describes, after all tree lines have
        # been drawn. A small white backing keeps crossings from obscuring it.
        for parent, node in supported_branches:
            support = bootstrap_support(node)
            x = (dx[parent] + dx[node]) / 2
            ax.annotate(
                f"{support:g}",
                (x, ys[node]),
                xytext=(0, 2.5),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=5.5,
                color=muted,
                bbox=dict(facecolor="white", edgecolor="none", pad=0.4),
                zorder=5,
            )
    for t in tips:
        y = ys[t]
        genome_id = tip_to_genome.get(t.name)
        g = genes.get(genome_id)
        ax.plot(
            [dx[t], label_x - 0.06],
            [y, y],
            color="#909090",
            lw=args.guide_width,
            ls=(0, (3, 2)),
            zorder=0,
        )
        label_text = labels[t]
        if "  " in label_text:
            species, accession = label_text.split("  ", 1)
            # Native italic text avoids mathtext interpreting species punctuation.
            ax.text(
                label_x,
                y,
                species,
                fontsize=7.3,
                va="center",
                color=ink,
                fontstyle="italic",
                parse_math=False,
            )
            species_width = (
                fig.canvas.get_renderer().get_text_width_height_descent(
                    species,
                    FontProperties(family="DejaVu Sans", size=7.3, style="italic"),
                    False,
                )[0]
                * 72
                / fig.dpi
            )
            ax.annotate(
                "  " + accession,
                (label_x, y),
                xytext=(species_width, 0),
                textcoords="offset points",
                fontsize=7.3,
                va="center",
                color=ink,
                parse_math=False,
            )
        else:
            ax.text(
                label_x,
                y,
                label_text,
                fontsize=7.3,
                va="center",
                color=ink,
                parse_math=False,
            )
        if g is None:
            ax.scatter(
                [gene_x],
                [y],
                marker="s",
                s=14,
                facecolors="#CDD4D8",
                edgecolors="#84939A",
                linewidths=0.45,
            )
        elif g["gene_present"]:
            ax.scatter([gene_x], [y], s=22, c=teal, linewidths=0)
        else:
            ax.scatter(
                [gene_x], [y], s=22, facecolors="white", edgecolors=teal, linewidths=0.9
            )
        for key, alias, title, x, L in panels:
            s = sites.get(genome_id, {}).get(key)
            if g is None:
                message = "No data available"
            elif not g["gene_present"]:
                message = "Gene not detected"
            elif s["status"] == "detected":
                message = None
            elif s["status"] == "unresolved_bases":
                message = "Ambiguous base(s)"
            elif s["status"] == "indel":
                message = "Indel in target region"
            elif s["status"] == "outside_extracted_region":
                message = "Target not present in gene extract"
            else:
                message = "Greater than {} mismatch{}".format(
                    s["threshold"], "es" if s["threshold"] != 1 else ""
                )
            if key == "PAM" and message:
                message = {
                    "Indel in target region": "Indel",
                    "Ambiguous base(s)": "N",
                    "Target not present in gene extract": "Missing",
                }.get(message, "—")
                if s and s.get("candidate_sequence"):
                    # Show observed ambiguous PAM bases rather than hiding them.
                    s = {**s, "sequence": s["candidate_sequence"]}
                    message = None
            if message:
                ax.text(x, y, message, fontsize=7, va="center", color=muted)
            else:
                pattern = comparison_sequence(key, oligos[alias])
                mismatch = set(mismatch_positions(pattern, s["sequence"]))
                # An unresolved assembly base is not a confirmed substitution.
                mismatch = {i for i in mismatch if s["sequence"][i - 1] in "ACGT"}
                for j, b in enumerate(s["sequence"]):
                    if j + 1 in mismatch:
                        ax.add_patch(
                            Rectangle(
                                (x + j * step - step * 0.48, y - 0.36),
                                step * 0.96,
                                0.72,
                                facecolor="#FBE7E8",
                                lw=0,
                            )
                        )
                    ax.text(
                        x + j * step,
                        y,
                        b,
                        fontfamily="DejaVu Sans Mono",
                        fontsize=7.4,
                        ha="center",
                        va="center",
                        color=muted if b not in "ACGT" else (red if j + 1 in mismatch else ink),
                        fontweight="bold" if j + 1 in mismatch else "normal",
                    )
                if key != "PAM":
                    ax.text(
                        x + L * step + 0.13,
                        y,
                        str(s["mismatches"]),
                        fontsize=7.5,
                        fontweight="bold",
                        ha="center",
                        va="center",
                        color=red if s["mismatches"] else teal,
                    )
    hy = n - 1
    header_y = hy + 0.48 / row_h
    heading = dict(
        fontsize=args.header_size, fontweight="bold", color="black", va="baseline"
    )
    ax.text(left, header_y, "Phylogeny", **heading)
    ax.text(label_x, header_y, "Isolate Species and Accession Number", **heading)
    ax.text(gene_x, header_y, "Gene Presence", ha="center", **heading)
    for key, alias, title, x, L in panels:
        ax.text(x, header_y, title, **heading)
        for j, b in enumerate(comparison_sequence(key, oligos[alias])):
            ax.text(
                x + j * step,
                hy + 0.23 / row_h,
                b,
                fontfamily="DejaVu Sans Mono",
                fontsize=7.4,
                ha="center",
                color=ink,
                va="center",
                fontweight="bold",
            )
        if key != "PAM":
            ax.text(
                x + L * step + 0.13,
                hy + 0.23 / row_h,
                "MM",
                fontsize=7,
                ha="center",
                color=muted,
                va="center",
            )
        for pos in [1] + list(range(5, L + 1, 5)):
            ax.text(
                x + (pos - 1) * step,
                hy + 0.075 / row_h,
                str(pos),
                fontsize=5,
                ha="center",
                color=muted,
            )
        ax.plot(
            [x - 0.055, x + L * step + 0.33],
            [hy + 0.13 / row_h] * 2,
            color="#777777",
            lw=0.6,
        )
    # Root placement changes, but branch-length units and unrooted splits do not.
    scale = 10 ** math.floor(math.log10(mx / 5))
    scale *= max(1, int((mx / 5) / scale))
    bar = tree_width * scale / mx
    y = -0.28 / row_h
    ax.plot([left, left + bar], [y, y], color=ink, lw=1)
    ax.text(
        left,
        y - 0.14 / row_h,
        f"{scale:g} branch-length units",
        fontsize=7,
        ha="left",
        color=muted,
    )
    y = -0.62 / row_h
    legend_items = [
        (left, "o", teal, teal, "Gene detected"),
        (left + 1.85, "o", "white", teal, "Gene not detected"),
    ]
    has_unmatched_tips = any(t.name not in tip_to_genome for t in tips)
    if has_unmatched_tips:
        legend_items.append(
            (left + 4.45, "s", "#CDD4D8", "#84939A", "No data available")
        )
    for x, marker, fc, ec, text_ in legend_items:
        ax.scatter(
            [x + 0.04],
            [y],
            s=22,
            marker=marker,
            facecolors=fc,
            edgecolors=ec,
            linewidths=0.7,
        )
        ax.text(x + 0.17, y, text_, fontsize=8, va="center", color=ink)
    mismatch_x = left + (7.2 if has_unmatched_tips else 4.45)
    ax.text(mismatch_x, y, "A", fontsize=9, fontweight="bold", color=red, va="center")
    ax.text(mismatch_x + 0.19, y, "Mismatch", fontsize=8, va="center", color=ink)
    ax.text(
        mismatch_x + 1.6, y, "MM: mismatch count", fontsize=8, va="center", color=ink
    )
    fig.savefig(output / "assay_phylogeny_figure.pdf")
    fig.savefig(output / "assay_phylogeny_figure.png", dpi=args.dpi)
    plt.close(fig)


# Command-line interface and workflow


def parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    for name in ["genomes", "gene-reference", "assay-oligos", "tree", "output-dir"]:
        p.add_argument("--" + name, required=True, type=Path)
    p.add_argument(
        "--rooting",
        choices=["minvar", "midpoint", "input"],
        default="minvar",
        help="Root placement for the figure. minvar minimizes root-to-tip distance variance; input preserves the supplied root.",
    )
    p.add_argument(
        "--min-identity",
        type=float,
        default=60,
        help="Minimum BLAST nucleotide identity, percent.",
    )
    p.add_argument(
        "--min-coverage",
        type=float,
        default=50,
        help="Minimum reference coverage across nonoverlapping chained HSPs, percent.",
    )
    p.add_argument(
        "--min-length-factor",
        type=float,
        default=0.25,
        help="Minimum extracted length / reference length.",
    )
    p.add_argument(
        "--max-length-factor",
        type=float,
        default=1.5,
        help="Maximum extracted length / reference length.",
    )
    p.add_argument(
        "--max-hsp-gap",
        type=int,
        default=1000,
        help="Maximum unaligned genomic gap between collinear HSPs.",
    )
    p.add_argument(
        "--max-query-gap",
        type=int,
        default=200,
        help="Maximum unaligned reference gap between collinear HSPs.",
    )
    p.add_argument(
        "--max-hsp-overlap",
        type=int,
        default=10,
        help="Maximum query/subject HSP overlap allowed when chaining; overlapping columns are counted once.",
    )
    p.add_argument(
        "--alignment-mode",
        choices=["pairwise", "joint"],
        default="joint",
        help="One joint MAFFT alignment of the reference and all candidates (default, as in v1.4), or independent reference/candidate alignments.",
    )
    p.add_argument(
        "--max-mismatches",
        type=int,
        default=6,
        help="Default maximum substitutions per full-length ungapped site; a sequence-display threshold, not an activity prediction. Override each component as needed.",
    )
    p.add_argument(
        "--min-amplicon-len",
        type=int,
        default=100,
        help="Minimum length of a complete RPA product, including both primers.",
    )
    p.add_argument(
        "--max-amplicon-len",
        type=int,
        default=250,
        help="Maximum length of a complete RPA product, including both primers.",
    )
    for key in ["fwd", "rev", "grna"]:
        p.add_argument("--" + key + "-max-mismatches", type=int)
    p.add_argument(
        "--reference-max-mismatches",
        type=int,
        default=0,
        help="Allowed oligo mismatches when anchoring to reference gene.",
    )
    p.add_argument(
        "--reverse-input",
        choices=["oligo", "target"],
        default="oligo",
        help="oligo: preserve supplied reverse sequence; target: reverse-complement it before comparison.",
    )
    p.add_argument(
        "--guide-input",
        choices=["pam-guide", "guide-only"],
        default="pam-guide",
        help="gRNA/crRNA record contains PAM plus guide (default), or guide only. No repeat/scaffold should be supplied.",
    )
    p.add_argument(
        "--pam-length",
        type=int,
        default=4,
        help="Number of PAM bases to separate or extract adjacent to the guide.",
    )
    p.add_argument(
        "--pam-side",
        choices=["5", "3"],
        default="5",
        help="PAM side relative to the supplied guide/target sequence, read 5-prime to 3-prime. PAM bases are reported, not motif-filtered.",
    )
    p.add_argument(
        "--threads",
        type=int,
        help="Total CPU thread budget. If omitted, prompt in an interactive terminal.",
    )
    p.add_argument(
        "--genome-jobs",
        type=int,
        help="Maximum concurrent genome searches (default: the full --threads budget). Reduce for large assemblies or limited RAM.",
    )
    p.add_argument("--blastn", default="blastn")
    p.add_argument("--mafft", default="mafft")
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument(
        "--branch-width", type=float, default=0.9, help="Tree branch width in points."
    )
    p.add_argument(
        "--guide-width",
        type=float,
        default=0.65,
        help="Dashed tip-to-label guide width in points.",
    )
    p.add_argument(
        "--header-size",
        type=float,
        default=11,
        help="Uniform black bold column-heading size in points.",
    )
    p.add_argument(
        "--row-height", type=float, default=0.20, help="Leaf row spacing in inches."
    )
    p.add_argument(
        "--show-support",
        dest="show_support",
        action="store_true",
        default=True,
        help="Show bootstrap support as numbers out of 100 (default).",
    )
    p.add_argument(
        "--hide-support",
        dest="show_support",
        action="store_false",
        help="Hide branch support labels.",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Replace the output directory only if it belongs to this workflow.",
    )
    return p


def discover_genomes(root, excluded=()):
    """Skip known workflow outputs and explicit reference/oligo inputs."""
    excluded = {Path(p).resolve() for p in excluded}
    files = []
    for directory, subdirs, names in os.walk(root):
        folder = Path(directory)
        if (folder / ".assay_phylogeny_output").is_file():
            subdirs[:] = []
            continue
        for name in names:
            path = folder / name
            if (
                path.is_file()
                and path.resolve() not in excluded
                and re.search(r"\.(fna|fa|fasta)(\.gz)?$", name, re.I)
            ):
                files.append(path)
    return sorted(files)


def main(argv=None):
    args = parser().parse_args(argv)
    if args.pam_length < 1:
        raise ValueError("--pam-length must be a positive integer.")
    args.threads = resolve_threads(args.threads)
    if args.genome_jobs is not None and args.genome_jobs < 1:
        raise ValueError("--genome-jobs must be a positive integer.")
    if not 0 < args.min_identity <= 100 or not 0 < args.min_coverage <= 100:
        raise ValueError("Identity and coverage must be percentages in (0,100].")
    if not 0 < args.min_length_factor <= args.max_length_factor:
        raise ValueError("Invalid length-factor limits.")
    for name in [
        "max_mismatches",
        "reference_max_mismatches",
        "fwd_max_mismatches",
        "rev_max_mismatches",
        "grna_max_mismatches",
    ]:
        value = getattr(args, name)
        if value is not None and value < 0:
            raise ValueError("Mismatch limits must be nonnegative.")
    if min(args.max_hsp_gap, args.max_query_gap, args.max_hsp_overlap) < 0:
        raise ValueError("HSP gap/overlap limits must be nonnegative.")
    if not 1 <= args.min_amplicon_len <= args.max_amplicon_len:
        raise ValueError("Invalid amplicon length limits.")
    if min(args.branch_width, args.guide_width, args.header_size) <= 0:
        raise ValueError("Plot line widths and heading size must be positive.")
    if args.threads < 1 or args.dpi < 72 or args.row_height < 0.13 or args.evalue <= 0:
        raise ValueError("Invalid threads, dpi, row height or evalue.")
    log("[1/6] Validating inputs and mapping oligos to the reference gene")
    if not args.genomes.is_dir():
        raise ValueError(
            "--genomes must be a directory containing unpacked FASTA assemblies."
        )
    out = args.output_dir.resolve()
    inputs = [
        args.genomes.resolve(),
        args.tree.resolve(),
        args.gene_reference.resolve(),
        args.assay_oligos.resolve(),
    ]
    if any(x == out or out in x.parents for x in inputs):
        raise ValueError("Output directory must not contain any input path.")
    if args.genomes.resolve() in out.parents:
        raise ValueError(
            "Output directory must be outside --genomes, so generated FASTAs cannot be searched as assemblies on reruns."
        )
    files = discover_genomes(args.genomes, [args.gene_reference, args.assay_oligos])
    if not files:
        raise ValueError("No .fna, .fa, .fasta or corresponding .gz genomes found.")
    names = [p.name for p in files]
    if len(names) != len(set(names)):
        raise ValueError(
            "Duplicate genome basenames; identifiers must be unique across subdirectories."
        )
    if any(re.search(r"\s", n) for n in names):
        raise ValueError("Genome filenames must not contain whitespace.")
    for executable in [args.blastn, args.mafft]:
        if not shutil.which(executable):
            raise ValueError(
                f"Missing executable: {executable}. Install dependencies as described in README.md."
            )
    references = read_fasta(args.gene_reference)
    if len(references) != 1:
        raise ValueError("Reference gene FASTA must contain exactly one record.")
    reference = next(iter(references.values()))
    if set(reference) - set("ACGT"):
        raise ValueError(
            "Reference gene must contain only A/C/G/T for unambiguous site anchoring."
        )
    raw = read_fasta(args.assay_oligos)
    aliases = {
        "fwdprimer": "fwdPrimer",
        "forwardrpaprimer": "fwdPrimer",
        "revprimer": "revPrimer",
        "reverserpaprimer": "revPrimer",
        "grna": "gRNA",
        "crrna": "gRNA",
    }
    oligos = {}
    for name, seq in raw.items():
        canon = aliases.get(name.lower())
        if canon is None or canon in oligos:
            raise ValueError(
                "Oligo IDs must uniquely identify fwdPrimer, revPrimer and gRNA/crRNA (ForwardRPAprimer/ReverseRPAprimer also accepted)."
            )
        oligos[canon] = seq
    if set(oligos) != {"fwdPrimer", "revPrimer", "gRNA"}:
        raise ValueError("Supply all three oligos.")
    original = oligos.copy()
    if args.reverse_input == "target":
        oligos["revPrimer"] = rc(oligos["revPrimer"])
    anchors = {}
    for key, alias, _, _ in COMPONENTS:
        try:
            anchors[key] = locate_oligo(
                reference, oligos[alias], args.reference_max_mismatches
            )
        except ValueError as e:
            raise ValueError(f"{alias}: {e}") from e
    separate_pam(reference, oligos, anchors, args)
    for key, alias, title, _ in REPORT_COMPONENTS:
        a = anchors[key]
        log(
            f'  {title}: {len(oligos[alias])} nt; reference {a["start"]+1}-{a["end"]} ({a["strand"]})'
        )
    warnings = []

    def warn(msg):
        warnings.append(msg)
        log("WARNING: " + msg)

    _, reference_assessment, _ = evaluate_candidate(
        dict(
            start=1,
            end=len(reference),
            strand="+",
            bitscore=0,
            coverage=100,
            identity=100,
        ),
        reference,
        reference,
        oligos,
        anchors,
        args,
    )
    if not reference_assessment["assay_compatible"]:
        warn(
            "The reference assay is not complete under the configured primer orientation, guide placement, "
            f'amplicon length and mismatch criteria ({reference_assessment["site_geometry"]}); '
            "check oligo orientation and --min-amplicon-len/--max-amplicon-len."
        )
    if anchors["fwd"]["strand"] == anchors["rev"]["strand"]:
        warn(
            "Both primer inputs map to the same reference strand. Sequences are preserved as supplied. If the reverse entry is the forward-strand target, rerun with --reverse-input target. This run is not evidence of a valid primer pair."
        )
    tree = inspect_tree(args.tree, args.rooting)
    if out.exists() and any(out.iterdir()):
        if not args.force or not (out / ".assay_phylogeny_output").is_file():
            raise ValueError(
                "Output directory is nonempty; choose a new directory or use --force for a previous workflow output."
            )
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / ".assay_phylogeny_output").write_text(VERSION)
    work = out / "work"
    work.mkdir()
    logs = out / "logs"
    logs.mkdir()
    # Export the plotted rooting, with support on its original taxon splits.
    Phylo.write(
        tree,
        out / "figure_rooted_tree.nwk",
        "newick",
        format_branch_length="%1.12g",
        format_confidence="%1.8g",
    )
    fasta(work / "reference.fasta", {REF: reference})
    fasta(out / "effective_assay_oligos.fasta", oligos)
    fasta(out / "input_assay_oligos.fasta", original)
    table(
        out / "oligo_reference_mapping.csv",
        [
            dict(
                component=key,
                input_sequence=original.get(alias, ""),
                effective_sequence=oligos[alias],
                reference_start=anchors[key]["start"] + 1,
                reference_end=anchors[key]["end"],
                reference_strand=anchors[key]["strand"],
                reference_mismatches=anchors[key]["reference_mismatches"],
            )
            for key, alias, _, _ in REPORT_COMPONENTS
        ],
        [
            "component",
            "input_sequence",
            "effective_sequence",
            "reference_start",
            "reference_end",
            "reference_strand",
            "reference_mismatches",
        ],
    )
    genes = {}
    candidate_hits = {}
    candidates_by_genome = {}
    extracted = {}
    log(
        f"[2/6] Searching {len(files)} genome(s) with BLASTn; retaining all qualifying hits"
    )
    # Each BLAST -subject search uses a single core, so by default every genome
    # runs concurrently up to your full --threads budget. Pass --genome-jobs to
    # cap this lower if peak memory becomes a problem on very large assemblies.
    workers = min(len(files), args.threads, args.genome_jobs or args.threads)
    log(
        f"  {workers} concurrent genome search(es); total MAFFT budget: {args.threads} thread(s); alignment mode: {args.alignment_mode}"
    )
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(
                search_genome,
                index,
                path,
                work / "reference.fasta",
                len(reference),
                work,
                logs,
                args,
            )
            for index, path in enumerate(files, 1)
        ]
        for index, (path, future) in enumerate(zip(files, futures), 1):
            try:
                row, candidates = future.result()
            except Exception as exc:
                raise RuntimeError(
                    f"Genome search failed for {path.name}: {exc}"
                ) from exc
            genes[path.name] = row
            candidates_by_genome[path.name] = candidates
            candidate_hits.update({hit["candidate_id"]: hit for hit in candidates})
            log(
                f'  [{index}/{len(files)}] {path.name}: {row["valid_candidate_count"]} putative locus/loci from {len(candidates)} candidate chain(s)'
            )
    log("[3/6] Aligning every candidate hit with the gene reference")
    alignment_input = {
        REF: reference,
        **{cid: h["sequence"] for cid, h in candidate_hits.items()},
    }
    fasta(work / "alignment_input.fasta", alignment_input)
    fasta(
        out / "all_candidate_chains.fasta",
        {cid: h["sequence"] for cid, h in candidate_hits.items()},
    )
    all_aligned_path = out / "all_candidate_chains_aligned.fasta"
    if candidate_hits and args.alignment_mode == "pairwise":
        aligned = align_candidate_pairs(reference, candidate_hits, work, logs, args)
        fasta(all_aligned_path, aligned)
    elif candidate_hits:
        run(
            [
                args.mafft,
                "--auto",
                "--thread",
                str(args.threads),
                "--inputorder",
                work / "alignment_input.fasta",
            ],
            all_aligned_path,
            logs / "mafft.log",
        )
        aligned = read_fasta(all_aligned_path, aligned=True)
    else:
        aligned = alignment_input
        fasta(all_aligned_path, aligned)
        log(
            "  No qualifying gene hits; continuing with empty gene and assay-site outputs."
        )
    if set(aligned) != set(alignment_input):
        raise ValueError("MAFFT changed or lost candidate identifiers.")
    for name, seq in aligned.items():
        if seq.replace("-", "") != alignment_input[name]:
            raise ValueError(f"MAFFT sequence integrity check failed: {name}")
    log(
        "[4/6] Comparing all candidate loci and selecting the best assay match per genome"
    )
    evaluations = {}
    candidate_rows = []
    selected_ids = {}
    for cid, hit in candidate_hits.items():
        sites, assessment, rank = evaluate_candidate(
            hit, aligned[REF], aligned[cid], oligos, anchors, args
        )
        evaluations[cid] = dict(sites=sites, assessment=assessment, rank=rank)
    locus_representatives = set()
    for name, candidates in candidates_by_genome.items():
        if not candidates:
            continue
        groups = {}
        for hit in candidates:
            groups.setdefault(hit["locus_id"], []).append(hit)
        representatives = [
            max(group, key=lambda h: evaluations[h["candidate_id"]]["rank"])
            for group in groups.values()
        ]
        locus_representatives.update(h["candidate_id"] for h in representatives)
        best = max(
            representatives, key=lambda h: evaluations[h["candidate_id"]]["rank"]
        )
        selected_ids[name] = best["candidate_id"]
        extracted[name] = best["sequence"]
        genes[name].update(
            {
                k: v
                for k, v in best.items()
                if k not in ("sequence", "candidate_id", "genome_id")
                and not k.startswith("_")
            }
        )
        genes[name]["selected_candidate_id"] = best["candidate_id"]
        genes[name]["assay_compatible_candidate_count"] = sum(
            evaluations[h["candidate_id"]]["assessment"]["assay_compatible"]
            for h in representatives
        )
        chosen = evaluations[best["candidate_id"]]["assessment"]
        genes[name].update(chosen)
        log(
            f'  {name}: selected {best["candidate_id"]} on {best["contig"]}:{best["start"]}-{best["end"]}; '
            f'{chosen["sites_detected"]}/3 sites within mismatch limits; geometry={chosen["site_geometry"]}'
        )
    for cid, hit in candidate_hits.items():
        evaluation = evaluations[cid]
        row = {
            k: v for k, v in hit.items() if k != "sequence" and not k.startswith("_")
        }
        row.update(evaluation["assessment"])
        row["selected"] = int(selected_ids.get(hit["genome_id"]) == cid)
        row["locus_representative"] = int(cid in locus_representatives)
        for key, _, _, _ in REPORT_COMPONENTS:
            site = evaluation["sites"][key]
            for field in (
                "status",
                "present",
                "candidate_mismatches",
                "mismatch_positions",
                "candidate_sequence",
                "genomic_start",
                "genomic_end",
                "genomic_strand",
            ):
                row[f"{key}_{field}"] = site[field]
        candidate_rows.append(row)
    candidate_fields = [
        "genome_id",
        "locus_id",
        "candidate_id",
        "selected",
        "locus_representative",
        "contig",
        "start",
        "end",
        "strand",
        "identity",
        "coverage",
        "bitscore",
        "evalue",
        "extracted_length",
        "hsp_count",
        "assay_compatible",
        "site_geometry",
        "amplicon_length",
        "sites_detected",
        "total_mismatches",
    ] + [
        f"{key}_{field}"
        for key, _, _, _ in REPORT_COMPONENTS
        for field in (
            "status",
            "present",
            "candidate_mismatches",
            "mismatch_positions",
            "candidate_sequence",
            "genomic_start",
            "genomic_end",
            "genomic_strand",
        )
    ]
    table(out / "all_candidate_chains.csv", candidate_rows, candidate_fields)
    table(
        out / "all_candidate_loci.csv",
        [r for r in candidate_rows if r["locus_representative"]],
        candidate_fields,
    )
    fasta(
        out / "all_candidate_loci.fasta",
        {
            cid: candidate_hits[cid]["sequence"]
            for cid in candidate_hits
            if cid in locus_representatives
        },
    )
    fasta(
        out / "all_candidate_loci_aligned.fasta",
        {cid: aligned[cid] for cid in candidate_hits if cid in locus_representatives},
    )
    fasta(out / "extracted_genes.fasta", extracted)
    gene_fields = [
        "genome_id",
        "gene_present",
        "identity",
        "coverage",
        "contig",
        "start",
        "end",
        "strand",
        "bitscore",
        "evalue",
        "query_start",
        "query_end",
        "extracted_length",
        "hsp_count",
        "valid_candidate_count",
        "candidate_chain_count",
        "selected_candidate_id",
        "assay_compatible_candidate_count",
        "assay_compatible",
        "site_geometry",
        "amplicon_length",
    ]
    table(out / "gene_presence.csv", genes.values(), gene_fields)
    selected_alignment = {
        REF: aligned[REF],
        **{name: aligned[cid] for name, cid in selected_ids.items()},
    }
    fasta(out / "reference_anchored_alignment.fasta", selected_alignment)
    fasta(
        out / "extracted_genes_aligned.fasta",
        {name: aligned[cid] for name, cid in selected_ids.items()},
    )
    sites = {}
    summary = []
    detail = []
    for name, gene in genes.items():
        sites[name] = {}
        row = dict(
            genome_id=name,
            gene_present=gene["gene_present"],
            selected_candidate_id=gene["selected_candidate_id"],
            assay_compatible=gene["assay_compatible"],
            site_geometry=gene["site_geometry"],
            amplicon_length=gene["amplicon_length"],
        )
        for key, alias, _, _ in REPORT_COMPONENTS:
            threshold = component_threshold(args, key, oligos)
            if name in selected_ids:
                site = evaluations[selected_ids[name]]["sites"][key]
            else:
                site = dict(
                    present=0,
                    status="gene_not_detected",
                    mismatches="",
                    mismatch_positions="",
                    sequence="",
                    candidate_sequence="",
                    candidate_mismatches="",
                    insertions=0,
                    deletions=0,
                    start="",
                    end="",
                    strand=anchors[key]["strand"],
                    genomic_start="",
                    genomic_end="",
                    genomic_strand="",
                    threshold=threshold,
                )
            sites[name][key] = site
            for field in ("present", "mismatches", "status", "mismatch_positions"):
                row[f"{key}_{field}"] = site[field]
            detail.append(
                dict(genome_id=name, component=key, max_mismatches=threshold, **site)
            )
        summary.append(row)
    fields = [
        "genome_id",
        "gene_present",
        "selected_candidate_id",
        "assay_compatible",
        "site_geometry",
        "amplicon_length",
    ] + [
        f"{key}_{field}"
        for key, _, _, _ in REPORT_COMPONENTS
        for field in ("present", "mismatches", "status", "mismatch_positions")
    ]
    table(out / "assay_site_summary.csv", summary, fields)
    table(
        out / "assay_site_details.csv",
        detail,
        [
            "genome_id",
            "component",
            "max_mismatches",
            "present",
            "status",
            "mismatches",
            "mismatch_positions",
            "sequence",
            "candidate_sequence",
            "candidate_mismatches",
            "insertions",
            "deletions",
            "start",
            "end",
            "strand",
            "genomic_start",
            "genomic_end",
            "genomic_strand",
        ],
    )
    for key, alias, _, filename in REPORT_COMPONENTS:
        rec = {f"REFERENCE_{alias}": oligos[alias]}
        rec.update(
            {
                name: sites[name][key]["sequence"]
                for name in genes
                if sites[name][key]["present"] == 1
            }
        )
        if len({len(s) for s in rec.values()}) != 1:
            raise ValueError(f"Inconsistent short-site lengths for {key}")
        fasta(out / filename, rec)
        if key == "gRNA":
            fasta(out / "crRNA_aligned.fasta", rec)
    log("[5/6] Reconciling complete filenames with Newick tip labels")
    tipnames = [t.name for t in tree.get_terminals()]
    tip_to_genome, match_methods, missing, extra, ambiguous = reconcile_identifiers(
        tipnames, list(genes)
    )
    sections = [
        ("Tree tips with no matching gene-extraction record", missing),
        ("Gene-extraction records with no matching tree tip", extra),
    ]
    with (out / "identifier_mismatch_report.txt").open("w") as handle:
        handle.write(
            "Matching: exact filename, then recognized FASTA extension and whitespace differences, then unique species/accession when only one side lacks a version. Two explicit different versions never match.\n\n"
        )
        handle.write("Matched tree tip -> original genome filename (method)\n")
        for tip in tipnames:
            if tip in tip_to_genome:
                handle.write(f"{tip} -> {tip_to_genome[tip]} ({match_methods[tip]})\n")
        handle.write("\n")
        for title, items in sections:
            handle.write(
                title + f" ({len(items)})\n" + "\n".join(items or ["(none)"]) + "\n\n"
            )
            if items:
                preview = ", ".join(items[:5])
                suffix = (
                    f", ... ({len(items)-5} more in report)" if len(items) > 5 else ""
                )
                warn(f"{title} ({len(items)}): {preview}{suffix}")
        handle.write(f"Ambiguous unmatched tree tips ({len(ambiguous)})\n")
        for method, tip, candidates in ambiguous:
            handle.write(f'{tip} ({method}) -> {", ".join(candidates)}\n')
    if ambiguous:
        warn(
            f"{len(ambiguous)} ambiguous tip match(es) were left unmatched; see identifier_mismatch_report.txt"
        )
    normalized = sum(method != "exact" for method in match_methods.values())
    if normalized:
        log(f"  Matched {normalized} tip(s) after conservative label normalization.")
    table(
        out / "plotted_tip_order.csv",
        [
            dict(
                row=i + 1,
                tree_tip=t,
                genome_id=tip_to_genome.get(t, ""),
                match_method=match_methods.get(t, ""),
                has_genome_result=int(t in tip_to_genome),
            )
            for i, t in enumerate(tipnames)
        ],
        ["row", "tree_tip", "genome_id", "match_method", "has_genome_result"],
    )
    log(
        "[6/6] Drawing the shared-row phylogeny and assay panels; exporting PDF and PNG"
    )
    plot_figure(tree, genes, sites, oligos, out, args, tip_to_genome)
    root_distances = [tree.distance(t) for t in tree.get_terminals()]
    root_mean = math.fsum(root_distances) / len(root_distances)
    stats = dict(
        genomes_searched=len(genes),
        genes_detected=len(extracted),
        tree_rooting=args.rooting,
        root_to_tip_variance=math.fsum((d - root_mean) ** 2 for d in root_distances)
        / len(root_distances),
        candidate_hits_evaluated=len(locus_representatives),
        candidate_chains_evaluated=len(candidate_hits),
        assay_compatible_candidates=sum(
            evaluations[cid]["assessment"]["assay_compatible"]
            for cid in locus_representatives
        ),
        genomes_with_complete_assay=sum(
            g["assay_compatible_candidate_count"] > 0 for g in genes.values()
        ),
        tree_tips_plotted=len(tipnames),
        tree_tips_with_results=len(tip_to_genome),
        tree_tips_without_results=len(missing),
        records_without_tree_tip=len(extra),
        ambiguous_tip_matches=len(ambiguous),
        sites_detected={
            key: sum(s[key]["present"] == 1 for s in sites.values())
            for key, _, _, _ in COMPONENTS
        },
    )
    versions = {}
    for key, cmd in [
        ("blastn", [args.blastn, "-version"]),
        ("mafft", [args.mafft, "--version"]),
    ]:
        r = subprocess.run(cmd, capture_output=True, text=True)
        versions[key] = (r.stdout + r.stderr).strip()
    import Bio, matplotlib

    versions.update(
        python=sys.version, biopython=Bio.__version__, matplotlib=matplotlib.__version__
    )
    hashes = {
        str(p): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [args.gene_reference, args.assay_oligos, args.tree]
    }
    manifest = dict(
        workflow_version=VERSION,
        utc=datetime.now(timezone.utc).isoformat(),
        command=sys.argv,
        parameters={
            k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
        },
        versions=versions,
        input_sha256=hashes,
        genomes=[
            dict(path=str(p.resolve()), size_bytes=p.stat().st_size) for p in files
        ],
        summary=stats,
        warnings=warnings,
    )
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "run_summary.json").write_text(json.dumps(stats, indent=2) + "\n")
    (out / "warnings.txt").write_text("\n".join(warnings) + "\n")
    log("COMPLETE\n" + json.dumps(stats, indent=2) + "\nOutputs: " + str(out))
    return stats


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
