#!/usr/bin/env python3
"""
combinatorial_assembly.py — Real-Data BGC Chimera Pipeline (v3, MIBiG 4.0)

Reads real protein sequences from data/assembly_proteins.json (parsed from
local MIBiG 4.0 FASTA + JSON by fetch_mibig_data.py) and evaluates
combinatorial chimeric assembly line architectures.

Biological rationale:
  In modular PKS/NRPS assembly lines, individual mega-proteins are arranged
  in a linear order. Inter-protein communication is mediated by short
  N-terminal and C-terminal "docking domains" (~30 residues). We score the
  hydropathy compatibility at these real inter-protein boundaries using a
  Gaussian-decay kernel (DJCS), calibrated against natural same-cluster
  boundaries.

Data: 972 PKS/NRP BGCs, 6,523 assembly-line proteins from MIBiG 4.0.
Zero external dependencies. Deterministic (seed=42). Fully offline.
"""

import json
import random
import math
import csv
import os
import sys

# ── 0. DETERMINISTIC EXECUTION ──────────────────────────────────────
RANDOM_SEED = 42
random.seed(RANDOM_SEED)

# ── 1. HYDROPATHY DATA (Kyte & Doolittle 1982, J. Mol. Biol. 157:105) ─
HYDROPATHY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# E. coli K-12 most-frequent codons (for output annotation)
ECOLI_CODONS = {
    'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAC', 'C': 'TGC',
    'Q': 'CAG', 'E': 'GAA', 'G': 'GGC', 'H': 'CAC', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTC', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG'
}

# Docking domain window size (residues)
DOCK_WINDOW = 30

# How many BGCs to sample for chimera generation (full dataset = 972 BGCs
# which would produce millions of chimeras; we sample for tractability)
MAX_BGC_SAMPLE = 50

# Maximum chimeras to generate
MAX_CHIMERAS = 10000

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "assembly_proteins.json")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


# ── 2. DATA LOADING ────────────────────────────────────────────────
def load_bgc_data():
    """Load real BGC data from the MIBiG-parsed cache file."""
    if not os.path.exists(DATA_PATH):
        print(f"[ERROR] Data file not found: {DATA_PATH}")
        print("       Run 'python fetch_mibig_data.py' first to parse MIBiG 4.0 data.")
        sys.exit(1)
    with open(DATA_PATH) as f:
        return json.load(f)


# ── 3. SCORING ──────────────────────────────────────────────────────
def mean_hydropathy(seq):
    """Compute mean Kyte-Doolittle hydropathy for a peptide sequence."""
    vals = [HYDROPATHY.get(aa, 0.0) for aa in seq if aa in HYDROPATHY]
    return sum(vals) / max(1, len(vals))


def boundary_hydropathy_diff(seq_up, seq_down, window=DOCK_WINDOW):
    """
    Compute the hydropathy differential at the inter-protein boundary.
    Uses the C-terminal 'window' residues of the upstream protein and
    the N-terminal 'window' residues of the downstream protein.
    """
    c_term = seq_up[-window:] if len(seq_up) >= window else seq_up
    n_term = seq_down[:window] if len(seq_down) >= window else seq_down
    return abs(mean_hydropathy(c_term) - mean_hydropathy(n_term))


def junction_compatibility_score(diff, gamma):
    """Gaussian-decay DJCS: penalizes hydropathy mismatch at boundaries."""
    return round(100.0 * math.exp(-(diff ** 2) * gamma), 4)


def calibrate_gamma(bgc_entries):
    """
    Auto-calibrate gamma from natural intra-cluster consecutive protein pairs.
    These are known-functional boundaries; we set gamma so that their mean
    DJCS equals ~90 (more stringent than before to improve discrimination).
    """
    natural_diffs = []
    for entry in bgc_entries:
        proteins = entry.get("core_proteins", [])
        if len(proteins) < 2:
            continue
        # Proteins are ordered by genomic coordinate in the FASTA
        for i in range(len(proteins) - 1):
            d = boundary_hydropathy_diff(
                proteins[i]["sequence"],
                proteins[i + 1]["sequence"]
            )
            natural_diffs.append(d)

    if not natural_diffs:
        print("[WARN] No natural boundaries found for calibration, using default gamma=0.5")
        return 0.5, 0.0, []

    mean_diff = sum(natural_diffs) / len(natural_diffs)
    # Force natural boundaries to score ~90: 90 = 100 * exp(-gamma * diff^2)
    # gamma = -ln(0.90) / diff^2
    diff_sq = max(0.001, mean_diff ** 2)
    gamma = -math.log(0.90) / diff_sq

    return round(gamma, 6), round(mean_diff, 4), natural_diffs


# ── 4. LINKER RESCUE ────────────────────────────────────────────────
FLEXIBLE_LINKERS = [
    ("GGSGG", "Glycine-Serine flex"),
    ("GGSGGSGG", "Extended GS flex"),
    ("EAAAK", "Alpha-helix rigid"),
]

def try_rescue(seq_up, seq_down, gamma):
    """
    Attempt to improve a junction by computationally inserting a flexible
    linker peptide. NOTE: This is a computational heuristic; linker insertion
    between megasynthase subunits requires experimental validation.
    """
    original_diff = boundary_hydropathy_diff(seq_up, seq_down)
    original_djcs = junction_compatibility_score(original_diff, gamma)
    best = (None, original_djcs)
    for linker_seq, linker_name in FLEXIBLE_LINKERS:
        new_up = seq_up + linker_seq
        new_diff = boundary_hydropathy_diff(new_up, seq_down)
        new_djcs = junction_compatibility_score(new_diff, gamma)
        if new_djcs > best[1]:
            best = (linker_name, new_djcs)
    return best


# ── 5. CHIMERA GENERATION ──────────────────────────────────────────
def generate_chimeras(sampled_entries, gamma):
    """
    Generate chimeric assembly lines by swapping mega-proteins between BGCs.
    For each pair of BGCs, try replacing one protein in BGC-A's assembly line
    with a protein from BGC-B. Limit total output for tractability.
    """
    chimeras = []
    chimera_count = 0

    for i, entry_a in enumerate(sampled_entries):
        prots_a = entry_a.get("core_proteins", [])
        if len(prots_a) < 2:
            continue
        for j, entry_b in enumerate(sampled_entries):
            if i == j:
                continue
            prots_b = entry_b.get("core_proteins", [])
            if not prots_b:
                continue

            # For each position in BGC-A, try one random donor from BGC-B
            for pos in range(min(len(prots_a), 5)):  # Cap positions to avoid explosion
                donor = random.choice(prots_b)

                # Build chimeric assembly line
                chimeric_line = list(prots_a)
                chimeric_line[pos] = donor

                # Score all boundaries
                boundary_scores = []
                boundary_diffs = []
                for k in range(len(chimeric_line) - 1):
                    d = boundary_hydropathy_diff(
                        chimeric_line[k]["sequence"],
                        chimeric_line[k + 1]["sequence"]
                    )
                    s = junction_compatibility_score(d, gamma)
                    boundary_scores.append(s)
                    boundary_diffs.append(d)

                mean_djcs = round(sum(boundary_scores) / max(1, len(boundary_scores)), 4)
                mean_diff = round(sum(boundary_diffs) / max(1, len(boundary_diffs)), 4)

                # Try rescue on the worst boundary
                rescued = False
                rescued_linker = ""
                if boundary_scores:
                    worst_idx = boundary_scores.index(min(boundary_scores))
                    worst_djcs = boundary_scores[worst_idx]
                    if worst_djcs < 70.0:
                        linker_name, new_djcs = try_rescue(
                            chimeric_line[worst_idx]["sequence"],
                            chimeric_line[worst_idx + 1]["sequence"],
                            gamma
                        )
                        if linker_name and new_djcs > worst_djcs + 1.0:
                            rescued = True
                            rescued_linker = linker_name
                            boundary_scores[worst_idx] = new_djcs
                            mean_djcs = round(sum(boundary_scores) / len(boundary_scores), 4)

                chimeras.append({
                    "chimera_id": f"CHM_{chimera_count:05d}",
                    "host_bgc": entry_a["bgc_id"],
                    "host_organism": entry_a["organism"],
                    "host_compound": entry_a["compounds"][0] if entry_a["compounds"] else "",
                    "donor_bgc": entry_b["bgc_id"],
                    "donor_organism": entry_b["organism"],
                    "swap_position": pos,
                    "donor_protein": donor["protein_id"],
                    "donor_product": donor["product"],
                    "donor_length": donor["length"],
                    "assembly_size": len(chimeric_line),
                    "num_boundaries": len(boundary_scores),
                    "mean_djcs": mean_djcs,
                    "mean_delta_h": mean_diff,
                    "min_djcs": round(min(boundary_scores), 4) if boundary_scores else 0,
                    "max_djcs": round(max(boundary_scores), 4) if boundary_scores else 0,
                    "rescued": rescued,
                    "rescue_linker": rescued_linker,
                    "protein_ids": [p["protein_id"] for p in chimeric_line],
                })
                chimera_count += 1

                if chimera_count >= MAX_CHIMERAS:
                    chimeras.sort(key=lambda x: x["mean_djcs"], reverse=True)
                    return chimeras

    chimeras.sort(key=lambda x: x["mean_djcs"], reverse=True)
    return chimeras


# ── 6. PERMUTATION TEST ────────────────────────────────────────────
def permutation_test(chimeras, all_proteins, gamma, n_perm=1000):
    """
    Monte Carlo permutation test: compare top chimera DJCS against random
    protein boundary pairings drawn from the entire protein pool.
    """
    if not chimeras:
        return {"error": "No chimeras generated"}

    top_djcs = chimeras[0]["mean_djcs"]
    null_scores = []
    better_count = 0

    for _ in range(n_perm):
        p1, p2 = random.sample(all_proteins, 2)
        d = boundary_hydropathy_diff(p1["sequence"], p2["sequence"])
        s = junction_compatibility_score(d, gamma)
        null_scores.append(s)
        if s >= top_djcs:
            better_count += 1

    null_mean = round(sum(null_scores) / len(null_scores), 4)
    null_std = round(math.sqrt(sum((x - null_mean) ** 2 for x in null_scores) / len(null_scores)), 4)
    p_value = better_count / n_perm

    return {
        "observed_top_djcs": top_djcs,
        "null_mean_djcs": null_mean,
        "null_std_djcs": null_std,
        "n_permutations": n_perm,
        "p_value": p_value if p_value > 0 else f"< {1.0 / n_perm}",
        "significance": "Significant (p < 0.05)" if p_value < 0.05 else "Not Significant",
        "empirical_gamma": gamma,
        "n_chimeras_evaluated": len(chimeras),
        "n_proteins_in_pool": len(all_proteins),
    }


# ── 7. MAIN PIPELINE ───────────────────────────────────────────────
def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("=" * 72)
    print("ClawRxiv Generative BGC Forge v3 — MIBiG 4.0 Local Pipeline")
    print("=" * 72)

    # Load real data
    all_entries = load_bgc_data()
    print(f"\n[DATA] Loaded {len(all_entries)} PKS/NRP BGC entries from local MIBiG 4.0")

    # Collect all proteins for the permutation test pool
    all_proteins = []
    for entry in all_entries:
        all_proteins.extend(entry.get("core_proteins", []))
    print(f"[POOL] Total assembly-line proteins: {len(all_proteins)}")
    lengths = [p["length"] for p in all_proteins]
    print(f"       Size range: {min(lengths)}–{max(lengths)} aa")
    print(f"       Mean size: {sum(lengths) / len(lengths):.0f} aa")

    # Calibrate gamma from ALL natural boundaries
    print(f"\n[CALIBRATION] Computing gamma from natural boundaries...")
    gamma, natural_baseline, natural_diffs = calibrate_gamma(all_entries)
    print(f"  Natural boundary ΔH mean = {natural_baseline}")
    print(f"  Natural boundaries measured: {len(natural_diffs)}")
    print(f"  Empirical γ = {gamma}")

    # Sample BGCs for chimera generation (full dataset too large)
    sampled = random.sample(all_entries, min(MAX_BGC_SAMPLE, len(all_entries)))
    n_sample_prots = sum(len(e["proteins"]) for e in sampled)
    print(f"\n[SAMPLE] Selected {len(sampled)} BGCs ({n_sample_prots} proteins) for chimera generation")

    # Generate chimeras
    print(f"[GENERATION] Building chimeric assembly lines (max {MAX_CHIMERAS})...")
    chimeras = generate_chimeras(sampled, gamma)
    print(f"  Generated {len(chimeras)} chimeric architectures")

    rescued_count = sum(1 for c in chimeras if c["rescued"])
    print(f"  Rescued via linker insertion: {rescued_count}")

    if not chimeras:
        print("[ERROR] No chimeras generated.")
        sys.exit(1)

    # Statistics
    djcs_values = [c["mean_djcs"] for c in chimeras]
    print(f"\n[SCORES] DJCS range: {min(djcs_values):.2f} – {max(djcs_values):.2f}")
    print(f"         Mean DJCS: {sum(djcs_values)/len(djcs_values):.2f}")
    print(f"         Top 10 chimeras:")
    for c in chimeras[:10]:
        print(f"           {c['chimera_id']}: DJCS={c['mean_djcs']:.2f} "
              f"({c['host_bgc']}[{c['host_compound'][:15]}] ← "
              f"{c['donor_protein']} from {c['donor_bgc']}[{c['donor_organism'][:25]}])"
              f"{' [RESCUED]' if c['rescued'] else ''}")

    # Permutation test — use the FULL protein pool
    print(f"\n[PERMUTATION TEST] Running 1000 MC iterations against {len(all_proteins)}-protein pool...")
    stats = permutation_test(chimeras, all_proteins, gamma, 1000)
    print(f"  Observed top DJCS: {stats['observed_top_djcs']}")
    print(f"  Null mean DJCS:    {stats['null_mean_djcs']} ± {stats['null_std_djcs']}")
    print(f"  P-value:           {stats['p_value']}")
    print(f"  Significance:      {stats['significance']}")

    # Write outputs
    stats["data_source"] = "MIBiG 4.0 (local)"
    stats["total_bgcs_in_database"] = len(all_entries)
    stats["bgcs_sampled_for_chimeras"] = len(sampled)
    stats["natural_boundaries_measured"] = len(natural_diffs)
    stats["natural_mean_delta_h"] = natural_baseline

    stats_path = os.path.join(RESULTS_DIR, "permutation_stats.json")
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"\n[OUTPUT] {stats_path}")

    # Write chimeras TSV
    top_n = min(500, len(chimeras))
    tsv_path = os.path.join(RESULTS_DIR, "chimera_scores.tsv")
    with open(tsv_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "ChimeraID", "HostBGC", "HostOrganism", "HostCompound",
            "DonorBGC", "DonorOrganism", "DonorProtein", "DonorProduct",
            "DonorLength", "SwapPosition", "AssemblySize",
            "MeanDJCS", "MinDJCS", "MaxDJCS", "MeanDeltaH",
            "Rescued", "RescueLinker", "ProteinIDs"
        ])
        for c in chimeras[:top_n]:
            writer.writerow([
                c["chimera_id"], c["host_bgc"], c["host_organism"],
                c["host_compound"], c["donor_bgc"], c["donor_organism"],
                c["donor_protein"], c["donor_product"], c["donor_length"],
                c["swap_position"], c["assembly_size"],
                c["mean_djcs"], c["min_djcs"], c["max_djcs"], c["mean_delta_h"],
                c["rescued"], c["rescue_linker"], ";".join(c["protein_ids"])
            ])
    print(f"[OUTPUT] {tsv_path} (top {top_n} chimeras)")

    # Write detailed JSON for top 50
    detail_path = os.path.join(RESULTS_DIR, "chimera_details.json")
    with open(detail_path, "w") as f:
        json.dump(chimeras[:50], f, indent=2)
    print(f"[OUTPUT] {detail_path} (top 50 detailed)")

    # Write database summary
    summary_path = os.path.join(RESULTS_DIR, "database_summary.json")
    biosyn_class_counts = {}
    organism_counts = {}
    for e in all_entries:
        for cls in e.get("biosyn_class", []):
            biosyn_class_counts[cls] = biosyn_class_counts.get(cls, 0) + 1
        org = e.get("organism", "Unknown")
        organism_counts[org] = organism_counts.get(org, 0) + 1
    top_organisms = sorted(organism_counts.items(), key=lambda x: -x[1])[:20]
    with open(summary_path, "w") as f:
        json.dump({
            "mibig_version": "4.0",
            "total_bgcs": len(all_entries),
            "total_proteins": len(all_proteins),
            "protein_length_range": [min(lengths), max(lengths)],
            "protein_length_mean": round(sum(lengths)/len(lengths), 1),
            "biosynthetic_class_distribution": biosyn_class_counts,
            "top_organisms": dict(top_organisms),
        }, f, indent=2)
    print(f"[OUTPUT] {summary_path}")

    print(f"\n{'=' * 72}")
    print("[COMPLETE] All results in results/ — data source: MIBiG 4.0 (local)")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
