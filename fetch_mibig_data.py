#!/usr/bin/env python3
"""
fetch_mibig_data.py — Parse Local MIBiG 4.0 Data (Zero Network Access)

Reads from the locally downloaded MIBiG 4.0 files:
  - data/mibig_prot_seqs_4.0.fasta  (all protein sequences)
  - data/mibig_json_4.0/*.json       (BGC metadata)

Extracts PKS and NRPS mega-enzyme proteins (>500 aa) with their BGC
context. Outputs a curated cache at data/assembly_proteins.json.

Zero external dependencies. No internet required.
"""

import json
import os
import sys
import hashlib

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FASTA_PATH = os.path.join(SCRIPT_DIR, "data", "mibig_prot_seqs_4.0.fasta")
JSON_DIR = os.path.join(SCRIPT_DIR, "data", "mibig_json_4.0")
OUTPUT_PATH = os.path.join(SCRIPT_DIR, "data", "assembly_proteins.json")

# Minimum protein length to be considered an assembly-line mega-enzyme
MIN_LENGTH = 500


def parse_fasta(fasta_path):
    """
    Parse the MIBiG FASTA file into a dict keyed by BGC ID.

    Header format:
    >BGC0000001.5|c1|32476-43413|+|AEK75503.1|type_1_polyketide_synthase|AEK75503.1
    Fields: bgc_version | contig | coordinates | strand | protein_id | product | locus_tag
    """
    proteins_by_bgc = {}
    current_header = None
    current_seq_parts = []

    def flush():
        if current_header and current_seq_parts:
            seq = "".join(current_seq_parts)
            parts = current_header.split("|")
            if len(parts) >= 6:
                bgc_full = parts[0]  # e.g. "BGC0000001.5"
                bgc_id = bgc_full.split(".")[0]  # e.g. "BGC0000001"
                protein_id = parts[4] if parts[4] else parts[6]
                product = parts[5].replace("_", " ")
                coords = parts[2]

                if bgc_id not in proteins_by_bgc:
                    proteins_by_bgc[bgc_id] = []
                proteins_by_bgc[bgc_id].append({
                    "protein_id": protein_id,
                    "product": product,
                    "coordinates": coords,
                    "strand": parts[3],
                    "length": len(seq),
                    "sequence": seq,
                })

    print(f"[PARSE] Reading {fasta_path} ...")
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                flush()
                current_header = line[1:]  # Remove '>'
                current_seq_parts = []
            else:
                current_seq_parts.append(line)
    flush()  # Last entry

    return proteins_by_bgc


def load_bgc_metadata(json_dir, bgc_ids):
    """Load metadata for specific BGC IDs from MIBiG 4.0 JSON files."""
    metadata = {}
    for bgc_id in bgc_ids:
        json_path = os.path.join(json_dir, f"{bgc_id}.json")
        if not os.path.exists(json_path):
            continue
        with open(json_path) as f:
            data = json.load(f)
        # MIBiG 4.0 schema: top-level keys, not wrapped in 'cluster'
        biosyn_classes = []
        for bc in data.get("biosynthesis", {}).get("classes", []):
            biosyn_classes.append(bc.get("class", ""))
        taxonomy = data.get("taxonomy", {})
        compounds = data.get("compounds", [])
        loci = data.get("loci", {})
        
        # Extract unique bioactivities
        unique_acts = set()
        for cpd in compounds:
            for bio in cpd.get("bioactivities", []):
                name = bio.get("name", "")
                if isinstance(name, dict):
                    name = name.get("activity", "")
                if name:
                    unique_acts.add(name.lower())

        metadata[bgc_id] = {
            "organism": taxonomy.get("name", "Unknown"),
            "biosyn_class": biosyn_classes,
            "compounds": [c.get("name", "") for c in compounds],
            "activities": list(unique_acts),
            "loci_accession": loci.get("accession", "") if isinstance(loci, dict) else "",
        }
    return metadata


def main():
    if os.path.exists(OUTPUT_PATH):
        print(f"[CACHE] {OUTPUT_PATH} already exists. Delete to re-parse.")
        with open(OUTPUT_PATH) as f:
            data = json.load(f)
        n_prots = sum(len(e.get("core_proteins", [])) for e in data)
        print(f"  → {len(data)} BGCs, {n_prots} assembly-line proteins")
        return

    if not os.path.exists(FASTA_PATH):
        print(f"[ERROR] FASTA not found: {FASTA_PATH}")
        print("  Download MIBiG 4.0 from https://mibig.secondarymetabolites.org/download")
        sys.exit(1)

    # Cryptographic Reproducibility Enforcement
    EXPECTED_SHA256 = "4b196343ed82d49f2b8f3babe89d512b92f0e026040fa50b15eeb97950717e8b"
    print(f"[VERIFY] Checking SHA256 checksum of {FASTA_PATH}...")
    sha256_hash = hashlib.sha256()
    with open(FASTA_PATH, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    actual_hash = sha256_hash.hexdigest()
    
    if actual_hash != EXPECTED_SHA256:
        print(f"[WARNING] SHA256 mismatch!\n  Expected: {EXPECTED_SHA256}\n  Actual:   {actual_hash}")
        print("  This may compromise deterministic execution and reproducibility.")
    else:
        print("  [SUCCESS] Checksum verified. Dataset is cryptographically pristine.")

    # Step 1: Parse all proteins from FASTA
    proteins_by_bgc = parse_fasta(FASTA_PATH)
    total_prots = sum(len(v) for v in proteins_by_bgc.values())
    print(f"  Found {total_prots} structural and auxiliary proteins across {len(proteins_by_bgc)} BGCs")

    # Step 2: Load metadata for BGCs that have large proteins
    print(f"[META] Loading BGC metadata from {JSON_DIR} ...")
    metadata = load_bgc_metadata(JSON_DIR, proteins_by_bgc.keys())
    print(f"  Loaded metadata for {len(metadata)} BGCs")

    # Step 3: Filter to PKS / NRP / hybrid BGCs only
    assembly_data = []
    for bgc_id, proteins in sorted(proteins_by_bgc.items()):
        meta = metadata.get(bgc_id, {})
        biosyn_class = meta.get("biosyn_class", [])

        # Keep only BGCs with PKS or NRP classes (MIBiG 4.0 uses "PKS", "NRP")
        has_pks_nrps = any(
            cls in ["Polyketide", "NRP", "PKS"]
            for cls in biosyn_class
        )
        if not has_pks_nrps:
            continue

        # Keep only BGCs with at least 2 large proteins (assembly line)
        core_proteins = [p for p in proteins if p["length"] >= MIN_LENGTH]
        aux_proteins = [p for p in proteins if p["length"] < MIN_LENGTH]
        
        if len(core_proteins) < 2:
            continue

        assembly_data.append({
            "bgc_id": bgc_id,
            "organism": meta.get("organism", "Unknown"),
            "biosyn_class": biosyn_class,
            "compounds": meta.get("compounds", []),
            "activities": meta.get("activities", []),
            "loci_accession": meta.get("loci_accession", ""),
            "core_proteins": core_proteins,
            "aux_proteins": aux_proteins,
        })

    # Write output
    with open(OUTPUT_PATH, "w") as f:
        json.dump(assembly_data, f, indent=2)

    n_prots = sum(len(e["core_proteins"]) for e in assembly_data)
    n_aux = sum(len(e["aux_proteins"]) for e in assembly_data)
    print(f"\n{'=' * 70}")
    print(f"[DONE] Wrote {len(assembly_data)} PKS/NRP BGCs to {OUTPUT_PATH}")
    print(f"       Total core assembly-line proteins: {n_prots}")
    print(f"       Total auxiliary/tailoring proteins: {n_aux}")
    print(f"       File size: {os.path.getsize(OUTPUT_PATH) / (1024*1024):.1f} MB")
    print(f"{'=' * 70}")

    # Print summary table
    print(f"\n{'BGC ID':14s} | {'Organism':45s} | {'Class':20s} | {'Core':>4s} | {'Aux':>3s} | Compound")
    print("-" * 115)
    for e in assembly_data[:30]:
        cls = ", ".join(e["biosyn_class"])
        cpd = e["compounds"][0] if e["compounds"] else ""
        print(f"{e['bgc_id']:14s} | {e['organism'][:45]:45s} | {cls[:20]:20s} | {len(e['core_proteins']):4d} | {len(e['aux_proteins']):3d} | {cpd[:30]}")
    if len(assembly_data) > 30:
        print(f"  ... and {len(assembly_data) - 30} more BGCs")


if __name__ == "__main__":
    main()
