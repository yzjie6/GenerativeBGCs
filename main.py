#!/usr/bin/env python3
"""
main.py — Interactive CLI for GenerativeBGCs v6 (AI-Enhanced)

Includes three zero-dependency lightweight AI models:
1. UCB1 (Multi-Armed Bandit Reinforcement Learning) for donor discovery.
2. Simulated Annealing for thermodynamic linker optimization.
3. TF-IDF Cosine Similarity for NLP-based tailoring gene substitution.
"""

import json
import os
import sys
import math
import random
import shutil
import subprocess
from collections import Counter
from combinatorial_assembly import (
    load_bgc_data,
    calibrate_gamma,
    boundary_hydropathy_diff,
    junction_compatibility_score,
    FLEXIBLE_LINKERS,
    RESULTS_DIR
)
from gbk_writer import export_chimeras_to_gbk

def prompt_menu(title, options):
    print(f"\n{'=' * 60}")
    print(f" {title}")
    print(f"{'=' * 60}")
    for i, (name, count) in enumerate(options, 1):
        if count is not None:
            print(f"  [{i:2d}] {name:35s} ({count} allowed options)")
        else:
            print(f"  [{i:2d}] {name:35s}")
    print(f"  [{len(options)+1:2d}] Any / Skip filter")
    while True:
        try:
            choice = input(f"\nEnter your choice (1-{len(options)+1}): ").strip()
            idx = int(choice) - 1
            if 0 <= idx < len(options):
                return options[idx][0]
            elif idx == len(options):
                return None
            else:
                print("Invalid selection.")
        except ValueError:
            print("Please enter a valid number.")

def prompt_yes_no(prompt):
    while True:
        ans = input(prompt + " [y/n]: ").strip().lower()
        if ans in ('y', 'yes'): return True
        if ans in ('n', 'no'): return False

# --- AI CORE 1: NLP TF-IDF Vectorizer ---
def compute_tfidf_similarity(query_text, documents):
    """Zero-dependency TF-IDF engine returning cosine similarities."""
    if not documents: return []
    texts = [query_text] + documents
    tokenized = [t.lower().replace('-', ' ').replace(',', '').split() for t in texts]
    vocab = list(set([word for doc in tokenized for word in doc if len(word) > 2]))
    
    tfs = []
    for doc in tokenized:
        tf_dict = Counter(doc)
        doc_len = max(len(doc), 1)
        tfs.append({word: count / doc_len for word, count in tf_dict.items()})
        
    N = len(tokenized)
    idfs = {}
    for word in vocab:
        df = sum(1 for doc in tokenized if word in doc)
        idfs[word] = math.log((1 + N) / (1 + df)) + 1
        
    vectors = []
    for tf_dict in tfs:
        vec = []
        norm = 0
        for word in vocab:
            val = tf_dict.get(word, 0) * idfs[word]
            vec.append(val)
            norm += val * val
        norm = math.sqrt(norm) if norm > 0 else 1
        vectors.append([v / norm for v in vec])
        
    query_vec = vectors[0]
    sims = []
    for doc_vec in vectors[1:]:
        dot = sum(q * v for q, v in zip(query_vec, doc_vec))
        sims.append(dot)
    return sims

# --- AI CORE 2: Simulated Annealing ---
def simulated_annealing_linker(seq_up, seq_down, gamma, initial_djcs):
    """Thermodynamic optimization of inter-protein boundaries."""
    best_linker = ""
    current_djcs = initial_djcs
    best_djcs = initial_djcs
    
    T = 5.0 # Initial Temperature
    cooling_rate = 0.85
    
    for linker_seq, linker_name in FLEXIBLE_LINKERS:
        new_up = seq_up + linker_seq
        new_diff = boundary_hydropathy_diff(new_up, seq_down)
        new_djcs = junction_compatibility_score(new_diff, gamma)
        
        delta_E = new_djcs - current_djcs # DJCS is like negative energy; higher is better
        
        if delta_E > 0 or math.exp(delta_E / T) > random.random():
            current_djcs = new_djcs
            if new_djcs > best_djcs:
                best_djcs = new_djcs
                best_linker = linker_name
        T *= cooling_rate
        
    return best_linker, best_djcs


# --- AI CORE 3: Reinforcement Learning (UCB1) Integration ---
def generate_targeted_chimeras(host_candidates, donor_candidates, gamma, do_tailoring=False, max_chimeras=2000):
    chimeras = []
    if not host_candidates or not donor_candidates: return []
    
    # UCB1 State Initialization
    counts = [0] * len(donor_candidates)
    values = [0.0] * len(donor_candidates)
    C_EXPLORE = 25.0
    
    for t in range(max_chimeras):
        # 1. RL Agent: Select Donor BGC Arm
        if t < len(donor_candidates):
            d_idx = t
        else:
            ucb_values = [
                values[i] + C_EXPLORE * math.sqrt(math.log(t) / counts[i])
                for i in range(len(donor_candidates))
            ]
            d_idx = ucb_values.index(max(ucb_values))
            
        donor_entry = donor_candidates[d_idx]
        host_entry = random.choice(host_candidates)
        
        prots_host = host_entry.get("core_proteins", [])
        host_aux = host_entry.get("aux_proteins", [])
        prots_donor = donor_entry.get("core_proteins", [])
        donor_aux = donor_entry.get("aux_proteins", [])
        
        if len(prots_host) < 2 or not prots_donor or host_entry["bgc_id"] == donor_entry["bgc_id"]:
            counts[d_idx] += 1
            if counts[d_idx] > 1: values[d_idx] = (values[d_idx] * (counts[d_idx]-1)) / counts[d_idx]
            continue
            
        # Chimera Scaffold Construction
        pos = random.randint(0, len(prots_host)-1)
        donor_core = random.choice(prots_donor)
        
        chimeric_line = list(prots_host)
        chimeric_line[pos] = donor_core
        
        # Physicochemical Scoring
        boundary_scores = []
        for k in range(len(chimeric_line) - 1):
            d = boundary_hydropathy_diff(chimeric_line[k]["sequence"], chimeric_line[k + 1]["sequence"])
            boundary_scores.append(junction_compatibility_score(d, gamma))
            
        mean_djcs = sum(boundary_scores) / max(1, len(boundary_scores))
        
        # Simulated Annealing Trigger
        rescued = False
        rescued_linker = ""
        if boundary_scores:
            worst_idx = boundary_scores.index(min(boundary_scores))
            worst_djcs = boundary_scores[worst_idx]
            if worst_djcs < 70.0:
                linker_name, new_djcs = simulated_annealing_linker(
                    chimeric_line[worst_idx]["sequence"],
                    chimeric_line[worst_idx + 1]["sequence"],
                    gamma, worst_djcs
                )
                if linker_name and new_djcs > worst_djcs + 1.0:
                    rescued = True
                    rescued_linker = linker_name
                    boundary_scores[worst_idx] = new_djcs
                    mean_djcs = sum(boundary_scores) / len(boundary_scores)

        # RL Reward Update
        reward = mean_djcs
        counts[d_idx] += 1
        n = counts[d_idx]
        values[d_idx] = ((n - 1) * values[d_idx] + reward) / n

        # NLP Tailoring Substitution
        final_aux = list(host_aux)
        swapped_aux_ids = []
        if do_tailoring and host_aux and donor_aux:
            valid_host_aux = [a for a in host_aux if 'hypothetical' not in a['product'].lower() and len(a['product']) > 5]
            if valid_host_aux:
                h_target = random.choice(valid_host_aux)
                h_idx = host_aux.index(h_target)
                
                d_docs = [a["product"] for a in donor_aux]
                sims = compute_tfidf_similarity(h_target["product"], d_docs)
                
                best_sim = max(sims) if sims else 0
                if best_sim > 0.4:  # High confidence NLP match
                    best_match_idx = sims.index(best_sim)
                    found_match = donor_aux[best_match_idx]
                    
                    final_aux[h_idx] = found_match
                    swapped_aux_ids.append(found_match["protein_id"])
                    
        chimeras.append({
            "chimera_id": f"TGT_{len(chimeras):05d}",
            "host_bgc": host_entry["bgc_id"],
            "host_organism": host_entry["organism"],
            "host_compound": host_entry["compounds"][0] if host_entry["compounds"] else "",
            "donor_bgc": donor_entry["bgc_id"],
            "donor_organism": donor_entry["organism"],
            "swap_position": pos,
            "donor_protein": donor_core["protein_id"],
            "donor_product": donor_core["product"],
            "donor_length": donor_core["length"],
            "assembly_size": len(chimeric_line),
            "mean_djcs": round(mean_djcs, 4),
            "rescued": rescued,
            "rescue_linker": rescued_linker,
            "protein_ids": [p["protein_id"] for p in chimeric_line],
            "aux_proteins": final_aux,
            "donor_aux": swapped_aux_ids
        })
        
    chimeras.sort(key=lambda x: x["mean_djcs"], reverse=True)
    return chimeras

def apply_deepbgc_scoring(chimeras, all_entries, top_n=20):
    """Uses external DeepBGC docker to recalculate ranking for the top candidates."""
    print(f"\n[ML EVALUATOR] Transferring Top {top_n} candidates to Deep Learning neural network...")
    
    top_candidates = chimeras[:top_n]
    tmp_dir = os.path.join(RESULTS_DIR, ".tmp_scoring")
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)
    
    # 1. Export temporary sequences
    from gbk_writer import export_chimeras_to_gbk
    export_chimeras_to_gbk(top_candidates, all_entries, output_dir=tmp_dir, top_n=top_n)
    
    # 2. Run sequential Docker validation
    for i, chimera in enumerate(top_candidates):
        gbk_file = f"{chimera['chimera_id']}.gbk"
        gbk_path = os.path.join(tmp_dir, gbk_file)
        out_path = os.path.join(tmp_dir, f"{chimera['chimera_id']}_out")
        
        final_score = 0.0
        if os.path.exists(gbk_path):
            print(f"               Evaluating {chimera['chimera_id']} [{i+1}/{top_n}] via DeepBGC...")
            cmd = [
                "docker", "run", "--rm",
                "-v", f"{os.getcwd()}:/data",
                "antibioti/deepbgc", "pipeline",
                f"/data/results/.tmp_scoring/{gbk_file}",
                "-o", f"/data/results/.tmp_scoring/{chimera['chimera_id']}_out"
            ]
            try:
                subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
                bgc_csv = os.path.join(out_path, f"{chimera['chimera_id']}.bgc.csv")
                if os.path.exists(bgc_csv):
                    with open(bgc_csv, 'r') as f:
                        lines = f.readlines()
                        if len(lines) > 1:
                            parts = lines[1].split(',')
                            scores = []
                            for p in parts:
                                try: scores.append(float(p))
                                except: pass
                            if scores: final_score = max(scores)
            except Exception:
                pass
        chimera["deepbgc_score"] = float(f"{final_score:.4f}")
        
    shutil.rmtree(tmp_dir, ignore_errors=True)
    
    # 3. Re-sort
    chimeras[:top_n] = sorted(top_candidates, key=lambda x: (x.get("deepbgc_score", 0.0), x["mean_djcs"]), reverse=True)
    print("               -> Structural Verification Complete.")
    return chimeras

def print_chimera_summary(chimeras, top_n=5):
    print(f"\n[SCORES] RL Agent generated & explored {len(chimeras)} customized chimeric architectures.")
    if not chimeras: return
        
    djcs_values = [c["mean_djcs"] for c in chimeras]
    print(f"         DJCS range: {min(djcs_values):.2f} – {max(djcs_values):.2f}")
    print(f"         Mean DJCS: {sum(djcs_values)/len(djcs_values):.2f}")
    print(f"         Top {top_n} RL-optimized designs:")
    
    for c in chimeras[:top_n]:
        host_name = f"{c['host_bgc']}[{c['host_compound'][:15]}]"
        donor_name = f"{c['donor_protein']} from {c['donor_bgc']}[{c['donor_organism'][:20]}]"
        
        tags = []
        if c['rescued']: tags.append(f"SA_Linker({c['rescue_linker']})")
        if c.get("donor_aux"): tags.append("NLP_TailorSwap")
        if "deepbgc_score" in c: tags.append(f"DeepBGC={c['deepbgc_score']:.2f}")
        
        tag_str = f" [+ {', '.join(tags)}]" if tags else ""
        print(f"           {c['chimera_id']}: DJCS={c['mean_djcs']:.2f} "
              f"({host_name} ← {donor_name}){tag_str}")

def main():
    print(r"""
   _______  _______  _______  _______  _______  _______  _______  _______  ______   _______  _______ 
  (  ____ \(  ____ \(  ____ \(  ____ \(  ____ )(  ___  )(  ____ \(  ____ \/ ___  \ (  ___  )(  ____ \
  | (    \/| (    \/| (    \/| (    \/| (    )|| (   ) || (    \/| (    \/\/   )  )| (   ) || (    \/
  | |      | (__    | |      | (__    | (____)|| (___) || |      | (__      /  /  | (___) || (_____ 
  | | ____ |  __)   | |      |  __)   |     __)|  ___  || +---.  |  __)    /  /   |  ___  |(_____  )
  | | \_  )| (      | |      | (      | (\ (   | (   ) || |      | (      /  /    | (   ) |      ) |
  | (___) || (____/\| (____/\| (____/\| ) \ \__| )   ( || (____/\| )     /  /____/\ )   ( |/\____) |
  (_______)(_______/(_______/(_______/|/   \__/|/     \|(_______/|/      \_______/|/     \|\_______)
     Generative BGC Platform v6 — AI-Enhanced Targeted Chimera Design
    """)

    print("Loading MIBiG database cache...")
    try:
        all_entries = load_bgc_data()
    except SystemExit:
        print("\nPlease run 'python fetch_mibig_data.py' to parse the database first.")
        return

    print("Calibrating structural scoring kernel (DJCS)...")
    gamma, _, _ = calibrate_gamma(all_entries)
    print(f"  → γ = {gamma}\n")
    
    # 0. Biosynthetic Class Filter
    class_opts = [
        ("PKS (Polyketide Synthase)", sum(1 for e in all_entries if "PKS" in str(e["biosyn_class"]).upper() or "POLYKETIDE" in str(e["biosyn_class"]).upper())),
        ("NRPS (Non-Ribosomal Peptide)", sum(1 for e in all_entries if "NRP" in str(e["biosyn_class"]).upper())),
        ("PKS-NRPS Hybrids", sum(1 for e in all_entries if ("PKS" in str(e["biosyn_class"]).upper() or "POLYKETIDE" in str(e["biosyn_class"]).upper()) and "NRP" in str(e["biosyn_class"]).upper()))
    ]
    sel_class = prompt_menu("Select BIOSYNTHETIC CLASS", class_opts)
    if sel_class:
        if "Hybrids" in sel_class:
            all_entries = [e for e in all_entries if ("PKS" in str(e["biosyn_class"]).upper() or "POLYKETIDE" in str(e["biosyn_class"]).upper()) and "NRP" in str(e["biosyn_class"]).upper()]
        elif "PKS" in sel_class:
            all_entries = [e for e in all_entries if "PKS" in str(e["biosyn_class"]).upper() or "POLYKETIDE" in str(e["biosyn_class"]).upper()]
        elif "NRPS" in sel_class:
            all_entries = [e for e in all_entries if "NRP" in str(e["biosyn_class"]).upper()]
            
    # 1. Activities
    activity_counts = Counter()
    genus_counts = Counter()
    for entry in all_entries:
        for act in entry.get("activities", []): activity_counts[act] += 1
        org = entry.get("organism", "").strip()
        if org and org.lower() != "unknown":
            parts = org.split()
            if parts: genus_counts[parts[0]] += 1
            
    top_activities = activity_counts.most_common(12)
    selected_activity = prompt_menu("Select TARGET BIOACTIVITY for synthetic chimera", top_activities)
    
    donor_candidates = all_entries
    if selected_activity:
        donor_candidates = [e for e in all_entries if selected_activity in e.get("activities", [])]
        print(f"\n[FILTER] Found {len(donor_candidates)} BGCs encoding '{selected_activity}' pathways (Donors).")
        
    top_genera = genus_counts.most_common(10)
    selected_chassis = prompt_menu("Select HOST CHASSIS GENUS for expression", top_genera)
    
    host_candidates = all_entries
    if selected_chassis:
        host_candidates = [e for e in all_entries if e.get("organism", "").startswith(selected_chassis)]
        print(f"\n[FILTER] Found {len(host_candidates)} matching backbones.")
        
    if host_candidates:
        molecule_opts = []
        for e in host_candidates:
            cpd = e["compounds"][0] if e["compounds"] else "Unknown"
            n_name = f"{e['bgc_id']} [{cpd}]"
            molecule_opts.append((n_name, None))
        
        molecule_opts.sort(key=lambda x: x[0])
        sel_molecule = prompt_menu("Select EXACT HOST BACKBONE (Molecule) to modify", molecule_opts[:25])
        if sel_molecule:
            target_id = sel_molecule.split()[0]
            host_candidates = [e for e in host_candidates if e["bgc_id"] == target_id]
            print(f"\n[TARGET] Locked Host Backbone: {target_id}")

    if not host_candidates or not donor_candidates:
        print("\n[ERROR] Filter conditions too strict. No combinations possible.")
        return

    do_tailoring = prompt_yes_no("\nEnable NLP-based Tailoring Gene Integration (TF-IDF Match)?")
    do_deepbgc   = prompt_yes_no("Enable DeepBGC Neural Network Verification & ML Sorting?")

    print(f"\n[AI AGENT] Running Multi-Armed Bandit (UCB1) RL exploration over {len(donor_candidates)} candidates...")
    chimeras = generate_targeted_chimeras(host_candidates, donor_candidates, gamma, do_tailoring=do_tailoring)
    
    if do_deepbgc:
        chimeras = apply_deepbgc_scoring(chimeras, all_entries, top_n=20)
        
    print_chimera_summary(chimeras)
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_file = os.path.join(RESULTS_DIR, "targeted_chimeras.json")
    
    clean_chimeras = []
    for c in chimeras[:50]:
        cc = dict(c)
        if "aux_proteins" in cc: del cc["aux_proteins"]
        clean_chimeras.append(cc)
        
    with open(out_file, "w") as f:
        json.dump(clean_chimeras, f, indent=2)
        
    print(f"\n[COMPLETE] Detailed JSON designs saved to: {out_file}")
    
    print("\n[EXPORT] Reverse-translating top 10 designs into synthetic full-cluster GenBank DNA (.gbk)...")
    gbk_dir = os.path.join(RESULTS_DIR, "gbk")
    count, _ = export_chimeras_to_gbk(chimeras, all_entries, output_dir=gbk_dir, top_n=10)
    print(f"         Successfully exported {count} synthesis-ready GenBank files to: {gbk_dir}/")

if __name__ == "__main__":
    main()
