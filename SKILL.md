---
name: generative-bgc-forge
description: Native AI Chimeric PKS/NRPS assembly line generative model using MIBiG 4.0. Implements Sequential Markov Decision Processes (MDP), Simulated Annealing, zero-dependency K-mer Markov log-likelihood scoring, and ESMFold. Exclusively uses Structural Interface Score (SIS). Use when the user mentions BGCs, combinatorial biosynthesis, AI-driven synthetic biology, structural compatibility (SIS), or automated chimera generation.
---

# GenerativeBGCs from MIBiG 4.0 Secondary Metabolites (MDP + SA + Markov)

Fully native pipeline for generating synthetic biosynthetic assembly lines. Parses the MIBiG 4.0 database, uses an intelligent Markov Decision Process (MDP) to structurally grow compatible sequences left-to-right, refines boundaries using thermodynamic Simulated Annealing on multivariate biophysics (SIS), and substitutes tailoring genes using compute-light TF-IDF. Evaluated chimeras undergo native offline K-mer Markov plausibility checking and ESMFold 3D API boundary pulls before final synthesis to `.gbk`. Repo: https://github.com/yzjie6/GenerativeBGCs

## When the user asks about this pipeline, Claude should:

- **Always confirm** the required data is present: MIBiG 4.0 raw json and fasta in the `data/` directory.
- **Always ask** what target backbone (e.g., specific PKS/NRPS class) they wish to anchor the assembly on.
- **Flag immediately** if user wants to use complex non-modular pathways, as the junction compatibility logic is optimized for modular synthase clusters.
- **Recommend preflight checks first**: run `fetch_mibig_data.py` to cryptographically verify data via SHA-256 before inference.
- **Distinguish** stochastic MC baseline scores from post-Annealing AI optimal scores — users frequently conflate raw matching with thermodynamic stabilization.
- **Ask about downstream expression** — the framework outputs sequence (.gbk) files theoretically optimized for E. coli chassis, but physical expression may require codon optimization.

## Required Tools

| Tool | Version | Purpose |
|------|---------|---------|
| Python | >= 3.7 | Main generative orchestration (Zero-dependency core) |
| Standard Library | all | json, math, random, os, hashlib, csv, collections |
| None | N/A | Strictly 100% native Python standard library. DeepBGC Docker completely excised. |

Quick install: The core platform is mathematically zero-dependency. Clone repository and verify data by running `python fetch_mibig_data.py`. No virtual environments needed.

**Critical env vars after install:**
None required. Structural validation is entirely native and ESMFold is hit dynamically via urllib.

## Pipeline Structure

The pipeline is split into cryptographically secure data parsing, autonomous RL assembly, and statistical validation.

- **Parser path:** Read MIBiG JSON/FASTA → compute multivariate SIS (Hydropathy + Charge) → SHA-256 lock
- **Assembly path:** Host template selection → MDP sequence growth → SA linker insertion → NLP TF-IDF tailoring
- **Validation path:** Offline native K-mer Markov transition filter → ESMFold 3D API boundary fetch → final sequence GenBank format generation

The `ablation_and_statistics.py` suite independently verifies the RL/SA logic via paired permutation tests against Monte Carlo models.

## Execution Order
```bash
python fetch_mibig_data.py         # parse data and hash lock
python ablation_and_statistics.py  # verify statistical integrity 
python main.py                     # execute AI generation and emit .gbk
```

## AI Generative Logic

Instead of stochastic random walks, the framework uses explicit algorithms to solve the combinatorial explosion of the 46,000 constituent proteins mapping:

1. **Sequential Markov Decision Process (MDP)** — Builds structurally left-to-right via an epsilon-greedy algorithm, incrementally appending heterologous domains based upon contextual evaluations of preceding C-termini.
2. **Multivariate Structural Interface Score (SIS)** — Calculates mathematically robust string similarity via Kyte-Doolittle hydropathy differentials combined with deterministic electrostatic point repulsions.
3. **Simulated Annealing (SA)** — When SIS < 70, SA inserts a flexible linker and accepts suboptimal boundary states strictly under the Boltzmann probability.
4. **TF-IDF NLP substitution** — Employs a text-based semantic alignment to preserve absolute offline compute independence for tailoring enzymes.
5. **K-mer Markov & ESMFold** — Validates generation internally via custom standard-library Di-peptide log-likelihood models, followed by a zero-dependency remote Meta ESM Atlas `.pdb` pull to verify 3D limits.

## Parameter Sensitivity Analysis

Users can evaluate different hyperparameter combinations in the ablation study to observe distributional robustness under varied physical constraints. 

Empirical results from MIBiG statistical tests (RL + SA):

| T (Cooling Temp) | C (Exploration) | Mean DJCS | Recommended use |
|------------------|----------------|-----------|-----------------|
| 0.80 | 5.0 | ~90.5 | Rapid high-exploitation convergence |
| 0.85 | 25.0 | ~90.8 | General-purpose stochastic assembly (Default) |
| 0.90 | 50.0 | ~90.6 | High exploration, aggressive thermodynamic search |

The algorithmic outputs empirically deviate by less than 0.7 units across all extreme tested parameters, signifying strict thermodynamic buffering. Actual values depend on cluster families chosen for combination.

## Key Diagnostic Decisions

**DeepBGC vs Pure Inference — which is better?**
- The core algorithm is highly conservative and analytically explicit. DeepBGC provides an independent Random Forest assessment. Leaving DeepBGC on ensures double validation, but if Docker fails or limits execution, the deterministic output alone remains statistically sound.

**DJCS score lower than expected?**
- Check the input classes. Bridging wildly distinct families (e.g., pure NRPS with pure PKS without linker regions) causes sharp hydropathy clash. The algorithm will aggressively try to Anneal linkers, but fundamental bio-physics dictates a ceiling on poorly-matched boundaries.

**NLP TF-IDF finds no tailoring matches?**
- Ensure MIBiG annotations are present. Uncharacterized cryptic clusters lack the functional text data the TF-IDF vectorizer requires to match evolutionary analogues. 

## Common Issues

| Symptom | Cause | Fix |
|---------|-------|-----|
| SHA-256 error on startup | Incomplete MIBiG file | Verify raw fasta exists in `data/` directory |
| DeepBGC evaluation returns 0.00 | Docker daemon not running | Start docker daemon (`systemctl start docker` or `open /Applications/Docker.app`) |
| Ablation test too slow | Large dataset parsing | Subset constraints automatically applied in `ablation.py` |
| TF-IDF throws division by zero | Empty token arrays | Handled gracefully via try-rescue logic in main parser |

## References

- `ablation_and_statistics.py` — Complete parameter sweeps, empirical robustness checks, and p-value generation
- Auer, P., Cesa-Bianchi, N., & Fischer, P. (2002). Finite-time analysis of the multiarmed bandit problem. *Machine learning*, 47(2), 235-256.
- Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983). Optimization by simulated annealing. *Science*, 220(4598), 671-680.
- Medema, M. H., et al. (2015). Minimum information about a biosynthetic gene cluster. *Nature chemical biology*, 11(9), 625-631.
- Hannigan, G. D., et al. (2019). A deep learning genome-mining strategy for biosynthetic gene cluster prediction. *Nucleic acids research*, 47(18), e110-e110.
