---
name: generative-bgc-forge
description: AI-Enhanced Chimeric PKS/NRPS assembly line screening pipeline using MIBiG 4.0. Implements Reinforcement Learning (UCB1), Simulated Annealing, and NLP TF-IDF semantic matching. Includes rigorous statistical validation (Bootstrap, Permutation Testing) to demonstrate distributional robustness, secured by SHA-256 data integrity. Use when the user mentions BGCs, combinatorial biosynthesis, AI-driven synthetic biology, structural compatibility (DJCS), or automated chimera generation.
---

# GenerativeBGCs from MIBiG 4.0 Secondary Metabolites (UCB1 + SA + TF-IDF)

Fully offline pipeline for generating synthetic biosynthetic assembly lines. Parses the complete MIBiG 4.0 database, uses Multi-Armed Bandit (RL) to discover compatible components, refines inter-protein boundaries using thermodynamic Simulated Annealing on hydropathy scores (DJCS), and substitutes tailoring genes using semantic NLP matching (TF-IDF). Evaluated chimeras undergo DeepBGC orthogonal neural network check before final output as `.gbk`. Repo: https://github.com/yzjie6/GenerativeBGCs

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
| Docker | optional | Only required for orthogonal `antibioti/deepbgc` validation triage |

Quick install: The core platform is zero-dependency. Clone repository and verify data by running `python fetch_mibig_data.py`. No virtual environments needed.

**Critical env vars after install:**
None required for pure-Python execution. Docker daemon must be running in the background for DeepBGC verification stage.

## Pipeline Structure

The pipeline is split into cryptographically secure data parsing, autonomous RL assembly, and statistical validation.

- **Parser path:** Read MIBiG JSON/FASTA → compute structural hydropathy boundaries → SHA-256 lock
- **Assembly path:** Host template selection → UCB1 RL donor exploration → SA linker insertion → NLP TF-IDF tailoring gene transplantation
- **Validation path:** Native DeepBGC neural network processing → final sequence GenBank format generation

The `ablation_and_statistics.py` suite independently verifies the RL/SA logic via paired permutation tests against Monte Carlo models.

## Execution Order
```bash
python fetch_mibig_data.py         # parse data and hash lock
python ablation_and_statistics.py  # verify statistical integrity 
python main.py                     # execute AI generation and emit .gbk
```

## AI Generative Logic

Instead of stochastic random walks, the framework uses explicit algorithms to solve the combinatorial explosion of the 46,000 constituent proteins mapping:

1. **UCB1 (Exploration vs Exploitation)** — Treats potential donor sequences as independent arms in a stochastic bandit problem to optimize discovery of highly compatible parts, escaping the random assembly trap.
2. **Domain Junction Compatibility Score (DJCS)** — Evaluates hydropathy differentials at protein boundaries to predict folding disruption. 
3. **Simulated Annealing (SA)** — When DJCS < 70, SA inserts a flexible linker and accepts suboptimal boundary states strictly under the Boltzmann probability $e^{-\Delta E / T}$. This prevents the sequence from collapsing into local minima.
4. **TF-IDF NLP substitution** — Selects valid active downstream tailoring enzymes using natural language processing (tokenized dropping fragments ≤ 2 chars and applying a cosine-similarity activation threshold of 0.40) over evolutionary functional notes.
5. **DeepBGC ML Validation** — As an orthogonal truth check, sequences are piped through a pre-trained offline DeepBGC docker container to evaluate mean Random Forest active-probability.

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
