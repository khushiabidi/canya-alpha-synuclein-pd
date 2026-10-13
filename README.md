Overview

This project asks: does a single amino acid substitution in α-synuclein alter CANYA's predicted propensity to aggregate? Seven Parkinson's disease-associated missense mutations were tested against the wild type using the CANYA model on the Vecura platform. Window-level Δ scores were calculated for each variant, ranked by magnitude, and compared with published experimental aggregation data.
Background
α-Synuclein (140 aa, UniProt P37840) is the primary protein implicated in Parkinson's disease. Its aggregation into Lewy bodies is the pathological hallmark of PD. Familial PD mutations cluster in the N-terminal/NAC region (residues 1–95).
CANYA (Convolution Attention Network for amYloid Aggregation) is a hybrid CNN-attention neural network trained on >100,000 random peptide sequences. It predicts amyloid nucleation propensity using a 20-residue sliding window, outputting a score per window. Higher score = greater predicted propensity for aggregation.
Vecura (NYB.AI) is the platform that runs CANYA. It hosts 200+ computational biology models and charges per token.

Repository Structure

canya-alpha-synuclein-pd/
├── README.md
├── data/
│   ├── a53t_scores.txt
│   ├── a53v_scores.txt
│   ├── a53e_scores.txt
│   ├── a30p_scores.txt
│   ├── e46k_scores.txt
│   ├── h50q_scores.txt
│   └── g51d_scores.txt
├── scripts/
│   └── canya_visualization.R
└── figures/
    └── canya_pd_variants.png

Each data/ file has three columns: window, WT_score, mutant_score. The R script reads from these files and produces the figures.

Methods & How to Run

Sequences tested: Wild-type α-synuclein (UniProt P37840) and 7 single-residue mutants: A30P, E46K, H50Q, G51D, A53T, A53V, A53E.
CANYA run: Each sequence submitted to Vecura's CANYA model. Window size = 20 aa, sliding by 1 residue. Output = one score per window.
Δ score: mutant_score − WT_score per window. Positive = mutation increases predicted aggregation propensity in that window.
Four metrics per variant:
 • Δ Mean — average change across all 20 windows
 • Max Δ — peak positive shift
 • % Windows Increased — fraction with Δ > 0
 • Peak Position — window start with largest absolute Δ
To reproduce the figures:

# Install once
install.packages(c("ggplot2", "patchwork", "dplyr", "scales"))
# Run
source("scripts/canya_visualization.R")
# Output saved to figures/canya_pd_variants.png

Results

CANYA correctly predicted the direction of change in aggregation for 5 of 7 variants (71%).

| Variant | Δ Mean  | Max Δ   | % Windows ↑ | Peak Position | Literature    | Match |
| A53V    | +0.0126 | +0.0234 |     95%     |      46       | ↑ accelerates |  Yes  |
| A53T    | +0.0020 | +0.0046 |     85%     |      35       | ↑ accelerates |  Yes  |
| H50Q    | −0.0021 | +0.0026 |     15%     |      49       | ↑ accelerates |   No  |
| A30P    | −0.0039 | −0.0011 |     0%      |      29       | ↓ slows       |  Yes  |
| E46K    | −0.0051 | +0.0017 |     5%      |      44       | ↑ accelerates |   No  |
| G51D    | −0.0075 | −0.0036 |     0%      |      48       | ↓ slows       |  Yes  |
| A53E    | −0.0145 | −0.0022 |     0%      |      46       | ↓ slows       |  Yes  |

Key Findings & Limitations

A53V produces the strongest signal (+0.0234 peak Δ), consistent with valine's higher hydrophobicity vs alanine. A53E produces the strongest negative signal (−0.0207 peak Δ), consistent with glutamate's charge disrupting the hydrophobic core.
CANYA misses E46K and H50Q. Both experimentally accelerate aggregation, but CANYA predicts a net decrease. These mutations act through electrostatic mechanisms — charge changes that alter membrane binding and intermolecular interactions — rather than direct hydrophobic changes. CANYA is trained on hydrophobic sequence features and does not capture charge-driven aggregation.
Limitation: CANYA's 20-residue window scores reflect only sequence-level hydrophobicity and β-sheet propensity. It does not model:
• Electrostatic interactions
• Membrane-mediated aggregation
• Post-translational modifications
• Cellular context

References

• Flagmeier et al. (2016). Mutations associated with familial Parkinson's disease alter the initiation and amplification steps of α-synuclein aggregation. PNAS. PubMed
• Ghosh et al. (2013). The newly discovered Parkinson's disease associated Finnish mutation (A53E) attenuates α-synuclein aggregation. J Mol Biol.
• Lesage et al. (2013). G51D α-synuclein mutation causes a novel Parkinsonian-pyramidal syndrome. Ann Neurol. PubMed
• Khalaf et al. (2014). The H50Q mutation enhances α-synuclein aggregation, secretion, and toxicity. J Biol Chem. PubMed
• Vecura / CANYA platform: vecura.com



