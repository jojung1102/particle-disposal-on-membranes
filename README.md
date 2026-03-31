# Membrane bending, bifurcations, and flow

This folder contains MATLAB scripts for bifurcation figures and ODE continuation, Jupyter notebooks for figures and analysis, and Mathematica notebooks for membrane geometry illustrations. 

---

## Folder structure

| Path | Role |
|------|------|
| **`output/`** | Numerical outputs (CSV). Created when you run the MATLAB scripts or when notebooks run their setup cell (`os.makedirs`). Typical subfolders: `Bifurcations`, `trajectories_flow`, `score`, `score_trajectories`. |
| **`figures/`** | Exported figures. Notebooks write under **`figures/Paper/`** (SVG/PDF as configured in each notebook). |
| **`mathematica figures/`** | Mathematica (`.nb`) sources for membrane shapes and related graphics. |

---

## Jupyter notebooks (`.ipynb`)

| Notebook | Content |
|----------|--------|
| **`bifurcation_flow_paper.ipynb`** | Phase plane and bifurcation figures for the **2D** reduced model (no particles): shape space of $(A_M, A_C)$, continuation branches, shape flow diagrams and trajectories. |
| **`bifurcation_flow_particles_paper.ipynb`** | Same style of analysis with **particles**. Trajectories are also available for particle numbers.  |
| **`Score_paper.ipynb`** | **Score** plots: trajectory ensembles from `output/score/` and `output/score_trajectories/`. |
| **`number_equilib.ipynb`** | Equilibrium **particle fractions** $N_{i,x}/N_i$ as functions of **$A_M/A$** (no flow fields). |

---

## Workflow

1. Run the relevant **MATLAB** continuation or trajectory scripts to populate **`output/`**.
2. Open the matching **notebook**, run the setup cell, then run cells or “Run All” to regenerate **`figures/Paper/`**.