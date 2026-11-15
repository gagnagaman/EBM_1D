# 1-D Energy Balance Model (Python Implementation)

This directory contains a complete Python implementation of the 1-D diffusive energy balance model (EBM) derived in **[`Docs/EBM_eng3.pdf`](Docs/EBM_eng3.pdf)** (*"Derivation of an 1-D energy balance model with diffusive heat transports"*), translating and modernizing the R scripts in `../R_EBM/`.

---

## 1. Directory Structure

```
EBM_1D/
├── Docs/
│   └── EBM_eng3.pdf                   # Theoretical notes and derivation by Halldór Björnsson
├── figures/
│   ├── ebm_assignment1_results.png    # 4-panel comparison: baseline, D sensitivity, albedo(x), Lindzen-Farrell
│   └── ebm_sici_bifurcation.png       # SICI hysteresis loop: global mean T vs. S_0
├── Matlab_code/                       # MATLAB / GNU Octave implementation
│   ├── ebm.m                          # Baseline 1-D EBM solver
│   ├── ebmAlbedo.m                    # Temperature-dependent albedo solver
│   ├── IterTalb.m                     # Equilibrium relaxation iterator
│   ├── SICI.m                         # SICI hysteresis experiment driver
│   └── README.md                      # MATLAB documentation and usage guide
├── R_code/                            # R implementation
│   ├── ebm.R                          # Baseline 1-D EBM solver
│   ├── analyT.R                       # Exact analytical Legendre polynomial solution
│   ├── ebmAlbedo.R                    # Temperature-dependent albedo solver
│   ├── iterate_ebmAlbedo.R            # Equilibrium relaxation iterator
│   ├── SICI.R                         # SICI hysteresis experiment driver
│   ├── readncep.R                     # Symmetrized NCEP data loader
│   ├── test_ebm_ncep.R                # Assignment 1 test script
│   ├── ncep_temp.txt                  # Local observational NCEP dataset
│   └── README.md                      # R documentation and usage guide
├── ncep_temp.txt                      # Monthly longterm mean air temperature (NCEP Reanalysis)
├── EBM_1D_Analysis.ipynb              # Interactive Jupyter Notebook mirroring all derivations and experiments
├── ebm_model.py                       # Core model module and exact analytical Legendre solver
├── run_ebm.py                         # Driver script for Assignment 1 experiments
├── ebm_sici.py                        # Driver script for SICI hysteresis (Assignment 2)
└── README.md                          # Documentation (this file)
```

---

## 2. Description of Files

| File | Purpose |
|------|---------|
| **`EBM_1D_Analysis.ipynb`** | **Interactive Jupyter Notebook.** Comprehensive, cell-by-cell walkthrough containing baseline model runs, radiative flux budgets, poleward heat transport diagnostics, Assignment 1 sensitivity experiments, and the Assignment 2 SICI hysteresis sweep. |
| **`ebm_model.py`** | **Core model module.** Implements the `EnergyBalanceModel` class: staggered grid in $x=\sin\theta$, finite-volume discretization of diffusive transport $Q_T = \frac{d}{dx}[(1-x^2)D\frac{dT}{dx}]$, linearized $Q_{lw}=A+BT$, Legendre-approximated insolation $Q_{sw}$, matrix assembly $(\mathbf{K}\mathbf{T} = \mathbf{A} - \mathbf{Q}_{sw})$, and the exact 2-term analytical Legendre solution $T(x) = T_0 + T_2 P_2(x)$. |
| **`run_ebm.py`** | **Model runner and Assignment 1 experiments.** Runs the model against symmetric NCEP reanalysis observations (`ncep_temp.txt`), verifies numerical accuracy against the analytical solution, conducts sensitivity tests across diffusivity $D$, tests the Fanning & Weaver latitude-dependent albedo $\alpha(x) = (0.315 + 0.3615x^2)/f$, and evaluates the Lindzen & Farrell (1977) Hadley-cell enhanced diffusivity $D(x)$. Saves plots to `figures/ebm_assignment1_results.png`. |
| **`ncep_temp.txt`** | **Observational data.** Monthly long-term mean air temperatures from the NCEP Reanalysis used for model validation. |
| **`Matlab_code/`** | **MATLAB / GNU Octave code.** Standalone implementations (`ebm.m`, `ebmAlbedo.m`, `IterTalb.m`, `SICI.m`) with their own `README.md`. |
| **`R_code/`** | **R code.** Standalone implementations (`ebm.R`, `analyT.R`, `ebmAlbedo.R`, `iterate_ebmAlbedo.R`, `SICI.R`, `test_ebm_ncep.R`, and local `ncep_temp.txt`) with their own `README.md`. |
| **`ebm_sici.py`** | **Small Ice Cap Instability (SICI) and hysteresis experiment.** Implements temperature-dependent albedo feedback ($\alpha = 0.60$ for $T < -10^\circ\mathrm{C}$, else $0.33$), relaxation iterations to equilibrium, and sweeps $S_0$ through warm, cooling (snowball collapse), and warming (deglaciation) branches. Calculates the deglaciation radiative forcing $\Delta Q_{lw}$ and required $\mathrm{CO}_2$ multiplier $K$. Saves plots to `figures/ebm_sici_bifurcation.png`. |
| **`Docs/EBM_eng3.pdf`** | **Theoretical reference document.** Complete course notes describing the 1-D EBM derivation, Legendre polynomial analytical solution, cell-averaged finite-volume scheme, Small Ice Cap Instability, and student assignments. |

---

## 3. Mathematical Formulation

The model solves the steady-state thermal balance across latitude $\theta$ (using $x = \sin\theta \in [0, 1]$):

$$Q_{sw}(x) + Q_T(x) - Q_{lw}(x) = 0$$

1. **Shortwave radiation ($Q_{sw}$):**
   $$Q_{sw}(x) = \frac{S_0}{4} \left(1 - 0.492 P_2(x)\right) (1 - \alpha(x)), \qquad P_2(x) = \frac{3x^2 - 1}{2}$$

2. **Longwave radiation ($Q_{lw}$, linearized gray-body):**
   $$Q_{lw}(x) = A + B T(x), \qquad A = 204\ \mathrm{W/m^2}, \quad B = 1.9\ \mathrm{W/(m^2\,^\circ C)}$$

3. **Diffusive poleward heat transport ($Q_T$):**
   $$Q_T(x) = \frac{d}{dx} \left[ D(x) (1 - x^2) \frac{dT}{dx} \right], \qquad F(0) = 0, \quad F(1) = 0$$

4. **Exact analytical solution (EBM_eng3.pdf, Section 2):**
   $$T(x) = T_0 + T_2 P_2(x) = \frac{a Q - A}{B} + \frac{\beta a Q}{6D + B} P_2(x)$$
   where $a = 1 - \alpha = 0.67$, $Q = S_0/4 = 342.5\ \mathrm{W/m^2}$, $\beta = -0.492$. The numerical and analytical solutions agree within $1.1 \times 10^{-4}\ ^\circ\mathrm{C}$.

5. **Cell-averaged discretization (staggered grid, Section 3.1):**
   $$\mathbf{K} \mathbf{T} = \mathbf{A} - \mathbf{Q}_{sw}, \qquad \mathbf{K} = \mathbf{Q}_F \mathbf{D}_{mat} \mathbf{F}_T - \mathbf{B}$$
   where $\mathbf{F}_T$ computes temperature gradients on internal cell edges, $\mathbf{D}_{mat}$ holds $(1 - x_e^2) D(x_e)$, $\mathbf{Q}_F$ calculates flux divergence, and $\mathbf{B} = B \mathbf{I}$.

---

## 4. How to Run

### Interactive Jupyter Notebook:

```bash
cd /home/halldor/Work/ELJ/MyEBM/EBM_1D
jupyter notebook EBM_1D_Analysis.ipynb
```

### Standalone Python Scripts:

From this directory:

```bash
cd /home/halldor/Work/ELJ/MyEBM/EBM_1D

# 1. Run baseline model, analytical verification, and Assignment 1 experiments
python3 run_ebm.py

# 2. Run the SICI hysteresis simulation (Assignment 2)
python3 ebm_sici.py
```

Output figures are automatically saved to `figures/`:
- `figures/ebm_assignment1_results.png`
- `figures/ebm_sici_bifurcation.png`

Or import in Python / Jupyter:

```python
from ebm_model import EnergyBalanceModel

# Initialize and solve
ebm = EnergyBalanceModel(N=200, Diff=0.42)
sol = ebm.solve()

# Access results
lat = sol["lat"]          # Latitude in degrees
temp = sol["Temperature"] # Temperature in deg C
qsw = sol["Qsw"]          # Absorbed solar flux (W/m^2)
qlw = sol["Qlw"]          # Outgoing thermal flux (W/m^2)
```
