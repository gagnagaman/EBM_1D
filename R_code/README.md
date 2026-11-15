# R Implementation of 1-D Energy Balance Model

This directory contains R versions of the 1-D diffusive energy balance model (EBM) and the Small Ice Cap Instability (SICI) experiment.

---

## Files

| File | Purpose |
|------|---------|
| **`ebm.R`** | Core steady-state 1-D diffusive EBM solver on cell averages. Run as `out = ebm()` or with custom parameters `out = ebm(N=200, A=204, B=1.9, Diff=0.42)`. Returns list with `sinlat`, `lat`, `Temperature`, `Qsw`, and `Qlw`. |
| **`analyT.R`** | Exact analytical solution using Legendre polynomials: $T(x) = T_0 + T_2 P_2(x)$. |
| **`readncep.R`** | Reads the local `ncep_temp.txt` dataset and computes the hemispherically symmetrized observational temperature profile. |
| **`test_ebm_ncep.R`** | Demonstration script running `ebm.R` across different diffusivity ($Diff$), $A$, and $B$ values and comparing with NCEP reanalysis data (Assignment 1). |
| **`ebmAlbedo.R`** | Computes temperature-dependent albedo: $\alpha = 0.60$ where $T < -10^\circ\mathrm{C}$, else $0.33$, and solves for the resulting temperature profile given $S_0$. |
| **`iterate_ebmAlbedo.R`** | Iterates `ebmAlbedo.R` until convergence to an equilibrium climate state given a solar constant $S_0$. |
| **`SICI.R`** | Small Ice Cap Instability (SICI) simulation: sweeps $S_0$ across the warm, snowball collapse, and deglaciation branches to trace the hysteresis loop (Assignment 2). |
| **`ncep_temp.txt`** | Local copy of the NCEP reanalysis monthly longterm mean air temperature dataset. |

---

## Dependencies

The scripts make use of the `matlab` R package (for functions such as `linspace`, `zeros`, `ones`, `eye`, `diff`, and `flipud`).

To install in R:
```R
install.packages("matlab")
```

---

## How to Run

From this directory in R or RStudio:

```R
# 1. Run baseline model and compare with observations
source("test_ebm_ncep.R")

# 2. Run the SICI hysteresis simulation
source("SICI.R")
```
