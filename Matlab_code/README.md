# MATLAB / GNU Octave Implementation of 1-D Energy Balance Model

This directory contains MATLAB (and GNU Octave compatible) versions of the 1-D diffusive energy balance model (EBM) and the Small Ice Cap Instability (SICI) experiment.

---

## Files

| File | Purpose |
|------|---------|
| **`ebm.m`** | Core steady-state 1-D diffusive EBM solver on cell averages. Can be run directly as a script (`ebm`) to produce variables `T` and `xm` in the workspace, or uncomment the function header for functional calls `[T, xm] = ebm(N, A, B, Diff)`. |
| **`ebmAlbedo.m`** | Function `[T, xm] = ebmAlbedo(T, So)` that adjusts the albedo according to local temperature ($\alpha = 0.60$ if $T < -10^\circ\mathrm{C}$, else $0.33$) and solves for the updated temperature profile. |
| **`IterTalb.m`** | Function `Talb = IterTalb(T, So)` that iteratively relaxes `ebmAlbedo` until $\|\mathbf{T}^{(k+1)} - \mathbf{T}^{(k)}\|_1 < 10^{-5}$, finding the equilibrium climate for a given solar constant $S_0$. |
| **`SICI.m`** | Script that conducts the Small Ice Cap Instability (SICI) experiment (Assignment 2): sweeps $S_0$ across the warm, snowball collapse, and deglaciation branches and plots the resulting hysteresis loop. |

---

## How to Run

### In MATLAB or GNU Octave:

```matlab
% 1. Run baseline model
ebm

% Plot temperature vs. sine of latitude
plot(xm, T)
xlabel('sin(latitude)')
ylabel('Temperature (^\circ C)')
title('Baseline 1-D EBM Temperature Profile')

% 2. Run the SICI hysteresis simulation (Assignment 2)
SICI
```
