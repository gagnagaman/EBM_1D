"""Run and explore the 1-D Energy Balance Model (Assignment 1 from EBM_eng3.pdf).

This script performs:
  1. Default run vs. NCEP observations (from ncep_temp.txt)
  2. Verification against the exact analytical Legendre solution
  3. Sensitivity to the diffusivity parameter D (D in {0.1, 0.42, 0.5, 1.0, 10.0})
  4. Latitude-dependent albedo: alpha(x) = (0.315 + 0.3615*x^2) / f
  5. Lindzen & Farrell (1977) variable diffusivity: D(x) = 0.45 * [1 + 9*exp(-(x/sin 30deg)^6)]
  6. Saves summary plots to 'figures/ebm_assignment1_results.png'
"""

from pathlib import Path
from typing import Optional, Tuple
import matplotlib.pyplot as plt
import numpy as np

from ebm_model import EnergyBalanceModel


def load_ncep_data(filepath: Optional[Path] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load NCEP longterm mean temperature data and compute symmetric profile."""
    if filepath is None:
        current_dir = Path(__file__).resolve().parent
        local_path = current_dir / "ncep_temp.txt"
        if local_path.exists():
            filepath = local_path
        else:
            filepath = current_dir.parent / "R_EBM" / "ncep_temp.txt"

    # Header is 6 lines (5 comment lines + 1 column header)
    data = np.genfromtxt(filepath, skip_header=6)
    lat = data[:, 0]
    temp = data[:, 2]

    # Symmetrize between Northern and Southern hemispheres: (T(lat) + T(-lat)) / 2
    temp_sym = 0.5 * (temp + temp[::-1])
    # Restrict to Northern hemisphere (lat >= 0)
    mask = lat >= 0.0
    lat_nh = lat[mask]
    temp_nh = temp_sym[mask]
    sinlat_nh = np.sin(np.deg2rad(lat_nh))

    # Sort from equator to pole
    order = np.argsort(lat_nh)
    return lat_nh[order], sinlat_nh[order], temp_nh[order]


def run_assignment1() -> None:
    """Execute all experiments for Assignment 1."""
    print("=" * 65)
    print("1-D Diffusive Energy Balance Model - Assignment 1 Experiments")
    print("=" * 65)

    model = EnergyBalanceModel(N=200, A=204.0, B=1.9, Diff=0.42, So=1370.0, alpha=0.33)
    sol_default = model.solve()
    sol_analytic = model.analytical_solution(slat=sol_default["sinlat"])

    max_diff = np.max(np.abs(sol_default["Temperature"] - sol_analytic["Temperature"]))
    print(f"Max absolute difference vs. exact analytical solution: {max_diff:.3e} °C")
    print(f"Global mean temperature (numerical):  {np.mean(sol_default['Temperature']):.2f} °C")
    print(f"Global mean temperature (analytical): {sol_analytic['T_mean']:.2f} °C")

    # Load NCEP data
    lat_ncep, sinlat_ncep, temp_ncep = load_ncep_data()

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # -------------------------------------------------------------
    # Panel 1: Baseline run vs NCEP observations
    # -------------------------------------------------------------
    ax = axes[0, 0]
    # Note: The numerical and exact analytical solutions match to within ~1e-4 °C.
    # We plot the analytical solution as a distinct dashed line and the numerical
    # solution with markers/staggered style so both curves are clearly visible.
    ax.plot(sol_analytic["lat"], sol_analytic["Temperature"], color="#d95f02", ls="-", lw=3.0, alpha=0.9, label="Analytical (P2 Legendre)")
    ax.plot(sol_default["lat"], sol_default["Temperature"], color="#1f78b4", ls="--", lw=2.0, dashes=(4, 3), label="Numerical EBM (D=0.42, alpha=0.33)")
    ax.plot(lat_ncep, temp_ncep, "ko", ms=4, alpha=0.75, label="NCEP Reanalysis (symmetrized)")
    ax.set_xlabel("Latitude (°N)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("1. Baseline Model vs. Analytical & NCEP")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=8.5)

    # -------------------------------------------------------------
    # Panel 2: Varying diffusivity Diff
    # -------------------------------------------------------------
    ax = axes[0, 1]
    diff_values = [0.0, 0.1, 0.42, 0.5, 1.0, 10.0]
    colors = plt.cm.viridis(np.linspace(0, 1, len(diff_values)))
    for d_val, c in zip(diff_values, colors):
        sol_d = model.solve(diffusion=d_val)
        label = f"D = {d_val}" if d_val > 0 else "D = 0 (no transport)"
        ax.plot(sol_d["lat"], sol_d["Temperature"], label=label, color=c, lw=1.7)
    ax.plot(lat_ncep, temp_ncep, "k.", ms=3, label="NCEP")
    ax.set_xlabel("Latitude (°N)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("2. Sensitivity to Diffusivity D")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=8)

    # -------------------------------------------------------------
    # Panel 3: Latitude-dependent albedo alpha(x) = (0.315 + 0.3615*x^2) / f
    # -------------------------------------------------------------
    ax = axes[1, 0]
    f_values = [1.0, 1.1, 1.2, 1.3]
    for f in f_values:
        def alb_func(x, factor=f):
            return (0.315 + 0.3615 * x**2) / factor

        sol_alb = model.solve(albedo=alb_func)
        ax.plot(sol_alb["lat"], sol_alb["Temperature"], lw=1.7, label=f"f = {f}")
    ax.plot(lat_ncep, temp_ncep, "ro", ms=4, alpha=0.8, label="NCEP")
    ax.set_xlabel("Latitude (°N)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title(r"3. Variable Albedo: $\alpha(x) = (0.315 + 0.3615 x^2) / f$")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=9)

    # -------------------------------------------------------------
    # Panel 4: Lindzen & Farrell (1977) Variable Diffusivity
    # -------------------------------------------------------------
    ax = axes[1, 1]
    def d_lindzen(x_inner):
        sin30 = np.sin(np.deg2rad(30.0))
        return 0.45 * (1.0 + 9.0 * np.exp(-((x_inner / sin30) ** 6)))

    sol_lf_const_alb = model.solve(diffusion=d_lindzen)
    # Combined: Lindzen-Farrell D + Fanning-Weaver albedo with tuned f
    def alb_fw(x):
        return (0.315 + 0.3615 * x**2) / 1.15
    sol_lf_tuned = model.solve(diffusion=d_lindzen, albedo=alb_fw)

    ax.plot(sol_default["lat"], sol_default["Temperature"], "k--", lw=1.5, label="Constant D=0.42")
    ax.plot(sol_lf_const_alb["lat"], sol_lf_const_alb["Temperature"], "b-", lw=1.8, label="Lindzen-Farrell D(x)")
    ax.plot(sol_lf_tuned["lat"], sol_lf_tuned["Temperature"], "g-", lw=2, label="Lindzen-Farrell D(x) + alb(x)")
    ax.plot(lat_ncep, temp_ncep, "ro", ms=4, alpha=0.8, label="NCEP")
    ax.set_xlabel("Latitude (°N)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("4. Lindzen & Farrell (1977) Variable D(x)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=8.5)

    plt.tight_layout()
    fig_dir = Path(__file__).resolve().parent / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    out_png = fig_dir / "ebm_assignment1_results.png"
    plt.savefig(out_png, dpi=180)
    plt.close()
    print(f"Saved figure to: {out_png}")
    print("Assignment 1 runs completed successfully.")


if __name__ == "__main__":
    run_assignment1()
