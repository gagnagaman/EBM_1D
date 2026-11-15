"""Small Ice Cap Instability (SICI) and Temperature-Albedo Feedback Experiment.

Implements the SICI experiment described in:
  - EBM_eng3.pdf, Appendix A (Small Ice Cap Instability) & Assignment 2
  - R_EBM/ebmAlbedo.R, iterate_ebmAlbedo.R, and SICI.R

Mechanism:
  The planetary albedo depends on local temperature:
      alpha(x) = 0.60   if T(x) < -10 °C  (ice/snow cover)
      alpha(x) = 0.33   otherwise         (ice-free)

Iterative Relaxation:
  Given an initial temperature profile T^(0) and solar constant S_0,
  we update alpha^(k) = alpha(T^(k)), solve for T^(k+1), and repeat until
      ||T^(k+1) - T^(k)||_1 < 1e-5.

Hysteresis Sweep:
  Branch 1: Warm branch swept upwards (S_0 from 1370 to 1370 * 1.30)
  Branch 2: Warm branch swept downwards to trigger collapse (S_0 from 1370 down to 1370 * 0.85)
  Branch 3: Snowball branch warmed upwards to trigger deglaciation (S_0 from 1150 up to 1370 * 1.30)

Outputs:
  - Bifurcation diagram: Global mean temperature vs. Solar Constant S_0
  - Transition values of S_0 for snowball collapse and deglaciation
  - Diagnostic estimates of delta Q_lw and required CO_2 multiplier for deglaciation
  - Saves plot to 'ebm_sici_bifurcation.png'
"""

from pathlib import Path
from typing import Dict, Optional, Tuple
import matplotlib.pyplot as plt
import numpy as np

from ebm_model import EnergyBalanceModel


class SICIExperiment:
    """Solves temperature-albedo feedback iterations and traces hysteresis curves."""

    def __init__(
        self,
        N: int = 200,
        A: float = 204.0,
        B: float = 1.9,
        Diff: float = 0.42,
        alpha_warm: float = 0.33,
        alpha_ice: float = 0.60,
        T_crit: float = -10.0,
    ) -> None:
        self.model = EnergyBalanceModel(N=N, A=A, B=B, Diff=Diff)
        self.alpha_warm = alpha_warm
        self.alpha_ice = alpha_ice
        self.T_crit = T_crit

    def ebm_albedo_step(self, Temp: np.ndarray, So: float) -> Dict[str, np.ndarray]:
        """One step: update albedo based on temperature, then solve EBM."""
        alb = np.full(self.model.N - 1, self.alpha_warm)
        alb[Temp < self.T_crit] = self.alpha_ice
        sol = self.model.solve(albedo=alb, So=So)
        return sol

    def iterate_equilibrium(
        self,
        T_init: np.ndarray,
        So: float,
        tol: float = 1e-5,
        max_iter: int = 60,
        verbose: bool = False,
    ) -> Dict[str, np.ndarray]:
        """Iterate ebm_albedo_step until convergence to an equilibrium state."""
        T_current = np.copy(T_init)
        sol = self.ebm_albedo_step(T_current, So)

        for counter in range(1, max_iter + 1):
            sol = self.ebm_albedo_step(T_current, So)
            T_new = sol["Temperature"]
            diff = np.sum(np.abs(T_new - T_current))
            T_current = T_new

            if verbose:
                mean_alb = np.mean(sol["albedo"])
                print(
                    f"So={So:.1f} | Iter {counter:2d} | diff={diff:.2e} | mean alb={mean_alb:.3f}"
                )

            if diff < tol:
                break

        return sol

    def run_sweep(
        self,
        n_up: int = 40,
        n_down: int = 71,
        n_warmup: int = 250,
    ) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """Trace the three hysteresis branches matching SICI.R."""
        So_start = 1370.0

        # Branch 1: Warm branch going up
        print("Tracing Branch 1: Warm branch upwards...")
        So_up = np.linspace(So_start, So_start * 1.30, n_up)
        T_mean_up = np.zeros(n_up)
        sol_init = self.model.solve(albedo=self.alpha_warm, So=So_start)
        T_state = sol_init["Temperature"]

        for i, so_val in enumerate(So_up):
            sol = self.iterate_equilibrium(T_state, so_val)
            T_state = sol["Temperature"]
            T_mean_up[i] = np.mean(T_state)

        # Branch 2: Warm branch going down to snowball collapse
        print("Tracing Branch 2: Warm branch downwards to snowball collapse...")
        So_down = np.linspace(So_start, So_start * (1.0 - 0.15), n_down)
        T_mean_down = np.zeros(n_down)
        T_state = sol_init["Temperature"]

        collapse_So = None
        for i, so_val in enumerate(So_down):
            sol = self.iterate_equilibrium(T_state, so_val)
            T_state = sol["Temperature"]
            T_mean_down[i] = np.mean(T_state)
            mean_alb = np.mean(sol["albedo"])
            if collapse_So is None and np.isclose(mean_alb, self.alpha_ice, atol=1e-3):
                collapse_So = so_val

        # Branch 3: Snowball branch warmed up to deglaciation
        print("Tracing Branch 3: Snowball branch warmed up to deglaciation...")
        So_cold = 1150.0
        # Initialize cold snowball state
        sol_cold_init = self.iterate_equilibrium(
            np.full(self.model.N - 1, -50.0), So_cold
        )
        T_state = sol_cold_init["Temperature"]

        So_up2 = np.linspace(So_cold, So_start * 1.30, n_warmup)
        T_mean_up2 = np.zeros(n_warmup)
        deglaciation_So = None

        for i, so_val in enumerate(So_up2):
            sol = self.iterate_equilibrium(T_state, so_val)
            T_state = sol["Temperature"]
            T_mean_up2[i] = np.mean(T_state)
            mean_alb = np.mean(sol["albedo"])
            if deglaciation_So is None and mean_alb < 0.55:
                deglaciation_So = so_val

        print("-" * 65)
        print("SICI Bifurcation Summary (Assignment 2 in EBM_eng3.pdf):")
        if collapse_So is not None:
            print(f"  Snowball collapse threshold:   S_0 = {collapse_So:.2f} W/m^2")
        if deglaciation_So is not None:
            print(f"  Deglaciation threshold:        S_0 = {deglaciation_So:.2f} W/m^2")
        if collapse_So is not None and deglaciation_So is not None:
            delta_So = deglaciation_So - collapse_So
            delta_Qlw = (delta_So / 4.0) * (1.0 - self.alpha_ice)
            # CO2 forcing: delta Q = 5.35 * ln(K) -> K = exp(delta Q / 5.35)
            K_mult = np.exp(delta_Qlw / 5.35)
            pCO2_modern = 0.0004  # bars
            pCO2_required = pCO2_modern * K_mult

            print(f"  Delta S_0 (hysteresis width):  {delta_So:.2f} W/m^2")
            print(f"  Delta Q_lw (radiative forcing): {delta_Qlw:.2f} W/m^2")
            print(f"  Required CO_2 multiplier K:    {K_mult:.1f}x modern")
            print(f"  Estimated pCO_2 for deglaciation: {pCO2_required:.4f} bars")
            print("  (Compare with Pierrehumbert / Higgins & Schrag: ~0.1 - 0.4 bars)")
        print("-" * 65)

        return {
            "branch_up": (So_up, T_mean_up),
            "branch_down": (So_down, T_mean_down),
            "branch_warmup": (So_up2, T_mean_up2),
        }

    def plot_bifurcation(
        self,
        results: Dict[str, Tuple[np.ndarray, np.ndarray]],
        out_png: Optional[Path] = None,
    ) -> None:
        """Plot the hysteresis diagram matching Figure 4 in EBM_eng3.pdf."""
        So_up, T_up = results["branch_up"]
        So_down, T_down = results["branch_down"]
        So_up2, T_up2 = results["branch_warmup"]

        fig, ax = plt.subplots(figsize=(9, 6.5))

        ax.plot(So_up2, T_up2, "k-", lw=1.8, label="Snowball branch (warming up)")
        ax.plot(So_up, T_up, "r-", lw=2.0, label="Warm branch (warming up)")
        ax.plot(So_down, T_down, "g-", lw=2.0, label="Warm branch (cooling down)")

        # Highlight modern solar constant
        ax.axvline(1370.0, color="gray", linestyle=":", lw=1.2, label=r"Modern $S_0 = 1370\ \mathrm{W/m^2}$")
        ax.axhline(self.T_crit, color="blue", linestyle="--", lw=1.0, alpha=0.6, label=r"$T_c = -10\ ^\circ\mathrm{C}$")

        ax.set_xlabel(r"Solar Constant $S_0\ (\mathrm{W/m^2})$", fontsize=11)
        ax.set_ylabel(r"Global Mean Temperature $(^\circ\mathrm{C})$", fontsize=11)
        ax.set_title("Small Ice Cap Instability (SICI) Hysteresis Loop\n(1-D EBM with Temperature-Albedo Feedback)", fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=9.5)

        plt.tight_layout()
        if out_png is None:
            fig_dir = Path(__file__).resolve().parent / "figures"
            fig_dir.mkdir(parents=True, exist_ok=True)
            out_png = fig_dir / "ebm_sici_bifurcation.png"
        plt.savefig(out_png, dpi=180)
        plt.close()
        print(f"Saved bifurcation diagram to: {out_png}")


def main() -> None:
    exp = SICIExperiment(N=200, A=204.0, B=1.9, Diff=0.42)
    results = exp.run_sweep()
    exp.plot_bifurcation(results)


if __name__ == "__main__":
    main()
