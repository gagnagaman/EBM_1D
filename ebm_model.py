"""1-D Energy Balance Model (EBM) with Diffusive Heat Transport.

Implements the 1-D steady-state zonally averaged energy balance model derived
in EBM_eng3.pdf (North 1975, JAS 32), formulated on cell averages with a
staggered finite-volume discretization.

The governing equation is:
    Q_sw(x) + Q_T(x) - Q_lw(x) = 0
where:
    x = sin(latitude), x in [0, 1] (symmetric northern hemisphere)
    Q_sw = (S_0 / 4) * (1 - 0.492 * P_2(x)) * (1 - alpha(x))
    Q_lw = A + B * T(x)
    Q_T  = d/dx [ D(x) * (1 - x^2) * dT/dx ]

Boundary conditions:
    F(0) = 0 (symmetry at equator)
    F(1) = 0 (no flux across pole)
"""

from typing import Callable, Dict, Optional, Union
import numpy as np


def legendre2(x: np.ndarray) -> np.ndarray:
    """Second Legendre polynomial: P_2(x) = (3*x^2 - 1) / 2."""
    return (3.0 * x**2 - 1.0) / 2.0


class EnergyBalanceModel:
    """1-D Diffusive Energy Balance Model."""

    def __init__(
        self,
        N: int = 200,
        A: float = 204.0,
        B: float = 1.9,
        Diff: float = 0.42,
        So: float = 1370.0,
        alpha: float = 0.33,
    ) -> None:
        """Initialize standard parameters.

        Args:
            N: Number of cell edges (N-1 grid cells).
            A: Outgoing longwave radiation intercept (W/m^2).
            B: Outgoing longwave radiation slope (W/(m^2 degC)).
            Diff: Thermal diffusivity constant (dimensionless / scaled).
            So: Solar constant (W/m^2).
            alpha: Baseline planetary albedo (fraction).
        """
        self.N = N
        self.A = A
        self.B = B
        self.Diff = Diff
        self.So = So
        self.alpha = alpha

    def setup_grid(self) -> Dict[str, np.ndarray]:
        """Create staggered grid in x = sin(latitude).

        Returns:
            Dictionary with:
                xe: cell edges of length N (x in [0, 1])
                xm: cell midpoints of length N-1
                dx: cell width
                lat: latitudes at midpoints in degrees
        """
        xe = np.linspace(0.0, 1.0, self.N)
        xm = xe[:-1] + np.diff(xe) / 2.0
        dx = float(xe[1] - xe[0])
        lat = np.rad2deg(np.arcsin(xm))
        return {"xe": xe, "xm": xm, "dx": dx, "lat": lat}

    def solve(
        self,
        albedo: Optional[Union[float, np.ndarray, Callable[[np.ndarray], np.ndarray]]] = None,
        diffusion: Optional[Union[float, np.ndarray, Callable[[np.ndarray], np.ndarray]]] = None,
        So: Optional[float] = None,
        A: Optional[float] = None,
        B: Optional[float] = None,
    ) -> Dict[str, np.ndarray]:
        """Solve for steady-state temperature profile.

        Args:
            albedo: Scalar, vector of length N-1, or function alpha(xm).
            diffusion: Scalar, vector of length N-2 (at inner edges xe[1:-1]),
                       or function D(xe_inner).
            So: Solar constant override (W/m^2).
            A: OLR intercept override (W/m^2).
            B: OLR slope override (W/(m^2 degC)).

        Returns:
            Dictionary with:
                sinlat: midpoints xm (length N-1)
                lat: latitudes in degrees
                Temperature: steady-state T(x) in deg C
                Qsw: absorbed shortwave radiation (W/m^2)
                Qlw: outgoing longwave radiation (W/m^2)
                albedo: planetary albedo profile
                heat_transport: poleward flux F at cell edges (length N)
        """
        grid = self.setup_grid()
        xe, xm, dx = grid["xe"], grid["xm"], grid["dx"]
        N = self.N

        # Resolve parameters
        So_val = self.So if So is None else float(So)
        A_val = self.A if A is None else float(A)
        B_val = self.B if B is None else float(B)

        # 1. Insolation profile S(x)
        p = legendre2(xm)
        S = (So_val / 4.0) * (1.0 - 0.492 * p)

        # 2. Albedo profile
        if albedo is None:
            alb = np.full(N - 1, self.alpha)
        elif callable(albedo):
            alb = np.asarray(albedo(xm), dtype=float)
        elif np.isscalar(albedo):
            alb = np.full(N - 1, float(albedo))
        else:
            alb = np.asarray(albedo, dtype=float)

        Qsw = (1.0 - alb) * S

        # 3. Diffusivity profile D(x) at internal cell edges xe[1:-1]
        xe_inner = xe[1:-1]
        if diffusion is None:
            D = np.full(N - 2, self.Diff)
        elif callable(diffusion):
            D = np.asarray(diffusion(xe_inner), dtype=float)
        elif np.isscalar(diffusion):
            D = np.full(N - 2, float(diffusion))
        else:
            D = np.asarray(diffusion, dtype=float)

        # 4. Matrix operator setup (matches Section 3.1 in EBM_eng3.pdf)
        # Dmat = diag((1 - x_i^2) * D(x_i)) for i = 2 .. N-1 (size: (N-2) x (N-2))
        diagline = D * (1.0 - xe_inner**2)
        Dmat = np.diag(diagline)

        # FT = gradient matrix dT/dx at internal cell edges (size: (N-2) x (N-1))
        # row i: [-1, 1] / dx
        FT = np.zeros((N - 2, N - 1))
        for i in range(N - 2):
            FT[i, i] = -1.0
            FT[i, i + 1] = 1.0
        FT /= dx

        # QF = divergence matrix dF/dx from edges to cell centers (size: (N-1) x (N-2))
        # incorporates boundary conditions F(0) = 0 and F(1) = 0
        QF = np.zeros((N - 1, N - 2))
        for i in range(N - 2):
            QF[i, i] = 1.0
            QF[i + 1, i] = -1.0
        QF /= dx

        # Bmat = B * I (size: (N-1) x (N-1))
        Bmat = np.eye(N - 1) * B_val

        # Full system matrix: K = QF @ Dmat @ FT - Bmat
        # Governing equation: K @ T = A - Qsw
        K = QF @ Dmat @ FT - Bmat
        rhs = A_val - Qsw

        Temperature = np.linalg.solve(K, rhs)
        Qlw = A_val + B_val * Temperature

        # Diagnosed poleward heat flux at edges: F_edges = [0, Dmat @ FT @ T, 0]
        F_inner = Dmat @ (FT @ Temperature)
        heat_transport = np.zeros(N)
        heat_transport[1:-1] = F_inner

        return {
            "sinlat": xm,
            "lat": grid["lat"],
            "Temperature": Temperature,
            "Qsw": Qsw,
            "Qlw": Qlw,
            "albedo": alb,
            "heat_transport": heat_transport,
        }

    @staticmethod
    def analytical_solution(
        slat: Optional[np.ndarray] = None,
        A: float = 204.0,
        B: float = 1.9,
        Diff: float = 0.42,
        So: float = 1370.0,
        alpha: float = 0.33,
    ) -> Dict[str, np.ndarray]:
        """Exact 2-term Legendre polynomial solution (EBM_eng3.pdf Section 2).

        T(x) = (a*Q - A)/B + [beta * a * Q / (6*D + B)] * P_2(x)
        """
        if slat is None:
            slat = np.linspace(0.0, 1.0, 1000)
        Q = So / 4.0
        a = 1.0 - alpha
        beta = -0.492
        p = legendre2(slat)

        T_0 = (a * Q - A) / B
        T_2 = (beta * a * Q) / (6.0 * Diff + B)
        Temperature = T_0 + T_2 * p
        Qsw = Q * a * (1.0 + beta * p)
        Qlw = A + B * Temperature

        return {
            "sinlat": slat,
            "lat": np.rad2deg(np.arcsin(slat)),
            "Temperature": Temperature,
            "Qsw": Qsw,
            "Qlw": Qlw,
            "T_mean": T_0,
            "T_contrast": T_2,
        }
