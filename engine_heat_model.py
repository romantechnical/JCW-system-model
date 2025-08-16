from dataclasses import dataclass, asdict
from typing import Optional, Tuple

# Physical constants
CP_WATER = 4180.0  # J/(kg*K)

@dataclass
class EngineThermalParams:
    """
    Parameters for the two-capacitance engine thermal model (metal ↔ jacket water).
    All units are SI.
    """
    C_m: float = 5.24e6      # J/K   - thermal capacitance of metal (effective)
    C_w: float = 2.65e6      # J/K   - thermal capacitance of jacket water volume
    UA_mw: float = 80e3      # W/K   - effective metal↔water heat transfer
    cp_w: float = CP_WATER   # J/(kg*K) - water specific heat (can be tuned slightly)
    # Defaults derived from your data:
    #   m_dot ≈ 45.83 kg/s (from 165 m^3/h), Qgen75 ≈ 2.30 MW, Tin≈72°C, Tout≈84°C at 75%

class EngineHeatModel:
    """
    Two-capacitance thermal model of MAN B&W 6L70MC-C8 engine jacket cooling.
    States:
      T_m: effective metal temperature [°C]
      T_w: jacket water bulk temperature (well-mixed in block) [°C] ~ outlet temperature
    Dynamics:
      C_m * dT_m/dt = Q_gen(t) - UA*(T_m - T_w)
      C_w * dT_w/dt = m_dot*cp*(T_in - T_w) + UA*(T_m - T_w)
    Use step(T_in, Qgen, m_dot, dt) to advance the model.
    """
    def __init__(self, params: Optional[EngineThermalParams] = None,
                 Tm0: float = 110.0, Tw0: float = 84.0):
        self.p = params or EngineThermalParams()
        self.T_m = float(Tm0)
        self.T_w = float(Tw0)

    def reset(self, Tm0: float, Tw0: float) -> None:
        """Reset internal temperatures (°C)."""
        self.T_m = float(Tm0)
        self.T_w = float(Tw0)

    def get_state(self) -> Tuple[float, float]:
        """Return (T_m, T_w) in °C."""
        return self.T_m, self.T_w

    def set_params(self, **kwargs) -> None:
        """Update parameters (C_m, C_w, UA_mw, cp_w)."""
        for k, v in kwargs.items():
            if hasattr(self.p, k):
                setattr(self.p, k, float(v))

    # ---- core integration (Heun's method / RK2 for stability) ----
    def _rhs(self, Tm: float, Tw: float, T_in: float, Qgen: float, m_dot: float) -> Tuple[float, float]:
        p = self.p
        # Heat flows: metal->water positive when Tm>Tw
        qm_w = p.UA_mw * (Tm - Tw)              # W
        q_in_mix = m_dot * p.cp_w * (T_in - Tw)  # W (convective exchange with inlet flow)
        dTm = (Qgen - qm_w) / max(1e-9, p.C_m)   # K/s
        dTw = (q_in_mix + qm_w) / max(1e-9, p.C_w)
        return dTm, dTw

    def step(self, T_in: float, Qgen: float, m_dot: float, dt: float) -> float:
        """
        Advance model by dt seconds.
        Inputs:
          T_in  [°C] : water temperature at engine inlet (from cooler)
          Qgen  [W]  : heat flow from engine to jacket (use ~2.30e6 W at 75% as nominal)
          m_dot [kg/s]: mass flow rate of water through jacket
          dt    [s]  : time step
        Returns:
          T_out [°C]  : outlet water temperature (≈ T_w)
        """
        dt = float(dt)
        # RK2 (Heun)
        k1_m, k1_w = self._rhs(self.T_m, self.T_w, T_in, Qgen, m_dot)
        Tm_e = self.T_m + dt * k1_m
        Tw_e = self.T_w + dt * k1_w
        k2_m, k2_w = self._rhs(Tm_e, Tw_e, T_in, Qgen, m_dot)

        self.T_m += dt * 0.5 * (k1_m + k2_m)
        self.T_w += dt * 0.5 * (k1_w + k2_w)

        return self.T_w  # outlet ≈ bulk water temperature in jacket

# Tiny demo when run directly (no plotting, just prints)
if __name__ == "__main__":
    # Defaults from your 75% condition
    m_dot = 45.83         # kg/s (165 m^3/h)
    Qgen = 2.30e6         # W (≈ m_dot*cp*12K)
    model = EngineHeatModel(Tm0=109.0, Tw0=84.0)
    Tin = 72.0
    dt = 0.5
    for k in range(4):
        Tout = model.step(Tin, Qgen, m_dot, dt)
        print(f"step={k} Tin={Tin:.1f}C -> Tout={Tout:.2f}C, Tm={model.T_m:.2f}C")
```0