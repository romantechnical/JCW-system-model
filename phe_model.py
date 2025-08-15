
import argparse
import numpy as np
import matplotlib.pyplot as plt

# --- Constants ---
RHO_WATER = 1000.0      # kg/m^3
CP_WATER  = 4180.0      # J/(kg*K)

def lmtd(dT1, dT2):
    # Log Mean Temperature Difference with safe handling
    if dT1 <= 0 or dT2 <= 0:
        return 0.0
    if abs(dT1 - dT2) < 1e-9:
        return dT1  # limit as dT1 -> dT2
    return (dT1 - dT2) / np.log(dT1 / dT2)

def ua_from_spec(U_W_m2K=6017.0, A_m2=19.44):
    # Returns UA in W/K
    return U_W_m2K * A_m2

def solve_hot_out(Th_in, Tc_in, m_hot_kg_s, m_cold_kg_s, UA_W_K, tol=1e-6, maxit=100):
    """
    Solve for hot outlet temperature for a counter-current plate heat exchanger at steady state:
    Q = m_h*cp*(Th_in - Th_out) = UA * LMTD(Th_in, Th_out, Tc_in, Tc_out)
    and Tc_out = Tc_in + Q/(m_c*cp)
    """
    cp = CP_WATER

    # Bisection bounds for Th_out: it must be between Tc_in+eps and Th_in-eps
    lo = Tc_in + 1e-6
    hi = Th_in - 1e-6
    lo = min(lo, hi - 1e-6)  # safety

    def f(Th_out):
        Qh = m_hot_kg_s * cp * (Th_in - Th_out)
        Tc_out = Tc_in + Qh / (m_cold_kg_s * cp)
        dT1 = Th_in - Tc_out
        dT2 = Th_out - Tc_in
        L = lmtd(dT1, dT2)
        Q_UA = UA_W_K * L
        return Qh - Q_UA  # root when energy balance equals UA*LMTD

    # If UA is huge/small, f may be monotonic; use bisection
    flo = f(lo)
    fhi = f(hi)
    # Ensure we have a sign change; if not, clamp result to bounds
    if flo == 0:
        return lo, Tc_in + m_hot_kg_s*cp*(Th_in-lo)/(m_cold_kg_s*cp)
    if fhi == 0:
        return hi, Tc_in + m_hot_kg_s*cp*(Th_in-hi)/(m_cold_kg_s*cp)
    # If no sign change, try to expand a bit
    for _ in range(5):
        if flo * fhi < 0:
            break
        # nudge bounds
        lo = (lo + hi) / 2.0 * 0.99
        hi = (lo + hi) / 2.0 * 1.01
        flo = f(lo)
        fhi = f(hi)
    if flo * fhi > 0:
        # fallback: return a physically plausible estimate (effectively pinch)
        Th_out = (Th_in + Tc_in) / 2.0
        Qh = m_hot_kg_s*cp*(Th_in - Th_out)
        Tc_out = Tc_in + Qh/(m_cold_kg_s*cp)
        return Th_out, Tc_out

    for _ in range(maxit):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if abs(fmid) < tol:
            Th_out = mid
            Qh = m_hot_kg_s * cp * (Th_in - Th_out)
            Tc_out = Tc_in + Qh / (m_cold_kg_s * cp)
            return Th_out, Tc_out
        if flo * fmid < 0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    # Max iterations: return mid
    Th_out = 0.5 * (lo + hi)
    Qh = m_hot_kg_s * cp * (Th_in - Th_out)
    Tc_out = Tc_in + Qh / (m_cold_kg_s * cp)
    return Th_out, Tc_out

def simulate_phe(duration1=200.0, duration2=200.0, dt=0.5,
                 Th_min1=80.0, Th_max1=95.0,
                 Th_min2=70.0, Th_max2=85.0,
                 Tc_in=32.0,
                 U=6017.0, A=19.44,
                 hot_flow_m3_h=165.0, cold_flow_m3_h=214.0):
    """
    Generate Th_in(t) as requested, then compute Th_out(t) using a steady-state PHE model at each step.
    """
    UA = ua_from_spec(U, A)
    # mass flows (kg/s)
    m_hot = hot_flow_m3_h * RHO_WATER / 3600.0
    m_cold = cold_flow_m3_h * RHO_WATER / 3600.0

    # time arrays
    n1 = int(duration1 / dt)
    n2 = int(duration2 / dt)
    t1 = np.arange(n1) * dt
    t2 = np.arange(n2) * dt + duration1

    # input waveforms
    Th_in1 = (Th_min1 + Th_max1)/2.0 + (Th_max1 - Th_min1)/2.0 * np.sin(2*np.pi * (1.0/duration1) * t1)
    Th_in2 = (Th_min2 + Th_max2)/2.0 + (Th_max2 - Th_min2)/2.0 * np.sin(2*np.pi * (1.0/duration2) * t2)
    t = np.concatenate([t1, t2])
    Th_in = np.concatenate([Th_in1, Th_in2])

    Th_out = np.zeros_like(Th_in)
    Tc_out = np.zeros_like(Th_in)

    for k in range(len(t)):
        Th_out[k], Tc_out[k] = solve_hot_out(Th_in[k], Tc_in, m_hot, m_cold, UA)

    return t, Th_in, Th_out, Tc_out

def plot_results(t, Th_in, Th_out, Tc_out, save_png=None):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(11,5))
    plt.plot(t, Th_in, label="Hot in (°C)")
    plt.plot(t, Th_out, label="Hot out (°C)")
    plt.plot(t, Tc_out, label="Cold out (°C)")
    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (°C)")
    plt.title("Plate Heat Exchanger (Counter-current) — Steady-state per step")
    plt.legend()
    plt.tight_layout()
    if save_png:
        plt.savefig(save_png, dpi=160)
    try:
        plt.show()
    except Exception:
        pass

def main():
    ap = argparse.ArgumentParser(description="Plate Heat Exchanger model with spec-based UA")
    ap.add_argument("--dt", type=float, default=0.5, help="Time step (s)")
    ap.add_argument("--save", type=str, default="/storage/emulated/0/Download/logs/phe_plot.png",
                    help="PNG path to save plot")
    ap.add_argument("--hot_flow", type=float, default=165.0, help="Hot side flow (m^3/h)")
    ap.add_argument("--cold_flow", type=float, default=214.0, help="Cold side flow (m^3/h)")
    args = ap.parse_args()

    t, Th_in, Th_out, Tc_out = simulate_phe(dt=args.dt,
                                            hot_flow_m3_h=args.hot_flow,
                                            cold_flow_m3_h=args.cold_flow)
    plot_results(t, Th_in, Th_out, Tc_out, save_png=args.save)
    print("Saved plot to:", args.save)

if __name__ == "__main__":
    main()
