import argparse
import csv
import math
import os
from typing import List, Tuple

RHO_WATER = 1000.0      # kg/m^3
CP_WATER  = 4180.0      # J/(kg*K)

def lmtd(dT1: float, dT2: float) -> float:
    if dT1 <= 0 or dT2 <= 0:
        return 0.0
    if abs(dT1 - dT2) < 1e-12:
        return dT1
    return (dT1 - dT2) / math.log(dT1 / dT2)

def ua_from_spec(U_W_m2K: float = 6017.0, A_m2: float = 19.44) -> float:
    return U_W_m2K * A_m2  # W/K

def solve_hot_out(Th_in: float, Tc_in: float, m_hot_kg_s: float, m_cold_kg_s: float,
                  UA_W_K: float, tol: float = 1e-6, maxit: int = 100) -> Tuple[float, float]:
    cp = CP_WATER
    lo = Tc_in + 1e-6
    hi = Th_in - 1e-6
    if lo >= hi:
        # physically impossible pinch; return tiny approach
        Th_out = (Th_in + Tc_in) / 2.0
        Qh = m_hot_kg_s * cp * (Th_in - Th_out)
        Tc_out = Tc_in + Qh / (m_cold_kg_s * cp)
        return Th_out, Tc_out

    def f(Th_out: float) -> float:
        Qh = m_hot_kg_s * cp * (Th_in - Th_out)
        Tc_out = Tc_in + Qh / (m_cold_kg_s * cp)
        dT1 = Th_in - Tc_out
        dT2 = Th_out - Tc_in
        L = lmtd(dT1, dT2)
        return Qh - UA_W_K * L

    flo = f(lo)
    fhi = f(hi)
    if flo == 0.0:
        Qh = m_hot_kg_s * cp * (Th_in - lo)
        Tc_out = Tc_in + Qh / (m_cold_kg_s * cp)
        return lo, Tc_out
    if fhi == 0.0:
        Qh = m_hot_kg_s * cp * (Th_in - hi)
        Tc_out = Tc_in + Qh / (m_cold_kg_s * cp)
        return hi, Tc_out

    # Ensure sign change
    if flo * fhi > 0:
        # fallback
        Th_out = (Th_in + Tc_in) / 2.0
        Qh = m_hot_kg_s * cp * (Th_in - Th_out)
        Tc_out = Tc_in + Qh / (m_cold_kg_s * cp)
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

    Th_out = 0.5 * (lo + hi)
    Qh = m_hot_kg_s * cp * (Th_in - Th_out)
    Tc_out = Tc_in + Qh / (m_cold_kg_s * cp)
    return Th_out, Tc_out

def read_th_in(csv_path: str, dt_default: float = 1.0) -> Tuple[List[float], List[float]]:
    """Read input CSV.
    Accepted formats:
      1) with header 't,Th_in' (seconds, degC)
      2) two columns without header -> interpreted as t, Th_in
      3) one column -> Th_in; time will be 0, dt, 2dt, ...
    """
    t: List[float] = []
    th: List[float] = []
    with open(csv_path, 'r', newline='') as f:
        r = csv.reader(f)
        rows = list(r)
    if not rows:
        return t, th
    # Detect header
    def is_float(s: str) -> bool:
        try:
            float(s)
            return True
        except:
            return False
    start_idx = 0
    if not all(is_float(x) for x in rows[0]):
        start_idx = 1  # skip header
    cols = len(rows[start_idx])
    if cols >= 2:
        for row in rows[start_idx:]:
            if len(row) < 2:
                continue
            try:
                t.append(float(row[0]))
                th.append(float(row[1]))
            except:
                continue
    else:
        # single column: synthesize time
        for i, row in enumerate(rows[start_idx:]):
            try:
                th.append(float(row[0]))
                t.append(i * dt_default)
            except:
                continue
    return t, th

def write_th_out(csv_path: str, t: List[float], th_out: List[float]) -> None:
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['t', 'Th_out'])
        for ti, tho in zip(t, th_out):
            w.writerow([ti, f\"{tho:.6f}\"])

def compute_series(t: List[float], th_in: List[float], Tc_in: float,
                   U: float, A: float, hot_flow_m3_h: float, cold_flow_m3_h: float) -> List[float]:
    UA = ua_from_spec(U, A)
    m_hot = hot_flow_m3_h * RHO_WATER / 3600.0
    m_cold = cold_flow_m3_h * RHO_WATER / 3600.0
    th_out: List[float] = []
    for Th_in in th_in:
        Th_out, _Tc_out = solve_hot_out(Th_in, Tc_in, m_hot, m_cold, UA)
        th_out.append(Th_out)
    return th_out

def main():
    ap = argparse.ArgumentParser(description='PHE: Th_out(t) from real Th_in(t) CSV (no plotting)')
    ap.add_argument('--input', required=True, help='CSV with t,Th_in or just Th_in per row')
    ap.add_argument('--output', default='/storage/emulated/0/Download/logs/phe_thout.csv', help='Output CSV path')
    ap.add_argument('--dt', type=float, default=1.0, help='dt if input CSV has single column')
    ap.add_argument('--Tc_in', type=float, default=32.0, help='Cold side inlet temperature (°C)')
    ap.add_argument('--U', type=float, default=6017.0, help='Overall heat transfer coeff (W/m^2·K)')
    ap.add_argument('--A', type=float, default=19.44, help='Heat transfer area (m^2)')
    ap.add_argument('--hot_flow', type=float, default=165.0, help='Hot-side flow (m^3/h)')
    ap.add_argument('--cold_flow', type=float, default=214.0, help='Cold-side flow (m^3/h)')
    args = ap.parse_args()

    t, th_in = read_th_in(args.input, dt_default=args.dt)
    if not t:
        raise SystemExit('Input CSV is empty or invalid. Expect header t,Th_in or single column of Th_in.')

    th_out = compute_series(t, th_in, args.Tc_in, args.U, args.A, args.hot_flow, args.cold_flow)
    write_th_out(args.output, t, th_out)
    print('Saved:', args.output)

if __name__ == '__main__':
    main()