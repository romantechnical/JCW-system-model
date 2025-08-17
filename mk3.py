import numpy as np
import matplotlib.pyplot as plt
from collections import deque

# =========================================================
# USER PARAMETERS (edit here)
# =========================================================
# ---- Controller TZ4W (PB/Ti/Td) ----
PB          = 3.0    # proportional band, %  (Kp=100/PB)
Ti          = 50.0       # integral time, s      (Ki=Kp/Ti; 0 => off)
Td          = 15.0       # derivative time, s    (Kd=Kp*Td)
SAMPLE_TIME = 0.5       # controller period, s
REVERSE     = True      # PV>SP => U↑ (охолодження)

SP          = 86.0      # target engine outlet (°C)
LEAK_PCT    = 2.0       # % завжди на кулер, навіть у байпасі

# ---- Simulation time ----
DT_SIM      = 0.1       # integration step, s
T_END       = 2500.0    # total time, s

# ---- Plug-flow transport delays (s) ----
# «яка вода» приходить у відповідну точку в даний момент
DELAY_SENSOR_S          = 0.5   # від виходу двигуна до сенсора PV
DELAY_POS_S             = 5.0   # від виходу двигуна до того, «що бачить» позиціонер
DELAY_PHE_IN_S          = 7.0  # від виходу двигуна до входу кулера
DELAY_COOLER_TO_MIX_S   = 3.0   # від виходу кулера до точки змішування
DELAY_BYPASS_TO_MIX_S   = 9.0   # від виходу двигуна по байпасу до змішування
DELAY_MIX_TO_ENGIN_S    = 10.0   # від змішування до входу у двигун

# ---- Hydraulics / PHE ----
M_TOTAL      = 45.83              # kg/s  (165 m^3/h)
T_COLD_IN    = 32.0               # °C
U_PHE        = 6017.0             # W/(m^2·K)
A_PHE        = 19.44              # m^2
UA_PHE       = U_PHE * A_PHE
M_COLD       = 214000/3600.0      # kg/s (F.W. to cooler)
CP_WATER     = 4180.0             # J/(kg·K)

# ---- Engine thermal model ----
# Водяна частина моделюється як CSTR_N ємностей у серії (plug-flow усередині блоку).
CSTR_N       = 4                  # 1..4 зазвичай достатньо
C_W_TOTAL    = 1.2e7              # J/K (ефективна водяна ємність у блоці)
C_M          = 8.0e6              # J/K (ефективна металева ємність)
UA_MW        = 8.0e4              # W/K (метал ↔ вода сумарно)
Q_GEN        = 2.30e6             # W (тепловиділення на 75% навантаженні)
T_ENG_OUT0   = 84.0               # °C (початковий вихід двигуна)
T_ENG_IN0    = 72.0               # °C (початковий вхід двигуна)

# =========================================================
# Utilities
# =========================================================
def clamp(x, lo, hi): return max(lo, min(hi, x))

class PlugFlowDelay:
    """Кільцевий буфер для plug-flow затримок."""
    def __init__(self, delay_s, dt, init_value):
        n = max(1, int(round(delay_s / dt)))
        self.buf = deque([init_value]*n, maxlen=n)
    def push_pop(self, value):
        self.buf.append(value)
        return self.buf[0]

def lmtd(dT1, dT2):
    if dT1 <= 0 or dT2 <= 0:
        return 0.0
    if abs(dT1-dT2) < 1e-12:
        return dT1
    return (dT1 - dT2) / np.log(dT1 / dT2)

def phe_hot_out(Th_in, Tc_in, m_hot, m_cold, UA):
    """Рівноважний вихід гарячого потоку пластинчастого теплообмінника."""
    if m_hot <= 1e-9:
        return Th_in
    lo = min(Th_in-1e-6, Tc_in + 1e-6)
    hi = Th_in - 1e-6
    cp = CP_WATER

    def f(Th_out):
        Qh = m_hot * cp * (Th_in - Th_out)
        Tc_out = Tc_in + Qh / (max(1e-9, m_cold) * cp)
        return Qh - UA * lmtd(Th_in - Tc_out, Th_out - Tc_in)

    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return 0.5*(Th_in + Tc_in)
    for _ in range(60):
        mid = 0.5*(lo+hi)
        fm = f(mid)
        if abs(fm) < 1e-6: return mid
        if flo*fm < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5*(lo+hi)

# =========================================================
# Controllers / Actuators
# =========================================================
class PID_TZ4W:
    """
    TZ4W-style PID (PIDS): PB[%], Ti[s], Td[s]; derivative on measurement.
    Output: 4..20 mA. Reverse=True для охолодження.
    """
    def __init__(self, PB, Ti, Td, Ts, reverse=True, out_min_ma=4.0, out_max_ma=20.0):
        self.Kp = 0.0 if PB == 0 else 100.0/float(PB)
        self.Ki = 0.0 if (Ti is None or Ti <= 0) else self.Kp/float(Ti)
        self.Kd = self.Kp*float(Td)
        self.Ts = float(Ts)
        self.reverse = bool(reverse)
        self.out_min = float(out_min_ma)
        self.out_max = float(out_max_ma)
        self.I = 0.0
        self.last_meas = None
        self.last_u_pct = 0.0
        self.last_u_ma = 4.0

    def compute(self, SP, PV):
        e = (PV - SP) if self.reverse else (SP - PV)
        if self.last_meas is None:
            d_meas = 0.0
        else:
            d_meas = (PV - self.last_meas)/max(1e-9, self.Ts)
        self.last_meas = PV

        P = self.Kp*e
        D = self.Kd*d_meas
        u_unsat = P + self.I + D
        u_pct = clamp(u_unsat, 0.0, 100.0)
        if 0.0 < u_pct < 100.0 and self.Ki != 0.0:
            self.I += self.Ki*e*self.Ts
        self.last_u_pct = u_pct

        u_mA = 4.0 + (u_pct/100.0)*(20.0-4.0)
        u_mA = clamp(u_mA, self.out_min, self.out_max)
        self.last_u_ma = u_mA
        return u_mA

class Positioner:
    """Е-пневмо позиціонер з гістерезисом 16↑/8↓, 1 с хід, з «протічкою»."""
    def __init__(self, dt, up_th=16.0, down_th=8.0, jump_time=1.0, leak_frac=0.02):
        self.dt = dt
        self.up = up_th
        self.down = down_th
        self.jump_time = jump_time
        self.pos = 0.0
        self.target = 0.0
        self.t_left = 0.0
        self.prev_i = None
        self.leak = float(leak_frac)

    def step(self, i_mA):
        if self.prev_i is None: self.prev_i = i_mA
        if self.prev_i < self.up <= i_mA:
            self.target, self.t_left = 50.0, self.jump_time
        elif self.prev_i > self.down >= i_mA:
            self.target, self.t_left = 0.0, self.jump_time

        if self.t_left > 0.0:
            frac = min(self.dt, self.t_left)/max(1e-12, self.t_left)
            self.pos += (self.target - self.pos)*frac
            self.t_left -= self.dt
            if self.t_left <= 0.0: self.pos = self.target
        self.prev_i = i_mA

        f_cool = self.leak + (1.0 - self.leak)*(self.pos/100.0)
        return self.pos, f_cool

# =========================================================
# Engine model: metal lump + N water CSTRs in series (Heun)
# =========================================================
class EngineCSTR:
    def __init__(self, C_m, C_w_total, UA_mw, n_series, Tm0=109.0, Tw_out0=84.0):
        self.C_m = float(C_m)
        self.N = max(1, int(n_series))
        self.Cw = float(C_w_total)/self.N
        self.UA_seg = float(UA_mw)/self.N
        # ініціалізуємо всі сегменти води як Tw_out0
        self.Tw = np.full(self.N, float(Tw_out0))
        self.Tm = float(Tm0)

    def rhs(self, Tin, Qgen, m_dot):
        # послідовний потік через CSTRи
        Tw = self.Tw
        N = self.N
        UA = self.UA_seg
        Cp = CP_WATER

        dTw = np.zeros_like(Tw)
        Tin_i = Tin
        q_mw_sum = 0.0
        for i in range(N):
            q_mw = UA*(self.Tm - Tw[i])          # метал → вода, позитивний при Tm>Tw
            q_flow = m_dot*Cp*(Tin_i - Tw[i])    # конвекція
            dTw[i] = (q_flow + q_mw)/max(1e-9, self.Cw)
            q_mw_sum += q_mw
            Tin_i = Tw[i]                        # вихід і-го стає входом (i+1)-го
        dTm = (Qgen - q_mw_sum)/max(1e-9, self.C_m)
        return dTm, dTw

    def step(self, Tin, Qgen, m_dot, dt):
        k1m, k1w = self.rhs(Tin, Qgen, m_dot)
        Tm_e = self.Tm + dt*k1m
        Tw_e = self.Tw + dt*k1w
        # друга оцінка з екстрапольованих станів
        self.Tm = self.Tm + 0.5*dt*(k1m + self.rhs(Tin, Qgen, m_dot)[0])  # dTm знову від базових станів ладно
        # Точніша друга оцінка: перерахуємо RHS на екстрапольованих
        # (для води треба тимчасово підмінити self.Tm/self.Tw)
        Tm_save, Tw_save = self.Tm, self.Tw.copy()
        self.Tm, self.Tw = Tm_e, Tw_e
        k2m, k2w = self.rhs(Tin, Qgen, m_dot)
        self.Tm, self.Tw = Tm_save, Tw_save
        self.Tm += 0.5*dt*(k1m + k2m)
        self.Tw += 0.5*dt*(k1w + k2w)
        return self.Tw[-1]  # Tw_out

# =========================================================
# Simulation
# =========================================================
def run_sim():
    n = int(T_END/DT_SIM) + 1
    t = np.arange(n)*DT_SIM

    # блоки
    pid = PID_TZ4W(PB, Ti, Td, Ts=SAMPLE_TIME, reverse=REVERSE)
    positioner = Positioner(dt=DT_SIM, leak_frac=LEAK_PCT/100.0)
    engine = EngineCSTR(C_M, C_W_TOTAL, UA_MW, CSTR_N, Tm0=109.0, Tw_out0=T_ENG_OUT0)

    # затримки
    buf_sensor  = PlugFlowDelay(DELAY_SENSOR_S,        DT_SIM, T_ENG_OUT0)
    buf_pos     = PlugFlowDelay(DELAY_POS_S,           DT_SIM, T_ENG_OUT0)
    buf_phe_in  = PlugFlowDelay(DELAY_PHE_IN_S,        DT_SIM, T_ENG_OUT0)
    buf_coolmix = PlugFlowDelay(DELAY_COOLER_TO_MIX_S, DT_SIM, T_ENG_OUT0)
    buf_bypass  = PlugFlowDelay(DELAY_BYPASS_TO_MIX_S, DT_SIM, T_ENG_OUT0)
    buf_eng_in  = PlugFlowDelay(DELAY_MIX_TO_ENGIN_S,  DT_SIM, T_ENG_IN0)

    # логи
    PV = np.zeros(n)
    U  = np.zeros(n)
    POS = np.zeros(n)

    T_eng_out = T_ENG_OUT0
    next_ctrl = 0.0

    for k in range(n):
        # 1) вимір (сенсор)
        T_meas = buf_sensor.push_pop(T_eng_out)

        # 2) PID раз у SAMPLE_TIME
        if t[k] >= next_ctrl - 1e-12:
            U[k] = pid.compute(SP=SP, PV=T_meas)
            next_ctrl += SAMPLE_TIME
        else:
            U[k] = pid.last_u_ma

        # 3) позиціонер
        _ = buf_pos.push_pop(T_eng_out)
        pos_pct, f_cool = positioner.step(U[k])
        POS[k] = pos_pct

        # 4) розгалуження потоку
        m_hot_cool = f_cool * M_TOTAL
        m_bypass   = (1.0 - f_cool) * M_TOTAL

        # 5) кулер «бачить» затриману гарячу воду
        Th_in = buf_phe_in.push_pop(T_eng_out)
        if m_hot_cool > 1e-9:
            Th_out_cool = phe_hot_out(Th_in, T_COLD_IN, m_hot_cool, M_COLD, UA_PHE)
        else:
            Th_out_cool = Th_in

        # 6) транспорт після кулера і по байпасу до міксу
        T_cool_at_mix   = buf_coolmix.push_pop(Th_out_cool)
        T_bypass_at_mix = buf_bypass.push_pop(T_eng_out)

        # 7) змішування та затримка до входу в двигун
        T_mix = (m_bypass*T_bypass_at_mix + m_hot_cool*T_cool_at_mix)/max(1e-9, M_TOTAL)
        T_eng_in = buf_eng_in.push_pop(T_mix)

        # 8) динаміка двигуна → новий T_eng_out
        T_eng_out = engine.step(T_eng_in, Q_GEN, M_TOTAL, DT_SIM)

        # 9) лог
        PV[k] = T_eng_out

        # --- sanity checks ---
        if not np.isfinite(T_eng_out):
            raise RuntimeError("Non-finite T_eng_out encountered")

    return t, PV, U

def plot_results(t, PV, U):
    fig, ax1 = plt.subplots(figsize=(12,4.8))
    ax1.plot(t, PV, label="PV = T_eng_out (°C)", color="tab:blue")
    ax1.axhline(SP, color="tab:gray", linestyle="--", label="SP (°C)")
    ax1.set_xlabel("Time (s)"); ax1.set_ylabel("Temperature (°C)")
    ax1.set_title(f"JCW Plug-Flow — TZ4W: PB={PB}%, Ti={Ti}s, Td={Td}s | leakage={LEAK_PCT}% | CSTR_N={CSTR_N}")
    ax2 = ax1.twinx()
    ax2.plot(t, U, label="U (mA)", color="tab:orange")
    ax2.set_ylabel("Controller output (mA)")
    # legend
    l1, lab1 = ax1.get_legend_handles_labels()
    l2, lab2 = ax2.get_legend_handles_labels()
    ax1.legend(l1+l2, lab1+lab2, loc="best")
    fig.tight_layout(); plt.show()

if __name__ == "__main__":
    t, PV, U = run_sim()
    plot_results(t, PV, U)