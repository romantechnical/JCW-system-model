# Principles for Modeling the JCW Cooling System with Transport Delays

## 1) Time and Discretization
- Choose an integration step **Δt** (e.g., 0.1–0.2 s) — it must be:
  - Smaller than the smallest time constant in the loop (cooler, sensor filter, PID),  
  - Significantly smaller than the **minimum transport delay** of any section (see below).
- The PID can have **its own sampling period** (in your case, 0.5 s). This means the **model runs every Δt**, while **the PID updates less frequently** (every 0.5 s) and holds the last `U` between updates.

## 2) Objects and “Water Age” (Transport Delays)
Each pipe/volume defines a **transport delay**:
\[
	au = \frac{V}{\dot V} = \frac{\text{section volume}}{\text{flow rate}}
\]
where \(V\) is the hydraulic volume of the JCW section, \(\dot V\) is the volumetric flow rate.

**How to model transport simply and stably:**
- For each section, create a **ring buffer** (queue) with length \(N=\lceil 	au/\Delta t ceil\).  
  Each step, shift the temperature forward by one element — this is the **plug-flow (ideal piston flow)** method.
- If the flow rate changes → \(	au\) changes → update the buffer length and/or interpolate between cells.

This matches your description: at second 100, different points in the loop “see” the temperature that left the engine at 95, 90, 85… seconds.

## 3) Model Objects (States and Interfaces)

### 3.1 MainEngine (Heat Source)
- **Input:** \(T_{in,ME}(t)\) (from mixer), flow \(\dot V_{JCW}\).
- **State:** Effective thermal capacity of the “jacket + engine metal” \(C_{th}\) (J/K).
- **Heat balance model (FO/FO+Dead-Time):**
  \[
  C_{th}\,\frac{dT_{out,ME}}{dt} = Q_{in}(t) - \dot m\,c_p\,[T_{out,ME}(t)-T_{in,ME}(t)]
  \]
  where \(\dot m=ho\,\dot V\).  
  \(Q_{in}(t)\) — heat input from the engine (may be a function of load).
- **Output:** \(T_{out,ME}(t)\) (also the **PV** for PID after sensor filtering).

### 3.2 Sensor/Filter (Optional)
- **Input:** \(T_{out,ME}\); **output:** \(PV\).
- Simple first-order filter or moving average:  
  \(	au_s \dot{PV} = T_{out,ME}-PV\).  
  This adds realism (with your 0.5 s sampling).

### 3.3 PIDController
- **Input:** \(PV\) and setpoint \(SV\).
- **Update:** every 0.5 s. Between updates, output **U** is held constant.
- Anti-windup, D-filter, saturation — same as your real controller.

### 3.4 Positioner (Discrete, 0% / 50%)
- **Input:** \(U\) (with hysteresis/delay if present).
- **Output:** cooler flow fraction \(lpha \in \{0,\,0.5\}\).  
  Then \(\dot V_{cool}=lpha\,\dot V_{tot}\), \(\dot V_{bypass}=(1-lpha)\,\dot V_{tot}\).
- **Optional:** switching time (mechanical inertia) as a small delay/filter.

### 3.5 BypassLine (Transport)
- **Input:** \(T_{out,ME}\).
- **State:** buffer length \(N_{bp}\).
- **Output:** \(T_{bypass}= 	ext{shift}(T_{out,ME})\).

### 3.6 Cooler (Heat Exchange + Transport)
- **Inputs:** \(T_{in,C} = 	ext{shift}(T_{out,ME})\) (own buffer), \(\dot V_{cool}\); external circuit **CW_in = 32 °C** (constant or with its own dynamics).
- **Heat exchange:** convenient via **effectiveness-NTU** or simplified first-order:
  - NTU method: \(arepsilon = f(\mathrm{NTU}, C_r)\), \(Q = arepsilon\,Q_{max}\),  
    \(T_{out,C} = T_{in,C} - Q/(\dot m_{cool} c_p)\).
  - Simplified:  
    \(	au_C \dot T_{out,C} = T_{in,C}^{*} - T_{out,C}\), where \(T_{in,C}^{*}\) is the pseudo-output of ideal heat exchange to CW_in.
- **Output:** \(T_{out,C}\).

### 3.7 Mixer (Enthalpy Mixing)
- **Inputs:** \(T_{bypass}\), \(T_{out,C}\), corresponding flows \(\dot m_{bp}\), \(\dot m_{cool}\).
- **Model:**  
  \[
  T_{in,ME}=\frac{\dot m_{bp} c_p\,T_{bypass}+\dot m_{cool} c_p\,T_{out,C}}
                  {(\dot m_{bp}+\dot m_{cool}) c_p}
  pprox \frac{\dot m_{bp}T_{bypass}+\dot m_{cool}T_{out,C}}{\dot m_{tot}}
  \]
- **Optional:** own mixing volume (small buffer), if noticeable.

## 4) Computation Order at Each Δt
1. **Update transport buffers** for all sections (bypass and cooler line).  
   Get \(T_{bypass}(t)\), \(T_{in,C}(t)\).
2. **Cooler:** compute \(T_{out,C}(t)\) (NTU or simplified ODE).
3. **Mixer:** compute \(T_{in,ME}(t)\) from flows and temperatures.
4. **Engine:** integrate ODE to find \(T_{out,ME}(t+\Delta t)\).
5. **Sensor/filter:** update \(PV\).
6. **PID (every 0.5 s):** update `U` with anti-windup.
7. **Positioner:** based on `U` (with hysteresis) set \(lpha \in \{0,0.5\}\) → \(\dot V_{cool}, \dot V_{bp}\).
8. **Write to buffers:** push new \(T_{out,ME}\) into both transport queues.
9. Move to next step.

## 5) User-Defined vs. Computed Variables
**User-defined:**
- Geometry/volumes of sections, cooler parameters (UA or NTU data), \(\dot V_{tot}\) (or load-dependent curve), \(ho, c_p\).
- PID: \(P, I, D\), sampling period, setpoint \(SV\), limits on `U`, positioner hysteresis.
- Initial conditions: \(T\) in all buffers, \(T_{out,ME}(0)\), \(PV(0)\).

**Computed:**
- \(T_{in,ME}, T_{out,ME}, T_{in,C}, T_{out,C}, T_{bypass}\).
- \(lpha\), \(\dot V_{cool}, \dot V_{bp}\) (via positioner).
- \(PV\), `U`, PID error, integral state.

## 6) Details with Major Impact
- **Positioner hysteresis** (and small delay) — prevents 0/50% switching chatter.
- **PID anti-windup** and D filtering.
- **Time alignment:** PID at 0.5 s vs Δt=0.1 s; keep `U` between updates.
- **Stability:** if numerical oscillations occur — reduce Δt or add small mixing volumes (pseudo-capacities) at nodes.

## 7) Model Health Checks
- Step test: setpoint or load step → observe \(T_{out,ME}\), check realistic delay/overshoot.
- Energy balance: \(\int Q_{in} pprox \int Q_{cool} + C_{th}\Delta T\).
- Compare with real system logs (settling time, phase shifts).

---

**Summary:**  
- **Transport** — via **queue buffers** (plug-flow) length = delay/Δt.  
- **Mixing** — by enthalpy.  
- **Cooler** — NTU/efficiency or 1st-order with 32 °C CW reference.  
- **PID/positioner** — different timing, discrete choice 0%/50% + hysteresis.
