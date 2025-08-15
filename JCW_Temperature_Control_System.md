
# Jacket Cooling Water (JCW) Temperature Control System

## 🔧 General Structure

This diagram represents the temperature control system of the Jacket Cooling Water (JCW) circuit on a marine main engine. It includes:

- **Control loop**
- **JCW thermal loop**
- **Main components**: Main Engine, Cooler, Bypass Line, Positioner, PID Controller

---

## 🧠 Operating Principle

- The temperature of the water exiting the main engine (`ToutME`) is the **main controlled variable (PV)**.
- The **PID controller** compares this temperature to a setpoint (SV), calculates the error, and generates a control signal `U`.
- The control signal `U` is sent to a **positioner** with only **two discrete states**:  
  - `0%` → Full bypass  
  - `50%` → Partial cooling
- Water flow either:
  - Goes through the **bypass line** directly (no cooling),
  - Or partially through the **JCW cooler**, which is cooled by the low-temperature cooling water (CW in, t = 32°C).
- The mixed water after the cooler and bypass forms the temperature `TinME`, which re-enters the main engine.
- `TinME` is a function of `ToutC` (cooler outlet) and `Tby-pass`.

---

## 🔍 System Components

### 🔹 Controller

- Standard PID formula:  
  `U = f(Error, P, I, D)`
- Process Variable (PV): `ToutME`
- Control loop marked with dashed lines.

### 🔹 Positioner

- Implements **discrete switching**:
  - `0%` → 100% bypass
  - `50%` → 50% cooler flow  
- Function defined as `f(U)`.

### 🔹 JCW Cooler

- Input: `Tinc = ToutME`
- Output: `ToutC`
- Cooling via external circuit (`CW in`, temp = 32°C)

### 🔹 Mixing Point

- `TinME = f(Tbypass, ToutC)`  
  Defines temperature returning to the engine as a mix of bypass and cooled water.

---

## 🧾 Comments from the Original Diagram

- `"PV is the principal controlled variable of the complete system"` – correct
- `"JCW Temp is a principal variable of the JCW contour"` – true for the engine’s thermal regime
- `"heat in / heat out"` – clearly labeled energy flows

---

## ❗ Special Note

Although unusual in practical marine systems, the positioner in this setup operates with only **two discrete positions**: `0%` and `50%`. This simplified model was intentionally used in the diagram.
