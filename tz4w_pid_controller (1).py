
class TZ4WPIDController:
    """PID controller (simple form) that outputs 4–20 mA.
    
    - Internals compute an output in percent [0,100].
    - Returned value is current in mA, where 0% -> 4 mA and 100% -> 20 mA.
    """
    def __init__(self, Kp=0.4, Ki=60.0, Kd=15.0, dt=1.0,
                 percent_limits=(0.0, 100.0), ma_range=(4.0, 20.0)):
        self.Kp = float(Kp)
        self.Ki = float(Ki)
        self.Kd = float(Kd)
        self.dt = float(dt)
        self.pmin, self.pmax = percent_limits
        self.mamin, self.mamax = ma_range
        self.integral = 0.0
        self.prev_error = 0.0

    # ---- helpers ----
    def _clamp(self, v, lo, hi):
        return max(lo, min(hi, v))

    def _percent_to_mA(self, pct):
        # Map [pmin,pmax] -> [mamin,mamax]
        span_p = (self.pmax - self.pmin) if (self.pmax != self.pmin) else 1.0
        frac = (pct - self.pmin) / span_p
        return self.mamin + frac * (self.mamax - self.mamin)

    def _mA_to_percent(self, ma):
        # Map [mamin,mamax] -> [pmin,pmax]
        span_m = (self.mamax - self.mamin) if (self.mamax != self.mamin) else 1.0
        frac = (ma - self.mamin) / span_m
        return self.pmin + frac * (self.pmax - self.pmin)

    # ---- main step ----
    def compute(self, SV, PV):
        error = SV - PV

        # PID terms (discrete, derivative on error)
        P = self.Kp * error

        self.integral += error * self.dt
        I = self.Ki * self.integral

        derivative = (error - self.prev_error) / self.dt if self.dt > 0 else 0.0
        D = self.Kd * derivative

        # Raw percent command
        u_pct = P + I + D

        # Saturate in percent first
        u_pct = self._clamp(u_pct, self.pmin, self.pmax)

        # Save for next call
        self.prev_error = error

        # Convert to mA and clamp to [4,20] mA
        u_mA = self._percent_to_mA(u_pct)
        u_mA = self._clamp(u_mA, self.mamin, self.mamax)
        return u_mA
