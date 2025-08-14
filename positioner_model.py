
from dataclasses import dataclass

@dataclass
class PositionerParams:
    up_threshold_ma: float = 16.0  # Rising threshold
    down_threshold_ma: float = 8.0 # Falling threshold
    low_pos_pct: float = 0.0       # "0% cooler"
    mid_pos_pct: float = 50.0      # "50% cooler"
    jump_time_s: float = 1.0       # time to complete jump (s)

class ElectropneumaticPositioner:
    """Electro-pneumatic positioner with hysteresis and jump transitions.
    Call step(i_ma) each dt to update position based on input current (mA).
    """
    def __init__(self, params: PositionerParams, dt: float):
        self.p = params
        self.dt = dt
        self.pos = self.p.low_pos_pct
        self._target = self.pos
        self._transition_left = 0.0
        self._prev_i_ma = None

    def step(self, i_ma: float) -> float:
        if self._prev_i_ma is None:
            self._prev_i_ma = i_ma

        # Hysteresis logic
        if self._prev_i_ma < self.p.up_threshold_ma <= i_ma:
            self._target = self.p.mid_pos_pct
            self._transition_left = self.p.jump_time_s
        elif self._prev_i_ma > self.p.down_threshold_ma >= i_ma:
            self._target = self.p.low_pos_pct
            self._transition_left = self.p.jump_time_s

        # Transition
        if self._transition_left > 0.0:
            step_frac = min(self.dt, self._transition_left) / max(1e-12, self._transition_left)
            self.pos = self.pos + (self._target - self.pos) * step_frac
            self._transition_left -= self.dt
            if self._transition_left <= 0.0:
                self.pos = self._target

        self._prev_i_ma = i_ma
        return self.pos

# Example usage
if __name__ == "__main__":
    # Demo: step through a sequence of input currents
    params = PositionerParams()
    pos_model = ElectropneumaticPositioner(params, dt=0.1)
    inputs = [4, 10, 16.1, 16.1, 19, 19, 15, 9, 8.0, 7.9, 7.9, 6]
    for i_ma in inputs:
        pos = pos_model.step(i_ma)
        print(f"Input: {i_ma:.1f} mA -> Position: {pos:.1f}% cooler")
