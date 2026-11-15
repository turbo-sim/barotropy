"""Saturation-based cases shared by regression tests and benchmarks."""

from dataclasses import dataclass
from scipy.integrate import solve_ivp
import CoolProp.CoolProp as CP

FLUIDS = ("CO2", "Water", "Nitrogen", "R134a", "Novec649")


@dataclass(frozen=True)
class Case:
    fluid: str
    family: str
    mode: str = "equilibrium"
    efficiency: float = 1.0

    @property
    def id(self):
        return f"{self.fluid}-{self.family}-{self.mode}-eta{self.efficiency:g}"

    def parameters(self):
        f, family = self.fluid, self.family
        tc = CP.PropsSI("Tcrit", f)
        sat = lambda t: CP.PropsSI("P", "T", t, "Q", 0, f)
        lo, hi = sat(0.8 * tc), sat(0.9 * tc)
        if family.startswith("gas_"):
            pin, tin = 0.2 * CP.PropsSI("pcrit", f), 1.1 * tc
            pout = pin * (0.7 if family.endswith("expansion") else 1.3)
            rho = CP.PropsSI("D", "P", pin, "T", tin, f)
        elif family.startswith("liquid_"):
            pin, tin = 1.5 * lo, 0.8 * tc
            pout = pin * (0.8 if family.endswith("expansion") else 1.2)
            rho = CP.PropsSI("D", "P", pin, "T", tin, f)
        elif family == "flashing":
            pin, tin, pout = 1.05 * hi, 0.9 * tc, lo
            rho = CP.PropsSI("D", "P", pin, "T", tin, f)
        elif family == "wet_liquid":
            pin, pout = lo, hi
            rho = CP.PropsSI("D", "P", pin, "Q", 0.05, f)
        elif family == "wet_inside":
            pin, pout = hi, 0.9 * hi
            rho = CP.PropsSI("D", "P", pin, "Q", 0.5, f)
        elif family in ("vapor_entry", "wet_vapor"):
            # Dry Novec crosses the vapor boundary in the opposite direction.
            fraction = 0.65 if f == "R134a" and self.efficiency < 1 else 0.85
            s = CP.PropsSI("S", "T", fraction * tc, "Q", 1, f)
            center = sat(fraction * tc)
            pin, pout = (
                (1.05 * center, 0.95 * center)
                if f != "Novec649"
                else (0.95 * center, 1.05 * center)
            )
            if f == "R134a" and self.efficiency < 1:
                pin, pout = 1.001 * center, 0.9 * center
            if family == "wet_vapor":
                pin, pout = pout, pin
            rho = CP.PropsSI("D", "P", pin, "S", s, f)
        elif family in ("meta_liquid", "meta_vapor"):
            q = 0 if family == "meta_liquid" else 1
            fraction = 0.65 if f == "R134a" and q == 1 and self.efficiency < 1 else 0.8
            center = sat(fraction * tc)
            s = CP.PropsSI("S", "T", fraction * tc, "Q", q, f)
            direction = -1 if q == 0 or f != "Novec649" else 1
            pin, pout = center * (1 - direction * 0.0001), center * (
                1 + direction * 0.005
            )
            rho = CP.PropsSI("D", "P", pin, "S", s, f)
        elif family in ("blend_liquid", "blend_vapor"):
            # Narrow, complete transition before homogeneous branch termination.
            q = 0 if family == "blend_liquid" else 1
            fraction = 0.65 if f == "R134a" and q == 1 and self.efficiency < 1 else 0.8
            center = sat(fraction * tc)
            s = CP.PropsSI("S", "T", fraction * tc, "Q", q, f)
            direction = -1 if q == 0 or f != "Novec649" else 1
            pin, pout = center * (1 - direction * 0.01), center * (1 + direction * 0.03)
            rho = CP.PropsSI("D", "P", pin, "S", s, f)
        else:
            raise ValueError(family)
        result = dict(
            fluid_name=f,
            p_in=pin,
            rho_in=rho,
            p_out=pout,
            efficiency=self.efficiency,
            calculation_type=self.mode,
            ODE_solver="RK45",
            ODE_tolerance=1e-8,
            HEOS_tolerance=1e-9,
            polynomial_degree=4,
        )
        if self.mode == "blending":
            vapor = (
                family.endswith("vapor")
                or family.startswith("gas")
                or family == "wet_inside"
            )
            if family.startswith("blend_"):
                oracle = CP.AbstractState("HEOS", f)
                oracle.update(CP.DmassP_INPUTS, rho, pin)
                factor = self.efficiency if pout < pin else 1 / self.efficiency

                def rhs(p, h):
                    oracle.update(CP.HmassP_INPUTS, h[0], p)
                    return [factor / oracle.rhomass()]

                ref = solve_ivp(
                    rhs, (pin, pout), [oracle.hmass()], rtol=1e-10, atol=1e-6
                )
                oracle.update(CP.HmassP_INPUTS, ref.y[0, -1], pout)
                distance = (1 - oracle.Q()) if vapor else oracle.Q()
                if not 0 < distance < 1:
                    raise ValueError(f"Invalid transition recipe: {self.id}")
                width = 0.3 * distance
                result.update(
                    blending_onset=1 - 0.1 * distance if vapor else 0.1 * distance,
                    blending_width=width,
                    ODE_max_step=abs(pout - pin) / 100,
                )
            else:
                result.update(
                    blending_onset=0.99 if vapor else 0.01, blending_width=0.01
                )
        return result


SINGLE = ("gas_expansion", "gas_compression", "liquid_expansion", "liquid_compression")
EQUILIBRIUM = SINGLE + (
    "flashing",
    "wet_liquid",
    "wet_inside",
    "vapor_entry",
    "wet_vapor",
)
CASES = tuple(
    Case(f, family, mode, eta)
    for f in FLUIDS
    for eta in (1.0, 0.8)
    for mode, families in (
        ("equilibrium", EQUILIBRIUM),
        ("metastable", SINGLE + ("meta_liquid", "meta_vapor")),
        ("blending", SINGLE + ("blend_liquid", "blend_vapor", "wet_inside")),
    )
    for family in families
)
