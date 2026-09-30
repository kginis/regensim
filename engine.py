"""Fast regen-cooling model: gas side (CEA, isentropic flow, film), coolant side
(pressure drop, thermal march) and structural margin, as callable functions.

main.py uses it for single runs, batch.py for sweeps.  The physics follows the
original main.py; speed comes from
  * gas side computed once per unique (MR, film, ...) set and cached,
  * vectorised isentropic lookup / film integral,
  * coolant properties from a pre-built (P, T) table instead of CoolProp calls,
  * loop-invariant work hoisted out of the wall-temperature iteration.
"""
import csv
import math
import tomllib
from pathlib import Path

import numpy as np

# Keys that live in the input file but are not model parameters.
RESERVED_SECTIONS = ("sweep", "batch", "constraints")

DEFAULTS = {
    "identifier": "Lynx",
    "nozzle_file": "nozzle.csv",
    # performance
    "P_Ambient": 101300.0, "Cstar_efficiency": 1.0,
    # chamber
    "throat_r_D": 1.5, "chamber_diameter": 0.0592, "total_mdot": 1.3,
    "Mass_Ratio": 2.0, "Expansion_Ratio": 3.5,
    # cooling settings
    "Two_Pass": False, "generatrix_angle": 0.0, "Constant_Rib": False,
    "Variable_Width": True, "Film_Cooling": True,
    # cooling inputs
    "Regen_Coolant": "Ethanol", "Film_Coolant": "Ethanol",
    "Surface_Roughness": 35e-6, "Film_Fraction": 0.25, "Film_Inlet_Temp": 375.0,
    "Coolant_Inlet_Temp": 298.15, "Coolant_Inlet_Pressure_Bar": 45.0,
    "Channel_Conductivity": 160.0,
    # channel geometry
    "Channel_Width": 0.0010, "Channel_Height": 0.0010, "Channel_Count": 46,
    "Channel_Wall": 0.0008, "Channel_Rib": 0.0010,
    "Downsteam_Pass_Channels": 18, "Upstream_Pass_Channels": 18,
    "pass_start_location": 0.0,
    "Channel_Width_Injector": 0.0024, "Channel_Width_Throat": 0.0010,
    "Constant_Rib_Len": 0.09, "Channel_Width_Manifold": 0.0014,
    # modifiers
    "x_pdms": 0.0, "pdms_modifier": 1.0, "performance_modifier": 1.0,
    "bartz_coeff": 0.85,
    # structural material
    "alfa": 2.0e-5, "poissons": 0.33, "youngsmodulus": 5.0e10, "Yield_S": 5.5e7,
    "shutdown_B": 0.5,
}

# Parameters that change the gas side (and therefore need a new CEA/film solve).
GAS_KEYS = (
    "Mass_Ratio", "Expansion_Ratio", "Cstar_efficiency", "total_mdot", "x_pdms",
    "throat_r_D", "chamber_diameter", "P_Ambient", "Film_Cooling", "Film_Fraction",
    "Film_Inlet_Temp", "Film_Coolant", "Regen_Coolant",
)


class ModelError(Exception):
    """A configuration could not be evaluated (bad geometry, non-convergence...)."""


# =============================================================================
# Input file
# =============================================================================
def load_inputs(path):
    """Read the single input file. Returns (cfg, raw) where cfg is the flat
    parameter dict and raw is the parsed TOML (incl. sweep/batch/constraints)."""
    path = Path(path)
    with path.open("rb") as f:
        raw = tomllib.load(f)
    cfg = dict(DEFAULTS)
    cfg["_dir"] = str(path.resolve().parent)
    for key, value in raw.items():
        if key in RESERVED_SECTIONS:
            continue
        items = value.items() if isinstance(value, dict) else [(key, value)]
        for k, v in items:
            if k not in DEFAULTS:
                raise ValueError(f"Unknown input '{k}' in {path.name}")
            cfg[k] = v
    return cfg, raw


def finalize(cfg):
    """Derived values + consistency checks. Returns a new dict."""
    c = dict(cfg)
    if c["Variable_Width"] and c["Constant_Rib"]:
        raise ValueError("You cannot have both a variable width and constant rib.")
    if c["Two_Pass"]:
        c["Channel_Count"] = int(c["Downsteam_Pass_Channels"] + c["Upstream_Pass_Channels"])
    c["Channel_Count"] = int(c["Channel_Count"])
    c["Coolant_Mdot"] = c["total_mdot"] / (1 + c["Mass_Ratio"])
    c["Film_Mdot"] = c["Coolant_Mdot"] * c["Film_Fraction"]
    c["Coolant_Inlet_Pressure_Pa"] = c["Coolant_Inlet_Pressure_Bar"] * 1e5
    c["generic_modifier"] = c["pdms_modifier"] * c["performance_modifier"]
    return c


def load_nozzle(path):
    x, r = [], []
    with open(path, newline="") as f:
        for row in csv.reader(f):
            try:
                xv, rv = float(row[0]), float(row[1])
            except (ValueError, IndexError):
                continue  # header
            x.append(xv)
            r.append(rv)
    return x, r


# =============================================================================
# Coolant property table
# =============================================================================
class CoolantTable:
    """Bilinear (P, T) table of cp, rho, mu, k built once from CoolProp (pyfluids)."""

    def __init__(self, fluid_name, p_max_bar, cache_dir=None,
                 p_min_bar=1.0, dp_bar=1.0, t_min=250.0, t_max=560.0, dt=1.0):
        from pyfluids import Fluid, FluidsList, Input

        self.fluid_name = fluid_name
        p_max_bar = math.ceil(p_max_bar * 1.02) + 1
        self.p0 = p_min_bar * 1e5
        self.dp = dp_bar * 1e5
        self.t0, self.dt = t_min, dt
        self.np = int(round((p_max_bar - p_min_bar) / dp_bar)) + 1
        self.nt = int(round((t_max - t_min) / dt)) + 1
        self.p_max = self.p0 + (self.np - 1) * self.dp
        self.t_max = self.t0 + (self.nt - 1) * self.dt
        self.inv_dp, self.inv_dt = 1.0 / self.dp, 1.0 / self.dt

        cache = None
        if cache_dir is not None:
            cache = Path(cache_dir) / f".coolant_{fluid_name}_{self.np}x{self.nt}_{t_min:g}_{dt:g}.npz"
        if cache is not None and cache.exists():
            data = np.load(cache)
            grid, tsat = data["grid"], data["tsat"]
        else:
            fl = getattr(FluidsList, fluid_name)
            grid = np.zeros((self.np, self.nt, 4))
            for ip in range(self.np):
                p = self.p0 + ip * self.dp
                for it in range(self.nt):
                    t = self.t0 + it * self.dt
                    try:
                        s = Fluid(fl).with_state(Input.pressure(p), Input.temperature(t - 273.15))
                        grid[ip, it] = (s.specific_heat, s.density, s.dynamic_viscosity, s.conductivity)
                    except Exception:
                        grid[ip, it] = np.nan
            # Two-phase / failed nodes: fill from the neighbouring node in T.
            for ip in range(self.np):
                for it in range(self.nt):
                    if np.isnan(grid[ip, it, 0]):
                        src = it - 1 if it > 0 else it + 1
                        grid[ip, it] = grid[ip, src]
            tsat = np.zeros(self.np)
            for ip in range(self.np):
                p = self.p0 + ip * self.dp
                try:
                    tsat[ip] = Fluid(fl).dew_point_at_pressure(p).temperature + 273.15
                except Exception:
                    tsat[ip] = 1e9  # supercritical: no boiling
            if cache is not None:
                try:
                    np.savez(cache, grid=grid, tsat=tsat)
                except OSError:
                    pass
        self.nodes = [tuple(map(float, grid[ip, it])) for ip in range(self.np) for it in range(self.nt)]
        self.tsat_nodes = [float(v) for v in tsat]

    def lookup(self, p, t):
        """Returns (cp, rho, mu, k)."""
        fp = (p - self.p0) * self.inv_dp
        ft = (t - self.t0) * self.inv_dt
        if fp < 0.0 or fp > self.np - 1 or ft < 0.0 or ft > self.nt - 1:
            raise ModelError(f"coolant state out of table range (P={p/1e5:.2f} bar, T={t:.1f} K)")
        ip = int(fp)
        if ip > self.np - 2:
            ip = self.np - 2
        it = int(ft)
        if it > self.nt - 2:
            it = self.nt - 2
        wp, wt = fp - ip, ft - it
        i = ip * self.nt + it
        a, b = self.nodes[i], self.nodes[i + 1]
        c, d = self.nodes[i + self.nt], self.nodes[i + self.nt + 1]
        w00, w01, w10, w11 = (1 - wp) * (1 - wt), (1 - wp) * wt, wp * (1 - wt), wp * wt
        return (
            a[0] * w00 + b[0] * w01 + c[0] * w10 + d[0] * w11,
            a[1] * w00 + b[1] * w01 + c[1] * w10 + d[1] * w11,
            a[2] * w00 + b[2] * w01 + c[2] * w10 + d[2] * w11,
            a[3] * w00 + b[3] * w01 + c[3] * w10 + d[3] * w11,
        )

    def tsat(self, p):
        fp = (p - self.p0) * self.inv_dp
        ip = min(max(int(fp), 0), self.np - 2)
        w = min(max(fp - ip, 0.0), 1.0)
        return self.tsat_nodes[ip] * (1 - w) + self.tsat_nodes[ip + 1] * w


# =============================================================================
# Gas side (cached per unique gas-side parameter set)
# =============================================================================
_CEA_CACHE = {}


def _get_cea(x_pdms):
    if x_pdms in _CEA_CACHE:
        return _CEA_CACHE[x_pdms]
    import os
    import tempfile
    import rocketcea.cea_obj as cea_mod
    from rocketcea.cea_obj import add_new_fuel, add_new_oxidizer
    from rocketcea.cea_obj_w_units import CEA_Obj

    # rocketcea keeps its Fortran scratch files (temp.dat, f.out...) in one shared
    # directory; give every process its own so parallel workers can't corrupt each other.
    cea_mod.ROCKETCEA_DATA_DIR = os.path.join(tempfile.gettempdir(), f"regensim_cea_{os.getpid()}")
    from rocketcea.units import add_user_units

    fuel_name = f"IPA_PDMS_{x_pdms:g}wt"
    add_new_fuel(fuel_name, f"""
fuel C3H8O(L)  C 3 H 8 O 1  wt%={100.0 - x_pdms:.6f}
h,cal=-65133  t(k)=298.15  rho.g/cc=0.803

fuel PDMS(L)  C 2 H 6 O 1 Si 1  wt%={x_pdms:.6f}
h,cal=-178712.0  t(k)=298.15  rho.g/cc=0.970
""")
    add_new_oxidizer("N2O_L", """
oxid N2O(L)  N 2 O 1  wt%=100.0
h,cal=19721.03  t(k)=298.15
""")
    add_user_units("millipoise", "Pa-s", 1e-4)
    obj = CEA_Obj(
        oxName="N2O_L", fuelName=fuel_name,
        cstar_units="m/s", pressure_units="Bar", temperature_units="K",
        specific_heat_units="J/kg-K", viscosity_units="Pa-s",
        density_units="kg/m^3", enthalpy_units="J/kg",
    )
    _CEA_CACHE[x_pdms] = obj
    return obj


def hg_bartz(radius, Cp, mu, Pr, Pc, C_star, Area_Ratio, throat_curvature_radius):
    D_star = float(radius * 2)
    return (
        (0.026 / D_star**0.2) * (((mu**0.2) * Cp) / Pr**0.6) * ((Pc / C_star) ** 0.8)
        * ((D_star / throat_curvature_radius) ** 0.1) * (Area_Ratio) ** -0.9
    )


def bartz_boundary_sigma(Tw, Tc, gamma, Mach, omega):
    M_gamma = 1 + ((gamma - 1) / 2) * Mach**2
    return 1 / ((0.5 * (Tw / Tc) * M_gamma + 0.5) ** (0.8 - omega / 5) * M_gamma ** (omega / 5))


class Gas:
    """Plain container for the cached gas-side solution."""


def build_gas(cfg, nozzle):
    """CEA + isentropic flow + film cooling for one gas-side parameter set."""
    from pyfluids import Fluid, FluidsList, Input

    x = np.array(nozzle[0])
    r = np.array(nozzle[1])
    n = len(x)
    g = Gas()
    g.x, g.r, g.n = x, r, n
    g.throat_radius = float(r.min())
    throat_area = math.pi * g.throat_radius**2
    MR, eps = cfg["Mass_Ratio"], cfg["Expansion_Ratio"]
    mdot = cfg["total_mdot"]
    cea = _get_cea(cfg["x_pdms"])

    # Chamber pressure iteration (Pc depends on c*, c* on Pc)
    Pc_bar = 35.0
    Pc_Pa = Pc_bar * 1e5
    for _ in range(200):
        chamber_transport = cea.get_Chamber_Transport(Pc=Pc_bar, MR=MR, frozen=1)
        mw, gamma = cea.get_Chamber_MolWt_gamma(Pc=Pc_bar, MR=MR, eps=eps)
        gas_enthalpy = cea.get_Chamber_H(Pc=Pc_bar, MR=MR, eps=eps)
        cstar = cea.get_Cstar(Pc=Pc_bar, MR=MR) * cfg["Cstar_efficiency"]
        chamber_temp = cea.get_Temperatures(Pc=Pc_bar, MR=MR)[0]
        chamber_rho = cea.get_Chamber_Density(Pc=Pc_bar, MR=MR)
        Pc_true = cstar * mdot / throat_area
        dev = abs(Pc_true - Pc_Pa)
        Pc_Pa = Pc_true
        Pc_bar = Pc_Pa / 1e5
        if dev <= 0.01:
            break
    else:
        raise ModelError("chamber pressure iteration did not converge")

    g.gas_Cp, g.gas_mu, g.gas_Pr = chamber_transport[0], chamber_transport[1], chamber_transport[3]
    g.gamma, g.mw, g.gas_enthalpy = gamma, mw, gas_enthalpy
    g.cstar, g.chamber_temp, g.chamber_rho = cstar, chamber_temp, chamber_rho
    g.Pc_Pa, g.Pc_bar = Pc_Pa, Pc_bar
    R_gas = 8314 / mw
    g.R_gas = R_gas

    # Isentropic flow, nearest-Mach lookup on the same 0.01 grid as before
    mach = np.arange(1, 601) * 0.01
    ar_tab = (1.0 / mach) * ((1 + 0.5 * (gamma - 1.0) * mach**2.0) / (1.0 + 0.5 * (gamma - 1.0))) ** (
        (gamma + 1.0) / (2.0 * (gamma - 1.0)))
    pr_tab = (1 + 0.5 * (gamma - 1) * mach**2) ** (-(gamma / (gamma - 1)))
    tr_tab = (1 + 0.5 * (gamma - 1) * mach**2) ** (-1)
    g.area_ratio = (r**2 * math.pi) / throat_area
    sub = mach < 1
    idx = np.empty(n, dtype=int)
    for mask, sel in ((x <= 0.0, np.where(sub)[0]), (x > 0.0, np.where(~sub)[0])):
        rows = np.where(mask)[0]
        if len(rows):
            d = np.abs(ar_tab[sel][None, :] - g.area_ratio[rows][:, None])
            idx[rows] = sel[d.argmin(axis=1)]
    g.mach = mach[idx]
    g.gas_temp = chamber_temp * tr_tab[idx]
    g.gas_pressure = Pc_Pa * pr_tab[idx]
    g.gas_velocity = g.mach * np.sqrt(gamma * R_gas * g.gas_temp)
    rec = g.gas_Pr ** (1 / 3)
    g.taw_uncooled = g.gas_temp * (1 + rec * ((gamma - 1) / 2) * g.mach**2)
    g.taw = g.taw_uncooled.copy()

    # Segment geometry (segment i runs from station i-1 to i)
    dx = np.abs(np.diff(x))
    dr = np.diff(r)
    g.ds = np.concatenate([[0.0], np.sqrt(dx**2 + dr**2)])
    g.gas_area = np.concatenate([[0.0], math.pi * (r[1:] + r[:-1]) * g.ds[1:]])

    # Gas-side Bartz term that does not depend on wall temperature
    curv = cfg["throat_r_D"] * (2 * g.throat_radius)
    g.hg_base = np.array([
        hg_bartz(g.throat_radius, g.gas_Cp, g.gas_mu, g.gas_Pr, Pc_Pa, cstar, ar, curv)
        for ar in g.area_ratio
    ])
    M_gamma = 1 + ((gamma - 1) / 2) * g.mach**2
    g.sigma_a = 0.5 * M_gamma / chamber_temp      # sigma = 1/((a*Tw+0.5)^0.68 * Mg^0.12)
    g.sigma_m = M_gamma**0.12

    # Performance
    g.exit_pressure = float(g.gas_pressure[-1])
    g.exit_area = math.pi * float(r[-1]) ** 2
    g.exit_velocity = float(g.gas_velocity[-1])
    g.thrust = g.exit_velocity * mdot + (g.exit_pressure - cfg["P_Ambient"]) * g.exit_area
    g.isp = g.thrust / (mdot * 9.80665)

    # Film cooling
    g.film_state = ["None"] * n
    g.film_data = [dict() for _ in range(n)]
    g.film_summary = {}
    g.f_liquid_len = 0.0
    if cfg["Film_Cooling"]:
        _film(cfg, g, Fluid, FluidsList, Input, curv)
    return g


def _film(cfg, g, Fluid, FluidsList, Input, curv):
    film = getattr(FluidsList, cfg["Film_Coolant"])
    x, r, n = g.x, g.r, g.n
    Film_Mdot = cfg["Film_Mdot"]
    total_mdot = cfg["total_mdot"]
    gamma, gas_Cp, Pr = g.gamma, g.gas_Cp, g.gas_Pr
    Pc_Pa = g.Pc_Pa

    hg_1 = g.hg_base[1] * bartz_boundary_sigma(
        g.taw_uncooled[1], g.chamber_temp, gamma, g.mach[1], 0.6)

    # Huzel & Huang liquid film length
    delta = 0.5
    liquid = Fluid(film).bubble_point_at_pressure(Pc_Pa)
    vapour = Fluid(film).dew_point_at_pressure(Pc_Pa)
    T_sat = vapour.temperature + 273.15
    latent = vapour.enthalpy - liquid.enthalpy
    T_re = g.chamber_temp * (1 + (Pr ** (1 / 3) * ((gamma - 1) / 2) * g.mach[1] ** 2))
    B = gas_Cp * (T_re - T_sat) / latent
    nu_fcl = B / (B + 1)
    surface_tens = Fluid(film).with_state(Input.quality(0), Input.temperature(T_sat - 273.15)).surface_tension
    Xe = delta * (g.chamber_rho**0.5) * g.gas_velocity[1] * ((g.chamber_temp / T_sat) ** 0.25) / surface_tens
    Xr = Xe * surface_tens
    alfa = (1 + 3 * Xr**-0.8) * (7 / (20 * 10**4) * Xe + 1)
    A_f_C = ((1.3 / 200000 * Xe) + 0.1) / 0.0254
    St = hg_1 / (g.chamber_rho * g.gas_velocity[1] * gas_Cp)
    V = (math.pi * cfg["chamber_diameter"] * (g.chamber_rho * g.gas_velocity[1])) * St * B * alfa
    f_len = (1 / A_f_C) * math.log(1 + ((A_f_C * Film_Mdot) / V))
    g.f_liquid_len = f_len
    g.film_summary.update({
        "Liquid film length": (f_len, "m"),
        "Liquid film end X": (x[0] + f_len, "m"),
        "Film saturation temperature": (T_sat, "K"),
        "Film heat-input parameter B": (B, "-"),
        "Film latent heat": (latent, "J/kg"),
        "Film nu at liquid-film end": (nu_fcl, "-"),
        "Film entrainment parameter Xe": (Xe, "correlation units"),
        "Film roughness parameter Xr": (Xr, "correlation units"),
        "Film heat-transfer multiplier alpha": (alfa, "-"),
        "Film entrainment coefficient A": (A_f_C, "1/m"),
        "Film Stanton number": (St, "-"),
        "Film recovery temperature used": (T_re, "K"),
        "Film surface tension used": (surface_tens, "N/m"),
        "Film vaporization rate per length V": (V, "kg/(m*s)"),
        "Film reference heat-transfer coefficient": (hg_1, "W/(m^2*K)"),
    })

    liquid_mask = np.abs(x[0] - x) <= f_len
    gas_idx = np.where(~liquid_mask)[0]
    for i in np.where(liquid_mask)[0]:
        g.film_state[i] = "Liquid"
        g.film_data[i]["film_temperature"] = T_sat
    g.taw[liquid_mask] = T_sat

    if len(gas_idx):
        # phi_m (entrainment multiplier) per station
        phi = np.where(
            x <= 0,
            1.75 + x * (3.5 - 1.75) / x[0],
            0.277 + 1.4961 * np.exp(-0.1199 * g.area_ratio) + 0.14612 / g.area_ratio,
        )
        rho_g = g.gas_pressure / (g.R_gas * g.gas_temp)
        flux = rho_g * g.gas_velocity
        integrand = (r[0] / r) * (flux / flux[0]) * phi
        x_start = x[0] + f_len
        xa, xb = x[:-1], x[1:]
        left = np.maximum(xa, x_start)
        frac = (left - xa) / (xb - xa)
        f_left = integrand[:-1] + frac * (integrand[1:] - integrand[:-1])
        seg = np.where(xb > x_start, 0.5 * (f_left + integrand[1:]) * (xb - left), 0.0)
        cum = np.concatenate([[0.0], np.cumsum(seg)])
        x_bar = np.where(x > x_start, cum, 0.0)

        core_mdot = total_mdot - Film_Mdot
        if core_mdot <= 0 or Film_Mdot <= 0 or nu_fcl <= 0:
            raise ModelError("Film flow, core flow, and nu must be positive")
        liquid_ratio = 1 / (0.6 * nu_fcl) - 1
        mass_term = Film_Mdot * liquid_ratio / core_mdot
        dist_term = 0.004 * x_bar / r
        radicand = 1 - mass_term - dist_term
        if np.any(radicand[gas_idx] < 0):
            raise ModelError("film entrainment radicand < 0 (film too long / too much liquid)")
        radicand = np.maximum(radicand, 0.0)
        entr = (core_mdot / Film_Mdot) * 2 * dist_term * np.sqrt(radicand) + liquid_ratio
        ent_liq = Film_Mdot * liquid_ratio
        ent_i = Film_Mdot * entr
        shape = np.where(ent_i > ent_liq + 0.6 * Film_Mdot,
                         0.6 + 0.263 * ((ent_i - ent_liq) / Film_Mdot), 0.758)
        eta = 1 / (shape * (1 + entr))

        # vapour film coolant properties
        p1 = g.gas_pressure[1]
        cp_ref = Fluid(film).with_state(Input.quality(100), Input.pressure(p1))
        R_co = 8.314462618 / cp_ref.molar_mass
        gamma_co = cp_ref.specific_heat / (cp_ref.specific_heat - R_co)
        T_co = cfg["Film_Inlet_Temp"] * (g.gas_pressure / p1) ** ((gamma_co - 1) / gamma_co)
        T_ref = 298.15
        for i in gas_idx:
            Cp_i = Fluid(film).with_state(Input.quality(100), Input.pressure(g.gas_pressure[i])).specific_heat
            hg_i = gas_Cp * (g.gas_temp[i] - T_ref)
            hco_i = Cp_i * (T_co[i] - T_ref)
            e = eta[i]
            dH = g.gas_enthalpy - hco_i
            g.taw[i] = ((hg_i - e * dH) - (1 - Pr ** (1 / 3)) * dH) / (e * Cp_i + (1 - e) * gas_Cp)
            g.film_state[i] = "Gas"
            g.film_data[i].update({
                "film_temperature": T_co[i], "cp_used": Cp_i, "eta": e,
                "cp_mix": e * Cp_i + (1 - e) * gas_Cp,
                "coolant_sensible_enthalpy": hco_i, "gas_sensible_enthalpy": hg_i,
                "phi": phi[i], "x_bar": x_bar[i], "entrainment_ratio": entr[i],
                "shape_factor": shape[i], "gamma_used": gamma_co,
            })
    g.film_summary["Film liquid-station count"] = (int(liquid_mask.sum()), "-")
    g.film_summary["Film gas-station count"] = (len(gas_idx), "-")


_GAS_CACHE = {}


def get_gas(cfg, nozzle):
    key = tuple(cfg[k] for k in GAS_KEYS)
    g = _GAS_CACHE.get(key)
    if g is None:
        if len(_GAS_CACHE) > 64:
            _GAS_CACHE.clear()
        g = _GAS_CACHE[key] = build_gas(cfg, nozzle)
    return g


# =============================================================================
# Channel geometry
# =============================================================================
def _profile_coeffs(cfg, gas):
    """Variable-width profile as width_i = a_i*W_inj + b_i*W_throat + c_i*W_manifold."""
    xs = [float(v) for v in gas.x]
    n = gas.n
    end = min(range(n), key=lambda i: abs((xs[i] - xs[0]) - cfg["Constant_Rib_Len"]))
    coeffs = []
    for i in range(n):
        if xs[i] < xs[end]:
            coeffs.append((1.0, 0.0, 0.0))
        elif xs[i] <= 0:
            frac = (xs[i] - xs[end]) / -xs[end] if xs[end] != 0 else 1.0
            coeffs.append((1 - frac, frac, 0.0))
        else:
            frac = xs[i] / xs[-1]
            coeffs.append((0.0, 1 - frac, frac))
    return coeffs


def resolve_widths(cfg, gas):
    """Replace "max" channel widths with the widest value that keeps every rib
    >= cfg['Min_Rib'].  Resolved in priority order throat -> injector -> manifold
    (for a constant width, the throat is the binding station)."""
    c = dict(cfg)
    N = c["Channel_Count"]
    rib = c.get("Min_Rib", 0.0012)
    circ = [2 * math.pi * float(ri) / N for ri in gas.r]
    for key in ("Channel_Width", "Channel_Width_Throat", "Channel_Width_Injector", "Channel_Width_Manifold"):
        if isinstance(c[key], str) and c[key].lower() != "max":
            raise ValueError(f"{key} must be a number or \"max\"")
    if isinstance(c["Channel_Width"], str):
        c["Channel_Width"] = min(circ) - rib
    if c["Variable_Width"]:
        names = ("Channel_Width_Injector", "Channel_Width_Throat", "Channel_Width_Manifold")
        if any(isinstance(c[k], str) for k in names):
            coeffs = _profile_coeffs(c, gas)
            cur = [0.0 if isinstance(c[k], str) else c[k] for k in names]
            for key in ("Channel_Width_Throat", "Channel_Width_Injector", "Channel_Width_Manifold"):
                if not isinstance(c[key], str):
                    continue
                j = names.index(key)
                best = math.inf
                for ci, co in zip(circ, coeffs):
                    if co[j] > 1e-12:
                        other = sum(co[m] * cur[m] for m in range(3) if m != j)
                        best = min(best, (ci - rib - other) / co[j])
                c[key] = cur[j] = best
    return c


def geometry(cfg, gas):
    x, r, n = gas.x, gas.r, gas.n
    N, H = cfg["Channel_Count"], cfg["Channel_Height"]
    circ = [2 * math.pi * float(ri) / N for ri in r]
    widths, ribs = [], []
    if cfg["Constant_Rib"]:
        for c in circ:
            widths.append(c - cfg["Channel_Rib"])
            ribs.append(cfg["Channel_Rib"])
    elif cfg["Variable_Width"]:
        w_inj, w_thr, w_man = (cfg["Channel_Width_Injector"], cfg["Channel_Width_Throat"],
                               cfg["Channel_Width_Manifold"])
        for c, (a, b, m) in zip(circ, _profile_coeffs(cfg, gas)):
            w = a * w_inj + b * w_thr + m * w_man
            widths.append(w)
            ribs.append(c - w)
    else:
        for c in circ:
            widths.append(cfg["Channel_Width"])
            ribs.append(c - cfg["Channel_Width"])
    if min(widths) <= 0 or min(ribs) <= 0 or H <= 0:
        raise ModelError("non-positive channel width/rib/height")
    area = [w * H for w in widths]
    perim = [2 * H + 2 * w for w in widths]
    hr = [a / p for a, p in zip(area, perim)]
    dh = [4 * v for v in hr]
    return dict(widths=widths, ribs=ribs, area=area, perim=perim, hr=hr, dh=dh)


# =============================================================================
# Coolant side + structural
# =============================================================================
def simulate(cfg, gas, tab, detail=False, margin_floor=None, min_rib=None):
    """Run the coolant-side thermal march and structural check.

    margin_floor: if set (single-pass only), stop as soon as any station's
    structural margin drops below it (result['aborted'] = True).
    """
    n = gas.n
    cfg = finalize(resolve_widths(cfg, gas))
    geo = geometry(cfg, gas)
    widths, ribs, area = geo["widths"], geo["ribs"], geo["area"]
    perim, hr, dh = geo["perim"], geo["hr"], geo["dh"]
    if min_rib is not None and min(ribs) < min_rib - 1e-9:
        return {"rejected": "rib", "min_rib_m": min(ribs), "min_width_m": min(widths),
                "max_width_m": max(widths), "min_channel_area_m2": min(area),
                "mean_channel_area_m2": sum(area) / n}

    two_pass = cfg["Two_Pass"]
    N = cfg["Channel_Count"]
    H = cfg["Channel_Height"]
    wall = cfg["Channel_Wall"]
    kc = cfg["Channel_Conductivity"]
    mdot = cfg["Coolant_Mdot"]
    p_in = cfg["Coolant_Inlet_Pressure_Pa"]
    t_in = cfg["Coolant_Inlet_Temp"]
    gm = cfg["generic_modifier"]
    bartz_coeff = cfg["bartz_coeff"]
    theta = math.radians(cfg["generatrix_angle"])
    cos_t = math.cos(theta)
    x = gas.x
    xs = [float(v) for v in x]

    cp_in, rho_in, mu_in, k_in = tab.lookup(p_in, t_in)
    nu_in = mu_in / rho_in

    # Friction factor (Colebrook), coolant velocity
    rough = cfg["Surface_Roughness"]
    f_list, v_ret, v_down = [], [], []
    n_up = cfg["Upstream_Pass_Channels"] if two_pass else N
    n_dn = cfg["Downsteam_Pass_Channels"] if two_pass else N
    for i in range(n):
        Re = 4 * (mdot / N) / (nu_in * perim[i] * rho_in)
        f = 0.02
        for _ in range(200):
            f_new = 1 / (-2 * math.log10(rough / (14.8 * hr[i]) + 2.51 / (Re * math.sqrt(f)))) ** 2
            d = abs(f_new - f)
            f = f_new
            if d < 1e-6:
                break
        f_list.append(f)
        v_ret.append(mdot / (rho_in * area[i] * n_up))
        if two_pass:
            v_down.append(mdot / (rho_in * area[i] * n_dn))
    coolant_velocity = v_ret                       # original naming
    coolant_velocity_downstream = v_down

    # pass start station
    pass_start = min(range(n), key=lambda i: abs((xs[i] - xs[0]) - cfg["pass_start_location"]))
    if two_pass and xs[-1] - xs[0] <= cfg["pass_start_location"]:
        raise ModelError("pass start location is larger than total nozzle length")

    def pass_pressures(stations, inlet, vel):
        p = [None] * n
        prev = stations[0]
        p[prev] = inlet
        for i in stations[1:]:
            seg = abs(xs[i] - xs[prev]) / cos_t
            f = 0.5 * (f_list[prev] + f_list[i])
            d = 0.5 * (dh[prev] + dh[i])
            v = 0.5 * (vel[prev] + vel[i])
            p[i] = p[prev] - f * (seg / d) * (rho_in * v * v / 2)
            prev = i
        return p

    if two_pass:
        p_first = pass_pressures(list(range(pass_start, n)), p_in, v_down)
        p_ret_in = p_first[-1]
    else:
        p_first, p_ret_in = None, p_in
    p_ret = pass_pressures(list(range(n - 1, -1, -1)), p_ret_in, v_ret)
    pressures = p_first if two_pass else p_ret

    if two_pass:
        thermal_indices = list(range(pass_start + 1, n))
        n_down, n_return = cfg["Downsteam_Pass_Channels"], cfg["Upstream_Pass_Channels"]
    else:
        thermal_indices = list(range(n - 1, 0, -1))
        n_down, n_return = N, 0

    ds, garea, taw = gas.ds, gas.gas_area, gas.taw
    hg_base, sig_a, sig_m = gas.hg_base, gas.sigma_a, gas.sigma_m
    lookup = tab.lookup

    # structural constants
    E, alfa, nu_p = cfg["youngsmodulus"], cfg["alfa"], cfg["poissons"]
    yield_s = cfg["Yield_S"]
    S_mech = cfg["shutdown_B"] * p_in   # times (a/h)^2 per station below
    tan_coeff = E * alfa * wall / (2 * (1 - nu_p) * kc)
    ax_coeff = E * alfa
    mech_h2 = 1.0 / (wall * wall)

    def hl_rpe(cp, rho, mu, k, qty, vel, i):
        return 0.023 * cp * (mdot / (area[i] * qty)) * ((dh[i] * vel * rho) / mu) ** -0.2 * ((mu * cp) / k) ** (-2 / 3)

    def fin_area(hl, i, qty, length):
        m_fin = math.sqrt(2 * hl / (kc * ribs[i]))
        arg = m_fin * H
        eff = math.tanh(arg) / arg if arg > 1e-12 else 1.0
        return (widths[i] + 2 * eff * H) * qty * length

    def wall_solve(i, Taw_i, coolant_T_eff, R_l, R_w, A):
        Tg = 0.5 * (Taw_i + coolant_T_eff)
        hgb = hg_base[i] * bartz_coeff
        a_i, m_i = sig_a[i], sig_m[i]
        for _ in range(150):
            hg = hgb / (((a_i * Tg + 0.5) ** 0.68) * m_i)
            R_g = 1.0 / (hg * A)
            R_base = R_g + R_w + R_l
            R_mod = R_base * (1.0 / gm - 1.0)
            Q = (Taw_i - coolant_T_eff) / (R_base + R_mod)
            Tg_calc = Taw_i - Q * (R_g + R_mod)
            if abs(Tg_calc - Tg) < 0.01:
                return hg, Q, Tg_calc, R_w
            Tg += 0.5 * (Tg_calc - Tg)
        raise ModelError("wall temperature iteration failed")

    # ---------------------------------------------------------------- march
    res_rows = []  # (i, hg, hl, hl2, Twg, Twl, Tc, Tc2, Q, flux, props)
    mos_by_station = {}
    vm_by_station = {}
    aborted = False

    initial_guess = t_in + 80
    for outer in range(100):
        t_cool = t_in
        t_ret = initial_guess
        res_rows = []
        aborted = False
        for i in thermal_indices:
            A_gas = garea[i]
            length = ds[i] / cos_t
            A_wall = A_gas
            p_i = pressures[i]
            cp, rho, mu, k = lookup(p_i, t_cool)
            hl = hl_rpe(cp, rho, mu, k, n_down, coolant_velocity[i], i)
            A_cool = fin_area(hl, i, n_down, length)
            cond_down = hl * A_cool
            if two_pass:
                pr = p_ret[i]
                cpr, rhor, mur, kr = lookup(pr, t_ret)
                hl2 = hl_rpe(cpr, rhor, mur, kr, n_return, coolant_velocity_downstream[i], i)
                A_cool2 = fin_area(hl2, i, n_return, length)
                cond_ret = hl2 * A_cool2
            else:
                hl2, A_cool2, cond_ret = 0.0, 0.0, 0.0
            cond = cond_down + cond_ret
            R_l = 1.0 / cond
            T_eff = (cond_down * t_cool + cond_ret * t_ret) / cond
            R_w = wall / (kc * A_wall)
            hg, Q, Twg, _ = wall_solve(i, taw[i], T_eff, R_l, R_w, A_gas)
            Twl = Twg - Q * R_w
            Q1 = hl * A_cool * (Twl - t_cool)
            Q2 = hl2 * A_cool2 * (Twl - t_ret) if two_pass else 0.0
            if two_pass:
                t_ret = t_ret - Q2 / (mdot * lookup(p_ret[i], t_ret)[0])
            t_cool = t_cool + Q1 / (mdot * lookup(p_i, t_cool)[0])
            if not math.isclose(Q1 + Q2, Q, rel_tol=1e-6, abs_tol=1e-6):
                raise ModelError(f"heat balance failed at station {i}")
            flux = Q / A_gas
            res_rows.append((i, hg, hl, hl2 if two_pass else None, Twg, Twl, t_cool, t_ret, Q, flux))

            if not two_pass:
                a_w = widths[i]
                s_tan = tan_coeff * flux + S_mech * (a_w * a_w) * mech_h2
                s_ax = ax_coeff * (Twg - Twl)
                vm = math.sqrt(s_tan * s_tan + s_ax * s_ax - s_tan * s_ax)
                mos_by_station[i] = yield_s / vm - 1
                if margin_floor is not None and mos_by_station[i] < margin_floor:
                    aborted = True
                    break
        if aborted or not two_pass:
            break
        error = res_rows[-1][6] - res_rows[-1][7]
        if abs(error) < 0.01:
            break
        initial_guess += 0.8 * error
    else:
        raise ModelError("two-pass temperature calculation did not converge")

    coolant_outlet = initial_guess
    # Return-only section (two-pass, stations pass_start .. 1)
    if two_pass:
        t_ret = initial_guess
        pressures = list(pressures)
        for i in range(pass_start, 0, -1):
            A_gas = garea[i]
            length = ds[i] / cos_t
            pr = p_ret[i]
            cp, rho, mu, k = lookup(pr, t_ret)
            hl = hl_rpe(cp, rho, mu, k, n_return, coolant_velocity[i], i)
            A_cool = fin_area(hl, i, n_return, length)
            R_l = 1.0 / (hl * A_cool)
            R_w = wall / (kc * A_gas)
            hg, Q, Twg, _ = wall_solve(i, taw[i], t_ret, R_l, R_w, A_gas)
            Twl = Twg - Q * R_w
            t_ret = t_ret + Q / (mdot * cp)
            pressures[i] = pr
            res_rows.append((i, hg, hl, None, Twg, Twl, t_ret, t_ret, Q, Q / A_gas))
        coolant_outlet = t_ret
    elif res_rows and not aborted:
        coolant_outlet = res_rows[-1][6]
    elif res_rows:
        coolant_outlet = res_rows[-1][6]

    # ------------------------------------------------------------ structural
    if two_pass:
        for row in res_rows:
            i, flux, Twg, Twl = row[0], row[9], row[4], row[5]
            a_w = widths[i]
            s_tan = tan_coeff * flux + S_mech * (a_w * a_w) * mech_h2
            s_ax = ax_coeff * (Twg - Twl)
            vm = math.sqrt(s_tan * s_tan + s_ax * s_ax - s_tan * s_ax)
            mos_by_station[i] = yield_s / vm - 1

    # ---------------------------------------------------------------- results
    st = min(mos_by_station, key=mos_by_station.get)
    out = {
        "thrust_N": gas.thrust, "isp_s": gas.isp, "Pc_bar": gas.Pc_bar,
        "cstar_m_s": gas.cstar, "gamma": gas.gamma,
        "min_margin": mos_by_station[st], "min_margin_x_m": xs[st],
        "min_rib_m": min(ribs), "min_width_m": min(widths), "max_width_m": max(widths),
        "min_channel_area_m2": min(area), "mean_channel_area_m2": sum(area) / n,
        "max_Twg_K": max(r[4] for r in res_rows), "max_Twl_K": max(r[5] for r in res_rows),
        "coolant_dT_K": coolant_outlet - t_in,
        "heat_W": sum(r[8] for r in res_rows),
        "dP_bar": (p_in - p_ret[0]) / 1e5,
        "max_coolant_velocity_m_s": max(coolant_velocity),
        "max_flux_W_m2": max(r[9] for r in res_rows),
        "aborted": aborted,
    }
    boils = False
    if not aborted:
        for r_ in res_rows:
            if r_[6] > tab.tsat(pressures[r_[0]]):
                boils = True
                break
    out["coolant_boils"] = boils

    if detail:
        out["_detail"] = dict(
            cfg=cfg, geo=geo, rows=res_rows, pressures=pressures, p_ret=p_ret,
            coolant_velocity=coolant_velocity, thermal_indices=[r[0] for r in res_rows],
            mos=mos_by_station, outlet=coolant_outlet, pass_start=pass_start, tab=tab,
            margin_parts=None,
        )
    return out
