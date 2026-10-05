import csv
import math
from pathlib import Path
from time import perf_counter

start_time = perf_counter()

# =============================================================================
# USER-EDITABLE INPUTS
# =============================================================================
identifier = "Lynx"
workbook_csv = Path.cwd() / f"ThermalOutput_{identifier}.csv"
rpa_file = Path.cwd() / "Lynx_C_RPA_Output.txt"

# Material properties (conductivity k is read from the CSV)
alfa = 2.00e-05          # 1/K, thermal expansion coefficient
poissons = 0.35          # Poisson's ratio
youngsmodulus = 5.00e10  # Pa
Yield_S = 6.00e07        # Pa

B = 0.5                  # plate-bending coefficient

# Column indices in the CSV station table (same as the original script)
CSV_COL_X = 0            # X Position (m)
CSV_COL_WIDTH = 18       # channel width (m)

# Shutdown load case: coolant pressure acting with no hot-gas back-pressure.
#   "inlet" -> Coolant_Inlet_Pressure_Bar from the CSV, applied everywhere
#              (the original script's behaviour)
#   "local" -> local coolant pressure from the RPA table at each station
SHUTDOWN_PRESSURE_MODE = "inlet"

# Heat flux column for the tangential thermal stress: "total" or "conv"
HEAT_FLUX_COLUMN = "total"
# =============================================================================
# END USER-EDITABLE INPUTS
# =============================================================================

RPA_COLUMNS = [
    "x_mm", "radius_mm", "h_conv_kW_m2K", "q_conv_kW_m2", "q_rad_kW_m2",
    "q_total_kW_m2", "Twg_K", "Twi_K", "Twc_K", "Tc_K", "pc_MPa",
    "wc_m_s", "rho_kg_m3",
]


def read_workbook_csv(csv_path):
    """Return (settings dict, sorted list of (x_m, width_m)) from the workbook CSV."""
    settings = {}
    widths = []
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        for row in reader:
            if row and row[0] == "X Position (m)":
                break
            if len(row) >= 2:
                settings[row[0]] = row[1]
        for row in reader:
            if not row or not row[0]:
                continue
            try:
                widths.append((float(row[CSV_COL_X]), float(row[CSV_COL_WIDTH])))
            except (ValueError, IndexError):
                continue
    if not widths:
        raise ValueError(f"No station rows found in {csv_path}")
    return settings, sorted(widths)


def read_rpa_thermal_table(path):
    """Return a list of dicts, one per station, from an RPA thermal output file."""
    stations = []
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tokens = stripped.split()
            if len(tokens) < len(RPA_COLUMNS):
                continue
            try:
                values = [float(t) for t in tokens[:len(RPA_COLUMNS)]]
            except ValueError:
                continue  # title / header / units line
            row = dict(zip(RPA_COLUMNS, values))
            row["comment"] = " ".join(tokens[len(RPA_COLUMNS):])
            stations.append(row)
    if not stations:
        raise ValueError(f"No thermal-analysis rows found in {path}")
    return stations


def interp(x, table):
    """Linear interpolation in a sorted [(x, y), ...] table, clamped at the ends."""
    if x <= table[0][0]:
        return table[0][1]
    if x >= table[-1][0]:
        return table[-1][1]
    for (x0, y0), (x1, y1) in zip(table, table[1:]):
        if x0 <= x <= x1:
            return y0 if x1 == x0 else y0 + (y1 - y0) * (x - x0) / (x1 - x0)


# --- Read inputs --------------------------------------------------------------
settings, width_table = read_workbook_csv(workbook_csv)
h = float(settings["Channel_Wall"])                                # m
k = float(settings["Channel_Conductivity"])                        # W/(m*K)
p_inlet_Pa = float(settings["Coolant_Inlet_Pressure_Bar"]) * 1e5   # Pa

stations = read_rpa_thermal_table(rpa_file)

# --- Per-station calculation -------------------------------------------------
results = []
for s in stations:
    x_m = s["x_mm"] * 1e-3
    a = interp(x_m, width_table)

    q_flux = (s["q_total_kW_m2"] if HEAT_FLUX_COLUMN == "total"
              else s["q_conv_kW_m2"]) * 1e3          # W/m^2
    Twg = s["Twg_K"]
    Twc = s["Twc_K"]                                  # coolant-side wall temp
    delta_t = Twg - Twc

    p_cool = p_inlet_Pa if SHUTDOWN_PRESSURE_MODE == "inlet" else s["pc_MPa"] * 1e6
    q_shutdown = -p_cool

    S_mech = -B * q_shutdown * (a / h) ** 2
    S_tan = youngsmodulus * alfa * q_flux * h / (2 * (1 - poissons) * k)
    S_ax = youngsmodulus * alfa * delta_t

    total_tan = S_tan + S_mech
    S_vm = math.sqrt(total_tan**2 + S_ax**2 - total_tan * S_ax)
    margin = Yield_S / S_vm - 1

    h_implied = k * delta_t / q_flux if q_flux > 0 else float("nan")

    results.append({
        "X Pos (m)": x_m,
        "Radius (m)": s["radius_mm"] * 1e-3,
        "Channel Width (m)": a,
        "Heat Flux (W/m2)": q_flux,
        "Twg (K)": Twg,
        "Twc (K)": Twc,
        "Coolant P (Pa)": s["pc_MPa"] * 1e6,
        "S_ax": S_ax,
        "S_tan": S_tan,
        "S_mech": S_mech,
        "S_vm": S_vm,
        "Margin": margin,
        "h implied (m)": h_implied,
    })

# --- Output ------------------------------------------------------------------
output_filename = Path.cwd() / f"StructuralOutput_{identifier}.csv"
with output_filename.open("w", newline="", encoding="utf-8") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=list(results[0].keys()),
                            quoting=csv.QUOTE_ALL)
    writer.writeheader()
    writer.writerows(results)

worst = min(results, key=lambda r: r["Margin"])
print(f"Inputs:              {workbook_csv.name} + {rpa_file.name} ({len(results)} stations)")
print(f"Wall h / k:          {h * 1e3:.3f} mm / {k:.1f} W/m-K")
print(f"Coolant inlet P:     {p_inlet_Pa / 1e6:.3f} MPa  (shutdown mode: {SHUTDOWN_PRESSURE_MODE})")
print(f"Max von Mises:       {max(r['S_vm'] for r in results) / 1e6:.2f} MPa")
print(f"Lowest margin:       {worst['Margin']:.3f}")
print(f"Lowest margin X pos: {worst['X Pos (m)']:.5f} m")
print(f"Written to:          {output_filename.name}")

# --- Consistency checks ------------------------------------------------------
x_rpa_min, x_rpa_max = stations[0]["x_mm"] * 1e-3, stations[-1]["x_mm"] * 1e-3
x_csv_min, x_csv_max = width_table[0][0], width_table[-1][0]
tol = 0.02 * (x_rpa_max - x_rpa_min)
if x_rpa_min < x_csv_min - tol or x_rpa_max > x_csv_max + tol:
    print(f"WARNING: RPA stations span {x_rpa_min:.4f}-{x_rpa_max:.4f} m but CSV widths "
          f"span {x_csv_min:.4f}-{x_csv_max:.4f} m; widths outside are held constant. "
          f"Check both files use the same x origin.")

heated = [r["h implied (m)"] for r in results if r["Heat Flux (W/m2)"] > 5e5]
if heated:
    h_rpa = sum(heated) / len(heated)
    if abs(h_rpa - h) / h > 0.10:
        print(f"WARNING: RPA wall dT implies h ~ {h_rpa * 1e3:.3f} mm at k = {k}, "
              f"but the CSV gives Channel_Wall = {h * 1e3:.3f} mm.")

print(f"Run time:            {perf_counter() - start_time:.3f} s")
