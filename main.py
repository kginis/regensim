"""regensim entry point.

    python main.py [inputs.toml]            single run, full CSV output
    python main.py [inputs.toml] --batch    sweep [sweep] in the input file, pick the optimum
"""
import argparse
import csv
import math
from pathlib import Path
from time import perf_counter

import engine as E


def print_range(label, values, unit):
    if len(values):
        print(f"{label}: {min(values):.6g} to {max(values):.6g} {unit}")
    else:
        print(f"{label}: N/A")


def single_run(cfg, out_dir):
    start = perf_counter()
    cfg = E.finalize(cfg)
    ident = cfg["identifier"]
    nozzle = E.load_nozzle(Path(cfg["_dir"]) / cfg["nozzle_file"])
    tab = E.CoolantTable(cfg["Regen_Coolant"], cfg["Coolant_Inlet_Pressure_Bar"], cache_dir=cfg["_dir"])
    gas = E.get_gas(cfg, nozzle)
    out = E.simulate(cfg, gas, tab, detail=True)
    d = out["_detail"]
    cfg = d["cfg"]
    geo, rows, pressures = d["geo"], d["rows"], d["pressures"]
    two_pass = cfg["Two_Pass"]
    film_on = cfg["Film_Cooling"]

    # per-station coolant properties for the output file
    props = {}
    for r in rows:
        i, t_c = r[0], r[6]
        cp, rho, mu, k = tab.lookup(pressures[i], t_c)
        props[i] = (cp, rho, mu, mu / rho, k)

    # ---- ThermalOutput
    thermal_path = out_dir / f"ThermalOutput_{ident}.csv"
    with thermal_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, quoting=csv.QUOTE_ALL)
        w.writerow(["Input", "Value", "Unit"])
        for label, key, unit in [
            ("throat_r_D", "throat_r_D", "-"), ("chamber_diameter", "chamber_diameter", "m"),
        ]:
            w.writerow([label, cfg[key], unit])
        w.writerow(["Chamber_Pressure_bar", gas.Pc_bar, "bar"])
        w.writerow(["Thrust", gas.thrust, "N"])
        for label, key, unit in [
            ("Mass_Ratio", "Mass_Ratio", "-"), ("Expansion_Ratio", "Expansion_Ratio", "-"),
            ("Two_Pass", "Two_Pass", "bool"), ("Film_Cooling", "Film_Cooling", "bool"),
            ("Regen_Coolant", "Regen_Coolant", "-"), ("Film_Coolant", "Film_Coolant", "-"),
            ("Surface_Roughness", "Surface_Roughness", "m"), ("Coolant_Mdot", "Coolant_Mdot", "kg/s"),
            ("Channel_Width", "Channel_Width", "m"), ("Channel_Height", "Channel_Height", "m"),
            ("Channel_Count", "Channel_Count", "-"), ("Channel_Wall", "Channel_Wall", "m"),
            ("Channel_Conductivity", "Channel_Conductivity", "W/(m·K)"),
            ("Coolant_Inlet_Temp", "Coolant_Inlet_Temp", "K"),
            ("Coolant_Inlet_Pressure_Bar", "Coolant_Inlet_Pressure_Bar", "bar"),
            ("total_mdot", "total_mdot", "kg/s"),
        ]:
            w.writerow([label, cfg[key], unit])
        w.writerow(["Film_Mdot (active)", cfg["Film_Mdot"] if film_on else 0.0, "kg/s"])
        w.writerow(["Film_Inlet_Temp", cfg["Film_Inlet_Temp"] if film_on else "", "K"])
        w.writerow(["Film / total mass flow", cfg["Film_Mdot"] / cfg["total_mdot"] if film_on else 0.0, "-"])
        w.writerow([])
        w.writerow(["Film summary", "Value", "Unit"])
        w.writerow(["Film cooling status", "Enabled" if film_on else "Disabled", "-"])
        if film_on:
            for label, (value, unit) in gas.film_summary.items():
                w.writerow([label, value, unit])
            w.writerow(["Gas-film Cp basis", "Saturated-liquid Cp (current model)", "-"])
        w.writerow([])
        w.writerow([])
        w.writerow(
            ["X Position (m)", "Chamber Radius (m)", "Area Ratio", "Mach Number", "Gas Temperature (K)",
             "Gas Velocity (m/s)", "Adiabatic Wall Temperature (K)", "Gas Pressure (Pa)", "",
             "Coolant Velocity (m/s)", "Coolant Pressure (Pa)", "Coolant Specific Heat (J/kg-K)",
             "Coolant Density (kg/m^3)", "Coolant Dynamic Viscosity (Pa-s)",
             "Coolant Kinematic Viscosity (m^2/s)", "Coolant Thermal Conductivity (W/m-K)", "",
             "Channel Rib Width (m)", "Channel Width (m)", "Coolant-Side Heat Transfer Coefficient (W/m^2-K)"]
            + (["Coolant-Side Heat Transfer Coefficient, Pass 2 (W/m^2-K)"] if two_pass else [])
            + ["Gas-Side Heat Transfer Coefficient (W/m^2-K)", "Gas-Side Wall Temperature (K)",
               "Coolant-Side Wall Temperature (K)", "Coolant Temperature (K)",
               "Coolant Temprature, Pass 2 (K)", "", "Heat Transferred (W)", "Heat Flux (W/m^2)", "",
               "Film Coolant State", "Distance from Film Injection (m)",
               "Uncooled Adiabatic Wall Temperature (K)", "Film Taw Reduction (K)",
               "Film Coolant Temperature Used (K)", "Film Coolant Cp Used (J/kg-K)",
               "Film Mixing Fraction eta (-)", "Film Mixture Cp Used (J/kg-K)",
               "Film Coolant Sensible Enthalpy Used (J/kg)", "Film Gas Sensible Enthalpy Used (J/kg)",
               "Film Flow Multiplier phi (-)", "Film Equivalent Distance x_bar (m)",
               "Film Entrainment Ratio (-)", "Film Shape Factor (-)",
               "Film Reference Vapor Gamma Used (-)"])
        for r in sorted(rows, key=lambda r: r[0]):
            i, hg, hl, hl2, Twg, Twl, Tc, Tc2, Q, flux = r
            fd = gas.film_data[i]
            cp, rho, mu, nu, k = props[i]
            w.writerow(
                [gas.x[i], gas.r[i], gas.area_ratio[i], gas.mach[i], gas.gas_temp[i], gas.gas_velocity[i],
                 gas.taw[i], gas.gas_pressure[i], "", d["coolant_velocity"][i], pressures[i],
                 cp, rho, mu, nu, k, "", geo["ribs"][i], geo["widths"][i], hl]
                + ([hl2 if hl2 is not None else ""] if two_pass else [])
                + [hg, Twg, Twl, Tc, Tc2, "", Q, flux, "", gas.film_state[i],
                   abs(gas.x[i] - gas.x[0]) if film_on else "", gas.taw_uncooled[i],
                   gas.taw_uncooled[i] - gas.taw[i] if film_on else "",
                   fd.get("film_temperature", ""), fd.get("cp_used", ""), fd.get("eta", ""),
                   fd.get("cp_mix", ""), fd.get("coolant_sensible_enthalpy", ""),
                   fd.get("gas_sensible_enthalpy", ""), fd.get("phi", ""), fd.get("x_bar", ""),
                   fd.get("entrainment_ratio", ""), fd.get("shape_factor", ""), fd.get("gamma_used", "")])

    # ---- StructuralOutput
    E_, al, nu_p, kc = cfg["youngsmodulus"], cfg["alfa"], cfg["poissons"], cfg["Channel_Conductivity"]
    h = cfg["Channel_Wall"]
    struct_path = out_dir / f"StructuralOutput_{ident}.csv"
    with struct_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, quoting=csv.QUOTE_ALL)
        w.writerow(["X Pos", "S_ax", "S_tan", "S_mech", "S_vm", "Margin"])
        for r in sorted(rows, key=lambda r: r[0]):
            i, Twg, Twl, flux = r[0], r[4], r[5], r[9]
            s_mech = cfg["shutdown_B"] * cfg["Coolant_Inlet_Pressure_Pa"] * (geo["widths"][i] / h) ** 2
            s_tan = E_ * al * flux * h / (2 * (1 - nu_p) * kc)
            s_ax = E_ * al * (Twg - Twl)
            tot = s_tan + s_mech
            vm = math.sqrt(tot**2 + s_ax**2 - tot * s_ax)
            w.writerow([gas.x[i], s_ax, s_tan, s_mech, vm, cfg["Yield_S"] / vm - 1])

    # ---- console summary
    Twg_l = [r[4] for r in rows]
    Twl_l = [r[5] for r in rows]
    flux_l = [r[9] for r in rows]
    qsum = sum(r[8] for r in rows)
    print(f"Output CSV created at: {thermal_path.resolve()}")
    print(f"Structural CSV created at: {struct_path.resolve()}")
    print("\nPerformance Results")
    print("Chamber Pressure:", gas.Pc_bar, "bar")
    print("Thrust:", gas.thrust, "N")
    print("Isp:", gas.isp, "s")
    print("Gamma:", gas.gamma)
    print("C star:", gas.cstar)
    print("Exit Pressure:", gas.exit_pressure, "Pa")
    print("\nThermal results")
    print_range("Effective adiabatic wall temperature", list(gas.taw), "K")
    print_range("Uncooled adiabatic wall temperature", list(gas.taw_uncooled), "K")
    print_range("Gas temperature", list(gas.gas_temp), "K")
    print_range("Gas velocity", list(gas.gas_velocity), "m/s")
    print_range("Regen coolant temperature (all modeled passes)", [r[6] for r in rows] + [r[7] for r in rows], "K")
    print(f"Regen Coolant Delta T: {out['coolant_dT_K']:.3f} K")
    print(f"Total heat absorbed: {qsum:.1f} W")
    print_range("Regen coolant velocity max:", d["coolant_velocity"], "m/s")
    print_range("Regen coolant channel widths", geo["widths"], "m")
    print_range("Regen coolant channel ribs", geo["ribs"], "m")
    print_range("Coolant-side wall temperature", Twl_l, "K")
    print("Coolant pressure drop", out["dP_bar"], "Bar")
    print_range("Gas-side wall temperature", Twg_l, "K")
    print_range("Heat flux", flux_l, "W/m^2")
    k_hot = max(range(len(Twg_l)), key=lambda k: Twg_l[k])
    i_hot = rows[k_hot][0]
    print("\nHottest segment:", i_hot - 1, "->", i_hot)
    print("Maximum wall temperature:", Twg_l[k_hot], "K")
    print("Cooling region:", "return only" if two_pass and i_hot <= d["pass_start"]
          else "both passes" if two_pass else "single pass")
    if out["coolant_boils"]:
        print("WARNING: bulk coolant exceeds saturation temperature somewhere; "
              "liquid property table is not valid there.")
    print("\nStructural results")
    print(f"Lowest Margin: {out['min_margin']:.4f} at x = {out['min_margin_x_m']:.5f} m")
    print("\nFilm cooling:", "Enabled" if film_on else "Disabled")
    if film_on:
        print(f"Film mass flow: {cfg['Film_Mdot']:.6g} kg/s")
        for label, (value, unit) in gas.film_summary.items():
            print(f"{label}: {value:.6g} {unit}")
    print(f"\nRun completed in {perf_counter() - start:.2f} seconds")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="?", default="inputs.toml")
    ap.add_argument("--batch", action="store_true", help="run the [sweep] and pick the optimum")
    ap.add_argument("--workers", type=int, help="worker processes (default: all cores)")
    ap.add_argument("--max-runs", type=int, help="cap number of runs (quick test)")
    ap.add_argument("--out", type=Path, help="output directory")
    a = ap.parse_args()

    cfg, raw = E.load_inputs(a.inputs)
    if a.batch:
        import batch
        batch.run_batch(cfg, raw, workers=a.workers, out_dir=a.out, max_runs=a.max_runs)
    else:
        out_dir = a.out or Path.cwd()
        out_dir.mkdir(exist_ok=True)
        single_run(cfg, out_dir)


if __name__ == "__main__":
    main()
