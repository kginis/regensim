import math
import csv
import time
from time import perf_counter
from pathlib import Path
import datetime

start_time = perf_counter()
# =============================================================================
# USER-EDITABLE INPUTS
# =============================================================================
identifier = 'Lynx'

#material properties
alfa=2.00E-05	#1/K
poissons=0.35 #poissons ratio	
youngsmodulus=5.00E+10	#Pa
Yield_S=6.00E+07  #Pa
k = 160.0 #W/M*K

# =============================================================================
# END USER-EDITABLE INPUTS
# =============================================================================

filename = Path.cwd() / f"ThermalOutput_{identifier}.csv"

import csv


def iter_workbook_inputs(csv_path):
    with open(csv_path, newline="", encoding="utf-8-sig") as file:
        reader = csv.reader(file)
        settings = {}

        for row in reader:
            if row and row[0] == "X Position (m)":
                break
            if len(row) >= 2:
                settings[row[0]] = row[1]

        B = 0.5
        h = float(settings["Channel_Wall"])  # m
        a_thermal = float(settings["Channel_Conductivity"])  # W/(m·K)
        q_shutdown = -float(settings["Coolant_Inlet_Pressure_Bar"]) * 100_000  # Pa

        for row in reader:
            if not row or not row[0]:
                continue

            x = float(row[0])                 # station position, m
            a = float(row[18])                # channel width, m
            q = float(row[7]) - float(row[10])  # gas - coolant pressure, Pa
            q_thermal = float(row[27])        # heat flux, W/m²
            Twg = float(row[22])              # K
            Twl = float(row[23])              # K
            delta_t = Twg - Twl           # K

            yield (
                x, B, a, h, q, q_shutdown,
                q_thermal, a_thermal, Twg, Twl, delta_t
            )

x_pos=[]

vm_stress=[]
thermal_axial=[]
thermal_tan=[]
mech_shutdown=[]
mos=[]

output_filename = Path.cwd() / f"StructuralOutput_{identifier}.csv"

for x, B, a, h, q, q_shutdown, q_thermal, a_thermal, Twg, Twl, delta_t in (
    iter_workbook_inputs("ThermalOutput_Lynx.csv")):
    local_mech_stress_shutdown=(0.5*youngsmodulus*a*2)/(h**2)
    local_mech_stress_shutdown = -B * q_shutdown * (a / h)**2

    local_thermal_tan = (
        youngsmodulus * alfa * q_thermal * h
        / (2 * (1 - poissons) * k)
    )
    local_thermal_axial = youngsmodulus * alfa * delta_t

    total_tan = local_thermal_tan + local_mech_stress_shutdown
    local_vm_stress = math.sqrt(
        total_tan**2
        + local_thermal_axial**2
        - total_tan * local_thermal_axial
    )
    local_mos = Yield_S / local_vm_stress - 1
    x_pos.append(x)
    vm_stress.append(local_vm_stress)
    thermal_axial.append(local_thermal_axial)
    thermal_tan.append(local_thermal_tan)
    mech_shutdown.append(local_mech_stress_shutdown)
    mos.append(local_mos)

savefilename = Path.cwd() / f"StructuralOutput_{identifier}.csv"

with savefilename.open("w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
    writer.writerow(["X Pos","S_ax", "S_tan", "S_mech", "S_vm", "Margin"])
    for i in range(len(vm_stress)):
        writer.writerow([x_pos[i],thermal_axial[i],thermal_tan[i],mech_shutdown[i],vm_stress[i],mos[i]])

print('Max Stress', max(vm_stress))
lowest_margin_station = min(range(len(mos)), key=lambda i: mos[i])
print("Lowest Margin:", mos[lowest_margin_station])
print("Lowest Margin X pos:", x_pos[lowest_margin_station])