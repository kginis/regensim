import datetime
import csv
import math
from pathlib import Path
from rocketcea.cea_obj import CEA_Obj, add_new_fuel

#Inputs
Chamber_Pressure_bar = 20
Mass_Ratio = 2.0
Expansion_Ratio = 4.0

#dont touch this :)
Chamber_Pressure_psi = Chamber_Pressure_bar*14.5038
#or this please :3
propanol_card = """
fuel C3H8O(L)  C 3 H 8 O 1  wt%=100.0
h,cal=-76052.0  t(k)=298.15  rho.g/cc=0.803
"""

add_new_fuel("1Propanol", propanol_card)

ispObj = CEA_Obj( oxName='N2O', fuelName='1Propanol')

#Transport Parameters:  heat capacity, viscosity, thermal conductivity, Prandtl number 
chamber_transport = ispObj.get_Chamber_Transport(Pc=Chamber_Pressure_psi, MR=Mass_Ratio)
throat_transport = ispObj.get_Throat_Transport(Pc=Chamber_Pressure_psi, MR=Mass_Ratio)
exit_transport= ispObj.get_Exit_Transport(Pc=Chamber_Pressure_psi, MR=Mass_Ratio, eps=Expansion_Ratio)

#for debug use
#print(chamber_transport)
#print(throat_transport)
#print(exit_transport)

#Gamma: molecular weight, gamma
moluecularweight_gamma = ispObj.get_Chamber_MolWt_gamma(Pc=Chamber_Pressure_psi, MR=Mass_Ratio, eps=Expansion_Ratio)
gamma=moluecularweight_gamma[1]
#print(gamma[1])

#Cstar
Cstar_fps = ispObj.get_Cstar(Pc=Chamber_Pressure_psi, MR=Mass_Ratio)
Cstar_meters_sec = (Cstar_fps*0.3048)
#print(Cstar_meters_sec)

#Temps: Chamber, Nozzle, Exit
temperatures = ispObj.get_Temperatures(Pc=Chamber_Pressure_psi, MR=Mass_Ratio)
print(temperatures[0]*0.555556) #Chamber temp in kelvin
#Isentropic Flow Calculations
time = datetime.datetime.now()
timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
filename = f"IsentropicFlow_{timestamp}.csv"


timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

filename = Path.cwd() / f"IsentropicFlow_{timestamp}.csv"

with filename.open("w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
    mach = 0.01
    for index in range(1, 601):
        mach = index * 0.01

        area_ratio = (
            (1.0 / mach) * ((1 + 0.5 * (gamma - 1.0) * mach**2.0) / (1.0 + 0.5 * (gamma - 1.0)))** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
            )
        pressure_ratio = (
            (1+0.5*(gamma-1)*mach**2)**(-(gamma/(gamma-1)))
        )
        temp_ratio = (
            (1+0.5*(gamma-1)*mach**2)**(-1)
        )
        writer.writerow([mach, area_ratio, pressure_ratio, temp_ratio])

print(f"Isentropic Flow CSV created at: {filename.resolve()}")

CEAout = ispObj.get_full_cea_output(Pc=Chamber_Pressure_psi, MR=Mass_Ratio, eps=Expansion_Ratio)
print(CEAout)