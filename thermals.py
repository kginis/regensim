import math
import csv
import datetime
from pathlib import Path
from pyfluids import Fluid, FluidsList, Input
from rocketcea.cea_obj import add_new_fuel
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.units import add_user_units

#Chamber_Inputs
throat_r_D=1.5
chamber_diameter=0.05 #meters

#CEA Inputs
Chamber_Pressure_bar = 20
Mass_Ratio = 2.0
Expansion_Ratio = 5.5

#Cooling Inputs
Regen_Coolant = FluidsList.Ethanol
Channel_Count = 55
Coolant_Inlet_Temp = 298.15 #kelvin
Coolant_Viscocity = 0.00
Coolant_Inlet_Pressure_Bar = 40 


#conversions
Chamber_Pressure_Pa = Chamber_Pressure_bar*100000
Coolant_Inlet_Pressure_PA = Coolant_Inlet_Pressure_Bar*100000

#IPA properties (dont touch)
propanol_card = """
fuel C3H8O(L)  C 3 H 8 O 1  wt%=100.0
h,cal=-76052.0  t(k)=298.15  rho.g/cc=0.803
"""

add_new_fuel("1Propanol", propanol_card)

add_user_units('millipoise', 'Pa-s', 1e-4)

ispObj = CEA_Obj(
    oxName='N2O',
    fuelName='1Propanol',
    cstar_units='m/s',
    pressure_units='Bar',
    temperature_units='K',
    specific_heat_units='J/kg-K',
    viscosity_units='Pa-s',
)

#Transport Parameters:  heat capacity, viscosity, thermal conductivity, Prandtl number 
chamber_transport = ispObj.get_Chamber_Transport(Pc=Chamber_Pressure_bar, MR=Mass_Ratio)
throat_transport = ispObj.get_Throat_Transport(Pc=Chamber_Pressure_bar, MR=Mass_Ratio)
exit_transport= ispObj.get_Exit_Transport(Pc=Chamber_Pressure_bar, MR=Mass_Ratio, eps=Expansion_Ratio)

gas_Cp=chamber_transport[0]
gas_viscocity=chamber_transport[1]
gas_prandtl=chamber_transport[3]

#Gamma: molecular weight, gamma
moluecularweight_gamma = ispObj.get_Chamber_MolWt_gamma(Pc=Chamber_Pressure_bar, MR=Mass_Ratio, eps=Expansion_Ratio)
gamma=moluecularweight_gamma[1]
#print(gamma[1])

#Cstar
Cstar_meters_sec = ispObj.get_Cstar(Pc=Chamber_Pressure_bar, MR=Mass_Ratio)
#print(Cstar_meters_sec)

#Temps: Chamber, Nozzle, Exit
temperatures = ispObj.get_Temperatures(Pc=Chamber_Pressure_bar, MR=Mass_Ratio)
chamber_temp = temperatures[0] #in kelvin

#Isentropic Flow Calculations
time = datetime.datetime.now()
timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
filename = f"IsentropicFlow_{timestamp}.csv"


timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

filename = Path.cwd() / f"IsentropicFlow_{timestamp}.csv"
#Saving it to a sheet
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

#CEAout = ispObj.get_full_cea_output(Pc=Chamber_Pressure_bar, MR=Mass_Ratio, eps=Expansion_Ratio)
#print(CEAout)

def find_throat(list):
    minimum = None
    for row in nozzlegeoemtry:
            value=(row[1])
            if minimum == None or value<minimum:
                minimum=value
    return minimum    

def coolant_rho(pressure, temperature):
    coolant = Fluid(Regen_Coolant).with_state(
        Input.pressure(pressure),     
        Input.temperature(temperature-273.15),   
    )

    return coolant.density

def coolant_velocity(pressure, temperature):
    coolant = Fluid(Regen_Coolant).with_state(
        Input.pressure(pressure),     
        Input.temperature(temperature-273.15),   
    )

    return coolant.dynamic_viscosity

def coolant_specific_heat(pressure, temperature):
    coolant = Fluid(Regen_Coolant).with_state(
        Input.pressure(pressure),     
        Input.temperature(temperature-273.15),   
    )

    return coolant.specific_heat

def coolant_specific_heat(pressure, temperature):
    coolant = Fluid(Regen_Coolant).with_state(
        Input.pressure(pressure),     
        Input.temperature(temperature-273.15),   
    )

    return coolant.conductivity

def hg_bartz(radius,Cp,mu,Pr,Pc,C_star,Area_Ratio,throat_curvature_radius): #note that sigma is not yet applied
    D_star=float(radius*2)

    bartz=(
        (0.026/D_star**0.2)*(((mu**0.2)*Cp)/Pr**0.6)*((Pc/C_star)**0.8)*((D_star/throat_curvature_radius)**0.1)*(Area_Ratio)**-0.9
        )

    return bartz
    
def hl_RPE(c_cp,c_mdot,c_rho,c_mu,c_hyddiameter,thermalconductivity):
    print('hi :3')

def hl_hezel_huang(c_cp,c_mdot,c_rho,c_mu,c_hyddiameter,thermalconductivity):
    print('hi :3')
with open('nozzle.csv', 'r') as nozzle:

    csv_reader = csv.reader(nozzle)
    nozzlegeoemtry = list(csv_reader)
    
    throat_radius=float((find_throat(nozzlegeoemtry)))
    throat_area=((throat_radius**2)*math.pi)
    #print("Throat Area:") 
    #print(Athroat)

#Isentropic Flow lookup table
x_positions=[]
nozzle_radii=[]
area_ratios=[]
mach_number=[]
gas_temp=[]
adiabatic_wall_temp=[]

for radius in nozzlegeoemtry:
    with filename.open("r", newline="", encoding="utf-8") as isenflow:

        isentropicflow = csv.reader(isenflow)
        isentropicflowlookup = list(isentropicflow)

        try:
            x_pos=float(radius[0])
            area=(((float(radius[1]))**2)*math.pi)
            local_area_ratio = (area / throat_area)

            lookuprange=[]

            for mach in isentropicflowlookup:

                if x_pos<=0.0:
                    lookupcondition = float(mach[0]) < 1
                else:
                    lookupcondition = float(mach[0]) >= 1

                if lookupcondition:
                    lookuprange.append(mach)
                    #print(mach)
                    

            #print(lookuprange[1])
            #exit(000)
            #print(local_area_ratio)
            closest_row = min(
                lookuprange,
            key=lambda row: abs(float(row[1]) - local_area_ratio),
            )
            Mach = float(closest_row[0])
            Pres_Ratio_Flow = float(closest_row[2])
            Temp_Ratio_Flow = float(closest_row[3])

            Gas_Temp=(chamber_temp*Temp_Ratio_Flow)
            Pressure=(Chamber_Pressure_bar*Pres_Ratio_Flow)

            r = gas_prandtl ** (1 / 3)

            Adiabatic_Wall = Gas_Temp * (
                1
                + r * ((gamma - 1) / 2) * Mach**2
                )
            x_positions.append(float(x_pos))

            try:
                nozzle_radii.append(float(radius[1]))
            except ValueError:
                continue

            area_ratios.append(float(local_area_ratio))
            mach_number.append(float(Mach))
            gas_temp.append(float(Gas_Temp))
            adiabatic_wall_temp.append(float(Adiabatic_Wall))
        except ValueError:
            continue

hl=[]
hg=[]
# Main Convergent solver
for i in range(len(area_ratios)):
    hg.append(hg_bartz(throat_radius, 
                   gas_Cp,
                   gas_viscocity,
                   gas_prandtl,
                   Chamber_Pressure_Pa,
                   Cstar_meters_sec,
                   area_ratios[i],
                   throat_r_D*chamber_diameter))
