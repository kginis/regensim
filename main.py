import math
import csv
import time
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

#Cooling Settings
Two_Pass = False #flowing from top to bottom, to top again

#Cooling Inputs
Regen_Coolant = FluidsList.Ethanol
Surface_Roughness = 0.000025 #25 Ra
Coolant_Mdot = 0.35 #kg/s
Channel_Width = 0.001 #note this is in meters
Channel_Height = 0.001
Channel_Count = 55.0
Channel_Wall = 0.001 #1mm
Channel_Conductivity = 135 
Coolant_Inlet_Temp = 298.15 #kelvin
Coolant_Inlet_Pressure_Bar = 40 

#Cooling Values
Hydraulic_Diameter=4*(Channel_Width*Channel_Height)/(2*Channel_Width+2*Channel_Height)
Channel_Area=(Channel_Width*Channel_Height)

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

def coolant_abs_viscocity(pressure, temperature):
    coolant = Fluid(Regen_Coolant).with_state(
        Input.pressure(pressure),     
        Input.temperature(temperature-273.15),   
    )
    return coolant.dynamic_viscosity

def coolant_kin_viscocity(pressure, temperature):
    coolant = Fluid(Regen_Coolant).with_state(
        Input.pressure(pressure),     
        Input.temperature(temperature-273.15),   
    )
    return coolant.kinematic_viscosity
    

def coolant_conductivity(pressure, temperature):
    coolant = Fluid(Regen_Coolant).with_state(
        Input.pressure(pressure),     
        Input.temperature(temperature-273.15),   #lets just say there are more efficient ways to do this but all code should have charachter
    )

    return coolant.conductivity

def coolant_specific_heat(pressure, temperature):
    coolant = Fluid(Regen_Coolant).with_state(
        Input.pressure(pressure),     
        Input.temperature(temperature-273.15),   
    )

    return coolant.specific_heat

def coolant_connductivity(pressure, temperature):
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

def hl_RPE(c_cp,c_mdot,c_rho,c_mu,c_conductivity, channelqty, c_velocity):

    hl = 0.023*c_cp*(c_mdot/(Channel_Area*channelqty))*((Hydraulic_Diameter*c_velocity*c_rho)/c_mu)**-0.2*((c_mu*c_cp)/c_conductivity)**(-2/3)

    return hl
#def hl_hezel_huang(coolant):
    #maybe later

def deltaP(f,L,V,rho):
    deltaP=(
          (f*L*(V)**2*rho) /
          (2*Hydraulic_Diameter)
     )
    return deltaP

def bartz_boundary_sigma(Tw, Tc, gamma, Mach, omega):
    M_gamma = 1 + ((gamma - 1) / 2) * Mach**2

    sigma = 1 / (
        (0.5 * (Tw / Tc) * M_gamma + 0.5)**(0.8 - omega / 5)
        * M_gamma**(omega / 5)
    )

    return sigma

with open('nozzle.csv', 'r') as nozzle:

    csv_reader = csv.reader(nozzle)
    nozzlegeoemtry = list(csv_reader)
    
    throat_radius=float((find_throat(nozzlegeoemtry)))
    throat_area=((throat_radius**2)*math.pi)
    #print("Throat Area:") 
    #print(Athroat)

#Isentropic Flow lookup table
x_positions=[]
chamber_radii=[]
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
                chamber_radii.append(float(radius[1]))
            except ValueError:
                continue

            area_ratios.append(float(local_area_ratio))
            mach_number.append(float(Mach))
            gas_temp.append(float(Gas_Temp))
            adiabatic_wall_temp.append(float(Adiabatic_Wall))
        except ValueError:
            continue
#coolant velocity math
#assuming constant density for this, b/c lazy
coolant_pressures=[]
#friction factor
dev=1.5
hydraulicpermiter=2*Channel_Height+2*Channel_Width #maybe i will refractor the variable names to be less cooked
hydraulicradius=(Channel_Width*Channel_Height)/(hydraulicpermiter)

Reynolds = (
    4 * (Coolant_Mdot / Channel_Count)
    / (
        coolant_kin_viscocity(
            Coolant_Inlet_Pressure_PA,
            Coolant_Inlet_Temp,
        )
        * hydraulicpermiter
        * coolant_rho(
            Coolant_Inlet_Pressure_PA,
            Coolant_Inlet_Temp,
        )
    )
)

f = 0.02
dev = float("inf")

while dev >= 1e-6:
    log_argument = (
        Surface_Roughness / (14.8 * hydraulicradius)
        + 2.51 / (Reynolds * math.sqrt(f))
    )

    inverse_sqrt_f = -2 * math.log10(log_argument)
    f_new = 1 / inverse_sqrt_f**2

    dev = abs(f_new - f)
    f = f_new

#print("Friction factor:", f)

coolant_velocity=Coolant_Mdot/(coolant_rho(Coolant_Inlet_Pressure_PA,Coolant_Inlet_Temp)*Channel_Area*Channel_Count)  #m/s
first_pass_dist=0
second_pass_dist=0
j=0
#print("Coolant Velocity:", coolant_velocity)

if Two_Pass:
    #print(inletpos)
    for i in range(len(area_ratios)):
        first_pass_dist=first_pass_dist+abs(x_positions[i]-x_positions[j])
        #print(first_pass_dist)
        coolant_pressures.append(
            Coolant_Inlet_Pressure_PA
            - deltaP(f,first_pass_dist,coolant_velocity,coolant_rho(Coolant_Inlet_Pressure_PA,Coolant_Inlet_Temp))
        )
        j=i

for i in range(len(area_ratios) - 1, -1, -1):
    second_pass_dist=second_pass_dist+abs(x_positions[i]-x_positions[j])
    #print(second_pass_dist+first_pass_dist)
    coolant_pressures.append(
                Coolant_Inlet_Pressure_PA
                - deltaP(f,first_pass_dist+second_pass_dist,coolant_velocity,coolant_rho(Coolant_Inlet_Pressure_PA,Coolant_Inlet_Temp))
            )
    j=i

hl=[]
hg=[]
Twg_list=[]
Twl_list=[]
Tc_list=[]
Q_list=[]
q_heatflux=[] #heat flux

coolant_temp = Coolant_Inlet_Temp

for i in range(len(area_ratios) - 1, 0, -1):
    dx = abs(x_positions[i] - x_positions[i - 1])
    dr = chamber_radii[i] - chamber_radii[i - 1]
    ds = math.sqrt(dx**2 + dr**2)   

    gas_area = (
        math.pi
        * (chamber_radii[i] + chamber_radii[i - 1])
        * ds
    )

    coolant_area = (
        hydraulicpermiter
        * Channel_Count
        * ds
    )

    wall_area = gas_area  # thin-wall approximation

    Twg_guess = 0.5 * (
        adiabatic_wall_temp[i] + coolant_temp
    )

    tolerance = 0.01       # K
    relaxation = 0.5
    max_iterations = 150

    for iteration in range(max_iterations):
        hg_local = (
            hg_bartz(
                throat_radius,
                gas_Cp,
                gas_viscocity,
                gas_prandtl,
                Chamber_Pressure_Pa,
                Cstar_meters_sec,
                area_ratios[i],
                throat_r_D * (2*throat_radius)
            )
            * bartz_boundary_sigma(
                Twg_guess,
                chamber_temp,
                gamma,
                mach_number[i],
                0.6,
            )
        )

        hl_local = hl_RPE(
            coolant_specific_heat(
                coolant_pressures[i], coolant_temp
            ),
            Coolant_Mdot,
            coolant_rho(
                coolant_pressures[i], coolant_temp
            ),
            coolant_abs_viscocity(
                coolant_pressures[i], coolant_temp
            ),
            coolant_conductivity(
                coolant_pressures[i], coolant_temp
            ),
            Channel_Count,
            coolant_velocity,
        )

        R_g = 1.0 / (hg_local * gas_area)
        R_w = Channel_Wall / (
            Channel_Conductivity * wall_area
        )
        R_l = 1.0 / (hl_local * coolant_area)

        Q = (
            adiabatic_wall_temp[i] - coolant_temp
        ) / (R_g + R_w + R_l)

        Twg_calculated = (
            adiabatic_wall_temp[i] - Q * R_g
        )

        dev = abs(Twg_calculated - Twg_guess)

        if dev < tolerance:
            Twg = Twg_calculated
            break

        Twg_guess += relaxation * (
            Twg_calculated - Twg_guess
        )
    else:
        raise RuntimeError(
            f"Wall-temperature iteration failed at index {i}"
        )

    Twl = Twg - Q * R_w

    coolant_temp=coolant_temp+Q/(Coolant_Mdot*coolant_specific_heat(coolant_pressures[i],coolant_temp))

    hg.append(hg_local)
    hl.append(hl_local)
    Twg_list.append(Twg_calculated)
    Twl_list.append(Twl)
    Tc_list.append(coolant_temp)
    Q_list.append(Q)
    q_heatflux.append(Q/gas_area)

    
    coolant_temperature_check = Twl - Q * R_l

savefilename = Path.cwd() / f"ThermalOutput_{timestamp}.csv"

with savefilename.open("w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile, quoting=csv.QUOTE_ALL)

    # Input settings
    writer.writerow(["Input", "Value", "Unit"])
    writer.writerow(["throat_r_D", throat_r_D, "-"])
    writer.writerow(["chamber_diameter", chamber_diameter, "m"])
    writer.writerow(["Chamber_Pressure_bar", Chamber_Pressure_bar, "bar"])
    writer.writerow(["Mass_Ratio", Mass_Ratio, "-"])
    writer.writerow(["Expansion_Ratio", Expansion_Ratio, "-"])
    writer.writerow(["Two_Pass", Two_Pass, "bool"])
    writer.writerow(["Regen_Coolant", str(Regen_Coolant), "-"])
    writer.writerow(["Surface_Roughness", Surface_Roughness, "m"])
    writer.writerow(["Coolant_Mdot", Coolant_Mdot, "kg/s"])
    writer.writerow(["Channel_Width", Channel_Width, "m"])
    writer.writerow(["Channel_Height", Channel_Height, "m"])
    writer.writerow(["Channel_Count", Channel_Count, "-"])
    writer.writerow(["Channel_Wall", Channel_Wall, "m"])
    writer.writerow([
        "Channel_Conductivity",
        Channel_Conductivity,
        "W/(m·K)",
    ])
    writer.writerow(["Coolant_Inlet_Temp", Coolant_Inlet_Temp, "K"])
    writer.writerow([
        "Coolant_Inlet_Pressure_Bar",
        Coolant_Inlet_Pressure_Bar,
        "bar",
    ])

    # Blank rows separating inputs from results
    writer.writerow([])
    writer.writerow([])

    # Thermal-output headings
    writer.writerow([
        "X Pos (m)",
        "Radius (m)",
        "Area Ratio",
        "Mach Number",
        "Gas Temp (K)",
        "Adiabatic Wall Temp (K)",
        "",
        "Coolant Velocity (m/s)",
        "Coolant Pressure (Pa)",
        "",
        "Coolant h (W/m^2-K)",
        "Gas h (W/m^2-K)",
        "Gas-Side Wall Temp (K)",
        "Coolant-Side Wall Temp (K)",
        "Coolant Temp (K)",
        "",
        "Heat Transferred (W)",
        "Heat Flux (W/m^2)",
    ])

    for i in range(len(Twg_list)):
        writer.writerow([
            x_positions[i],
            chamber_radii[i],
            area_ratios[i],
            mach_number[i],
            gas_temp[i],
            adiabatic_wall_temp[i],
            "",
            coolant_velocity,
            coolant_pressures[i],
            "",
            hl[i],
            hg[i],
            Twg_list[i],
            Twl_list[i],
            Tc_list[i],
            "",
            Q_list[i],
            q_heatflux[i],
            ])

print(f"Output CSV created at: {savefilename.resolve()}")

print()
print()
print("Maximum Adiabatic Wall Temp:", max(adiabatic_wall_temp))
print("Maximum Gas Temp:", max(gas_temp))
print()
print("Maximum Coolant Temp:", max(Tc_list))
print()
print("Maximum Coolant-Side Wall Temp:", max(Twl_list))
print("Maximum Gas-Side Wall Temp:", max(Twg_list))
print()
print("Maximum Heat Flux:", max(q_heatflux))
