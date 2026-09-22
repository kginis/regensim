import math
import csv
import time
import datetime
from pathlib import Path
from pyfluids import Fluid, FluidsList, Input
from rocketcea.cea_obj import add_new_fuel
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.units import add_user_units

# =============================================================================
# USER-EDITABLE INPUTS
# =============================================================================
identifier='Lynx'

#Chamber_Inputs
throat_r_D=1.5
chamber_diameter=0.03*2 #meters
total_mdot = 1.75

#CEA Inputs
Chamber_Pressure_bar = 25 #this is a 'guess'
Mass_Ratio = 1.75
Expansion_Ratio = 5.444
Cstar_efficiency = 0.95 #80% Cstar

#Cooling Settings
Two_Pass = False #flowing from top to bottom, to top again
coolant_parameter=0.35 #fraction of total hydraulic perimiter that is being cooled. conservative estimate = 0.25. one of many parameters to tune
Film_Cooling=True 

#Cooling Inputs
Regen_Coolant = FluidsList.Ethanol
Film_Coolant = FluidsList.Ethanol
Surface_Roughness = 0.000025 #25 Ra
Coolant_Mdot = total_mdot*1/(1+Mass_Ratio) #kg/s
Film_Mdot = Coolant_Mdot*0.15
Film_Inlet_Temp = 340
Channel_Width = 0.0015 #meters
Channel_Height = 0.0015
Channel_Count = 40.0
Channel_Wall = 0.001 #1mm
Channel_Conductivity = 130
Coolant_Inlet_Temp = 298.15 #kelvin
Coolant_Inlet_Pressure_Bar = 50

# =============================================================================
# END USER-EDITABLE INPUTS
# =============================================================================

#Cooling Values
Hydraulic_Diameter=4*(Channel_Width*Channel_Height)/(2*Channel_Width+2*Channel_Height)
Channel_Area=(Channel_Width*Channel_Height)

#conversions
Chamber_Pressure_Pa = Chamber_Pressure_bar*100000
Coolant_Inlet_Pressure_PA = Coolant_Inlet_Pressure_Bar*100000

if Two_Pass:
    Channel_Count=Channel_Count/2

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
    density_units="kg/m^3",
    enthalpy_units="J/kg",
)


def find_throat(list):
    minimum = None
    for row in nozzlegeoemtry:
            value=(row[1])
            if minimum == None or value<minimum:
                minimum=value
    return minimum  

with open('nozzle.csv', 'r') as nozzle:

    csv_reader = csv.reader(nozzle)
    nozzlegeoemtry = list(csv_reader)
    
    throat_radius=float((find_throat(nozzlegeoemtry)))
    throat_area=((throat_radius**2)*math.pi)
    #print("Throat Area:") 
    #print(Athroat)

dev = float("inf")
while dev > 0.01:

    #Transport Parameters:  heat capacity, viscosity, thermal conductivity, Prandtl number 
    chamber_transport = ispObj.get_Chamber_Transport(Pc=Chamber_Pressure_bar, MR=Mass_Ratio, frozen=1)
    throat_transport = ispObj.get_Throat_Transport(Pc=Chamber_Pressure_bar, MR=Mass_Ratio)
    exit_transport= ispObj.get_Exit_Transport(Pc=Chamber_Pressure_bar, MR=Mass_Ratio, eps=Expansion_Ratio)

    gas_Cp=chamber_transport[0]
    gas_viscocity=chamber_transport[1]
    gas_prandtl=chamber_transport[3]

    #Gamma: molecular weight, gamma
    moluecularweight_gamma = ispObj.get_Chamber_MolWt_gamma(Pc=Chamber_Pressure_bar, MR=Mass_Ratio, eps=Expansion_Ratio)
    gamma=moluecularweight_gamma[1]
    moluecularweight=moluecularweight_gamma[0]
    #print(gamma[1])

    #Gas Constant
    R=8314
    R_gas=(8314/moluecularweight)

    #Enthalpy
    gas_enthalpy = ispObj.get_Chamber_H(Pc=Chamber_Pressure_bar, MR=Mass_Ratio, eps=Expansion_Ratio)

    #Cstar
    Cstar_meters_sec = Cstar_efficiency*(ispObj.get_Cstar(Pc=Chamber_Pressure_bar, MR=Mass_Ratio))
    #print(Cstar_meters_sec)

    #Temps: Chamber, Nozzle, Exit
    temperatures = ispObj.get_Temperatures(Pc=Chamber_Pressure_bar, MR=Mass_Ratio)
    chamber_temp = Cstar_efficiency*(temperatures[0]) #in kelvin

    #Chamber Rho
    chamber_rho = ispObj.get_Chamber_Density(Pc=Chamber_Pressure_bar, MR = Mass_Ratio)

    chamber_pressure_true_PA = (Cstar_meters_sec*total_mdot)/throat_area
    dev=chamber_pressure_true_PA-Chamber_Pressure_Pa
    Chamber_Pressure_Pa=chamber_pressure_true_PA
    Chamber_Pressure_bar=Chamber_Pressure_Pa/100000

#Isentropic Flow Calculations
time = datetime.datetime.now()
timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
filename = f"IsentropicFlow_{timestamp}.csv"

timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

filename = Path.cwd() / f"IsentropicFlow_{identifier}.csv"
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
        Input.temperature(temperature-273.15),   #maybe this will be fixed eventually
    )

    return coolant.specific_heat

def coolant_s_tens(pressure, temperature):
    coolant = Fluid(Regen_Coolant).with_state(
            Input.quality(0),   
            Input.temperature(temperature-273.15)   #maybe this will be fixed eventually
        )
    
    return coolant.surface_tension 

def coolant_enthalpy_evap(pressure):
    liquid = Fluid(Film_Coolant).bubble_point_at_pressure(pressure)
    vapour = Fluid(Film_Coolant).dew_point_at_pressure(pressure)

    return vapour.enthalpy - liquid.enthalpy

def coolant_sat_temp(pressure):
    coolant = Fluid(Film_Coolant).dew_point_at_pressure(pressure)
    return coolant.temperature+273.15

def film_coolant_cp(pressure):
    filmcoolant = Fluid(Film_Coolant).with_state(
            Input.quality(0),   
            Input.pressure(pressure)   #maybe this will be fixed eventually
        )
    
    return filmcoolant.specific_heat

def film_coolant_gamma_vapor(pressure):
    filmcoolant = Fluid(Film_Coolant).with_state(
            Input.quality(100),   
            Input.pressure(pressure)   #maybe this will be fixed eventually
        )
    cp = filmcoolant.specific_heat
    R_co = 8.314462618 / filmcoolant.molar_mass  # molar mass: kg/mol
    cv = cp - R_co
    gamma_co = cp / cv
    return gamma_co

def film_enthalpy(pressure):
    filmcoolant = Fluid(Film_Coolant).with_state(
            Input.quality(0),   
            Input.pressure(pressure)   #maybe this will be fixed eventually
        )
    return filmcoolant.enthalpy

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

def film_liquid_length(hg,v_gas,mach_no): #Huzel & Huang empirical method, how on earth did they come up with this 
    delta=1.3 #Empircal value for film which is injected paralel to combustion gas
    T_sat=coolant_sat_temp(Chamber_Pressure_Pa) #it's intresting to me that entry temp isnt added into this. also seem to be yielding shorter liquid length than RPA.
    T_re=chamber_temp*(1+(gas_prandtl**(1/3)*((gamma-1)/2)*mach_no**2))
    latent_heat = coolant_enthalpy_evap(Chamber_Pressure_Pa)
    B=( #near wall heat input parameter
        (gas_Cp*(T_re-T_sat))
        / latent_heat       
    )
    nu_fcl=(B/(B+1))

    surface_tens=coolant_s_tens(Chamber_Pressure_Pa,T_sat)
    Xe=( #entrainment parameter, entrainment means when the film becomes mixed with the hot gas (and becomes mostly inneffective)
        delta*((chamber_rho)**0.5)*(v_gas)*((chamber_temp/T_sat)**0.25)/surface_tens
    )
    Xr = Xe*surface_tens
    alfa_f_c = (1+3*Xr**-0.8)*(7/(20*10**4)*Xe+1)
    A_f_c_in = ((1.3/(200000)*Xe)+0.1)
    A_f_C = A_f_c_in / 0.0254

    St=(
        hg
        / (chamber_rho*v_gas*gas_Cp)
        ) #stanton number? note to self: understand this
    V = (
        ((math.pi*chamber_diameter*(chamber_rho*v_gas))) 
        * St*B*alfa_f_c
    )
    Film_Liquid_Length = (
        (1/A_f_C)*
        math.log(1+((A_f_C*Film_Mdot)/V))
    )
    # Capture existing model values for reporting; do not recalculate them later.
    film_summary.update({
        "Liquid film length": (Film_Liquid_Length, "m"),
        "Liquid film end X": (x_positions[0] + Film_Liquid_Length, "m"),
        "Film saturation temperature": (T_sat, "K"),
        "Film heat-input parameter B": (B, "-"),
        "Film latent heat": (latent_heat, "J/kg"),
        "Film nu at liquid-film end": (nu_fcl, "-"),
        "Film entrainment parameter Xe": (Xe, "correlation units"),
        "Film roughness parameter Xr": (Xr, "correlation units"),
        "Film heat-transfer multiplier alpha": (alfa_f_c, "-"),
        "Film entrainment coefficient A": (A_f_C, "1/m"),
        "Film Stanton number": (St, "-"),
        "Film recovery temperature used": (T_re, "K"),
        "Film surface tension used": (surface_tens, "N/m"),
        "Film vaporization rate per length V": (V, "kg/(m*s)"),
        "Film reference heat-transfer coefficient": (hg, "W/(m^2*K)"),
    })
    return Film_Liquid_Length, T_sat,nu_fcl

def phi_m(station):
    x = x_positions[station]

    if x <= 0:
        Lc = -x_positions[0]
        return 1.75 + x * (3.5 - 1.75) / (-Lc)

    Ar = area_ratios[station]
    return (
        0.277
        + 1.4961 * math.exp(-0.1199 * Ar)
        + 0.14612 / Ar
    )

def shape_factor(film_entrainment,nu):
    liquid_film_entrainment_ratio= (1 / (0.6*nu)) -1
    entrained_liquid=Film_Mdot*liquid_film_entrainment_ratio
    entrained_i=Film_Mdot*film_entrainment
    if (entrained_i > entrained_liquid + 0.6*Film_Mdot):
        return 0.6 + 0.263*((entrained_i-entrained_liquid)/Film_Mdot)
    else:
        return 0.758

def film_entrainment_ratio(station, nu):
    phi_L = 0.004
    core_mdot = total_mdot - Film_Mdot

    if core_mdot <= 0 or Film_Mdot <= 0 or nu <= 0:
        raise ValueError("Film flow, core flow, and nu must be positive")

    phi = phi_m(station)
    x_bar_i = x_bar(station)
    radius = chamber_radii[station]

    liquid_ratio = 1 / (0.6 * nu) - 1
    entrained_liquid = Film_Mdot * liquid_ratio

    mass_term = entrained_liquid / core_mdot
    distance_term = phi_L * x_bar_i / radius
    radicand = 1 - mass_term - distance_term

    if not math.isfinite(radicand) or radicand < 0:
        raise ValueError(
            f"Invalid entrainment calculation at station {station}: "
            f"sqrt_argument={radicand:.6g}, "
            f"mass_term={mass_term:.6g}, "
            f"distance_term={distance_term:.6g}, "
            f"phi={phi:.6g}, x_bar={x_bar_i:.6g}, "
            f"film_length={f_liquid_len:.6g}, nu={nu:.6g}"
        )

    film_station_data[station].update({
        "phi": phi, "x_bar": x_bar_i,
    })
    return (
        (core_mdot / Film_Mdot)
        * 2 * distance_term
        * math.sqrt(radicand)
        + liquid_ratio
    )

def nu_i(station,nu):
    entrainment=film_entrainment_ratio(station,nu)
    film_station_data[station]["entrainment_ratio"] = entrainment
    film_station_data[station]["shape_factor"] = shape_factor(entrainment, nu)
    nu_i=(
        1/(shape_factor(entrainment,nu)*(1+film_entrainment_ratio(station,nu)))
    )
    return nu_i

def film_vapor_temperature(station, T_co_1):
    gamma_co = film_coolant_gamma_vapor(chamber_pressure[1])
    P_i = chamber_pressure[station]
    film_station_data[station]["gamma_used"] = gamma_co

    return T_co_1 * (P_i / chamber_pressure[1])**(
        (gamma_co - 1) / gamma_co
    )

def film_Taw_effective(station,nu, T_co_1):
    T_co_i=film_vapor_temperature(station,T_co_1)

    T_ref = 298.15  # Use the reference temperature specified by the source.
    Cp_i=film_coolant_cp(chamber_pressure[station])

    hg_i= gas_Cp * (gas_temp[station] - T_ref)
    hco_i = Cp_i * (T_co_i - T_ref)

    eta = nu_i(station, nu)
    cp_mix = eta * Cp_i + (1 - eta) * gas_Cp

    Taw_eff=(
        ((hg_i-nu_i(station,nu)*(gas_enthalpy-hco_i))-(1-gas_prandtl**(1/3))*(gas_enthalpy-hco_i))
        / (nu_i(station, nu) * Cp_i+ (1 - nu_i(station, nu)) * gas_Cp)
    )
    n = nu_i(station, nu)
    r = gas_prandtl**(1 / 3)

    numerator = hco_i + (r - n) * (gas_enthalpy - hco_i)
    denominator = n * Cp_i
    # These are the quantities actually used by the current model.
    # Cp_i currently comes from the saturated-liquid Cp helper.
    film_station_data[station].update({
        "film_temperature": T_co_i,
        "cp_used": Cp_i,
        "eta": eta,
        "cp_mix": cp_mix,
        "coolant_sensible_enthalpy": hco_i,
        "gas_sensible_enthalpy": hg_i,
    })
    return(Taw_eff)

#Isentropic Flow lookup table
x_positions=[]
chamber_radii=[]
area_ratios=[]
mach_number=[]
gas_temp=[]
chamber_pressure=[]
adiabatic_wall_temp=[]
gas_velocity=[]

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
            Pressure=(Chamber_Pressure_Pa*Pres_Ratio_Flow)
            Velocity=(Mach*(gamma*R_gas*Gas_Temp)**0.5)

            r = gas_prandtl ** (1 / 3) #note to self to remember what R is, I remember that prandtl**1/3 is turbulent and prandtl**1/2 is laminar, but otherwise lost lol

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
            gas_velocity.append(float(Velocity))
            chamber_pressure.append(float(Pressure))
        except ValueError:
            continue

#cstar*mdot/throatarea
chamber_pressure_check = (Cstar_meters_sec*total_mdot)/throat_area

if not math.isclose(Chamber_Pressure_Pa,chamber_pressure_check, rel_tol=0.001):
    print('FATAL ERROR: Input Chamber pressure & Calculated Chamber Pressure mismatch')
    print('Calculated:',chamber_pressure_check, ' PA')
    print('Input:',Chamber_Pressure_Pa,' PA')
    exit(12345)


def x_bar(endstation):
    reference_station = 0
    x_start = x_positions[reference_station] + f_liquid_len
    x_end = x_positions[endstation]

    if x_end <= x_start:
        return 0.0

    r_c = chamber_radii[reference_station]

    def gas_mass_flux(j):
        # Ideal-gas density, consistent with your constant-R flow model.
        # chamber_pressure MUST be in Pa; gas_temp in K.
        rho_g = chamber_pressure[j] / (R_gas * gas_temp[j])
        return rho_g * gas_velocity[j]

    reference_flux = gas_mass_flux(reference_station)

    if reference_flux <= 0:
        raise ValueError("Reference gas mass flux must be positive")

    def integrand(j):
        return (
            (r_c / chamber_radii[j])
            * (gas_mass_flux(j) / reference_flux)
            * phi_m(j)
        )

    integral = 0.0

    for j in range(reference_station, endstation):
        xa = x_positions[j]
        xb = x_positions[j + 1]

        if xb <= xa:
            raise ValueError("x_positions must increase downstream")

        if xb <= x_start:
            continue

        left = max(xa, x_start)

        fa = integrand(j)
        fb = integrand(j + 1)

        # Interpolate the integrand at a film end inside this segment.
        fraction = (left - xa) / (xb - xa)
        f_left = fa + fraction * (fb - fa)

        integral += 0.5 * (f_left + fb) * (xb - left)

    return integral

# Station-sized storage is valid with film cooling enabled or disabled.
uncooled_adiabatic_wall_temp = adiabatic_wall_temp.copy()
film_state = ["None"] * len(x_positions)
film_station_data = [{} for _ in x_positions]
film_summary = {}

if Film_Cooling:
    hg_1 = (
                hg_bartz(
                    throat_radius,
                    gas_Cp,
                    gas_viscocity,
                    gas_prandtl,
                    Chamber_Pressure_Pa,
                    Cstar_meters_sec,
                    area_ratios[1],
                    throat_r_D * (2*throat_radius)
                )
    )
    f_len=film_liquid_length(hg_1,gas_velocity[1],mach_number[1]) 
    f_liquid_len=f_len[0]
    t_sat=f_len[1]
    nu_fcl=f_len[2]
    #print(f_liquid_len)

    #print(f_liquid_len)
    for i in range(len(adiabatic_wall_temp)):
        x_1=x_positions[0]
        dist=abs(x_1-x_positions[i])
        #print(dist)
        if dist<=f_liquid_len:
            film_state[i] = "Liquid"
            film_station_data[i]["film_temperature"] = t_sat
            adiabatic_wall_temp[i]=t_sat
            #print(adiabatic_wall_temp[i])
            #exit(000)
        else:
            film_state[i] = "Gas"
            adiabatic_wall_temp[i]=film_Taw_effective(i,nu_fcl,Film_Inlet_Temp)
    film_summary["Film liquid-station count"] = (film_state.count("Liquid"), "-")
    film_summary["Film gas-station count"] = (film_state.count("Gas"), "-")



#coolant velocity math
#assuming constant density for this, b/c lazy

coolant_pressures=[]
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

coolant_velocity=Coolant_Mdot/(coolant_rho(Coolant_Inlet_Pressure_PA,Coolant_Inlet_Temp)*Channel_Area*Channel_Count)  #m/s
first_pass_dist=0
second_pass_dist=0
j=0
#print("Coolant Velocity:", coolant_velocity)

if Two_Pass:
    #print(inletpos)
    for i in range(len(area_ratios)):
        first_pass_dist=first_pass_dist+abs((x_positions[i]-x_positions[j]))
        #print(first_pass_dist)
        coolant_pressures.append(
            Coolant_Inlet_Pressure_PA
            - deltaP(f,first_pass_dist,coolant_velocity,coolant_rho(Coolant_Inlet_Pressure_PA,Coolant_Inlet_Temp))
        )
        j=i

for i in range(len(area_ratios) - 1, -1, -1):
    second_pass_dist=second_pass_dist+abs((x_positions[i]-x_positions[j]))
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
Tc_2_list=[]
Q_list=[]
q_heatflux=[] #heat flux
coolant_specific_heat_list=[]
coolant_rho_list=[]
coolant_abs_viscocity_list=[]
coolant_kin_viscocity_list=[]
coolant_conductivity_list=[]

coolant_temp = Coolant_Inlet_Temp

if Two_Pass:
    startingvalue=1
else:
    startingvalue=-1

for i in range(1, len(area_ratios)):
    dx = abs(x_positions[i] - x_positions[i - 1])
    dr = chamber_radii[i] - chamber_radii[i - 1]
    ds = math.sqrt(dx**2 + dr**2)

    gas_area = (
        math.pi
        * (chamber_radii[i] + chamber_radii[i - 1])
        * ds
    )

    coolant_area = (
            hydraulicpermiter*coolant_parameter #arbitrary :P, to be one of the values to tune
            * Channel_Count
            * ds
        )
    if Two_Pass:
        coolant_area=coolant_area*2

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
    if Two_Pass:
        coolant_temp=coolant_temp+(Q*0.5)/(Coolant_Mdot*coolant_specific_heat(coolant_pressures[i],coolant_temp))
    else:
        coolant_temp=coolant_temp+Q/(Coolant_Mdot*coolant_specific_heat(coolant_pressures[i],coolant_temp))


    hg.append(hg_local)
    hl.append(hl_local)

    coolant_specific_heat_list.append(coolant_specific_heat(coolant_pressures[i],coolant_temp))
    coolant_rho_list.append(coolant_rho(coolant_pressures[i],coolant_temp))
    coolant_abs_viscocity_list.append(coolant_abs_viscocity(coolant_pressures[i],coolant_temp))
    coolant_kin_viscocity_list.append(coolant_kin_viscocity(coolant_pressures[i],coolant_temp))
    coolant_conductivity_list.append(coolant_conductivity(coolant_pressures[i],coolant_temp))

    Twg_list.append(Twg_calculated)
    Twl_list.append(Twl)
    Tc_list.append(coolant_temp)

    Q_list.append(Q)
    q_heatflux.append(Q/gas_area)
    coolant_temperature_check = Twl - Q * R_l
    #print(coolant_sat_temp(coolant_pressures[i]))

if Two_Pass:
    i=0
    for i in range(len(area_ratios)-1, 0, startingvalue*-1):
        coolant_temp=coolant_temp+(Q_list[i-1]*0.5)/(Coolant_Mdot*coolant_specific_heat(coolant_pressures[i+len(area_ratios)],coolant_temp))
        Tc_2_list.append(coolant_temp)

savefilename = Path.cwd() / f"ThermalOutput_{identifier}.csv"

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
    writer.writerow(["Film_Cooling", Film_Cooling, "bool"])
    writer.writerow(["Regen_Coolant", str(Regen_Coolant), "-"])
    writer.writerow(["Film_Coolant", str(Film_Coolant), "-"])
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

    writer.writerow(["total_mdot", total_mdot, "kg/s"])
    writer.writerow(["Film_Mdot (active)", Film_Mdot if Film_Cooling else 0.0, "kg/s"])
    writer.writerow(["Film_Inlet_Temp", Film_Inlet_Temp if Film_Cooling else "", "K"])
    writer.writerow(["Film / total mass flow", Film_Mdot / total_mdot if Film_Cooling else 0.0, "-"])
    writer.writerow([])
    writer.writerow(["Film summary", "Value", "Unit"])
    writer.writerow(["Film cooling status", "Enabled" if Film_Cooling else "Disabled", "-"])
    if Film_Cooling:
        for label, (value, unit) in film_summary.items():
            writer.writerow([label, value, unit])
        writer.writerow(["Gas-film Cp basis", "Saturated-liquid Cp (current model)", "-"])

    # Blank rows separating inputs from results
    writer.writerow([])
    writer.writerow([])

    # Thermal-output headings
    writer.writerow([
        "X Position (m)",
        "Chamber Radius (m)",
        "Area Ratio",
        "Mach Number",
        "Gas Temperature (K)",
        "Gas Velocity (m/s)",
        "Adiabatic Wall Temperature (K)",
        "Gas Pressure (Pa)",
        "",
        "Coolant Velocity (m/s)",
        "Coolant Pressure (Pa)",
        "Coolant Specific Heat (J/kg-K)",
        "Coolant Density (kg/m^3)",
        "Coolant Dynamic Viscosity (Pa-s)",
        "Coolant Kinematic Viscosity (m^2/s)",
        "Coolant Thermal Conductivity (W/m-K)",
        "",
        "Coolant-Side Heat Transfer Coefficient (W/m^2-K)",
        "Gas-Side Heat Transfer Coefficient (W/m^2-K)",
        "Gas-Side Wall Temperature (K)",
        "Coolant-Side Wall Temperature (K)",
        "Coolant Temperature (K)",
        "",
        "Heat Transferred (W)",
        "Heat Flux (W/m^2)",
        "",
        "Film Coolant State",
        "Distance from Film Injection (m)",
        "Uncooled Adiabatic Wall Temperature (K)",
        "Film Taw Reduction (K)",
        "Film Coolant Temperature Used (K)",
        "Film Coolant Cp Used (J/kg-K)",
        "Film Mixing Fraction eta (-)",
        "Film Mixture Cp Used (J/kg-K)",
        "Film Coolant Sensible Enthalpy Used (J/kg)",
        "Film Gas Sensible Enthalpy Used (J/kg)",
        "Film Flow Multiplier phi (-)",
        "Film Equivalent Distance x_bar (m)",
        "Film Entrainment Ratio (-)",
        "Film Shape Factor (-)",
        "Film Reference Vapor Gamma Used (-)"
    ])

    for k in range(len(Twg_list)):
            i=k+1
            
            writer.writerow([
            x_positions[i],
            chamber_radii[i],
            area_ratios[i],
            mach_number[i],
            gas_temp[i],
            gas_velocity[i],
            adiabatic_wall_temp[i],
            chamber_pressure[i],
            "",
            coolant_velocity,
            coolant_pressures[i],
            coolant_specific_heat_list[k],
            coolant_rho_list[k],
            coolant_abs_viscocity_list[k],
            coolant_kin_viscocity_list[k],
            coolant_conductivity_list[k],
            "",
            hl[k],
            hg[k],
            Twg_list[k],
            Twl_list[k],
            Tc_list[k],
            "",
            Q_list[k],
            q_heatflux[k],
            "",
            film_state[i],
            abs(x_positions[i] - x_positions[0]) if Film_Cooling else "",
            uncooled_adiabatic_wall_temp[i],
            uncooled_adiabatic_wall_temp[i] - adiabatic_wall_temp[i] if Film_Cooling else "",
            film_station_data[i].get("film_temperature", ""),
            film_station_data[i].get("cp_used", ""),
            film_station_data[i].get("eta", ""),
            film_station_data[i].get("cp_mix", ""),
            film_station_data[i].get("coolant_sensible_enthalpy", ""),
            film_station_data[i].get("gas_sensible_enthalpy", ""),
            film_station_data[i].get("phi", ""),
            film_station_data[i].get("x_bar", ""),
            film_station_data[i].get("entrainment_ratio", ""),
            film_station_data[i].get("shape_factor", ""),
            film_station_data[i].get("gamma_used", "")
        ])

print(f"Output CSV created at: {savefilename.resolve()}")


def print_range(label, values, unit):
    if values:
        print(f"{label}: {min(values):.6g} to {max(values):.6g} {unit}")
    else:
        print(f"{label}: N/A")

print("\nPerformance Results")
print("Chamber Pressure", Chamber_Pressure_bar,'bar')
print("\nThermal results")
print_range("Effective adiabatic wall temperature", adiabatic_wall_temp, "K")
print_range("Uncooled adiabatic wall temperature", uncooled_adiabatic_wall_temp, "K")
print_range("Gas temperature", gas_temp, "K")
print_range("Gas velocity", gas_velocity, "m/s")
print_range("Regen coolant temperature (all modeled passes)", Tc_list + Tc_2_list, "K")
print(f"Regen coolant velocity: {coolant_velocity:.6g} m/s")
print_range("Coolant-side wall temperature", Twl_list, "K")
print_range("Coolant temperature", Tc_list, "K")
print_range("Gas-side wall temperature", Twg_list, "K")
print_range("Heat flux", q_heatflux, "W/m^2")
print(f"Sum of segment heat-transfer values: {sum(Q_list):.6g} W")

print("\nFilm cooling:", "Enabled" if Film_Cooling else "Disabled")
if Film_Cooling:
    print(f"Film coolant: {Film_Coolant}")
    print(f"Film mass flow: {Film_Mdot:.6g} kg/s")
    print(f"Film inlet temperature: {Film_Inlet_Temp:.6g} K")
    print(f"Film fraction of total flow: {100 * Film_Mdot / total_mdot:.6g} %")
    for label, (value, unit) in film_summary.items():
        print(f"{label}: {value:.6g} {unit}")
    print_range("Film coolant temperature used", [d["film_temperature"] for d in film_station_data if "film_temperature" in d], "K")
    print_range("Gas-film mixing fraction eta", [d["eta"] for d in film_station_data if "eta" in d], "")
    print_range("Film Taw reduction (uncooled minus effective)", [a - b for a, b in zip(uncooled_adiabatic_wall_temp, adiabatic_wall_temp)], "K")
    print("Gas-film Cp basis: saturated-liquid Cp (current model)")
#to be completley transparent chatgpt rewrote the csv writer and helped debug the film cooling but otherwise nothing else :)))