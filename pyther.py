import numpy as np
from eos_selecction import eos, convert_argument
#from cubic_parameters_1 import Parameter_eos, getdel1, compressibility_factor_cal, acentric_factor_cal
from cubic_parameters_1 import Parameter_eos

from componente import Control_arguments
from pure_data import Data_parse
from eos_pure import getdel1, compressibility_factor_cal, acentric_factor_cal

from properties_ecuations import Thermodynamic_correlations

# -----------------------------------------------------------------------------

# properties thermodynamics

Solid_Density = "Solid Density", "[kmol/m^3]", "A+B*T+C*T^2+D*T^3+E*T^4", 0
Liquid_Density = "Liquid Density", "[kmol/m^3]", "A/B^(1+(1-T/C)^D)", 1
Vapour_Pressure = "Vapour Pressure", "[Pa]", "exp(A+B/T+C*ln(T)+D*T^E)", 2
Heat_of_Vaporization = "Heat of Vaporization", "[J/kmol]", "A*(1-Tr)^(B+C*Tr+D*Tr^2)", 3
Solid_Heat_Capacity = "Solid Heat Capacity", "[J/(kmol*K)]", "A+B*T+C*T^2+D*T^3+E*T^4", 4
Liquid_Heat_Capacity = "Liquid Heat Capacity", "[J/(kmol*K)]", "A^2/(1-Tr)+B-2*A*C*(1-Tr)-A*D*(1-Tr)^2-C^2*(1-Tr)^3/3-C*D*(1-Tr)^4/2-D^2*(1-Tr)^5/5", 5
Ideal_Gas_Heat_Capacity = "Ideal Gas Heat Capacity" "[J/(kmol*K)]", "A+B*(C/T/sinh(C/T))^2+D*(E/T/cosh(E/T))^2", 6
Second_Virial_Coefficient = "Second Virial  Coefficient", "[m^3/kmol]", "A+B/T+C/T^3+D/T^8+E/T^9", 7
Liquid_Viscosity = "Liquid  Viscosity", "[kg/(m*s)]", "exp(A+B/T+C*ln(T)+D*T^E)", 8
Vapour_Viscosity = "Vapour  Viscosity", "[kg/(m*s)]", "A*T^B/(1+C/T+D/T^2)", 9
Liquid_Thermal_Conductivity = "Liquid Thermal Conductivity", "[J/(m*s*K)]", "A+B*T+C*T^2+D*T^3+E*T^4", 10
Vapour_Thermal_Conductivity = "Vapour Thermal Conductivity", "[J/(m*s*K)]", "A*T^B/(1+C/T+D/T^2)", 11
Surface_Tension = "Surface Tension", "[kg/s^2]", "A*(1-Tr)^(B+C*Tr+D*Tr^2)", 12 

# -----------------------------------------------------------------------------

RGAS = 0.08314472

# Definir el significado fisicoquímico
A0, B0, C0 = 0.0017, 1.9681, -2.7238

# Definir el significado fisicoquímico
A1, B1, C1 = -2.4407, 7.4513, 12.504

# Definir el significado fisicoquímico
D = np.array([0.428363, 18.496215, 0.338426, 0.660, 789.723105, 2.512392])


def initial_data(omega, delta_1, NMODEL, ICALC, Pc, dinputs):

    Zc, OMa, OMb = compressibility_factor_cal(delta_1)

    # initial guess for k parameter
    rk = (A1 * Zc + A0) * omega**2 + (B1 * Zc + B0) * omega + (C1 * Zc + C0)
    #rk = rk * 1.2 # 1.1 #5.2 #3.2
    if ICALC == 'constants_eps' or ICALC == 'parameters_eps' or ICALC == 'rk_param':
        rk *= 1.5
        Tr = 0.7
        Pvdat = Pc * 10 ** -(1.0 + omega)         
    elif ICALC == 'density':
        rk = rk * 1.5 #5.2 
        Tr_calculada = dinputs[4] / dinputs[0]
        Tr = Tr_calculada
        Pvdat = Pc * 10 ** -((1.0 / Tr - 1.0) * 7 * (1.0 + omega) / 3)

    return rk, Pvdat, Tr


def data_in(dinputs):

	if  ICALC == 'constants_eps':
		# CONSTANTS SPECIFICATION (Tc,Pc,OM,Vceos)
		Tc, Pc, OM, Vceos = dinputs[0], dinputs[1], dinputs[2], dinputs[3]

	if ICALC == 'parameters_eps':
		ac, b, del1, rk = dinputs[0], dinputs[1], dinputs[2], dinputs[3]

	if ICALC == 'rk_param':
		# dinputs = np.array([Tc, Pc, OM, dc, zrat, ac, b, d, rk])
		Tc, Pc, OM, Vceos, delta_1 = dinputs[0], dinputs[1], dinputs[2], dinputs[3], dinputs[7]

	if ICALC == 'density':

		Tc, Pc, omega, Vceos, delta_1 = dinputs[0], dinputs[1], dinputs[2], dinputs[3], dinputs[4]
		T_especific, RHOLSat_esp = dinputs[5], dinputs[6]


def require_ID (func):
    def wrapper (*arg):
        
        return Control_arguments(arg[0], arg[1])
 
    return wrapper


#@require_ID
def models_eos_cal(NMODEL, ICALC, dinputs):

    if NMODEL == 'SRK' or NMODEL == 'PR':
        # CONSTANTS SPECIFICATION READ [Tc, Pc, OM]
        if ICALC == 'constants_eps':
            Tc = dinputs[0]
            Pc = dinputs[1]
            OM = dinputs[2]

            if NMODEL == 'SRK':
                rm = 0.48 + 1.574 * OM - 0.175 * OM**2
                del1 = 1.0
            elif NMODEL == 'PR':
                rm = 0.37464 + 1.54226 * OM - 0.26992 * OM ** 2
                del1 = 1.0 + np.sqrt(2.0)

            Zc, OMa, OMb = compressibility_factor_cal(del1)
            
            Vceos = (Zc * RGAS * Tc) / Pc

            ac = OMa * (RGAS * Tc) ** 2 / Pc
            b = OMb * (RGAS * Tc) / Pc
            
            params = [ac, b, rm, del1]

        # PARAMETERS SPECIFICATION READ [ac, b, rm]
        if ICALC == 'parameters_eps':

        	ac = dinputs[0]
        	b = dinputs[1]
        	rm = dinputs[2]

        	Tc = (OMb * ac) / (OMa * RGAS * b)
        	Pc = OMb * RGAS * Tc / b
        	Vceos = Zc * RGAS * Tc / Pc

        	if NMODEL == 'SRK':
        		del1 = 1.0
        		al = -0.175
        		be = 1.574
        		ga = 0.48 - rm
        	elif NMODEL == 'PR':
        		del1 = 1.0 + np.sqrt(2.0)
        		al = -0.26992
        		be = 1.54226
        		ga = 0.37464 - rm

        	OM = acentric_factor_cal(al, be, ga)
        	constants = [Tc, Pc, OM, Vceos]

    elif NMODEL == 'RKPR':
        if ICALC == 'constants_eps':
            # CONSTANTS SPECIFICATION READ [Tc, Pc, OM, Vceos]
            Tc = dinputs[0]
            Pc = dinputs[1]
            OM = dinputs[2]
            Vceos = dinputs[3]
            
            Zc = Pc * Vceos / (RGAS * Tc)
            
            del1ini = D[0] + D[1] * (D[2] - Zc) ** D[3] + D[4] * (D[2] - Zc)** D[5]
            print('del1ini = {0}'.format(del1ini))

            delta_1 = getdel1(Zc, del1ini)[0]

            Zc, OMa, OMb = compressibility_factor_cal(delta_1)
            
            print('Zc = {0}'.format(Zc))

            ac = OMa * (RGAS * Tc) ** 2 / Pc
            b = OMb * (RGAS * Tc) / Pc

            # calcular rk
            rk, Pvdat, Tr = initial_data(OM, delta_1, NMODEL, ICALC, Pc, dinputs)

            eos_calculation = Parameter_eos()
            rk_cal = eos_calculation.resolver_rk_cal(rk, delta_1, Pvdat, Pc, Tc, Tr)

            # rk = 1

            params = [ac, b, rk, delta_1]

        elif ICALC == 'parameters_eps':
        	# PARAMETERS SPECIFICATION READ [ac, b, rk, del1]

            ac = dinputs[0]
            b = dinputs[1]
            del1 = dinputs[2]
            rk = dinputs[3]

            Zc, OMa, OMb = compressibility_factor_cal(del1)

            Tc = OMb * ac / (OMa * RGAS * b)
            Pc = OMb * RGAS * Tc / b
            Vceos = Zc * RGAS * Tc / Pc

            al = A1 * Zc + A0
            be = B1 * Zc + B0
            ga = C1 * Zc + C0 - rk

            OM = acentric_factor_cal(al, be, ga)

            constants = [Tc, Pc, OM, Vceos]
        elif ICALC == 'rk_param':
        	# CONSTANTS SPECIFICATION and del1 READ [Tc, Pc, OM, del1]

            Tc = dinputs[0]
            Pc = dinputs[1]
            OM = dinputs[2]
        
            delta_1 = dinputs[7]

            rk, Pvdat, Tr = initial_data(OM, delta_1, NMODEL, ICALC, Pc, dinputs)

            eos_calculation = Parameter_eos()
            rk_cal = eos_calculation.resolver_rk_cal(rk, delta_1, Pvdat, Pc, Tc, Tr)

        elif ICALC == 'density':
        	# CONSTANTS SPECIFICATION and (T, RhoLsat) READ [Tc, Pc, OM, del1, T, RHOLsat]
            # Trho = T / Tc,  read initial value of del1

            Tc = dinputs[0]
            Pc = dinputs[1]
            OM = dinputs[2]

            delta_1 = dinputs[3]
            
            T_especific = dinputs[4]
            RHOLSat_esp = dinputs[5]

            rk, Pvdat, Tr = initial_data(OM, delta_1, NMODEL, ICALC, Pc, dinputs)
            eos_calculation = Parameter_eos()
            delta_1_parameter = eos_calculation.resolver_delta_1_cal(delta_1, rk, Pvdat, RHOLSat_esp, Pc, Tc, Tr)


    print('The NMODEL is eos_{0} and method ICALC is {1}'.format(NMODEL, ICALC))


    if ICALC == 'constants_eps':
        print("params = [ac, b, rm, del1]")
        #ac, b, rm, del1 = params
        #print("ac = {0} b = {1} rm = {2} del1 = {3}".format(ac, b, rm, del1))
        return params
    elif ICALC == 'parameters_eps':
        print("constants = [Tc, Pc, OM, Vceos]")
        print(constants)
        return constants
    elif ICALC == 'rk_param':
        print('The parameter rk_cal is {0}'.format(rk_cal))
        return rk_cal
    elif ICALC == 'density':
        print('The parameter delta1(rho,T) = {0}'.format(delta_1_parameter))
        return delta_1_parameter


def print_properties_component(component, properties_component):
    print('Component = {0}'.format(component))
    print('Acentric_factor = {0}'.format(properties_component[1]['Omega']))
    print('Critical_Temperature = {0} K'.format(properties_component[1]['Tc']))
    print('Critical_Pressure = {0} Bar'.format(properties_component[1]['Pc']))
    print('Critical_Volume = {0} cm3/mol'.format(properties_component[1]['Vc']))
    print('Compressibility_factor_Z = {0}'.format(properties_component[1]['Zc']))
    print("\n")


def main():

    print("-" * 79)
    # component = 'METHANE'
    # component = "ETHANE"
    # component = "3-METHYLHEPTANE"
    # component = "n-PENTACOSANE"
    component = "ISOBUTANE"
    NMODEL = "RKPR"
    ICALC = "constants_eps"
    # ICALC = "density"

    properties_data = Data_parse()
    properties_component = properties_data.selec_component(component)
    print_properties_component(component, properties_component)
    dinputs = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
                        properties_component[1]['Omega'], properties_component[1]['Vc']])

    component_eos = models_eos_cal(NMODEL, ICALC, dinputs)
    print(component_eos[0])
    print('-' * 79)


if __name__ == '__main__':
    main()


def n_main():

    print("-" * 79)

    dppr_file = "PureFull_mod_properties.xls"
    #print(dppr_file)

    thermodynamic_correlations = Thermodynamic_correlations(dppr_file)

    #data = data_base.read_dppr()       
    #data_name = data_base.data_name_cal()

    component = 'METHANE'
    #component = "ETHANE"
    #component = "3-METHYLHEPTANE"
    #component = "n-PENTACOSANE"
    #component = "ISOBUTANE"
    #component = "n-TETRADECANE"

    #components = ["METHANE", "n-TETRACOSANE"]

    temp = [180.4, 181.4, 185.3, 210, 85]
    #temp = 180.4

    property_thermodynamics = thermodynamic_correlations.property_cal(component, Vapour_Pressure, temp)
    #property_thermodynamics = property_cal(components, Vapour_Pressure, temp)
    #property_thermodynamics = property_cal(component, Vapour_Pressure, [180.4, 181.4, 185.3, 210, 85])
    #property_thermodynamics = thermodynamic_correlations.property_cal(component, Vapour_Pressure)
    print(property_thermodynamics)

    print('-' * 79)

#n_main()


