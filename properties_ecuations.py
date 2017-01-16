import numpy as np
import pandas as pd
import pyther as pt

# properties thermodynamics

Solid_Density = "Solid Density", "[kmol/m^3]", "A+B*T+C*T^2+D*T^3+E*T^4", 0
Liquid_Density = "Liquid Density", "[kmol/m^3]", "A/B^(1+(1-T/C)^D)", 1
Vapour_Pressure = "Vapour Pressure", "[Pa]", "exp(A+B/T+C*ln(T)+D*T^E)", 2
Heat_of_Vaporization = "Heat of Vaporization", "[J/kmol]", "A*(1-Tr)^(B+C*Tr+D*Tr^2)", 3
Solid_Heat_Capacity = "Solid Heat Capacity", "[J/(kmol*K)]", "A+B*T+C*T^2+D*T^3+E*T^4", 4
Liquid_Heat_Capacity = "Liquid Heat Capacity", "[J/(kmol*K)]", "A^2/(1-Tr)+B-2*A*C*(1-Tr)-A*D*(1-Tr)^2-C^2*(1-Tr)^3/3-C*D*(1-Tr)^4/2-D^2*(1-Tr)^5/5", 5
Ideal_Gas_Heat_Capacity = "Ideal Gas Heat Capacity" "[J/(kmol*K)]", "A+B*(C/T/sinh(C/T))^2+D*(E/T/cosh(E/T))^2", 6
Second_Virial_Coefficient = "Second	Virial	Coefficient", "[m^3/kmol]", "A+B/T+C/T^3+D/T^8+E/T^9", 7
Liquid_Viscosity = "Liquid	Viscosity", "[kg/(m*s)]", "exp(A+B/T+C*ln(T)+D*T^E)", 8
Vapour_Viscosity = "Vapour	Viscosity", "[kg/(m*s)]", "A*T^B/(1+C/T+D/T^2)", 9
Liquid_Thermal_Conductivity = "Liquid Thermal Conductivity", "[J/(m*s*K)]", "A+B*T+C*T^2+D*T^3+E*T^4", 10
Vapour_Thermal_Conductivity = "Vapour Thermal Conductivity", "[J/(m*s*K)]", "A*T^B/(1+C/T+D/T^2)", 11
Surface_Tension = "Surface Tension", "[kg/s^2]", "A*(1-Tr)^(B+C*Tr+D*Tr^2)", 12	

print(Solid_Density)
#----------------------------------------------------------------------------------------------------------------

dppr_file = "PureFull_mod_properties.xls"

#print(dppr_file)

#component = 'METHANE'
#component = "ETHANE"
#component = "3-METHYLHEPTANE"
#component = "n-PENTACOSANE"
component = "ISOBUTANE"

#components = ["METHANE", "n-TETRACOSANE"]

properties_data = pt.Data_parse()

def read_dppr(dppr_file):
        #self.dppr_data = pd.read_excel(dppr_file).head().set_index('Name').ix[:, 1:12]

        #dppr_data = pd.read_excel(dppr_file).set_index("Name").ix[:, 1:5]
        dppr_data = pd.read_excel(dppr_file).ix[:, 0:8]
        
        # component_names = dppr_data.index.get_values()
        return dppr_data

data = read_dppr(dppr_file)

#print(data)

g = 3
components_labels = [x for x in range(0, 13*g, 13)]
data_name = data.ix[components_labels, 0].get_values()

def property_cal(component, property_thermodynamics):
	
	rho_liquido_constans = [x + property_thermodynamics[3]  for x in range(0, 13*g, 13)]
	datos_rho_liquido = data.ix[rho_liquido_constans, 1:8].get_values()

	print(datos_rho_liquido)
	datos_rho_liquido = pd.DataFrame(data=datos_rho_liquido,index=data_name,
						 columns=["A", "B", "C", "D", "E", "Min", "Max"])

	rho_liquido_constans_component = datos_rho_liquido.loc[component]

	# Liquid_density [kmol/m^3]
	# Temperature [K]
	# A / B ^ (1 + (1 - T / C) ^ D)
	A, B, C, D, E, Min, Max = rho_liquido_constans_component
	Temp_vector = np.array([Temp_vector for Temp_vector in np.arange(Min, Max)])


	if property_thermodynamics == Solid_Density:
		pass
	elif property_thermodynamics == Liquid_Density:
		rho_liquido = A / B ** (1 + (1 - Temp_vector / C) ** D)
		print(rho_liquido)
		return rho_liquido_constans_component

	elif property_thermodynamics == Vapour_Pressure:
		pass
	elif property_thermodynamics == Heat_of_Vaporization:
		pass
	elif property_thermodynamics == Solid_Heat_Capacity:
		pass
	elif property_thermodynamics == Liquid_Heat_Capacity:
		pass
	elif property_thermodynamics == Ideal_Gas_Heat_Capacity:
		pass
	elif property_thermodynamics == Second_Virial_Coefficient:
		pass
	elif property_thermodynamics == Liquid_Viscosity:
		pass
	elif property_thermodynamics == Vapour_Viscosity:
		pass
	elif property_thermodynamics == Liquid_Thermal_Conductivity:
		pass
	elif property_thermodynamics == Vapour_Thermal_Conductivity:
		pass
	elif property_thermodynamics == Surface_Tension:
		pass




property_cal("ETHANE", Liquid_Density)














# Liquid_Viscosity	[kg/(m*s)]
# Temperature [K]
# exp(A+B/T+C*ln(T)+D*T^E)

#A, B, C, D, E, Min, Max = mu_liq_cal("ETHANE")
#Temp_vector = np.array([Temp_vector for Temp_vector in np.arange(Min, Max)])
#mu_liquid = np.exp(A + B / Temp_vector + C * np.log(Temp_vector) + D * Temp_vector ** E)
#print(mu_liquid)