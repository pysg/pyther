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

#----------------------------------------------------------------------------------------------------------------

dppr_file = "PureFull_mod_properties.xls"

#print(dppr_file)


properties_data = pt.Data_parse()

def read_dppr(dppr_file):
        #self.dppr_data = pd.read_excel(dppr_file).head().set_index('Name').ix[:, 1:12]

        #dppr_data = pd.read_excel(dppr_file).set_index("Name").ix[:, 1:5]
        dppr_data = pd.read_excel(dppr_file).ix[:, 0:8]
        
        # component_names = dppr_data.index.get_values()
        return dppr_data

data = read_dppr(dppr_file)

#print(data)

g = 30
components_labels = [x for x in range(0, 13*g, 13)]
data_name = data.ix[components_labels, 0].get_values()

def property_cal(component, property_thermodynamics, temperature = None):

	"""
	properties thermodynamics

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
	"""

	select_constans = [x + property_thermodynamics[3] for x in range(0, 13*g, 13)]
	values_constans = data.ix[select_constans, 1:8].get_values()
	table_constans = pd.DataFrame(data=values_constans,index=data_name,
						 columns=["A", "B", "C", "D", "E", "T Min [K]", "T Max [K]"])

	print(table_constans)
	component_constans = table_constans.loc[component]

	A, B, C, D, E, Min, Max = component_constans

	if temperature == None:
		Temp_vector = np.array([Temp_vector for Temp_vector in np.arange(Min, Max)])
	elif (Min < np.array(temperature) < Max):
		Temp_vector = np.array(temperature)
	else:
		return print("temperature is not valid")

	print(component_constans)

	if property_thermodynamics == Solid_Density:

		solid_Density = A + B * Temp_vector + C * Temp_vector ** 2 + D * Temp_vector ** 3 + E * Temp_vector **4
		#print(solid_Density)
		return solid_Density

	elif property_thermodynamics == Liquid_Density:
	

		liquid_Density = A / B ** (1 + (1 - Temp_vector / C) ** D)
		#print("liquid_Density = {0}".format(liquid_Density))		
		return liquid_Density

	elif property_thermodynamics == Vapour_Pressure:

		vapour_Pressure = np.exp(A + B/Temp_vector + C * np.log(Temp_vector)+D*Temp_vector **E)
		#print(vapour_Pressure)
		return vapour_Pressure

	elif property_thermodynamics == Heat_of_Vaporization:

		heat_of_Vaporization = A*(1-Tr) ** (B+C*Tr+D*Tr**2)
		#print(heat_of_Vaporization)
		return heat_of_Vaporization

	elif property_thermodynamics == Solid_Heat_Capacity:
		solid_Heat_Capacity = A + B * Temp_vector + C * Temp_vector ** 2 + D * Temp_vector ** 3 + E * Temp_vector ** 4
		#print(solid_Heat_Capacity)
		return solid_Heat_Capacity
	elif property_thermodynamics == Liquid_Heat_Capacity:
		liquid_Heat_Capacity = A ** 2 / (1-Tr) + B-2*A*C*(1-Tr)-A*D*(1-Tr)**2-C**2*(1-Tr)**3/3-C*D*(1-Tr)**4/2-D**2*(1-Tr)**5/5
		#print(liquid_Heat_Capacity)
		return(liquid_Heat_Capacity)
	elif property_thermodynamics == Ideal_Gas_Heat_Capacity:
		ideal_Gas_Heat_Capacity = A+B*(C/Temp_vector/np.sinh(C/Temp_vector))**2+D*(E/Temp_vector/np.cosh(E/Temp_vector))**2
		#print(ideal_Gas_Heat_Capacity)
		return ideal_Gas_Heat_Capacity
	elif property_thermodynamics == Second_Virial_Coefficient:
		second_Virial_Coefficient = A+B/Temp_vector+C/Temp_vector**3+D/Temp_vector**8+E/Temp_vector**9
		#print(second_Virial_Coefficient)
		return second_Virial_Coefficient
	elif property_thermodynamics == Liquid_Viscosity:
		liquid_Viscosity = np.exp(A + B / Temp_vector + C * np.log(Temp_vector) + D * Temp_vector ** E)
		#print(liquid_Viscosity)
		return liquid_Viscosity
	elif property_thermodynamics == Vapour_Viscosity:
		vapour_Viscosity = A*Temp_vector**B/(1+C/Temp_vector+D/Temp_vector**2)
		#print(vapour_Viscosity)
		return vapour_Viscosity
	elif property_thermodynamics == Liquid_Thermal_Conductivity:
		liquid_Thermal_Conductivity = A+B*Temp_vector+C*Temp_vector**2+D*Temp_vector**3+E*Temp_vector**4
		#print(liquid_Thermal_Conductivity)
		return liquid_Thermal_Conductivity
	elif property_thermodynamics == Vapour_Thermal_Conductivity:
		vapour_Thermal_Conductivity = A*Temp_vector**B/(1+C/Temp_vector+D/Temp_vector**2)
		#print(vapour_Thermal_Conductivity)
		return vapour_Thermal_Conductivity
	elif property_thermodynamics == Surface_Tension:
		surface_Tension = A*(1-Tr)**(B+C*Tr+D*Tr*2)
		#print(surface_Tension)
		return surface_Tension


component = 'METHANE'
#component = "ETHANE"
#component = "3-METHYLHEPTANE"
#component = "n-PENTACOSANE"
#component = "ISOBUTANE"
#component = "n-TETRADECANE"

#components = ["METHANE", "n-TETRACOSANE"]

#property_thermodynamics = property_cal(component, Vapour_Pressure, [180.4, 181.4, 185.3])
property_thermodynamics = property_cal(component, Vapour_Pressure, 180.4)
#property_thermodynamics = property_cal(component, Vapour_Pressure)
print(property_thermodynamics)

Min = 30
Max = 90
t = np.array([40, 34, 56, 98, 25])

#myList = [ x if x%2 else x*100 for x in range(1, 10)]

Temperature_valid = [i if Min < np.array(i) < Max else False for i in t]

#condition = Min < np.array(t) < Max

print(Temperature_valid)
#print(myList)

# Tareas

# estado de referencias del estado triple para las propiedades del gas ideal

# graficar la correlación de la densidad de liquido para comparar y sobreponer
# con la parte de la campana calculada con la ecuación de estado para el líquido

# revisar las ecuaciones para utilizar el enfoque de fugacidad de solidó para determinar 
# las demás propiedades termodinámicas

# para los graficos de velocidad de sonido tener en cuenta "por ejemplo" las isotermas  
# para el caso de una sola fase

# criterio de parar el calculo es que sea hasta el 5% de la presión de vapor de la
# isoterma a la cual se está calculando. En caso de las isotermas super criticas, 
# se toma como parametro el 5% de la presión critica

# incluir el calculo del corrimiento de volumen

# dato experiental de líquido saturado para determinar un corrimiento C
# para adicionar dicho corrimiento


# escribir un borrador del artículo

# introduccion# 
# arquitectura diagrama 
# calculo termodinámicas
# referencias del libro de mollerup, conceptos de algoritmos por ejemmplo la cmapana
# - subsecciópn calculo de la cmapana de saturación,, parametros de inicio
# sección para calcular los otros diagramas, descripción general del calculo  
# resultados o ejemplos
# -sub sección pantallazos del notebook para mostrar capacidades de ussuarios
# conclusiones


