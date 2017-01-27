import numpy as np
import pandas as pd
#import pyther as pt

# constan to size of data base of components
size_data = 30


class Thermodynamic_correlations(object):
	
	"""
		Thermodynamic_correlations

		This class is used to calculated Thermodynamic_correlations like a function of temperature.



		1. Solid_Density = "Solid Density", "[kmol/m^3]", "A+B*T+C*T^2+D*T^3+E*T^4", 0
		2. Liquid_Density = "Liquid Density", "[kmol/m^3]", "A/B^(1+(1-T/C)^D)", 1
		3. Vapour_Pressure = "Vapour Pressure", "[Pa]", "exp(A+B/T+C*ln(T)+D*T^E)", 2
		4. Heat_of_Vaporization = "Heat of Vaporization", "[J/kmol]", "A*(1-Tr)^(B+C*Tr+D*Tr^2)", 3
		5. Solid_Heat_Capacity = "Solid Heat Capacity", "[J/(kmol*K)]", "A+B*T+C*T^2+D*T^3+E*T^4", 4
		6. Liquid_Heat_Capacity = "Liquid Heat Capacity", "[J/(kmol*K)]", "A^2/(1-Tr)+B-2*A*C*(1-Tr)-A*D*(1-Tr)^2-C^2*(1-Tr)^3/3-C*D*(1-Tr)^4/2-D^2*(1-Tr)^5/5", 5
		7. Ideal_Gas_Heat_Capacity = "Ideal Gas Heat Capacity" "[J/(kmol*K)]", "A+B*(C/T/sinh(C/T))^2+D*(E/T/cosh(E/T))^2", 6
		8. Second_Virial_Coefficient = "Second	Virial	Coefficient", "[m^3/kmol]", "A+B/T+C/T^3+D/T^8+E/T^9", 7
		9. Liquid_Viscosity = "Liquid	Viscosity", "[kg/(m*s)]", "exp(A+B/T+C*ln(T)+D*T^E)", 8
		10. Vapour_Viscosity = "Vapour	Viscosity", "[kg/(m*s)]", "A*T^B/(1+C/T+D/T^2)", 9
		11. Liquid_Thermal_Conductivity = "Liquid Thermal Conductivity", "[J/(m*s*K)]", "A+B*T+C*T^2+D*T^3+E*T^4", 10
		12. Vapour_Thermal_Conductivity = "Vapour Thermal Conductivity", "[J/(m*s*K)]", "A*T^B/(1+C/T+D/T^2)", 11
		13. Surface_Tension = "Surface Tension", "[kg/s^2]", "A*(1-Tr)^(B+C*Tr+D*Tr^2)", 12	
	"""
	def __init__(self, dppr_file):
		
		self.dppr_file = dppr_file
		

	def read_dppr(self):
		
		self.dppr_data = pd.read_excel(self.dppr_file).ix[:, 0:8]
		
		return self.dppr_data


	def data_name_cal(self):

		components_labels = [x for x in range(0, 13*size_data, 13)]
		#self.data_name = data.ix[components_labels, 0].get_values()
		self.data_name = self.read_dppr().ix[components_labels, 0].get_values()

		return self.data_name

	# properties thermodynamics
	def select_property(self, property_thermodynamics):

		


		if property_thermodynamics == "Solid_Density":
			Solid_Density = "Solid Density", "[kmol/m^3]", "A+B*T+C*T^2+D*T^3+E*T^4", 0
			return Solid_Density
		elif property_thermodynamics == "Liquid_Density":
			Liquid_Density = "Liquid Density", "[kmol/m^3]", "A/B^(1+(1-T/C)^D)", 1
			return Liquid_Density 
		elif property_thermodynamics == "Vapour_Pressure":
			Vapour_Pressure = "Vapour Pressure", "[Pa]", "exp(A+B/T+C*ln(T)+D*T^E)", 2
			return Vapour_Pressure
		elif property_thermodynamics == "Heat_of_Vaporization":
			Heat_of_Vaporization = "Heat of Vaporization", "[J/kmol]", "A*(1-Tr)^(B+C*Tr+D*Tr^2)", 3
			return Heat_of_Vaporization 
		elif property_thermodynamics == "solid_Heat_Capacity":
			Solid_Heat_Capacity = "Solid Heat Capacity", "[J/(kmol*K)]", "A+B*T+C*T^2+D*T^3+E*T^4", 4
			return Solid_Heat_Capacity
		elif property_thermodynamics == "Liquid_Heat_Capacity":
			
			Liquid_Heat_Capacity = "Liquid Heat Capacity", "[J/(kmol*K)]", "A^2/(1-Tr)+B-2*A*C*(1-Tr)-A*D*(1-Tr)^2-C^2*(1-Tr)^3/3-C*D*(1-Tr)^4/2-D^2*(1-Tr)^5/5", 5
			return Liquid_Heat_Capacity
		
		elif property_thermodynamics == "Ideal_Gas_Heat_Capacity":
		
			Ideal_Gas_Heat_Capacity = "Ideal Gas Heat Capacity", "[J/(kmol*K)]", "A+B*(C/T/sinh(C/T))^2+D*(E/T/cosh(E/T))^2", 6
			return Ideal_Gas_Heat_Capacity
		
		elif property_thermodynamics == "Second_Virial_Coefficient": 
		
			Second_Virial_Coefficient = "Second	Virial	Coefficient", "[m^3/kmol]", "A+B/T+C/T^3+D/T^8+E/T^9", 7
			return Second_Virial_Coefficient
		
		elif property_thermodynamics == "Liquid_Viscosity": 
			Liquid_Viscosity = "Liquid	Viscosity", "[kg/(m*s)]", "exp(A+B/T+C*ln(T)+D*T^E)", 8
			return Liquid_Viscosity
		elif property_thermodynamics == "Vapour_Viscosity": 
			Vapour_Viscosity = "Vapour	Viscosity", "[kg/(m*s)]", "A*T^B/(1+C/T+D/T^2)", 9
			return Vapour_Viscosity
		elif property_thermodynamics == "Liquid_Thermal_Conductivity": 
			Liquid_Thermal_Conductivity = "Liquid Thermal Conductivity", "[J/(m*s*K)]", "A+B*T+C*T^2+D*T^3+E*T^4", 10
			return liquid_Thermal_Conductivity
		elif property_thermodynamics == "Vapour_Thermal_Conductivity": 
			Vapour_Thermal_Conductivity = "Vapour Thermal Conductivity", "[J/(m*s*K)]", "A*T^B/(1+C/T+D/T^2)", 11
			return Vapour_Thermal_Conductivity
		elif property_thermodynamics == "surface_Tension": 
			Surface_Tension = "Surface Tension", "[kg/s^2]", "A*(1-Tr)^(B+C*Tr+D*Tr^2)", 12	
			return Surface_Tension

	def property_cal(self, component, property_thermodynamics, temperature = None):

		self.property_label = self.select_property(property_thermodynamics)

		self.units = self.property_label[1]

		select_constans = [x + self.property_label[3] for x in range(0, 13*size_data, 13)]	
		values_constans = self.read_dppr().ix[select_constans, 1:8].get_values()
		self.table_constans = pd.DataFrame(data=values_constans,index=self.data_name_cal(),
							 columns=["A", "B", "C", "D", "E", "T Min [K]", "T Max [K]"])

		#print(table_constans)
		self.component_constans = self.table_constans.loc[component]
		#print(component_constans)

		#Min, Max = np.zeros(2), np.array(2)

		A, B, C, D, E, Min, Max = self.component_constans.get_values()

		#print("sss = ",Min, Max)

		if temperature == None:
			Temp_vector = np.array([Temp_vector for Temp_vector in np.arange(Min, Max)])
		else:			
			if type(temperature) != list: temperature = [temperature]		
			
			Temperature_enter = [i if Min < np.array(i) < Max
			 else "{0} K is a temperature not valid".format(i) for i in temperature]
			Temperature_invalid = [i for i in Temperature_enter if type(i) == str]
			Temperature_valid = [i for i in Temperature_enter if type(i) != str]

			print("Temperature_enter = {0}".format(Temperature_enter))
			print("Temperature_invalid = {0}".format(Temperature_invalid))
			print("Temperature_valid = {0}".format(Temperature_valid))
			
			Temp_vector = np.array(Temperature_valid)

		self.temperature = Temp_vector

		if property_thermodynamics == "Solid_Density":
			solid_Density = A + B * Temp_vector + C * Temp_vector ** 2 + D * Temp_vector ** 3 + E * Temp_vector **4		
			return solid_Density
		elif property_thermodynamics == "Liquid_Density":
			liquid_Density = A / B ** (1 + (1 - Temp_vector / C) ** D)		
			return liquid_Density
		elif property_thermodynamics == "Vapour_Pressure":
			vapour_Pressure = np.exp(A + B/Temp_vector + C * np.log(Temp_vector)+D*Temp_vector **E)		
			return vapour_Pressure
		elif property_thermodynamics == "Heat_of_Vaporization":
			heat_of_Vaporization = A*(1-Tr) ** (B+C*Tr+D*Tr**2)		
			return heat_of_Vaporization
		elif property_thermodynamics == "Solid_Heat_Capacity":
			solid_Heat_Capacity = A + B * Temp_vector + C * Temp_vector ** 2 + D * Temp_vector ** 3 + E * Temp_vector ** 4
			return solid_Heat_Capacity
		elif property_thermodynamics == "Liquid_Heat_Capacity":
			liquid_Heat_Capacity = A ** 2 / (1-Tr) + B-2*A*C*(1-Tr)-A*D*(1-Tr)**2-C**2*(1-Tr)**3/3-C*D*(1-Tr)**4/2-D**2*(1-Tr)**5/5
			return(liquid_Heat_Capacity)
		elif property_thermodynamics == "Ideal_Gas_Heat_Capacity":
			ideal_Gas_Heat_Capacity = A+B*(C/Temp_vector/np.sinh(C/Temp_vector))**2+D*(E/Temp_vector/np.cosh(E/Temp_vector))**2
			return ideal_Gas_Heat_Capacity
		elif property_thermodynamics == "Second_Virial_Coefficient":
			second_Virial_Coefficient = A+B/Temp_vector+C/Temp_vector**3+D/Temp_vector**8+E/Temp_vector**9
			return second_Virial_Coefficient
		elif property_thermodynamics == "Liquid_Viscosity":
			liquid_Viscosity = np.exp(A + B / Temp_vector + C * np.log(Temp_vector) + D * Temp_vector ** E)
			return liquid_Viscosity
		elif property_thermodynamics == "Vapour_Viscosity":
			vapour_Viscosity = A*Temp_vector**B/(1+C/Temp_vector+D/Temp_vector**2)
			return vapour_Viscosity
		elif property_thermodynamics == "Liquid_Thermal_Conductivity":
			liquid_Thermal_Conductivity = A+B*Temp_vector+C*Temp_vector**2+D*Temp_vector**3+E*Temp_vector**4
			return liquid_Thermal_Conductivity
		elif property_thermodynamics == "Vapour_Thermal_Conductivity":
			vapour_Thermal_Conductivity = A*Temp_vector**B/(1+C/Temp_vector+D/Temp_vector**2)
			return vapour_Thermal_Conductivity
		elif property_thermodynamics == "Surface_Tension":
			surface_Tension = A*(1-Tr)**(B+C*Tr+D*Tr*2)
			return surface_Tension

	


def main():

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

	#property_thermodynamics = thermodynamic_correlations.property_cal(component, "Vapour_Pressure", temp)
	property_thermodynamics = thermodynamic_correlations.property_cal(component, "Vapour_Pressure", temp)
	#property_thermodynamics = thermodynamic_correlations.property_cal(component, "Ideal_Gas_Heat_Capacity", temp)
	#property_thermodynamics = property_cal(components, Vapour_Pressure, temp)
	#property_thermodynamics = property_cal(component, Vapour_Pressure, [180.4, 181.4, 185.3, 210, 85])
	#property_thermodynamics = thermodynamic_correlations.property_cal(component, Vapour_Pressure)
	print(property_thermodynamics)
	print(thermodynamic_correlations.table_constans)
	print(thermodynamic_correlations.units)
	print(thermodynamic_correlations.temperature)

	print(thermodynamic_correlations.__doc__)

	print('-' * 79)

if __name__ == '__main__':
	main()
