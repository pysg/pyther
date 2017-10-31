import numpy as np
import pandas as pd
# import math
import matplotlib.pyplot as plt
import os

# import pyther as pt

# constan to size of data base of components
# size_data = 30
size_data = 115
dfile = "PureFull_mod_properties.xls"
path_file = os.path.dirname(__file__)
dppr_file = os.path.join(path_file, dfile)


class Thermodynamic_correlations(object):

	def read_dppr(self):
		
		#-self.dppr_data = pd.read_excel(self.dppr_file).ix[:, 0:8]
		self.dppr_data = pd.read_excel(dppr_file).ix[:, 0:8]
		
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
			Vapour_Pressure = "Vapour Pressure", "[Bar]", "exp(A+B/T+C*ln(T)+D*T^E)* 1e-5", 2
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

	def select_constans_cal(self, component, property_thermodynamics):

		self.property_label = self.select_property(property_thermodynamics)		

		select_constans = [x + self.property_label[3] for x in range(0, 13*size_data, 13)]	
		values_constans = self.read_dppr().ix[select_constans, 1:8].get_values()
		self.table_constans = pd.DataFrame(data=values_constans,index=self.data_name_cal(),
							 columns=["A", "B", "C", "D", "E", "T Min [K]", "T Max [K]"])

		self.component_constans = self.table_constans.loc[component]

		A, B, C, D, E, Min, Max = self.component_constans.get_values()

		self.Tmin = Min
		self.Tmax = Max

		range_temp = np.array([self.Tmin, self.Tmax])

		return range_temp

	
	def graphical(self, t, p, pp, u):
	    plt.plot(t, p)
	    plt.title(pp)
	    plt.xlabel(u[0])
	    plt.ylabel(u[1])

	    return

	def multi_graphical(self, components, temperature, property_thermodynamics):
    
	    fig, axes = plt.subplots(1, len(components), figsize=(10, 4))
	    for pure, label in zip(range(len(components)), components):
	        axes[pure].plot(temperature[pure], property_thermodynamics[pure])
	        axes[pure].set_title(label)
	        axes[pure].set_xlabel("Temperature " + self.units[0])
	        axes[pure].set_ylabel("Pressure " + self.units[1])


	def data_temperature(self, components, temperature, property_thermodynamics, temperature_enter):
    
	    datas = np.ones([len(components),len(temperature)])
	    
	    for k in range(len(components)):
	        con1 = [type(temperature_enter[k][i])!= str for i in range(len(temperature))]
	        vap1 = property_thermodynamics[k]
	        
	        if len(vap1) == 1:
	            datas[k] = [vap1 if con1[i]== True else None for i in range(len(temperature))]
	        else:
	            datas[k] = [vap1[i] if con1[i]== True else None for i in range(len(temperature))]

	    t = [str(temp)+" K" for temp in temperature]
	    datas_table = pd.DataFrame(data=datas,index= components, columns=t)

	        
	    return datas_table


	def control_temperature(self, components, temperature, Min, Max):

		if temperature == None:

			if len(components) == 1:
				print("-"*70)
				print("Pure substance without temperature especific: {0}". format(self.components))
				print("-"*70)	
				Temp_vector = np.arange(Min, Max)
				#Temp_vector = np.array([Temp for Temp in np.arange(Min, Max)])			
			else:
				Temp_vector = np.array([ np.array([Temp_vector for Temp_vector in np.arange(Min[i], Max[i])]) for i in range(0, len(components))])				
		else:
			if len(components) == 1:
				print("-"*70)
				print("Pure substance with a temperature especific: {0}". format(self.components))
				print("-"*70)		
				if type(temperature) != list: temperature = [temperature]		
				
				self.temperature_enter = [i if Min < np.array(i) < Max
				 else "{0} K is a temperature not valid".format(i) for i in temperature]

				self.temperature_invalid = [i for i in self.temperature_enter if type(i) == str]
				self.temperature_valid = [i for i in self.temperature_enter if type(i) != str]

				print("Temperature_enter = {0}".format(self.temperature_enter))
				print("Temperature_invalid = {0}".format(self.temperature_invalid))
				print("Temperature_valid = {0}".format(self.temperature_valid))
				print("-"*70)
				Temp_vector = self.temperature_valid
			else:
				print("-"*70)
				print("Pure substances with a temperature especific: {0}". format(self.components))
				print("-"*70)	

				self.temperature_enter = [([i if Min[j] < i < Max[j]
				 else "{0} K is a temperature not valid".format(i) for i in temperature]) for j in range(0, len(components))]

				print(temperature)

				self.temperature_invalid = [([i for i in self.temperature_enter[j] if type(i) == str]) for j in range(0, len(components))]
				self.temperature_valid = [([ i for i in self.temperature_enter[j] if type(i) != str]) for j in range(0, len(components))]

				self.temperature_valid = [ np.array([ i for i in self.temperature_enter[j] if type(i) != str]) for j in range(0, len(components))]

				print("Temperature_enter = {0}".format(self.temperature_enter))
				print("Temperature_invalid = {0}".format(self.temperature_invalid))
				print("Temperature_valid = {0}".format(self.temperature_valid))
				
				Temp_vector = np.array(self.temperature_valid)
				Temp_vector = np.reshape(Temp_vector, len(components))

				#print("temperature = ", Temp_vector, len(Temp_vector), np.shape(Temp_vector) )


		self.temperature = Temp_vector

		#print("temperature = ", self.temperature, len(self.temperature) )

		return self.temperature	


	def property_cal(self, components, property_thermodynamics, temperature = None):

		self.property_label = self.select_property(property_thermodynamics)
		
		#print(self.property_label)
		
		self.units = ("K", self.property_label[1])
		
		self.components = components

		#print(self.components, type(components))

		select_constans = [x + self.property_label[3] for x in range(0, 13*size_data, 13)]	
		values_constans = self.read_dppr().ix[select_constans, 1:8].get_values()
		self.table_constans = pd.DataFrame(data=values_constans,index=self.data_name_cal(),
							 columns=["A", "B", "C", "D", "E", "T Min [K]", "T Max [K]"])

		self.component_constans = self.table_constans.loc[components]
	
		if len(self.components) == 1:
			#A, B, C, D, E, Min, Max = self.component_constans.get_values()
			A = np.array(self.component_constans["A"].get_values())
			B = np.array(self.component_constans["B"].get_values())
			C = np.array(self.component_constans["C"].get_values())
			D = np.array(self.component_constans["D"].get_values())
			E = np.array(self.component_constans["E"].get_values())
			Min = np.array(self.component_constans["T Min [K]"].get_values())
			Max = np.array(self.component_constans["T Max [K]"].get_values())
		else:

			A = np.array(self.component_constans["A"].get_values())
			B = np.array(self.component_constans["B"].get_values())
			C = np.array(self.component_constans["C"].get_values())
			D = np.array(self.component_constans["D"].get_values())
			E = np.array(self.component_constans["E"].get_values())
			Min = np.array(self.component_constans["T Min [K]"].get_values())
			Max = np.array(self.component_constans["T Max [K]"].get_values())

		self.Min = Min
		self.Max = Max

		Temp_vector = self.control_temperature(components, temperature, Min, Max)

		if property_thermodynamics == "Solid_Density":
			solid_Density = A + B * Temp_vector + C * Temp_vector ** 2 + D * Temp_vector ** 3 + E * Temp_vector **4
			return solid_Density
		elif property_thermodynamics == "Liquid_Density":
			liquid_Density = A / B ** (1 + (1 - Temp_vector / C) ** D)		
			return liquid_Density
		elif property_thermodynamics == "Vapour_Pressure":
			# used 1 Pa = 1e-5 Bar

			if temperature == None:
				# without temperature especific

				if len(self.components) == 1:
					# one component without one temperature especific
					
					log_tem = np.log(np.float64(Temp_vector))
					log_vapour_Pressure = A + B/Temp_vector + C * log_tem + D*Temp_vector **E
					vapour_Pressure = np.exp(np.float64(log_vapour_Pressure)) * 1e-5				
					
				else:
					# many components 
					log_tem = np.array([np.log(Temp_vector) for Temp_vector in Temp_vector])
					vapour = np.array(A + B/Temp_vector + C * log_tem + D*Temp_vector **E)
					vapour_Pressure = np.array([np.exp(vapour) for vapour in vapour]) * 1e-5				
					#print("vapour_Pressure = ", vapour_Pressure, np.shape(vapour_Pressure))
				
			else:

				if len(self.components) == 1:
					# one component with one temperature especific 
					
					log_tem = np.array([np.log(Temp_vector) for Temp_vector in Temp_vector])
					log_vapour_Pressure = A + B/Temp_vector + C * log_tem + D*Temp_vector **E
					#print("log_vapour-pressure = ",log_vapour_Pressure)
					vapour_Pressure = np.array([np.exp(vapour) for vapour in log_vapour_Pressure]) * 1e-5
				else:
					# many components

					if len(temperature) == 1:
						# many components with one temperature especific

						log_tem = np.array([np.log(Temp_vector) for Temp_vector in Temp_vector])
						vapour = np.array(A + B/Temp_vector.T + C * log_tem.T + D*Temp_vector.T **E)
						# used 1 Pa = 1e-5 Bar
						vapour_Pressure = np.array([np.exp(vapour) for vapour in vapour]) * 1e-5				
						print("vapour_Pressure = ", vapour_Pressure, np.shape(vapour_Pressure))
						
					else:
						# many components with many temperatures especific

						log_tem = [np.log(Temp_vector) for Temp_vector in Temp_vector]
						log_vapour_Pressure = A + B/Temp_vector + C * log_tem + D*Temp_vector **E
						#print("log_vapour-pressure = ",log_vapour_Pressure)
						vapour_Pressure = np.array([np.exp(vapour) for vapour in log_vapour_Pressure]) * 1e-5

						


					print(np.size(vapour_Pressure))

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

	def table_properties(self):
		table_components = pd.DataFrame(data=self.property_cal(self.components, self.property_label),index= self.components, 
							 columns=[str(self.temperature)+"K"])
		return table_components


def main():

	print("-" * 79)

	thermodynamic_correlations = Thermodynamic_correlations()
	components = ["METHANE", "n-TETRACOSANE", "n-PENTACOSANE", "ETHANE", "ISOBUTANE", "PROPANE", "3-METHYLHEPTANE"]
	temp = [180.4, 181.4, 185.3, 210, 800]

	property_thermodynamics = thermodynamic_correlations.property_cal(components, "Vapour_Pressure")
	table_components = pd.DataFrame(data=property_thermodynamics,index= components, 
							 columns=[str(temp)+"K"])
	
	print(table_components)
	print('-' * 79)

if __name__ == '__main__':
	main()






def solid_Density(T):
	solid_Density = A + B * T + C * T ** 2 + D * T ** 3 + E * T ** 4
	return solid_Density

def liquid_Density(T):
	liquid_Density = A / B ** (1 + (1 - T / C) ** D)
	return liquid_Density

def vapour_Pressure(T):
	# used 1 Pa = 1e-5 Bar
	vapour_Pressure = np.exp(A + B / T + C * log_tem + D * T ** E) * 1e-5
	return vapour_Pressure

def heat_of_Vaporization(T, Tc):
	Tr = T / Tc
	heat_of_Vaporization = A * (1 - Tr) ** (B + C * Tr + D * Tr ** 2)
	pass

def solid_Heat_Capacity(T):
	solid_Heat_Capacity = A + B * T + C * T ** 2 + D * T ** 3 + E * T ** 4
	return solid_Heat_Capacity

def liquid_Heat_Capacity(T, Tc):
	Tr = T / Tc

	variable1 = A ** 2 / (1 - Tr) + B - 2 * A * C * (1 - Tr)
	variable2 = - A * D * (1 - Tr) ** 2 -C ** 2 * (1 - Tr) ** 3 / 3
	variable3 = - C * D * (1-Tr) ** 4 / 2 - D ** 2 * (1 - Tr) ** 5 / 5

	liquid_Heat_Capacity = variable1 + variable2 + variable3
	return liquid_Heat_Capacity

def ideal_Gas_Heat_Capacity(T):
	ideal_Gas_Heat_Capacity = A + B * (C / T/np.sinh(C / T)) ** 2 + D * (E/T/np.cosh(E / T)) ** 2
	return ideal_Gas_Heat_Capacity

def second_Virial_Coefficient(T):
	second_Virial_Coefficient = A + B / T + C / T ** 3 + D / T ** 8 + E / T ** 9
	return second_Virial_Coefficient

def liquid_Viscosity(T):
	liquid_Viscosity = np.exp(A + B / T + C * np.log(T) + D * T ** E)
	return liquid_Viscosity

def vapour_Viscosity(T):
	vapour_Viscosity = A * T ** B / (1 + C / T + D / T ** 2)
	return vapour_Viscosity

def liquid_Thermal_Conductivity(T):
	liquid_Thermal_Conductivity = A + B * T + C * T ** 2 + D * T ** 3 + E * T ** 4
	return liquid_Thermal_Conductivity

def vapour_Thermal_Conductivity(T):
	vapour_Thermal_Conductivity = A * T ** B / (1 + C / T + D / T ** 2)
	return vapour_Thermal_Conductivity

def surface_Tension(T, Tc):
	Tr = T / Tc
	surface_Tension = A * (1 - Tr) ** (B + C * Tr + D * Tr * 2)
	return surface_Tension


# if property_thermodynamics == "Solid_Density":
# 	return solid_Density(T)
# elif property_thermodynamics == "Liquid_Density":
# 	return liquid_Density(T)
# elif property_thermodynamics == "Vapour_Pressure":
# 	return vapour_Pressure(T)
# elif property_thermodynamics == "Heat_of_Vaporization":
# 	return heat_of_Vaporization(T, Tc)
# elif property_thermodynamics == "Solid_Heat_Capacity":
# 	return solid_Heat_Capacity(T)
# elif property_thermodynamics == "Liquid_Heat_Capacity":
# 	return liquid_Heat_Capacity(T, Tc)
# elif property_thermodynamics == "Ideal_Gas_Heat_Capacity":
# 	return ideal_Gas_Heat_Capacity(T)
# elif property_thermodynamics == "Second_Virial_Coefficient":
# 	return second_Virial_Coefficient(T)
# elif property_thermodynamics == "Liquid_Viscosity":
# 	return liquid_Viscosity(T)
# elif property_thermodynamics == "Vapour_Viscosity":
# 	return vapour_Viscosity(T)
# elif property_thermodynamics == "Liquid_Thermal_Conductivity":
# 	return liquid_Thermal_Conductivity(T)
# elif property_thermodynamics == "Vapour_Thermal_Conductivity":
# 	return vapour_Thermal_Conductivity(T)
# elif property_thermodynamics == "Surface_Tension":
# 	return surface_Tension(T, Tc)
# else:
# 	print("No hay propiedad termodinámica")

# if Solid_Density:
# 	return solid_Density(T)
# elif Liquid_Density:
# 	return liquid_Density(T)
# elif Vapour_Pressure:
# 	return vapour_Pressure(T)
# elif Heat_of_Vaporization:
# 	return heat_of_Vaporization(T, Tc)
# elif Solid_Heat_Capacity:
# 	return solid_Heat_Capacity(T)
# elif Liquid_Heat_Capacity:
# 	return liquid_Heat_Capacity(T, Tc)
# elif Ideal_Gas_Heat_Capacity:
# 	return ideal_Gas_Heat_Capacity(T)
# elif Second_Virial_Coefficient:
# 	return second_Virial_Coefficient(T)
# elif Liquid_Viscosity:
# 	return liquid_Viscosity(T)
# elif Vapour_Viscosity:
# 	return vapour_Viscosity(T)
# elif Liquid_Thermal_Conductivity:
# 	return liquid_Thermal_Conductivity(T)
# elif Vapour_Thermal_Conductivity:
# 	return vapour_Thermal_Conductivity(T)
# elif Surface_Tension:
# 	return surface_Tension(T, Tc)
# else:
# 	print("No hay propiedad termodinámica")


if dppr_file:
	print(dppr_file)









