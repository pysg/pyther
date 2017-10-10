

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
