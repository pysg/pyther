import numpy as np



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


if property_thermodynamics == "Solid_Density":
	return solid_Density(T)
elif property_thermodynamics == "Liquid_Density":
	return liquid_Density(T)
elif property_thermodynamics == "Vapour_Pressure":
	return vapour_Pressure(T)
elif property_thermodynamics == "Heat_of_Vaporization":
	return heat_of_Vaporization(T, Tc)
elif property_thermodynamics == "Solid_Heat_Capacity":
	return solid_Heat_Capacity(T)
elif property_thermodynamics == "Liquid_Heat_Capacity":
	return liquid_Heat_Capacity(T, Tc)
elif property_thermodynamics == "Ideal_Gas_Heat_Capacity":
	return ideal_Gas_Heat_Capacity(T)
elif property_thermodynamics == "Second_Virial_Coefficient":
	return second_Virial_Coefficient(T)
elif property_thermodynamics == "Liquid_Viscosity":
	return liquid_Viscosity(T)
elif property_thermodynamics == "Vapour_Viscosity":
	return vapour_Viscosity(T)
elif property_thermodynamics == "Liquid_Thermal_Conductivity":
	return liquid_Thermal_Conductivity(T)
elif property_thermodynamics == "Vapour_Thermal_Conductivity":
	return vapour_Thermal_Conductivity(T)
elif property_thermodynamics == "Surface_Tension":
	return surface_Tension(T, Tc)
else:
	print("No hay propiedad termodinámica")

if Solid_Density:
	return solid_Density(T)
elif Liquid_Density:
	return liquid_Density(T)
elif Vapour_Pressure:
	return vapour_Pressure(T)
elif Heat_of_Vaporization:
	return heat_of_Vaporization(T, Tc)
elif Solid_Heat_Capacity:
	return solid_Heat_Capacity(T)
elif Liquid_Heat_Capacity:
	return liquid_Heat_Capacity(T, Tc)
elif Ideal_Gas_Heat_Capacity:
	return ideal_Gas_Heat_Capacity(T)
elif Second_Virial_Coefficient:
	return second_Virial_Coefficient(T)
elif Liquid_Viscosity:
	return liquid_Viscosity(T)
elif Vapour_Viscosity:
	return vapour_Viscosity(T)
elif Liquid_Thermal_Conductivity:
	return liquid_Thermal_Conductivity(T)
elif Vapour_Thermal_Conductivity:
	return vapour_Thermal_Conductivity(T)
elif Surface_Tension:
	return surface_Tension(T, Tc)
else:
	print("No hay propiedad termodinámica")
















