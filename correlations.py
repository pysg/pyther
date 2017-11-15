import numpy as np
import constantes

class Correlations(object):
	"""docstring for Correlations"""
	def __init__(self, constantes, T, Tc):		
		
		self.A = constantes[0]
		self.B = constantes[1]
		self.C = constantes[2]
		self.D = constantes[3]
		self.E = constantes[4]
		self.T = T
		self.Tc = Tc
		print(self.A)

	def solidDensity(self):

		VAR1 = self.A + self.B * self.T + self.C * self.T ** 2 
		VAR2 = self.D * self.T ** 3 + self.E * self.T ** 4
		
		#solidDensity = self.A + self.B * self.T + self.C * self.T ** 2 + self.D * self.T ** 3 + self.E * self.T ** 4
		solidDensity = VAR1 + VAR2

		return solidDensity

	def liquidDensity(self):

		liquidDensity = self.A / self.B ** (1 + (1 - self.T / self.C) ** self.D)

		return liquidDensity

	def presionVapor(self):

		VAR1 = self.A + self.B / self.T
		VAR2 = self.C * np.log(self.T) + self.D * self.T ** self.E

		vaporCalculado = np.exp(VAR1 + VAR2) * 1e-5

		#vaporCalculado = np.exp(self.A + self.B / self.T + self.C * np.log(self.T) + self.D * self.T ** self.E) * 1e-5

		return vaporCalculado

	def heat_of_Vaporization(self):

		Tr = self.T / self.Tc

		heat_of_Vaporization = self.A * (1 - Tr) ** (self.B + self.C * Tr + self.D * Tr ** 2)		

		return heat_of_Vaporization

	def solid_Heat_Capacity(self):

		VAR1 = self.A + self.B * self.T + self.C * self.T ** 2
		VAR2 = self.D * self.T ** 3 + self.E * self.T ** 4

		solid_Heat_Capacity = VAR1 + VAR2

		#solid_Heat_Capacity = self.A + self.B * self.T + self.C * self.T ** 2 + self.D * self.T ** 3 + self.E * self.T ** 4

		return solid_Heat_Capacity

	def liquid_Heat_Capacity(self):

		Tr = self.T / self.Tc

		VAR1 = self.A ** 2 / (1 - Tr)
		VAR2 = self.B - 2 * self.A * self.C * (1 - Tr)
		VAR2 = - self.A * self.D * (1 - Tr) ** 2
		VAR3 = - self.C ** 2 * (1 - Tr) ** 3 / 3
		VAR4 = - self.C * self.D * (1 - Tr) ** 4 / 2
		VAR5 = - self.D ** 2 * (1 - Tr) ** 5 / 5

		liquid_Heat_Capacity = VAR1 + VAR2 + VAR3 + VAR4 + VAR5
		
		return liquid_Heat_Capacity

	def ideal_Gas_Heat_Capacity(self):

		VAR1 = self.B * (self.C / self.T / np.sinh(self.C / self.T)) ** 2
		VAR2 = self.D * (self.E / self.T / np.cosh(self.E / self.T)) ** 2

		ideal_Gas_Heat_Capacity = self.A + VAR1 + VAR2

		return ideal_Gas_Heat_Capacity

	def second_Virial_Coefficient(self):

		VAR1 = self.B / self.T + self.C / self.T ** 3
		VAR2 = self.D / self.T ** 8
		VAR3 = self.E / self.T ** 9

		second_Virial_Coefficient = self.A + VAR1 + VAR2 + VAR3

		#second_Virial_Coefficient = self.A + self.B / self.T + self.C / self.T ** 3 + self.D / self.T ** 8 + self.E / self.T ** 9

		return second_Virial_Coefficient

	def liquid_Viscosity(self):

		liquid_Viscosity = np.exp(self.A + self.B / self.T + self.C * np.log(self.T) + self.D * self.T ** self.E)

		return liquid_Viscosity

	def vapour_Viscosity(self):

		vapour_Viscosity = self.A * self.T ** self.B / (1 + self.C / self.T + self.D / self.T ** 2)

		return vapour_Viscosity

	def liquid_Thermal_Conductivity(self):

		VAR1 = self.A + self.B * self.T
		VAR2 = self.C * self.T ** 2
		VAR3 = self.D * self.T ** 3
		VAR4 = self.E * self.T ** 4

		#liquid_Thermal_Conductivity = self.A + self.B * self.T + self.C * self.T ** 2 + self.D * self.T ** 3 + self.E * self.T ** 4

		liquid_Thermal_Conductivity = VAR1 + VAR2 + VAR3 + VAR4

		return liquid_Thermal_Conductivity

	def vapour_Thermal_Conductivity(self):

		vapour_Thermal_Conductivity = self.A * self.T ** self.B / (1 + self.C / self.T + self.D / self.T ** 2)

		return vapour_Thermal_Conductivity

	def surface_Tension(self):

		Tr = self.T / self.Tc

		surface_Tension = self.A * (1 - Tr) ** (self.B + self.C * Tr + self.D * Tr * 2)

		return surface_Tension


#vapour_Pressure = np.array(list(map(presonVapor, Temp_vector)))

correlations = Correlations()
presionVapor = correlations.presionVapor()

CONSTANTES = np.array([A, B, C, D, E]).T
print(CONSTANTES)
vapour_Pressure = np.array([np.array([presionVapor(constantes, T) for T in Temp_vector]) for constantes in CONSTANTES])

print(vapour_Pressure)
