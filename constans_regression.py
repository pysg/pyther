# encoding: UTF-8
import numpy as np

A0, B0, C0 = 0.0017, 1.9681, -2.7238
# Definir el significado fisicoquímico

A1, B1, C1 = -2.4407, 7.4513, 12.504
# Definir el significado fisicoquímico

D = np.array([0.428363, 18.496215, 0.338426, 0.660, 789.723105, 2.512392])

er = 2
if 0 < er < 5:
	print(42)

print(0 < er < 5)

NMODEL = "constants_eps"
ICALC = "RKPRa"

if ICALC in ["SRK", "PR", "RKPR"]:
	print(0)

if NMODEL in ["constants_eps", "parameters_eps", "rk_param", "density"]:
	print(99)


class Circulo(object):

	def __init__(self, radio):
		self.radio = radio

	@property
	def radio(self):
		return self._radio

	@radio.setter
	def radio(self, radio):
		if radio > 0:
			self._radio = radio
		else:
			self._radio = 0.5


c = Circulo(0)
print ("Radio:", c.radio)


class Eos(object):

	def __init__(self, NMODEL):
		self.NMODEL = NMODEL

	@property
	def NMODEL(self):
		return self._NMODEL

	@NMODEL.setter
	def radio(self, NMODEL):
		if NMODEL > 0:
			self._NMODEL = NMODEL
		self._NMODEL = 0.5





class CajaFuerte(object):

	def __init__(self, PIN):
		self.PIN = PIN

	@property
	def PIN(self):
		print ("Enviando copia a la NSA...")
		return self._PIN

	@PIN.setter
	def PIN(self, PIN):
		if len(str(PIN)) != 4:
			raise ValueError("’PIN’ ha de tener cuatro d´ıgitos")
		self._PIN = PIN

hucha = CajaFuerte(781)
print ("PIN:", hucha.PIN)
#hucha.PIN = 880

