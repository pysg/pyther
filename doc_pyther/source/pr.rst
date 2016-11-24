********************
Calculos flash (P,T)
********************

De esta forma, para trazar la envolvente de fase se utiliza la estrategia seguida por Mollerup y Mishelsen en su libro Thermodynamics Aspecs calations.



El balance de materia global es

.. math:: F = L + V
	:label:

El balance de materia por componente es

.. math:: Fz_i = Lx_i + Vy_i
	:label:

Las realciones de equilibrio de fases para C componentes son

.. math:: y_i = K_i x_i
	:label:

para

.. math:: i = 1, 2, ..., C

La primera clase es::

	import numpy as np
	class Termodinamica():
		def __ini__(P, T, V):
			self.P = P
			self.T = T
			self.V = V
			print "Temperatura = ", self.T

		return self.P ** 2


Modelo
******

El modelo del flash isotermico bifasico, corresponde al balance de materia global y por componente en el tanque separador que se muestra en la figura (), junto con la condición de equilibrio de fases líquido-vapor.

.. math:: Ki = \frac {yi} {xi}
	:label: