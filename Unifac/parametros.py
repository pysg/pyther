import numpy as np
import pandas as pd


def crearMezcla():

	# --------------------------------------------------------
	# Agua - Isoamil alcohol - ácido acético
	#     H2O CH3 CH2 CH  OH  COOH  COOCH3
	agua = {"agua" : np.array([1, 0, 0, 0, 0, 0, 0])}
	isoamil_acetato = {"isoamil_acetato" : np.array([0, 2, 2, 1, 0, 0, 1])}
	acido_acetico = {"acido_acetico" : np.array([0, 1, 0, 0, 0, 1, 0])}
	
	mezcla = np.array([agua, acido_acetico])

	return mezcla

mezcla = crearMezcla()
print("Mezcla = ", mezcla[0]["agua"])

# Agua - Isoamil acetato - ácido acético
#               H2O     CH3    CH2    CH     OH    COOH   COOCH3
R = np.array([0.9200, 0.9011, 0.6744, 0.4469, 1.0000, 1.3013, 1.9031])
Q = np.array([1.4000, 0.8480, 0.5400, 0.2280, 1.2000, 1.2240, 1.7280])

print("R = ", R)
filasRQ = ["R", "Q"]
labelsRQ = ["H2O", "CH3", "CH2", "CH", "OH", "COOH", "COOCH3"]

dfRQ = pd.DataFrame(data = [R, Q], index=filasRQ, columns=labelsRQ)
print(dfRQ)


# --------------------------------------------------------

# Agua - Isoamil alcohol - Ácido acético
# H2O     CH3     CH2     CH      OH      COOH    COOCH3
labels = ["H2O", "CH3", "CH2", "CH", "OH", "COOH", "COOCH3"]
a = np.array([[0, 300, 300, 300, -229.1, -14.09, 72.8700],
				[1318, 0, 0, 0, 986.5, 663.5, 232.100],
				[1318, 0, 0, 0, 986.5, 663.5, 232.100],
				[1318, 0, 0, 0, 986.5, 663.5, 232.100],
				[353.5, 156.4, 156.4, 156.4, 0, 199, 101.100],
				[-66.17, 315.3, 315.3, 315.3, -151, 0, -256.300],
				[200.800, 114.800, 114.800, 114.800, 245.400, 660.200, 0]])

df = pd.DataFrame(a, index=labels, columns=labels)
print(df)









