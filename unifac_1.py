import numpy as np

U  = 5

# Este es el número de moléculas en la mezcla
m = 2    

# Este es el número de grupos funcionales en la mezcla
#g = 3   
g = 7

#T = 331.15  # K
T = 328
#     Etanol - n-Hexano
#xj = [0.332 , 0.668]
#xj = [0.383 , 0.617]

#--------------------------------------------------------
# Agua - Isoamil alcohol - ácido acético
#     H2O CH3 CH2 CH  OH  COOH  COOCH3
v1 = np.array([1, 0, 0, 0, 0, 0, 0]) # Agua
v2 = np.array([0, 2, 2, 1, 0, 0, 1]) # Isoamil acetato
v3 = np.array([0, 1, 0, 0, 0, 1, 0]) # Ácido acético

#v = np.array([v1, v2, v3])

v = np.array([v1, v3])

print(v)

# Agua - Isoamil acetato - ácido acético
#     H2O     CH3    CH2    CH     OH    COOH   COOCH3
R = np.array([0.9200, 0.9011, 0.6744, 0.4469, 1.0000, 1.3013, 1.9031])
Q = np.array([1.4000, 0.8480, 0.5400, 0.2280, 1.2000, 1.2240, 1.7280])

###################################################

# Agua - Isoamil alcohol - Ácido acético
# H2O     CH3     CH2     CH      OH      COOH    COOCH3
a = np.array([[0, 300, 300, 300, -229.1, -14.09, 72.8700],
			  [1318, 0, 0, 0, 986.5, 663.5, 232.100],
			  [1318, 0, 0, 0, 986.5, 663.5, 232.100],
			  [1318, 0, 0, 0, 986.5, 663.5, 232.100],
			  [353.5, 156.4, 156.4, 156.4, 0, 199, 101.100],
			  [-66.17, 315.3, 315.3, 315.3, -151, 0, -256.300],
			  [200.800, 114.800, 114.800, 114.800, 245.400, 660.200, 0]
			  ])

print(a)

A = np.exp(-a / T)

print(A) 