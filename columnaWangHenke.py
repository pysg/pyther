import numpy as np
from scipy import optimize
from pyther import Data_parse
from sympy import *

# from .constants import RGAS
RGAS = 0.08314472
R = RGAS

NUMERO_COMPONENTES = 3
ETAPAS = 5
ETAPA_ALIMENTACION = 2

Z = np.zeros([NUMERO_COMPONENTES, ETAPAS])
Zf = np.array([0.3, 0.3, 0.4])
Z[:, ETAPA_ALIMENTACION] = Zf
print(Z)

Fm = np.array([0, 0, 100, 100, 100])
Fj = np.array([0, 0, 100, 0, 0])
Um = np.array([50, 50, 50, 50, 50])
Uj = np.array([50, 0, 0, 0, 0])
Wj = np.zeros(5)

L5 = sum(Fj) - sum(Uj)
# L1 = 2*Uj(1,1)
L1 = 2 * Uj[0]

Vj = np.array([0, 150.0, 150, 150, 150])
# Vj[0] = 0

Lj = np.array([L1, 50, 50, 50, L5])
Vj[1] = Lj[0] + Uj[0]

Te = np.array([350, 355, 360, 365, 370])

P = 1

components = ["ISOPROPANOL", "3-METHYL-1-BUTANOL", "1-BUTANOL"]

properties_data = Data_parse()
properties_component = properties_data.selec_component(components)

constans = properties_component[1].loc[:, ["Omega", "Tc", "Pc"]]

w = np.float64(constans["Omega"])
Tc = np.float64(constans["Tc"])
Pc = np.float64(constans["Pc"])


def Ki_wilson(T, P):
    """Equation of wilson for to calculate the Ki(T,P)"""
    variable_0 = 5.373 * (1 + w) * (1 - Tc / T)
    lnKi = np.log(Pc / P) + variable_0
    Ki = np.exp(lnKi)

    return Ki


def fugacity(T, P, yi, xi):
    m = 0.48 + (1.574 * w) - (0.176 * w ** 2)
    Tr = T / Tc
    alpha = (1 + m * (1 - (Tr ** 0.5))) ** 2
    ac = 0.42748 * (R * Tc) ** 2 / Pc
    a = ac * alpha
    b = 0.08664 * R * Tc / Pc

    # Vapor
    amv = np.sum(yi * a ** 0.5) ** 2
    bmv = np.sum(yi * b)
    Av = (amv * P) / ((R * T) ** 2)
    Bv = (bmv * P) / (R * T)
    Zv = np.max(np.roots([1, -1, (Av - Bv - Bv ** 2), (- Av * Bv)]))

    aav = (a / amv)
    bbv = (b / bmv)

    # Liquid
    aml = np.sum(xi * a ** 0.5) ** 2
    bml = np.sum(xi * b)
    Al = (aml * P) / ((R * T) ** 2)
    Bl = (bml * P) / (R * T)
    Zl = np.min(np.roots([1, -1, (Al - Bl - Bl ** 2), (- Al * Bl)]))

    aal = (a / aml)
    bbl = (b / bml)

    # Fugacity Coefficient
    factor_1 = (bbv - (2 * (aav ** 0.5))) * np.log((Zv + Bv) / Zv)
    ln_phi_v = bbv * (Zv - 1) - np.log(Zv - Bv) + (Av / Bv) * factor_1
    phi_v = np.exp(ln_phi_v)

    factor_2 = (bbl - (2 * (aal ** 0.5))) * np.log((Zl + Bl) / Zl)
    ln_phi_l = bbl * (Zl - 1) - np.log(Zl - Bl) + (Al / Bl) * factor_2
    phi_l = np.exp(ln_phi_l)

    return phi_l, phi_v


Tij = np.array([350, 355, 360, 365, 370])

Kij = np.array([Ki_wilson(Tj, P) for Tj in Tij])
Kij = np.transpose(Kij)
Kij_N = np.zeros([3, 1])
Kij = np.append(Kij, Kij_N, axis=1)

Aj = np.zeros([NUMERO_COMPONENTES, ETAPAS])
Bj = np.zeros([NUMERO_COMPONENTES, ETAPAS])
Cj = np.zeros([NUMERO_COMPONENTES, ETAPAS])
Dj = np.zeros([NUMERO_COMPONENTES, ETAPAS])
pj = np.zeros([NUMERO_COMPONENTES, ETAPAS])
qj = np.zeros([NUMERO_COMPONENTES, ETAPAS])
Xj = np.zeros([NUMERO_COMPONENTES, ETAPAS])


def fraccionesMolares(ETAPAS, NUMERO_COMPONENTES, Vj, Kij):
    Vj = np.append(Vj, 0)
    for j in range(ETAPAS):
        for i in range(NUMERO_COMPONENTES):
            # Calculo de los Coeficientes de la Matriz Tridiagonal
            Aj[i, j] = Vj[j] + (Fm[j - 1] - Um[i - 1]) - Vj[0]
            Bj[i, j] = -(Vj[j + 1] + (Fm[j] - Um[j]) + Uj[j] + Vj[j] * Kij[i, j])
            Cj[i, j] = Vj[j + 1] * Kij[i, j + 1]
            Dj[i, j] = -(Fj[j] * Z[i, j])

            # Calculo de las variable que se reemplazan en la Matriz Tridiagonal
            pj[i, j] = Cj[i, j] / (Bj[i, j] - Aj[i, j] * pj[i, j - 1])
            qj[i, j] = (Dj[i, j] - Aj[i, j] * qj[i, j - 1]) / (Bj[i, j] - Aj[i, j] * pj[i, j - 1])

    # Calculo de las fracciones molares del liquido
    Xj[:, -1] = qj[:, -1]
    for j in range(ETAPAS - 2, -1, -1):
        for i in range(NUMERO_COMPONENTES):
            Xj[i, j] = qj[i, j] - pj[i, j] * Xj[i, j + 1]

    return Xj


def normalizar(Xj):
    SXj = sum(Xj)
    Xj = Xj / SXj
    return Xj


fraccionesMolares(ETAPAS, NUMERO_COMPONENTES, Vj, Kij)
Xj = normalizar(fraccionesMolares(ETAPAS, NUMERO_COMPONENTES, Vj, Kij))


def equilibrio(variables, P, Xi):

    T = variables[0]
    Yi[0] = variables[1]
    Yi[1] = variables[2]
    Yi[2] = 1 - Yi[0] - Yi[1]
    Xi[2] = 1 - Xi[0] - Xi[1]

    Fi = fugacity(T, P, Yi, Xi)
    # Ki = Oil / Oiv
    Ki = Fi[0] / Fi[1]

    equilibrium = Yi - Ki * Xi

    return equilibrium


hFj = 2 * np.zeros(5)
hFj[2] = 348150
hFj


# ISOPROPANOL
# 106
CONSTANTES_VAP1 = np.array([6.31E07, 3.92E-01, 0.00E+00, 0.00E+00, 0.00E+00, 185.28, 508.3])
# 107
CONSTANTES_1 = np.array([5.72E+04, 1.91E+05, 1.42E+03, 1.22E+05, 6.26E+02, 150, 1500])

# 3-METHYL-1-BUTANOL
# 106
CONSTANTES_VAP2 = np.array([8.08E+07, 5.02E-01, 0.00E+00, 0.00E+00, 0.00E+00, 155.95, 577.2])
# 107
CONSTANTES_2 = np.array([1.11E+05, 2.21E+05, 8.76E+02, 1.22E+05, 2.94E+03, 298.15, 1200.15, 107])

# 1-BUTANOL
# 106
CONSTANTES_VAP3 = np.array([6.74E+07, 1.73E-01, 2.92E-01, 0.00E+00, 0.00E+00, 184.51, 563.05]) 
# 107
CONSTANTES_3 = np.array([7.45E+04, 2.59E+05, 1.61E+03, 1.73E+05, 7.12E+02, 200, 1500]) 

CONSTANTESVAP = np.array([CONSTANTES_VAP1, CONSTANTES_VAP2, CONSTANTES_VAP3])
CONSTANTES = np.array([CONSTANTES_1, CONSTANTES_2, CONSTANTES_3])


def entalpiaGasIdeal(constantes, T):
    # Cp = A + B * (C / T /sinh(C / T)) ** 2 + D * (E / T / cosh(E / T)) ** 2

    A = constantes[0]
    B = constantes[1]
    C = constantes[2]
    D = constantes[3]
    E = constantes[4]

    termino_1 = (1 / (2 * T)) * tanh(C / (2 * T))
    termino_2 = ((2 * tanh(E / (2 * T))) / (E * (tanh(E / (2 * T)) ** 2) + E))
    termino_3 = - D * E ** 2

    H = A * T + B * C ** 2 * termino_1 + termino_2 * termino_3
    return H


def entalpiaVaporizacion(constantes, Tc, T):

    A = constantes[0]
    B = constantes[1]
    C = constantes[2]
    D = constantes[3]

    Tr = T / Tc
    Hvap = A * (1 - Tr) ** (B + C * Tr + D * Tr ** 2)
    return Hvap


Tci = Tc


def entalpiaLiquidoIdeal(H, Hvap):
    h = H - Hvap
    return h


Qj = np.zeros(ETAPAS)
alfa = np.zeros(ETAPAS)
beta = np.zeros(ETAPAS)
gama = np.zeros(ETAPAS)


def primero(Hj):

    Hj = np.append(Hj, 0)

    for j in range(1, ETAPAS):

        alfa[j] = hj[j - 1] - Hj[j]
        beta[j] = Hj[j + 1] - hj[j]

        factor0 = sum(Fj[:j - 1] - Uj[:j - 1] - Wj[:j - 1] - Vj[0])
        factor1 = (hj[j] - hj[0])
        factor2 = Fj[j] * (hj[j] - hFj[j])
        factor3 = Wj[j] * (Hj[j] - hj[j]) + Qj[j]

        gama[j] = factor0 - factor1 + factor2 + factor3

    return alfa, beta, gama


def segundo(primer, Vj):

    alfa = primer[0]
    beta = primer[1]
    gama = primer[2]

    for j in range(2, ETAPAS):
        Vj[j] = (gama[j - 1] - alfa[j - 1] * Vj[j - 1]) / beta[j - 1]

    return Vj


def columna(Vj, Kij):

    i = 0
    TjInit1 = np.array([350, 355, 360, 365, 370])

    while True:
        i += 1
        print(i)

        fraccionesMolares(ETAPAS, NUMERO_COMPONENTES, Vj, Kij)
        Xj = normalizar(fraccionesMolares(ETAPAS, NUMERO_COMPONENTES, Vj, Kij))

        variables = np.array([378, 0.3, 0.3])
        Eqj = np.array([optimize.root(equilibrio, variables, args=(P, xj)) for xj in Xj.T])

        for n in range(ETAPAS):
            Teqj[n] = Eqj[n].x[0]
            yij[0, n] = Eqj[n].x[1]
            yij[1, n] = Eqj[n].x[2]

        print("TjInit1 = ", TjInit1)
        print("Teqj = ", Teqj)

        yij[2, :] = 1 - yij[0] - yij[1]
        Yj = yij
        Kij = Yj / Xj
        Kij_N = np.zeros([3, 1])
        Kij = np.append(Kij, Kij_N, axis=1)

        print("Kij = ", Kij)

        Hij = np.array([np.array([entalpiaGasIdeal(constantes, T) for T in Teqj]) for constantes in CONSTANTES])
        Hvapij = np.array([np.array([entalpiaVaporizacion(constantesVap, Tc, T) for T in Teqj]) for Tc, constantesVap in zip(Tci, CONSTANTESVAP)])
        hij = entalpiaLiquidoIdeal(Hij, Hvapij)

        hj = np.sum(Xj * hij, axis=0) / 1e3
        Hj = np.sum(Yj * Hij, axis=0) / 1e3

        Qj[0] = Vj[1] * Hj[1] - (Lj[0] + Uj[0]) * hj[0]
        Qj[-1] = np.sum(Fj * hFj - Uj * hj - Wj * Hj) - Lj[-1] * hj[-1]

        primer = primero(Hj)
        Vj = segundo(primer, Vj)

        print("Vj = ", Vj)

        TjInit2 = Teqj

        # tau = np.sum((TjInit1 - TjInit2) ** 2)

        TjInit1 = TjInit2

        print("Tempj2 = ", TjInit2)
        # print("tau = ", tau)

        criterio1 = 0.01 * ETAPAS
        print("criterio1 = ", criterio1)
        # if tau <= 0.01 * ETAPAS or i > 20:
        if i > 100:
            break

    return Vj


Te = np.array([350, 355, 360, 365, 370])
Eqj = np.zeros([NUMERO_COMPONENTES, ETAPAS])
yij = np.zeros([NUMERO_COMPONENTES, ETAPAS])
Teqj = np.zeros(ETAPAS)
# Teqj = np.zeros([10, ETAPAS])
Yi = np.zeros(NUMERO_COMPONENTES)
hj = np.zeros(ETAPAS)

columna(Vj, Kij)






