import numpy as np
import pandas as pd
import math
# import matplotlib as plt
import matplotlib.pyplot as plt
import os

from pure_data import Data_parse
# from .cubic_parameters_1 import Parameter_eos, getdel1, compressibility_factor_cal, acentric_factor_cal
from cubic_parameters_1 import Parameter_eos, getdel1, compressibility_factor_cal, acentric_factor_cal
from constans import RGAS, A0, B0, C0, A1, B1, C1, D

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------


def initial_data(omega, delta_1, NMODEL, ICALC, Pc, dinputs):

    Zc, OMa, OMb = compressibility_factor_cal(delta_1)

    # initial guess for k parameter
    rk = (A1 * Zc + A0) * omega ** 2 + (B1 * Zc + B0) * omega + (C1 * Zc + C0)
    # rk = rk * 1.2 # 1.1 #5.2 #3.2

    argumen_1 = "constants_eps"
    argumen_2 = "parameters_eps"
    argumen_3 = "rk_param"
    argumen_4 = "density"

    if ICALC == argumen_1 or ICALC == argumen_2 or ICALC == argumen_3:
        rk *= 1.5
        Tr = 0.7
        Pvdat = Pc * 10 ** -(1.0 + omega)
    elif ICALC == argumen_4:
        # 5.2 es otro valor que se puede usar en lugar de 1.5
        rk = rk * 1.5
        Tr_calculada = dinputs[4] / dinputs[0]
        Tr = Tr_calculada
        Pvdat = Pc * 10 ** -((1.0 / Tr - 1.0) * 7 * (1.0 + omega) / 3)

    return rk, Pvdat, Tr


def data_in(ICALC, dinputs):
	if ICALC == 'constants_eps':
		# CONSTANTS SPECIFICATION (Tc,Pc,OM,Vceos)
		Tc, Pc, OM, Vceos = dinputs[0], dinputs[1], dinputs[2], dinputs[3]

    if ICALC == 'parameters_eps':
        ac, b, del1, rk = dinputs[0], dinputs[1], dinputs[2], dinputs[3]

    if ICALC == 'rk_param':
        # dinputs = np.array([Tc, Pc, OM, dc, zrat, ac, b, d, rk])
        Tc, Pc, OM, Vceos, delta_1 = dinputs[0], dinputs[1], dinputs[2], dinputs[3], dinputs[7]

    if ICALC == 'density':

        Tc, Pc, omega, Vceos, delta_1 = dinputs[0], dinputs[1], dinputs[2], dinputs[3], dinputs[4]
        T_especific, RHOLSat_esp = dinputs[5], dinputs[6]

# --------------------------------------------------------------------------

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------

def func_ac_b(Tc, Pc, Zc, OMa, OMb):

    ac = OMa * (RGAS * Tc) ** 2 / Pc
    b = OMb * (RGAS * Tc) / Pc

    return ac, b

# --------------------------------------------------------------------------

def func_zc_ac_b(Tc, Pc, delta_1):
    
    Zc, OMa, OMb = compressibility_factor_cal(delta_1)
    
    # Vceos = (Zc * RGAS * Tc) / Pc

    ac, b = func_ac_b(Tc, Pc, Zc, OMa, OMb)

    return Zc, ac, b

def constans_criticals(NMODEL, ICALC, dinputs):

    # SPECIFICATION [Tc, Pc, OM]

    Tc = dinputs[0]
    Pc = dinputs[1]
    OM = dinputs[2]

    if NMODEL == 'SRK':
        rm = 0.48 + 1.574 * OM - 0.175 * OM**2
        delta_1 = 1.0
    elif NMODEL == 'PR':
        rm = 0.37464 + 1.54226 * OM - 0.26992 * OM ** 2
        delta_1 = 1.0 + np.sqrt(2.0)

    Zc, OMa, OMb = compressibility_factor_cal(delta_1)
    # Vceos = (Zc * RGAS * Tc) / Pc

    ac, b = func_ac_b(Tc, Pc, Zc, OMa, OMb)

    params = [ac, b, rm, delta_1]

    return params


def parameters_criticals(NMODEL, ICALC, dinputs):

    # PARAMETERS SPECIFICATION [ac, b, rm]

    ac = dinputs[0]
    b = dinputs[1]
    rm = dinputs[2]

    if NMODEL == 'SRK':
        del1, al, be, ga = 1.0, -0.175, 1.574, 0.48 - rm
    elif NMODEL == 'PR':
        del1, al, be, ga = 1.0 + np.sqrt(2.0), -0.26992, 1.54226, 0.37464 - rm

    Zc, OMa, OMb = compressibility_factor_cal(del1)
    Tc = (OMb * ac) / (OMa * RGAS * b)
    Pc = OMb * RGAS * Tc / b
    OM = acentric_factor_cal(al, be, ga)
    Vceos = Zc * RGAS * Tc / Pc

    constants = [Tc, Pc, OM, Vceos]

    return constants

# -------------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def call_rkpr_parameters(NMODEL, ICALC, dinputs):

    # SPECIFICATION [ac, b, rk, del1]

    ac = dinputs[0]
    b = dinputs[1]
    del1 = dinputs[2]
    rk = dinputs[3]

    Zc, OMa, OMb = compressibility_factor_cal(del1)

    Tc = OMb * ac / (OMa * RGAS * b)
    Pc = OMb * RGAS * Tc / b
    Vceos = Zc * RGAS * Tc / Pc

    al, be, ga = A1 * Zc + A0, B1 * Zc + B0, C1 * Zc + C0 - rk

    OM = acentric_factor_cal(al, be, ga)

    constants = [Tc, Pc, OM, Vceos]

    return constants




def call_rkpr_constans_v_critic(NMODEL, ICALC, dinputs):

    # CONSTANTS SPECIFICATION READ [Tc, Pc, OM, Vceos]

    Tc = dinputs[0]
    Pc = dinputs[1]
    OM = dinputs[2]
    Vceos = dinputs[3]

    Zc = Pc * Vceos / (RGAS * Tc)

    del1ini = D[0] + D[1] * (D[2] - Zc) ** D[3] + D[4] * (D[2] - Zc) ** D[5]
    print('del1ini = {0}'.format(del1ini))

    delta_1 = getdel1(Zc, del1ini)[0]

    Zc, OMa, OMb = compressibility_factor_cal(delta_1)

    print('Zc = {0}'.format(Zc))

    ac, b = func_ac_b(Tc, Pc, Zc, OMa, OMb)

    # calcular rk
    rk, Pvdat, Tr = initial_data(OM, delta_1, NMODEL, ICALC, Pc, dinputs)

    eos_calculation = Parameter_eos()
    rk_cal = eos_calculation.resolver_rk_cal(rk, delta_1, Pvdat, Pc, Tc, Tr)

    # rk = 1

    parameters = [ac, b, rk_cal, delta_1]

    return parameters


def call_rkpr_constans_delta_1(NMODEL, ICALC, dinputs):

    # SPECIFICATION [Tc, Pc, OM, del1]

    Tc = dinputs[0]
    Pc = dinputs[1]
    OM = dinputs[2]

    delta_1 = dinputs[7]

    Zc, OMa, OMb = compressibility_factor_cal(delta_1)

    print('Zc = {0}'.format(Zc))

    ac, b = func_ac_b(Tc, Pc, Zc, OMa, OMb)

    rk, Pvdat, Tr = initial_data(OM, delta_1, NMODEL, ICALC, Pc, dinputs)

    eos_calculation = Parameter_eos()
    rk_cal = eos_calculation.resolver_rk_cal(rk, delta_1, Pvdat, Pc, Tc, Tr)

    parameters = [ac, b, rk_cal, delta_1]

    return parameters


def call_rkpr_constans_density(NMODEL, ICALC, dinputs):

    # SPECIFICATION [Tc, Pc, OM, del1, Tspec, RHOLsat]
    # Trho = Tspec / Tc,  read initial value of del1

    Tc = dinputs[0]
    Pc = dinputs[1]
    OM = dinputs[2]

    delta_1 = dinputs[3]

    T_especific = dinputs[4]
    RHOLSat_esp = dinputs[5]

    Zc, OMa, OMb = compressibility_factor_cal(delta_1)

    print('Zc = {0}'.format(Zc))

    ac, b = func_ac_b(Tc, Pc, Zc, OMa, OMb)

    rk, Pvdat, Tr = initial_data(OM, delta_1, NMODEL, ICALC, Pc, dinputs)
    eos_calculation = Parameter_eos()
    list_args = [delta_1, rk, Pvdat, RHOLSat_esp, Pc, Tc, Tr]
    delta_1_parameter = eos_calculation.resolver_delta_1_cal(list_args)

    parameters = [ac, b, rk, delta_1_parameter]

    return parameters

# ------------------------------------------------------------------------


def call_eos(NMODEL, ICALC, dinputs):

    if NMODEL == 'SRK' or NMODEL == 'PR':

        if ICALC == 'constants_eps':
            # CONSTANTS SPECIFICATION READ [Tc, Pc, OM]
            constans_criticals(NMODEL, ICALC, dinputs)
        elif ICALC == 'parameters_eps':
            # PARAMETERS SPECIFICATION READ [ac, b, rm]
            parameters_criticals(NMODEL, ICALC, dinputs)

    elif ICALC == "RKPR":

        if ICALC == 'constants_eps':
            # CONSTANTS SPECIFICATION READ [Tc, Pc, OM, Vceos]
            call_rkpr_constans_v_critic(NMODEL, ICALC, dinputs)
        elif ICALC == 'rk_param':
            # CONSTANTS SPECIFICATION and del1 READ [Tc, Pc, OM, del1]
            call_rkpr_constans_delta_1()
        elif ICALC == 'density':
            # CONSTANTS SPECIFICATION and (T, RhoLsat)
            # READ [Tc, Pc, OM, del1, T, RHOLsat]
            # Trho = T / Tc,  read initial value of del1
            call_rkpr_constans_density()
        elif ICALC == 'parameters_eps':
            # PARAMETERS SPECIFICATION READ [ac, b, rk, del1]
            call_rkpr_parameters(NMODEL, ICALC, dinputs)


print(34)


