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


def constans_criticals(NMODEL, ICALC, dinputs):

    # CONSTANTS SPECIFICATION READ [Tc, Pc, OM]

    Tc = dinputs[0]
    Pc = dinputs[1]
    OM = dinputs[2]

    if NMODEL == 'SRK':
        rm = 0.48 + 1.574 * OM - 0.175 * OM**2
        del1 = 1.0
    elif NMODEL == 'PR':
        rm = 0.37464 + 1.54226 * OM - 0.26992 * OM ** 2
        del1 = 1.0 + np.sqrt(2.0)

    Zc, OMa, OMb = compressibility_factor_cal(del1)
    Vceos = (Zc * RGAS * Tc) / Pc

    ac = OMa * (RGAS * Tc) ** 2 / Pc
    b = OMb * (RGAS * Tc) / Pc

    params = [ac, b, rm, del1]

    return params


# -----------------------------------------------------------------------
# -----------------------------------------------------------------------


def parameters_criticals(NMODEL, ICALC, dinputs):

    # PARAMETERS SPECIFICATION READ [ac, b, rm]

    ac = dinputs[0]
    b = dinputs[1]
    rm = dinputs[2]

    if NMODEL == 'SRK':
        del1 = 1.0
        al = -0.175
        be = 1.574
        ga = 0.48 - rm
    elif NMODEL == 'PR':
        del1 = 1.0 + np.sqrt(2.0)
        al = -0.26992
        be = 1.54226
        ga = 0.37464 - rm

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

    # PARAMETERS SPECIFICATION READ [ac, b, rk, del1]

    ac = dinputs[0]
    b = dinputs[1]
    del1 = dinputs[2]
    rk = dinputs[3]

    Zc, OMa, OMb = compressibility_factor_cal(del1)

    Tc = OMb * ac / (OMa * RGAS * b)
    Pc = OMb * RGAS * Tc / b
    Vceos = Zc * RGAS * Tc / Pc

    al = A1 * Zc + A0
    be = B1 * Zc + B0
    ga = C1 * Zc + C0 - rk

    OM = acentric_factor_cal(al, be, ga)

    constants = [Tc, Pc, OM, Vceos]

    return constants

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------


def parameters_cal():
    ac = OMa * (RGAS * Tc) ** 2 / Pc
    b = OMb * (RGAS * Tc) / Pc

    return ac, b




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

    ac = OMa * (RGAS * Tc) ** 2 / Pc
    b = OMb * (RGAS * Tc) / Pc

    # calcular rk
    rk, Pvdat, Tr = initial_data(OM, delta_1, NMODEL, ICALC, Pc, dinputs)

    eos_calculation = Parameter_eos()
    rk_cal = eos_calculation.resolver_rk_cal(rk, delta_1, Pvdat, Pc, Tc, Tr)

    # rk = 1

    params = [ac, b, rk, delta_1]

    return params

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------


def call_rkpr_constans_delta_1(NMODEL, ICALC, dinputs):

    # CONSTANTS SPECIFICATION and del1 READ [Tc, Pc, OM, del1]

    Tc = dinputs[0]
    Pc = dinputs[1]
    OM = dinputs[2]

    delta_1 = dinputs[7]

    rk, Pvdat, Tr = initial_data(OM, delta_1, NMODEL, ICALC, Pc, dinputs)

    eos_calculation = Parameter_eos()
    rk_cal = eos_calculation.resolver_rk_cal(rk, delta_1, Pvdat, Pc, Tc, Tr)

    params = [ac, b, rk, delta_1]

    return params


def call_rkpr_constans_density(NMODEL, ICALC, dinputs):

    # CONSTANTS SPECIFICATION and (T, RhoLsat) READ [Tc, Pc, OM, del1, T, RHOLsat]
    # Trho = T / Tc,  read initial value of del1

    Tc = dinputs[0]
    Pc = dinputs[1]
    OM = dinputs[2]

    delta_1 = dinputs[3]

    T_especific = dinputs[4]
    RHOLSat_esp = dinputs[5]

    rk, Pvdat, Tr = initial_data(OM, delta_1, NMODEL, ICALC, Pc, dinputs)
    eos_calculation = Parameter_eos()
    delta_1_parameter = eos_calculation.resolver_delta_1_cal(delta_1, rk, Pvdat, RHOLSat_esp, Pc, Tc, Tr)

    params = [ac, b, rk, delta_1]

    return params

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------


def call_eos(NMODEL, ICALC, dinputs):
    """
    This function allow to calculate the eos paramaters for (SRK, PR, RKPR)
    """

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

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------


print(34)
