import numpy as np
from scipy.optimize import fsolve
# from .eos_selecction import eos, convert_argument
from eos_selecction import eos, convert_argument
# from .constants import RGAS
from constants import RGAS


'''
This module calculate parameters necessary to use the equations of state:
SRK, PR and RKPR
'''

# Constante universal de los gases R
# RGAS = 0.08314472

# Definir el significado fisicoquímico
# A0, B0, C0 = 0.0017, 1.9681, -2.7238

# Definir el significado fisicoquímico
# A1, B1, C1 = -2.4407, 7.4513, 12.504

# Definir el significado fisicoquímico
# D = np.array([0.428363, 18.496215, 0.338426, 0.660, 789.723105, 2.512392])


def getdel1(Zcin, del_1_init):
    del1 = del_1_init
    d1 = (1 + del1 ** 2) / (1 + del1)
    y = 1 + (2 * (1 + del1)) ** (1.0 / 3) + (4 / (1 + del1)) ** (1.0 / 3)
    Zc = y / (3 * y + d1 - 1.0)
    dold = del1

    if Zc > Zcin:
        del1 = 1.01 * del1
    else:
        del1 = 0.99 * del1

    error_Z_critico = abs(Zc - Zcin)

    while error_Z_critico >= 1e-6:

        d1 = (1 + del1 ** 2) / (1 + del1)
        y = 1 + (2 * (1 + del1)) ** (1.0 / 3) + (4 / (1 + del1)) ** (1.0 / 3)
        Zold = Zc
        Zc = y / (3 * y + d1 - 1.0)
        aux = del1
        del1 = del1 - (Zc - Zcin) * (del1 - dold) / (Zc - Zold)
        dold = aux
        error_Z_critico = abs(Zc - Zcin)

        if error_Z_critico <= 1e-6:
            break

    return del1, error_Z_critico


def get_del_1(Zcin, del_1_init):
    del1 = del_1_init
    d1 = (1 + del1 ** 2) / (1 + del1)
    y = 1 + (2 * (1 + del1)) ** (1.0 / 3) + (4 / (1 + del1)) ** (1.0 / 3)
    Zc = y / (3 * y + d1 - 1.0)
    dold = del1

    if Zc > Zcin:
        del1 = 1.01 * del1
    else:
        del1 = 0.99 * del1

    error_Z_critico = abs(Zc - Zcin)

    while error_Z_critico >= 1e-6:

        d1 = (1 + del1 ** 2) / (1 + del1)
        y = 1 + (2 * (1 + del1)) ** (1.0 / 3) + (4 / (1 + del1)) ** (1.0 / 3)
        Zold = Zc
        Zc = y / (3 * y + d1 - 1.0)
        aux = del1
        del1 = del1 - (Zc - Zcin) * (del1 - dold) / (Zc - Zold)
        dold = aux
        error_Z_critico = abs(Zc - Zcin)
    
    # else:
        # print("delta_1 = {0} with a error of = {1}".format(del1, error_Z_critico))

        if error_Z_critico <= 1e-6:
            break

    return del1, error_Z_critico


def compressibility_factor_cal(del1):
    d1 = (1 + del1 ** 2) / (1 + del1)
    y = 1 + (2 * (1 + del1)) ** (1.0 / 3) + (4 / (1 + del1)) ** (1.0 / 3)
    OMa = (3 * y * y + 3 * y * d1 + d1 ** 2 + d1 - 1.0) / (3 * y + d1 - 1.0) ** 2
    OMb = 1 / (3 * y + d1 - 1.0)
    Zc = y / (3 * y + d1 - 1.0)

    return Zc, OMa, OMb


def acentric_factor_cal(*arg):
    al = arg[0]
    be = arg[1]
    ga = arg[2]
    try:
        OM = 0.5 * (- be + np.sqrt(be ** 2 - 4 * al * ga)) / al
    except RuntimeWarning:
        raise RuntimeWarning
    else:
        OM = 0
    return OM


class Parameter_eos(object):
    '''
    Parameter_eos contains the methods to adjust the parameters
    delta 1 rk for the energy parameter to a function of temperature and
    setting the parameter to represent a point of density - temperature curve
    '''

    def energetic_parameter_cal(self, rk, delta_1_initial, Pc, Tc, Tr):

        self.ac = self.parameter_ab_cal(delta_1_initial, Pc, Tc)[0]
        self.a = self.ac * (3 / (2 + Tr)) ** rk

        return self.a

    def gvdW_Derivatives_cal(self, NDER, Volume, a, b, delta_1_initial, Tc, Tr):
        '''
        gvdW_Derivatives_cal: calculate the derivatives from Generalized van
        der Waals equation of state
        '''
        self.d = delta_1_initial

        c = (1 - self.d) / (1 + self.d)
        #aRT = self.a / (RGAS * T_especific)
        aRT = self.a / (RGAS * Tr * Tc)
        relation_covolume = 0.25 * b / Volume
        SUMC = c * b + Volume
        SUMD = self.d * b + Volume
        REP = -np.log(1 - 4 * relation_covolume)
        ATT = aRT * np.log(SUMD / SUMC) / (b * (c - self.d))
        ATTV = aRT / SUMC / SUMD
        REPV = 1 / (1 - 4 * relation_covolume) - 1
        REP2V = 1 / (1 - 4 * relation_covolume) ** 2 - 1
        ATT2V = aRT * Volume ** 2 * (1 / SUMD ** 2 - 1 / SUMC ** 2) \
        / (b * (c - self.d))
        F = REP + ATT
        F_V = (- REPV / Volume + ATTV)

        if NDER == 'Derivatives_of_V':
            F_2V = REP2V - ATT2V
            calculo_1 = 'F_V'
            calculo_2 = 'F_2V'
            return F, F_V, F_2V, calculo_1, calculo_2
        elif NDER == 'Derivatives_of_N':
            F_N = REP + ATT - Volume * F_V
            calculo = 'F_N'
            return F_N, calculo

    def volume_cal(self, ITYP, T, P, a, b, delta_1_initial, Tc, Tr):
        volume_iteration = 0
        VCP = b
        S3R = 1.0 / VCP
        relation_covolume_min, relation_covolume_max = 0.00, 0.99
        P_sur = P
        if ITYP >= 0.0:
            relation_covolume = 0.5
        else:
            # IDEAL GAS ESTIMATE
            relation_covolume = min(0.5, (VCP * P_sur) / (RGAS * T))
        while True:
            Volume = VCP / relation_covolume

            vdWg = self.gvdW_Derivatives_cal('Derivatives_of_V', Volume, a, b, delta_1_initial, Tc, Tr)
            F = vdWg[0]
            F_V = vdWg[1]
            F_2V = vdWg[2]
            pressure_cal_eos = RGAS * T * (1 / Volume - F_V)
            if pressure_cal_eos > P_sur:
                relation_covolume_max = relation_covolume
            else:
                relation_covolume_min = relation_covolume
            AT = F - np.log(Volume) + Volume * P_sur / (T * RGAS)
            DER = RGAS * T * (F_2V + 1.0) * S3R
            DEL = - (pressure_cal_eos - P_sur) / DER
            relation_covolume = relation_covolume + 1.0 * max(min(DEL, 0.1), -0.1)
            if relation_covolume > relation_covolume_max or relation_covolume < relation_covolume_min:
                relation_covolume = 0.5 * (relation_covolume_max + relation_covolume_min)
            if abs(DEL) < 1e-10 or volume_iteration >= 20:
                break
            volume_iteration += 1
        return Volume, pressure_cal_eos

    def fugacity_cal(self, T, P, Volume, delta_1_initial, Tc, Tr):

        Z = P * Volume / (RGAS * T)
        vdWg = self.gvdW_Derivatives_cal('Derivatives_of_N', Volume, self.a, self.b, delta_1_initial, Tc, Tr)
        F_N = vdWg[0]
        phi_fugacity = np.exp(F_N) / Z

        #print(T, T_especific)
        #print(T)
        return phi_fugacity

    def saturation_pressure_cal(self, PV_supuesta, rk_inicial, delta_1_initial, Pc, Tc, Tr):
        self.Tsat = Tr * Tc
        P_sur = PV_supuesta

        #print(self.Tsat)

        self.a = self.energetic_parameter_cal(rk_inicial, delta_1_initial, Pc, Tc, Tr)
        self.b = self.parameter_ab_cal(delta_1_initial, Pc, Tc)[1]

        self.volume_liquid = self.volume_cal(1, self.Tsat, P_sur, self.a,
        self.b, delta_1_initial, Tc, Tr)[0]

        phi_fugacity_liquid = self.fugacity_cal(self.Tsat, P_sur,
        self.volume_liquid, delta_1_initial, Tc, Tr)

        self.volume_vapor = self.volume_cal(-1, self.Tsat, P_sur, self.a,
        self.b, delta_1_initial, Tc, Tr)[0]

        phi_fugacity_vapor = self.fugacity_cal(self.Tsat, P_sur,
        self.volume_vapor, delta_1_initial, Tc, Tr)

        self.phase_equilibrium = abs(phi_fugacity_vapor - phi_fugacity_liquid)

        return self.phase_equilibrium

    def phase_equilibrium_cal(self, rk_inicial, delta_1_initial, Pvdat, Pc, Tc, Tr):
        PV_inicial = Pvdat

        rk_class = rk_inicial

        self.PV_calculed = fsolve(self.saturation_pressure_cal, PV_inicial,
        args=(rk_class, delta_1_initial, Pc, Tc, Tr), xtol=1e-4)

        return self.PV_calculed

    def funcion_saturacion_cal(self, rk_inicial, delta_1_initial, Pvdat, Pc, Tc, Tr):

        presion_saturada_modelo = self.phase_equilibrium_cal(rk_inicial, delta_1_initial, Pvdat, Pc, Tc, Tr)
        self.saturation_function = abs(presion_saturada_modelo - Pvdat) / Pvdat

        return self.saturation_function

    def resolver_rk_cal(self, rk_inicial, delta_1_initial, Pvdat, Pc, Tc, Tr):

        self.rk_calculated = fsolve(self.funcion_saturacion_cal, rk_inicial,
        args = (delta_1_initial, Pvdat, Pc, Tc, Tr), xtol=1e-4)

        return self.rk_calculated

    def parameter_ab_cal(self, delta_1_initial, Pc, Tc):

        Zc, OMa, OMb = compressibility_factor_cal(delta_1_initial)

        dc = Pc / Zc / (RGAS * Tc)
        Vceos = 1.0 / dc

        self.ac = OMa * (RGAS * Tc)**2 / Pc
        self.b = OMb * (RGAS * Tc) / Pc

        return self.ac, self.b

    def density_cal(self, delta_1_initial, rk_inicial, Pvdat, Pc, Tc, Tr):

        self.a = self.energetic_parameter_cal(rk_inicial, delta_1_initial, Pc, Tc, Tr)
        self.b = self.parameter_ab_cal(delta_1_initial, Pc, Tc)[1]
        
        self.rk_calculated = self.resolver_rk_cal(rk_inicial, delta_1_initial, Pvdat, Pc, Tc, Tr)

        self.volume_liquid_delta = self.volume_liquid
        self.density_liquid = 1 / self.volume_liquid_delta

        return self.density_liquid

    def density_function_cal(self, delta_1_initial, rk_inicial, Pvdat, RHOLSat_esp, Pc, Tc, Tr):

        RHOLSat_cal = self.density_cal(delta_1_initial, rk_inicial, Pvdat, Pc, Tc, Tr)
        self.density_function = abs(RHOLSat_cal - RHOLSat_esp)
        self.density_function = self.density_function / RHOLSat_esp

        #print(self.density_function)

        return self.density_function

    def resolver_delta_1_cal(self, delta_1_initial, rk_inicial, Pvdat, RHOLSat_esp, Pc, Tc, Tr):
        
        self.delta_1_calculated = fsolve(self.density_function_cal, delta_1_initial, args = (rk_inicial, Pvdat, RHOLSat_esp, Pc, Tc, Tr),
            xtol=1e-3)        

        return self.delta_1_calculated

    def abstract_cal(self):
        return self.a, self.b
