import numpy as np
from scipy.optimize import fsolve
from eos_selecction import eos, convert_argument
from models_eos import models_eos_cal

'''
This module calculate parameters necessary to use the equations of state:
SRK, PR and RKPR
'''

RGAS = 0.08314472
# Definir el significado fisicoquímico

A0, B0, C0 = 0.0017, 1.9681, -2.7238
# Definir el significado fisicoquímico

A1, B1, C1 = -2.4407, 7.4513, 12.504
# Definir el significado fisicoquímico

D = np.array([0.428363, 18.496215, 0.338426, 0.660, 789.723105, 2.512392])


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


def initial_data(omega, delta_1, NMODEL, ICALC, Pc):

    d1 = (1 + delta_1 ** 2) / (1 + delta_1)
    y = 1 + (2 * (1 + delta_1)) ** (1.0 / 3) + (4 / (1 + delta_1)) ** (1.0 / 3)
    OMa = (3 * y * y + 3 * y * d1 + d1 ** 2 + d1 - 1.0) \
    / (3 * y + d1 - 1.0) ** 2    
    OMb = 1 / (3 * y + d1 - 1.0)
    Zc = y / (3 * y + d1 - 1.0)

    # initial guess for k parameter
    rk = (A1 * Zc + A0) * omega**2 + (B1 * Zc + B0) * omega + (C1 * Zc + C0)
    #rk = rk * 1.2 # 1.1 #5.2 #3.2
    # if ICALC == 1 or ICALC == 2:
    if ICALC == 'parameters_eps' or ICALC == 'rk_param':
        rk = rk * 1.2
        Tr = 0.7
        Pvdat = Pc * 10 ** -(1.0 + omega)         
    elif ICALC == 'density':
        rk = rk * 5.2 
        Tr_calculada = dinputs[4] / dinputs[0]
        Tr = Tr_calculada
        Pvdat = Pc * 10 ** -((1.0 / Tr - 1.0) * 7 * (1.0 + omega) / 3)

    return rk, Pvdat, Tr


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

        #print(self.PV_calculed)

        return self.PV_calculed

    def funcion_saturacion_cal(self, rk_inicial, delta_1_initial, Pvdat, Pc, Tc, Tr):

        presion_saturada_modelo = self.phase_equilibrium_cal(rk_inicial, delta_1_initial, Pvdat, Pc, Tc, Tr)
        self.saturation_function = abs(presion_saturada_modelo - Pvdat) / Pvdat

        #print(self.saturation_function)

        return self.saturation_function

    def resolver_rk_cal(self, rk_inicial, delta_1_initial, Pvdat, Pc, Tc, Tr):

        self.rk_calculated = fsolve(self.funcion_saturacion_cal, rk_inicial,
        args = (delta_1_initial, Pvdat, Pc, Tc, Tr), xtol=1e-4)

        return self.rk_calculated

    def parameter_ab_cal(self, delta_1_initial, Pc, Tc):
        #RT = RGAS * T_especific
        #RT = RGAS * dinputs[0]
        d1 = (1 + delta_1_initial ** 2) / (1 + delta_1_initial)
        y = 1 + (2 * (1 + delta_1_initial)) ** (1.0 / 3) \
        + (4.0 / (1 + delta_1_initial)) ** (1.0 / 3.0)
        OMa = (3 * y * y + 3 * y * d1 + d1 ** 2 + d1 - 1.0) \
        / (3 * y + d1 - 1.0) ** 2
        OMb = 1 / (3 * y + d1 - 1.0)
        Zc = y / (3 * y + d1 - 1.0)
        dc = Pc / Zc / (RGAS * Tc)
        Vceos = 1.0 / dc

        #self.ac = OMa * RT**2 / Pc
        self.ac = OMa * (RGAS * Tc)**2 / Pc
        self.b = OMb * (RGAS * Tc) / Pc

        return self.ac, self.b

    def density_cal(self, delta_1_initial, rk_inicial, Pvdat, Pc, Tc, Tr):

        self.a = self.energetic_parameter_cal(rk_inicial, delta_1_initial, Pc, Tc, Tr)
        self.b = self.parameter_ab_cal(delta_1_initial, Pc, Tc)[1]
        
        self.rk_calculated = self.resolver_rk_cal(rk_inicial, delta_1_initial, Pvdat, Pc, Tc, Tr)

        #print(self.a, self.b)

        #print(self.rk_calculated)

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


def models_eos_cal(NMODEL, ICALC, dinputs):

    if NMODEL == 'SRK' or NMODEL == 'PR':
        # CONSTANTS SPECIFICATION READ(Tc,Pc,OM)
        if ICALC == 'constants_eps':
            Tc = dinputs[0]
            Pc = dinputs[1]
            OM = dinputs[2]
            RT = RGAS * Tc

            if NMODEL == 'SRK':
                rm = 0.48 + 1.574 * OM - 0.175 * OM**2
                del1 = 1.0
            elif NMODEL == 'PR':
                rm = 0.37464 + 1.54226 * OM - 0.26992 * OM ** 2
                del1 = 1.0 + np.sqrt(2.0)

            Zc, OMa, OMb = compressibility_factor_cal(del1)
            Vceos = Zc * RGAS * Tc / Pc
            ac = OMa * RT ** 2 / Pc
            b = OMb * RT / Pc
            params = [ac, b, rm]

        if ICALC == 'parameters_eps':
            Tc = OMb * ac / (OMa * RGAS * b)
            Pc = OMb * RGAS * Tc / b
            Vceos = Zc * RGAS * Tc / Pc

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

            OM = acentric_factor_cal(al, be, ga)

            constants = [Tc, Pc, OM, Vceos]

    elif NMODEL == 'RKPR':
        if ICALC == 'constants_eps':
            # CONSTANTS SPECIFICATION (Tc,Pc,OM,Vceos)
            Tc = dinputs[0]
            Pc = dinputs[1]
            OM = dinputs[2]
            Vceos = dinputs[3]
            RT = RGAS * Tc
            Zc = Pc * Vceos / RT
            del1ini = D[0] + D[1] * (D[2] - Zc) ** D[3] + D[4] * (D[2] - Zc)** D[5]
            print('del1ini = {0}'.format(del1ini))

            del_1 = getdel1(Zc, del1ini)[0]

            Zc, OMa, OMb = compressibility_factor_cal(del_1)
            print('Zc = {0}'.format(Zc))

            ac = OMa * RT ** 2 / Pc
            b = OMb * RT / Pc

            # calcular rk
            rk = 1

            params = [ac, b, del_1, rk]

        elif ICALC == 'parameters_eps':

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
        # RKPR EOS (Tc,Pc,OM)
        elif ICALC == 'rk_param':

            # dinputs = np.array([Tc, Pc, OM, dc, zrat, ac, b, d, rk])

            Tc = dinputs[0]
            Pc = dinputs[1]
            OM = dinputs[2]
            dc = dinputs[3]
            zrat = dinputs[4]
            ac = dinputs[5]
            b = dinputs[6]
            d = dinputs[7]
            delta_1 = d
            rk = dinputs[8]

            omega = OM

            #RT = RGAS * Tc

            #print(dinputs)

            # THEN  ! del1 SPECIFICATION together with Tc,Pc,OM
            # READ(NIN,*)del1  ! this line must be active when it is not a list

            #rk = dinputs[0]
            #delta_1 = dinputs[1]

            #print(rk, delta_1, omega)
            rk, Pvdat, Tr = initial_data(omega, delta_1, NMODEL, ICALC, Pc)

            #print(Pvdat)
            rk_cal = eos_calculation.resolver_rk_cal(rk, delta_1, Pvdat, Pc, Tc, Tr)
            #print('rk_cal = '.format(rk_cal))

        # RhoLsat SPECIFICATION together with Tc,Pc,OM
        elif ICALC == 'density':
            # (T, RhoLsat)
            # Trho = T / Tc
            # del1 = 2.0    #!  initial value
            # RHOld = 0.0

            Tc = dinputs[0]
            Pc = dinputs[1]

            omega = dinputs[2]
            delta_1 = dinputs[3]
            T_especific = dinputs[4]
            RHOLSat_esp = dinputs[5]

            rk, Pvdat, Tr = initial_data(omega, delta_1, NMODEL, ICALC, Pc)
            delta_1_parameter = eos_calculation.resolver_delta_1_cal(delta_1, rk, Pvdat, RHOLSat_esp, Pc, Tc, Tr)

            #print(delta_1_parameter)


    print('The NMODEL is eos_{0} and method ICALC is {1}'.format(NMODEL, ICALC))


    if ICALC == 'constants_eps':
        print(params)
        return params
    elif ICALC == 'parameters_eps':
        print(constants)
        return constants
    elif ICALC == 'rk_param':
        print('The parameter rk_cal is {0}'.format(rk_cal))
        return rk_cal
    elif ICALC == 'density':
        print('The parameter delta1(rho,T) is {0}'.format(delta_1_parameter))
        return delta_1_parameter


print("-" * 79)

eos_calculation = Parameter_eos()

# functions import from eos_selection
dinputs, NMODEL, ICALC = eos('RKPR_2')
NMODEL, ICALC = convert_argument(NMODEL, ICALC)

resultado = models_eos_cal(NMODEL, ICALC, dinputs)

print('-' * 79)












