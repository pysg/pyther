import numpy as np
from eos_selecction import eos, convert_argument
from cubic_parameters_1 import Parameter_eos, getdel1, compressibility_factor_cal, acentric_factor_cal

from componente import Control_arguments

RGAS = 0.08314472
# Definir el significado fisicoquímico

A0, B0, C0 = 0.0017, 1.9681, -2.7238
# Definir el significado fisicoquímico

A1, B1, C1 = -2.4407, 7.4513, 12.504
# Definir el significado fisicoquímico

D = np.array([0.428363, 18.496215, 0.338426, 0.660, 789.723105, 2.512392])

def initial_data(omega, delta_1, NMODEL, ICALC, Pc, dinputs):

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

        	ac = dinputs[0]
        	b = dinputs[1]
        	del1 = dinputs[2]
        	rk = dinputs[3]

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
            rk, Pvdat, Tr = initial_data(omega, delta_1, NMODEL, ICALC, Pc, dinputs)

            #print(Tr)
            eos_calculation = Parameter_eos()
            rk_cal = eos_calculation.resolver_rk_cal(rk, delta_1, Pvdat, Pc, Tc, Tr)
            #rk_cal = Parameter_eos.resolver_rk_cal(rk, delta_1, Pvdat, Pc, Tc, Tr)
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

            rk, Pvdat, Tr = initial_data(omega, delta_1, NMODEL, ICALC, Pc, dinputs)
            eos_calculation = Parameter_eos()
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

nm = Control_arguments("RKPR", "constants_eps")
print ("NMODEL: {0} and ICALC: {1}".format(nm.NMODEL, nm.ICALC))

def main():

	print("-" * 79)

	dinputs, NMODEL, ICALC = eos('PR_0')
	NMODEL, ICALC = convert_argument(NMODEL, ICALC)
	resultado = models_eos_cal(NMODEL, ICALC, dinputs)

	print('-' * 79)


if __name__ == '__main__':
	main()



