import numpy as np


RGAS = 0.08314472
D = np.array([0.428363, 18.496215, 0.338426, 0.660, 789.723105, 2.512392])

'''
ICALC = [0, 1, 2, 3]
NMODEL = [1, 2, 3]
Sin contar IVAP, puesto que se define actualmente siempre IVAP = 0,
se tienen 8 casos de posibles especificaciones
'''

# prueba = 'SRK_0' # ICALC, NMODEL, IVAP = 0, 1, 0 
# prueba = "SRK_1" # ICALC,NMODEL,IVAP = 1, 1, 0

# prueba = "PR_0" # ICALC, NMODEL, IVAP = 0, 2, 0
# prueba = "PR_1" # ICALC,NMODEL,IVAP = 1, 2, 0

prueba = 'RKPR_0' # ICALC, NMODEL, IVAP = 0, 3, 0
# prueba = "RKPR_1" # ICALC,NMODEL,IVAP = 1, 3, 0
# prueba = "RKPR_2" # ICALC,NMODEL,IVAP = 2, 3, 0
# prueba = "RKPR_3" # ICALC,NMODEL,IVAP = 3, 3, 0


def eos(prueba):

    if prueba == "RKPR_0":
        #305.32  48.72  0.09949  0.1724175  1.185   Tc, Pc, omega, Vceos(L/mol)  C2
        ICALC, NMODEL, IVAP = 0, 3, 0
        componente = "Etano_RKPR_0"
        Tc, Pc, omega, Vceos = 305.32, 48.72, 0.09949, 0.1724175 #  1.185
        print("componente = {0} \nTc = {1} \nPc = {2} \nomega = {3} \nVceos = {4}".format(componente, Tc, Pc, omega, Vceos))

        RT = RGAS * Tc
        Zc = Pc * Vceos / RT
        print("Zc = {0}".format(Zc))
        del1ini = D[0] + D[1] * (D[2] - Zc)**D[3] + D[4] * (D[2] - Zc)**D[5]
        print("del1ini = {0}".format(del1ini))

        dinputs = np.array([Tc, Pc, omega, Vceos])

    elif prueba == "RKPR_1":
        # 1  3  0                           ICALC,NMODEL,IVAP
        # 2422.44    3.90091    1.80323 8.60102     ac,b,del1,rk
        # Me sale este output:
        # 1377.9116    2.4027       NaN 15.22254
        # 2422.4400  3.900910  1.803230   8.60102
        ICALC,NMODEL,IVAP = 1, 3, 0
        componente = "Prueba_1_RKPR_1"
        ac, b, del1, rk = 2422.44, 3.90091, 1.80323, 8.60102
        print("componente = {0} \nac = {1} \nb = {2} \ndel1 = {3} \nrk = {4}".
            format(componente, ac,b,del1,rk))

        dinputs = np.array([ac, b, del1, rk])

    elif prueba == "RKPR_2":
        # 2  3  0                            ICALC,NMODEL,IVAP 
        # Tc    Pc(bar)     W   Zc  dc(mol/L)
        # CYCLOPROPANE       397.91  54.94956075 0.1269  0.27    6.142506143
        # del1 SPECIFICATION together with Tc,Pc,OM
        ICALC, NMODEL, IVAP = 2, 3, 0
        # componente = "CYCLOPROPANE_RKPR_2"
        # Tc, Pc, W, Zc, dc = 397.91, 54.94956075, 0.1269, 0.27, 6.142506143
        # print("componente = {0} \nTc = {1} \nPc = {2} \nW = {3} \nZc = {4} \ndc = {5}".format(componente, Tc, Pc, W, Zc, dc))

        # dinputs = np.array([Tc, Pc, W, Zc, dc])

        # METHANE (1)
        Tc, Pc, ohm, vc, zrat = 190.564, 45.99, 0.01155, 0.115165, 1.00000173664
        ac, b, d, rk = 2.3277, 0.029962, 0.932475, 1.49541
        # Tr = 0.7 # Change here to use another Pv than the one at Tr 0.7
        # Tr = 0.7
        delta_1 = d
        dc = 1 / vc
        OM = ohm

        dinputs = np.array([Tc, Pc, OM, dc, zrat, ac, b, d, rk])

    elif prueba == "RKPR_3":
        # ! RhoLsat SPECIFICATION together with Tc,Pc,OM
        # READ(NIN,*)T, RhoLsat
        # Trho = T/Tc
        # ! initial value
        del1 = 2.0
        RHOld = 0.0

        ICALC, NMODEL, IVAP = 3, 3, 0

        # Carbon Dioxide
        # Tc, Pc, omega
        Tc, Pc, omega = 304.21, 73.83, 0.2236
        # T(K), RhoLsat (L/mol)
        T_especific, RHOLSat_esp = 270.0, 21.4626

        # valor initial of delta_1
        delta_1 = 0.5

        dinputs = np.array([Tc, Pc, omega, delta_1, T_especific, RHOLSat_esp])
        # dinputs = np.array([omega, delta_1, T_especific, RHOLSat_esp])

    elif prueba == "PR_0":
        # 0  2  0                           ICALC,NMODEL,IVAP

        ICALC, NMODEL, IVAP = 0, 2, 0

        #                Tc,     Pc,   omega
        # Saturates     930     11.98   0.9     1.35
        # Arom+Resins   1074    10.85   1.5     1.35
        # Asphaltenes   1274    6.84    1.75    1.35
        componente = "Asphaltenes"
        Tc, Pc, OM = 1274.0, 6.84, 1.75
        print("componente = {0} \nTc = {1} \nPc = {2} \nOM = {3}".format(componente, Tc, Pc, OM))
        dinputs = np.array([Tc, Pc, OM])
        del1 = 1
    elif prueba == "PR_1":
        pass

    return dinputs, NMODEL, ICALC


def convert_argument(NMODEL, ICALC):

    if ICALC == 0:
        ICALC = "constants_eps"
    elif ICALC == 1:
        ICALC = "parameters_eps"
    elif ICALC == 2:
        ICALC = "rk_param"
    elif ICALC == 3:
        ICALC = "density"

    if NMODEL == 1:
        NMODEL = "SRK"
    elif NMODEL == 2:
        NMODEL = "PR"
    elif NMODEL == 3:
        NMODEL = "RKPR"

    print('NMODEL = {0} and ICALC = {1}'.format(NMODEL, ICALC))

    return NMODEL, ICALC


def main():

    dinputs, NMODEL, ICALC = eos(prueba)
    NMODEL, ICALC = convert_argument(NMODEL, ICALC)
    print((NMODEL))


if __name__ == "__main__":
    # execute only if run as a script
    main()
