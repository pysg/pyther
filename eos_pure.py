import numpy as np


def func_zc_d1(del1):

    d1 = (1 + del1 ** 2) / (1 + del1)
    y = 1 + (2 * (1 + del1)) ** (1.0 / 3) + (4 / (1 + del1)) ** (1.0 / 3)
    Zc = y / (3 * y + d1 - 1.0)

    return Zc


def get_del_1(Zcin, del1):

    # del1 = del_1_init
    # d1 = (1 + del1 ** 2) / (1 + del1)
    # y = 1 + (2 * (1 + del1)) ** (1.0 / 3) + (4 / (1 + del1)) ** (1.0 / 3)
    # Zc = y / (3 * y + d1 - 1.0)

    # dold = del_1_init
    # Zc = func_zc_d1(del_1_init)

    # if Zc > Zcin:
    #    del1 = 1.01 * del1
    # else:
    #    del1 = 0.99 * del1

    # error_Z_critico = abs(Zc - Zcin)

    while True:
    	Zc = func_zc_d1(del1)
    	aux = del1
    	del1 = del1 - (Zc - Zcin) * (del1 - dold) / (Zc - Zold)
    	dold = aux
    	error_Z_critico = abs(Zc - Zcin)
    	if error_Z_critico <= 1e-6:
    		break

    while True:

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


def getdel1(Zcin, del_1_init):

    del1 = del_1_init
    d1 = (1 + del1 ** 2) / (1 + del1)
    y = 1 + (2 * (1 + del1)) ** (1.0 / 3) + (4 / (1 + del1)) ** (1.0 / 3)
    Zc = y / (3 * y + d1 - 1.0)

    dold = del_1_init

    if Zc > Zcin:
        del1 = 1.01 * del1
    else:
        del1 = 0.99 * del1

    error_Z_critico = abs(Zc - Zcin)

    while True:

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


def acentric_factor_cal(*arg):

    al, be, ga = arg[0], arg[1], arg[2]
    try:
        OM = 0.5 * (- be + np.sqrt(be ** 2 - 4 * al * ga)) / al
    except RuntimeWarning:
        raise RuntimeWarning
    else:
        OM = 0
    return OM


def compressibility_factor_cal(del1):

    d1 = (1 + del1 ** 2) / (1 + del1)
    y = 1 + (2 * (1 + del1)) ** (1.0 / 3) + (4 / (1 + del1)) ** (1.0 / 3)
    # numerator_OMa = (3 * y * y + 3 * y * d1 + d1 ** 2 + d1 - 1.0)
    numerator_OMa = (3 * y **2 + 3 * y * d1 + d1 ** 2 + d1 - 1.0)
    denominator_OMa = (3 * y + d1 - 1.0) ** 2
    OMa = numerator_OMa / denominator_OMa
    OMb = 1 / (3 * y + d1 - 1.0)
    Zc = y / (3 * y + d1 - 1.0)

    return Zc, OMa, OMb
