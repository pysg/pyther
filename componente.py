import numpy as np
from pure_data import Data_parse


def argument_gpec(func):
    def box(*arg):
        try:
            if 1 <= arg[0] <= 3 and 0 <= arg[1] <= 3: return func(*arg)
        except ValueError:
            return ValueError("NMODEL = 1,2,3 and ICALC = 0,1,2,3")
        else:
            return func(3, 2)
    return box

@argument_gpec
def gpec2pyther_argument(NMODEL, ICALC):    

    if ICALC == 0: ICALC = "constants_eps"
    elif ICALC == 1: ICALC = "parameters_eps"
    elif ICALC == 2: ICALC = "rk_param"
    elif ICALC == 3: ICALC = "density"

    if NMODEL == 1: NMODEL = "SRK"
    elif NMODEL == 2: NMODEL = "PR"
    elif NMODEL == 3: NMODEL = "RKPR"

    return NMODEL, ICALC

def selec_eos_cal(NMODEL, ICALC):
    try:
        model = NMODEL in ["SRK", "PR", "RKPR"]
        calc = ICALC in ["constants_eps", "parameters_eps", "rk_param", "density"]
        if model and calc:
            return NMODEL, ICALC
    except ValueError:
        msg = """NMODEL is not ["SRK", "PR", "RKPR"] or 
                 ICALC is not ["constants_eps", "parameters_eps", "rk_param", "density"]"""             
        return ValueError(msg)
    else:
        NMODEL = "RKPR"
        ICALC = "rk_param"
        return NMODEL, ICALC

def main():
    dppr_file = "PureFull.xls"
    component = 'METHANE'
    #component = "ETHANE"

    component_eos = Data_parse()
    properties_component = component_eos.selec_component(dppr_file, component)

    print ('acentric_factor = {0}'.format(properties_component[1]['Omega']))
    print ('critical_Temperature = {0} K'.format(properties_component[1]['Tc']))
    print ('critical_Pressure = {0} Bar'.format(properties_component[1]['Pc']))
    print ('critical_Volume = {0} cm3/mol'.format(properties_component[1]['Vc']))
    print ('compressibility_factor_Zz = {0}'.format(properties_component[1]['Zc']))
    NMODEL = 1
    ICALC = -3    

    if type(NMODEL) == int or type(ICALC) == int:
        NMODEL, ICALC = gpec2pyther_argument(NMODEL, ICALC)
    else:#if type(NMODEL) == str or type(ICALC) == str:
        NMODEL, ICALC = selec_eos_cal(NMODEL, ICALC)
        
    print('NMODEL = {0} and ICALC = {1}'.format(NMODEL, ICALC))


if __name__ == "__main__":
    # execute only if run as a script
    main()
