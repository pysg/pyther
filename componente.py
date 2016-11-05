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


class Control_arguments(object):

    def __init__(self, NMODEL, ICALC):
        self.NMODEL = NMODEL
        self.ICALC = ICALC

    @property
    def NMODEL(self):
        print ("NMODEL is valid ...")
        return self._NMODEL

    @property
    def ICALC(self):
        print("ICALC is valid ...")
        return self._ICALC

    @NMODEL.setter
    def NMODEL(self, NMODEL):
        if NMODEL not in ["SRK", "PR", "RKPR"]:
            raise AttributeError("’NMODEL’ not in ['SRK', 'PR', 'RKPR']")
        self._NMODEL = NMODEL    

    @ICALC.setter
    def ICALC(self, ICALC):
        if ICALC not in ["constants_eps", "parameters_eps", "rk_param", "density"]:
            raise AttributeError("'ICALC' not in ['constants_eps', 'parameters_eps', 'rk_param', 'density']")
        self._ICALC = ICALC

nm = Control_arguments("PR", "constants_eps")
print ("NMODEL: {0} and ICALC: {1}".format(nm.NMODEL, nm.ICALC))


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
        
    print('NMODEL = {0} and ICALC = {1}'.format(NMODEL, ICALC))


if __name__ == "__main__":
    main()
