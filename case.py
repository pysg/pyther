import numpy as np
import pyther as pt

properties_data = pt.Data_parse()


dppr_file = "PureFull.xls"
#component = 'METHANE'
#component = "ETHANE"
#component = "3-METHYLHEPTANE"
#component = "n-PENTACOSANE"
component = "ISOBUTANE"

properties_component = properties_data.selec_component(dppr_file, component)

print ('Component = {0}'.format(component))
print ('Acentric_factor = {0}'.format(properties_component[1]['Omega']))
print ('Critical_Temperature = {0} K'.format(properties_component[1]['Tc']))
print ('Critical_Pressure = {0} Bar'.format(properties_component[1]['Pc']))
print ('Critical_Volume = {0} cm3/mol'.format(properties_component[1]['Vc']))
print ('Compressibility_factor_Z = {0}'.format(properties_component[1]['Zc']))


#dinputs = np.array[properties_component[1]['Tc']] #, properties_component[1]['Pc'], properties_component[1]['Omega']]

dinputs = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
 					properties_component[1]['Omega'], properties_component[1]['Vc']])

NMODEL = "RKPR"
ICALC = "constants_eps"

component_eos = pt.models_eos_cal(NMODEL, ICALC, dinputs)

print('-' * 79)

