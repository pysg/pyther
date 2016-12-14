import numpy as np
import pyther as pt

properties_data = pt.Data_parse()

dppr_file = "PureFull.xls"
#component = 'METHANE'
#component = "ETHANE"
#component = "3-METHYLHEPTANE"
#component = "n-PENTACOSANE"
component = "ISOBUTANE"
component = "CARBON DIOXIDE"

NMODEL = "RKPR"
NMODEL = "PR"
ICALC = "constants_eps"
ICALC = "density"

properties_component = properties_data.selec_component(dppr_file, component)
pt.print_properties_component(component, properties_component)
#dinputs = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
#                    properties_component[1]['Omega'], properties_component[1]['Vc']])

T_especific = 270.0
RHOLSat_esp = 21.4626
# valor initial of delta_1
delta_1 = 1.5

dinputs = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
                    properties_component[1]['Omega'], delta_1, T_especific, RHOLSat_esp])


component_eos = pt.models_eos_cal(NMODEL, ICALC, dinputs)

#ac = component_eos[0]
print(component_eos)
print('-' * 79)

