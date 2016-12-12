import numpy as np
import pyther as pt

properties_data = pt.Data_parse()

dppr_file = "PureFull.xls"
#component = 'METHANE'
#component = "ETHANE"
#component = "3-METHYLHEPTANE"
#component = "n-PENTACOSANE"
component = "ISOBUTANE"
NMODEL = "RKPR"
ICALC = "constants_eps"

properties_component = properties_data.selec_component(dppr_file, component)
pt.print_properties_component(component, properties_component)
dinputs = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
                    properties_component[1]['Omega'], properties_component[1]['Vc']])

component_eos = pt.models_eos_cal(NMODEL, ICALC, dinputs)

ac = component_eos[0]
print(ac)
print('-' * 79)

