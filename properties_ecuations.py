import numpy as np
import pandas as pd
import pyther as pt

dppr_file = "PureFull_mod_properties.xls"

print(dppr_file)

#component = 'METHANE'
#component = "ETHANE"
#component = "3-METHYLHEPTANE"
#component = "n-PENTACOSANE"
component = "ISOBUTANE"

#components = ["METHANE", "n-TETRACOSANE"]

properties_data = pt.Data_parse()

def read_dppr(dppr_file):
        #self.dppr_data = pd.read_excel(dppr_file).head().set_index('Name').ix[:, 1:12]

        #dppr_data = pd.read_excel(dppr_file).set_index("Name").ix[:, 1:5]
        dppr_data = pd.read_excel(dppr_file).ix[:, 0:8]
        
        # component_names = dppr_data.index.get_values()
        return dppr_data

data = read_dppr(dppr_file)

print(data)

g = 3
components_labels = [x for x in range(0, 13*g, 13)]
data_name = data.ix[components_labels, 0].get_values()

def rho_liq_cal(component):
	
	rho_liquido_constans = [x+1 for x in range(0, 13*g, 13)]
	datos_rho_liquido = data.ix[rho_liquido_constans, 1:8].get_values()

	print(datos_rho_liquido)
	datos_rho_liquido = pd.DataFrame(data=datos_rho_liquido,index=data_name,
						 columns=["A", "B", "C", "D", "E", "Min", "Max"])

	rho_liquido_constans_component = datos_rho_liquido.loc[component]


	return rho_liquido_constans_component

# Liquid_density [kmol/m^3]
# A / B ^ (1 + (1 - T / C) ^ D)
A, B, C, D, E, Min, Max = rho_liq_cal("ETHANE")

Temp_vector = np.array([Temp_vector for Temp_vector in np.arange(Min, Max)])
rho_liquido = A / B ** (1 + (1 - Temp_vector / C) ** D)
print(rho_liquido)