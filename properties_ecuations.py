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

        dppr_data = pd.read_excel(dppr_file).set_index("Name").ix[:, 1:5]
        # component_names = dppr_data.index.get_values()
        return dppr_data
    

print(read_dppr(dppr_file))


#print(properties_data)

def selec_component(self, dppr_file, component):
        self.name = str(component)
        self.properties = self.read_dppr(dppr_file).loc[self.name]
        self.label = {self.name : self.properties}

        #print(self.properties.index)
        #print ('name = {0}'.format(self.name))
        #print ('acentric_factor = {0}'.format(self.properties['Omega']))
        #print ('critical_Temperature = {0}'.format(self.properties['Tc']))
        #print ('critical_Pressure = {0}'.format(self.properties['Pc']))
        #print ('compressibility_factor_Z = {0}'.format(self.properties['Zc']))
        #print(self.label.keys())
        #print(self.label.values())

        return self.name, self.properties


#component_eos_list = np.zeros( (len(components),4) )
#dinputs = np.zeros( (len(components),4) )


#for index, component in enumerate(components):

#    properties_component = properties_data.selec_component(dppr_file, component)
#    pt.print_properties_component(component, properties_component)
#    dinputs[index] = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
#                        properties_component[1]['Omega'], properties_component[1]['Vc']])



#components_table = pd.DataFrame(component_eos_list, index=components, columns=['ac', 'b', 'rm', 'del1'])

#print(dinputs)
