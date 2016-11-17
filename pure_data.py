import pandas as pd


class Data_parse(object):
    """
    Data_parse read the dppr_file to get_values properties from a component
    """
    
    def read_dppr(self, dppr_file):
        #self.dppr_data = pd.read_excel(dppr_file).head().set_index('Name').ix[:, 1:12]

        self.dppr_data = pd.read_excel(dppr_file).set_index("Name").ix[:, 1:12]#.head().ix[:, 1:12]
        

        #component_names = dppr_data.index.get_values()
        return self.dppr_data
        
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


#dppr_file = "PureFull.xls"
#component = 'METHANE'

#component_eos = Data_parse()
#component_eos.selec_component(dppr_file, component)

#component = 'ETHANE'
#selec_component_cal(component)


def main():
    
    dppr_file = "PureFull.xls"
    #component = "METHANE"
    #component = "ISOBUTANE"
    component = "TRIPHENYLMETHANE"
    #component = "PYRENE"
    
    component_eos = Data_parse()

    #print(component_eos.read_dppr(dppr_file).index.get_values())
    #print(component_eos.read_dppr(dppr_file))
    # index.get_values())
    
    component_eos.selec_component(dppr_file, component)
    print(component_eos.selec_component(dppr_file, component))


if __name__ == "__main__":
    # execute only if run as a script
    main()
