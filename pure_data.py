import pandas as pd


class Data_parse(object):
	"""
	Data_parse read the dppr_fiel to get_values properties from a component
	"""

	def read_dppr(self, dppr_fiel):
		self.dppr_data = pd.read_excel(dppr_fiel).head().set_index('Name').ix[:, 1:12]
		#component_names = dppr_data.index.get_values()

		return self.dppr_data

	def selec_component(self, dppr_fiel, component):
		self.name = str(component)

		self.properties = self.read_dppr(dppr_fiel).loc[self.name]

		print(self.properties.index)

		print ('name = {0}'.format(self.name))
		print ('acentri_factor = {0}'.format(self.properties['Omega']))
		print ('critical_Temperature = {0}'.format(self.properties['Tc']))
		print ('critical_Pressure = {0}'.format(self.properties['Pc']))
		print ('compressibility_factor_Z = {0}'.format(self.properties['Zc']))

		return self.properties


dppr_fiel = "PureFull.xls"
component = 'METHANE'

component_eos = Data_parse()
component_eos.selec_component(dppr_fiel, component)

#component = 'ETHANE'
#selec_component_cal(component)


