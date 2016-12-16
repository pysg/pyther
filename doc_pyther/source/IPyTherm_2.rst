5b. PyTher: graphical interface to solid-fluid equilibrium
**********************************************************
**********************************************************


Importar las librerías
======================

Cargar la tabla de datos
========================

.. code:: python

    import scipy as sp
    from scipy import optimize
    from scipy.optimize import fsolve
    import numpy as np
    from matplotlib import pyplot
    %matplotlib inline
    import pandas as pd
    from numpy import linalg as LA
    from IPython.html import widgets
    from IPython.display import display
    from IPython.display import clear_output
    
    # encoding: utf-8
    
    from pandas import read_csv
    



.. parsed-literal::

    /home/andres-python/anaconda3/lib/python3.5/site-packages/IPython/html.py:14: ShimWarning: The `IPython.html` package has been deprecated. You should import from `notebook` instead. `IPython.html.widgets` has moved to `ipywidgets`.
      "`IPython.html.widgets` has moved to `ipywidgets`.", ShimWarning)


Cargar la tabla de datos
========================

.. code:: python

    f = pd.read_excel("PureFull.xls")
    f.head()
    data2 = pd.DataFrame(f)
    data2 = data2.set_index('Name')
    data2 = data2.ix[:, 1:12]
    Etiquetas = data2.index.get_values()
    Etiquetas




.. parsed-literal::

    array(['METHANE', 'ETHANE', 'PROPANE', ..., 'TITANIUM', 'PHOSPHORUS',
           'PHOSPHORUS'], dtype=object)



.. code:: python

    Componentes_1 = widgets.SelectMultiple(
                description="Component 1",
                options=list(Etiquetas))
    display(Componentes_1)

.. code:: python

    class Thermophysical_Properties():
        
        def __init__(self, nameData):
            self.nameData = nameData
            
        
        def cargar_Datos(self):
            f = pd.read_excel(self.nameData)
            f.head()
            data2 = pd.DataFrame(f)
            data2 = data2.set_index('Name')
            data2 = data2.ix[:, 1:12]
            self.Etiquetas = data2.index.get_values()
            
            print("Los datos del archivo: {0}, se han cargado correctamente !!!".format(self.nameData))
            
            return self.Etiquetas
        
        def seleccionar_Datos(self):
            Componentes_1 = widgets.SelectMultiple(
                description="Component 1",
                options=list(Etiquetas))
            display(Componentes_1)
            
            
        
        def mostrar_Datos(self):
            print ("Nombre componente: {0}".format(self.Etiquetas))
            
    
    
            
        
        def agregar_Datos(self):
            pass
        
        def borrar_Datos(self):
            pass
        
        def modificar_Datos(self):
            pass
        
        def crear_Datos(self):
            pass

.. code:: python

    nameData = "PureFull.xls"
    
    propiedades = Thermophysical_Properties(nameData)
    propiedades.cargar_Datos()
    propiedades.mostrar_Datos()
    propiedades.seleccionar_Datos()



.. parsed-literal::

    Los datos del archivo: PureFull.xls, se han cargado correctamente !!!
    Nombre componente: ['METHANE' 'ETHANE' 'PROPANE' ..., 'TITANIUM' 'PHOSPHORUS' 'PHOSPHORUS']


.. code:: python

    
    
        
        
        
    class System_Conditions():
        
        def __init__(self, Temperature, Pressure, Volume, Mole_fraction, Model_fluid, Model_solid):
            self.Temperature = Temperature
            self.Mole_fraction = Mole_fraction
            pass
        
        def normalizar(self):
            self.Mole_fraction_normal = Mole_fraction / sum(self.Mole_fraction)
            return self.Mole_fraction_normal
        
        def convertir(self):
            pass
            
            
        
    class Componentes(Thermophysical_Properties, System_Conditions):
        """
        Las variables aux_ se utilizan para presentar de forma más clara y acotada
        las expresiones necesarias en los calculos. Estas, se numeran de acuerdo al orden de
        aparición dentro de una clase.
        
        """
        
        def __init__(self):
            pass
        
        def cal_SRK_model(self):
            # Soave-Redlich-Kwong (SRK)
            self.s1, self.s2 = 1, 2
            self.m = 0.480 + 1.574 * self.w - 0.175 * self.w ** 2
            self.ac = 0.077796070 * self.R ** 2, self.Tc ** 2 / self.Pc
            self.bc = 0.086640 * self.R * self.Tc / self.Pc        
            
            return self.m, self.ac, self.bc
        
        def cal_PR_model(self):
            # Peng-Robinson (PR)
            self.s1, self.s2 = 1 + 2 ** 0.5, 1 - (2 ** 0.5)
            self.m = 0.37464 + 1.54226 * self.w - 0.26992 * self.w ** 2
            self.ac = 0.45723553 * self.R ** 2 * self.Tc ** 2 / self.Pc
            self.bc = 0.077796070 * self.R * self.Tc / self.Pc            
            
            self.alfa = (1 + self.m * (1 - (self.T / self.Tc) ** 0.5)) ** 2
            aux_1 = - (self.m / self.T) * (self.T / self.Tc) ** 0.5
            aux_2 = (self.m * (- (self.T / self.Tc) ** 0.5 + 1) + 1)
            self.dalfadT = aux_1 * aux_2
            
            aux_3 = 0.5 * self.m ** 2 * (self.T / self.Tc) ** 1.0 / self.T ** 2
            aux_4 = (self.m * (- (self.T / self.Tc) ** 0.5 + 1) + 1) / self.T ** 2
            aux_5 = 0.5 * self.m * (self.T / self.Tc) ** 0.5 * aux_4
            
            self.d2alfaT2 = aux_3 + aux_5
            
            self.a_ii = self.ac * self.alfa
            self.b_ii = self.bc
            self.da_iidT = self.ac * self.dalfadT
            d2adT2_puros = self.ac * self.d2alfaT2
            
            return self.m, self.a_ii, self.b_ii
        
        def cal_RKPR_model(self):
            pass
        
        def build_component(self):
            
            if self.eq == "SRK":
                # Soave-Redlich-Kwong (SRK)
                self.component = self.cal_SRK_model()
            elif self.eq == "PR":
                # Peng-Robinson (PR)
                self.component = self.cal_PR_model()
            elif self.eq == "RKPR":
                # (RKPR)
                #self.component = self.cal_RKPR_model()
                print ("No actualizada, intentalo de nuevo !!! ")
            else:
                print ("Che boludo... Modelo no valido, intentalo de nuevo !!! ")
                
                
                

.. code:: python

    Componentes_1 = widgets.SelectMultiple(
        description="Component 1",
        options=list(Etiquetas))
    
    Componentes_2 = widgets.SelectMultiple(
        description="Component 2",
        options=list(Etiquetas))
    
    button = widgets.Button(description="Upload Data")
    
    def cargarDatos(b):
        clear_output()
        print("Component 1: ", Componentes_1.value)
        Nombre = Componentes_1.value
        Propiedades = data2.loc[Nombre]
        Factor_Acentrico_1 = Propiedades[0]
        Temperatura_Critica_1 = Propiedades[1]
        Presion_Critica_1 = Propiedades[2]
        Z_Critico_1 = Propiedades[3]
    
        #print(Propiedades)
        print ("Acentric Factor = ", Factor_Acentrico_1)
        print ("Critical Temperature = ", Temperatura_Critica_1, "K")
        print ("Critical Pressure = ", Presion_Critica_1, "bar")
        print ("Z_Critical = ", Z_Critico_1, "\n")
        
        
        print("Component 2: ", Componentes_2.value)
        Nombre = Componentes_2.value
        Propiedades = data2.loc[Nombre]
        Factor_Acentrico_2 = Propiedades[0]
        Temperatura_Critica_2 = Propiedades[1]
        Presion_Critica_2 = Propiedades[2]
        Z_Critico_2 = Propiedades[3]
    
        #print(Propiedades)
        print ("Acentric Factor  = ", Factor_Acentrico_2)
        print ("Critical Temperature = ", Temperatura_Critica_2, "K")
        print ("Critical Pressure = ", Presion_Critica_2, "bar")
        print ("Z_Critical = ", Z_Critico_2)
        
        global TcDato, PcDato, wDato
        
        TcDato = np.array([Temperatura_Critica_1, Temperatura_Critica_2])
        PcDato = np.array([Presion_Critica_1, Presion_Critica_2])
        wDato = np.array([Factor_Acentrico_1, Factor_Acentrico_2])
    
    
    button.on_click(cargarDatos)

.. code:: python

    class Parameters_BD():
        
        def __init__(self):
            pass
        
        def cal_parameters_ij(self):           
    
            if self.nC > 1:
                self.aij = np.ones((len(self.ni), len(self.ni)))
                self.bij = np.ones((len(self.ni), len(self.ni)))
                self.daijdT = np.ones((len(self.ni), len(self.ni)))
    
                for j in range(self.nC):
                    for i in range(self.nC):
                        self.aij[i, j] = (self.a_ii[i] * self.a_ii[j]) ** 0.5
                        self.bij[i, j] = (self.b_ii[i] + self.b_ii[j]) / 2
                        self.bij[i, j] = self.bij[i, j]
                        self.daijdT[i, j] = (self.da_iidT[i] * self.da_iidT[j]) ** 0.5
    
                for i in range(self.nC):
                    for  j in range(self.nC):
                        if i == j:
                            self.aij[i, j] = self.a_ii[i] * (1 - self.kij[i, j])
                            self.daijdT[i, j] = self.da_iidT[i] * (1 - self.kij[i, j])
                        elif i != j:
                            self.aij[i, j] = self.aij[i, j] * (1 - self.kij[i, j])
                            self.daijdT[i, j] = self.daijdT[i, j] * (1 - self.kij[i, j])
                       
            if self.nC == 1:
                return self.a_ii, self.b_ii, self.da_iidT
            else:
                return self.aij, self.bij, self.daijdT
    
        def cal_parameter_D(self):
            if self.nC == 1:
                self.D = self.ni ** 2 * self.a_ii
                self.Di = 2 * self.ni * self.a_ii
            else:
                di = np.ones((len(self.ni), len(self.ni)))
                self.Di = np.ones((len(self.ni)))
                self.D = np.ones((len(self.ni)))
                for i in range(self.nC):
                    for j in range(self.nC):
                        di[i, j] = self.ni[j] * self.aij[i, j]
                        self.Di[i] = 2 * np.sum(di[i, :])
                self.D = 0.5 * np.sum(self.ni * self.Di)
    
            return self.D
        
        def cal_parameter_delta_1(self):
            
            if self.nC == 1:
                self.D1m = np.zeros((len(self.ni)-1))
                self.dD1i = np.ones((len(self.ni)))
                self.dD1ij = np.ones((len(self.ni), len(self.ni)))
                
                for i in range(self.nC):
                    self.D1m = self.D1m + self.ni[i] * self.delta_1[i]
                
                self.D1m = self.D1m / self.nT
                
            else:
                self.D1m = np.zeros((len(self.ni)-1))
                self.dD1i = np.ones((len(self.ni)))
                self.dD1ij = np.ones((len(self.ni), len(self.ni)))
                
                for i in range(self.nC):
                    self.D1m = self.D1m + self.ni[i] * self.delta_1[i]
                
                self.D1m = self.D1m / self.nT
                
                for i in range(self.nC):
                    self.dD1i[i] = (self.delta_1[i] - self.D1m) / self.nT
                    for j in range(self.nC):
                        self.dD1ij[i,j] = (2.0 * self.D1m - self.delta_1[i] - self.delta_1[j]) / self.nT ** 2
                        
            return self.D1m, self.dD1i, self.dD1ij
                        
        def cal_parameter_B(self):
            if self.nC == 1:
                self.B = self.ni * self.b_ii
            else:
                self.aux = np.zeros((len(self.ni)))
                for i in range(self.nC):
                    for j in range(self.nC):
                        self.aux[i] = self.aux[i] + self.ni[j] * self.bij[i, j]
    
                self.B = np.sum(self.ni * self.b_ii)
                #print("B = ", self.B)
    
            return self.B


.. code:: python

    class Fugacidad():
    
        def __init__(self, eq, w, Tc, Pc, Tr, R, ep, ni, nT, nC, V, T, P, kij, lij, delta_1, k, Avsl):
            self.eq = eq
            self.w = w
            self.Tc = Tc
            self.Pc = Pc
            self.Tr = Tr
            self.R = R
            self.ep = ep
            self.ni = ni
            self.nT = nT
            self.nC = nC
            self.V = V
            self.T = T
            self.P = P
            self.kij = kij
            self.lij = lij
            self.delta_1 = delta_1
            self.k = k
            self.Avsl = Avsl
            
            if self.eq == "SRK":
                # Soave-Redlich-Kwong (SRK)
                self.s1, self.s2 = 1, 2
                self.m = 0.480 + 1.574 * self.w - 0.175 * self.w ** 2
                self.ac = 0.077796070 * self.R ** 2, self.Tc ** 2 / self.Pc
                self.bc = 0.086640 * self.R * self.Tc / self.Pc
            elif self.eq == "PR":
                # Peng-Robinson (PR)
                self.s1, self.s2 = 1 + 2 ** 0.5, 1 - (2 ** 0.5)
                self.m = 0.37464 + 1.54226 * self.w - 0.26992 * self.w ** 2
                self.ac = 0.45723553 * self.R ** 2 * self.Tc ** 2 / self.Pc
                self.bc = 0.077796070 * self.R * self.Tc / self.Pc            
               
                self.alfa = (1 + self.m * (1 - (self.T / self.Tc) ** 0.5)) ** 2
                self.dalfadT = - (self.m / self.T) * (self.T / self.Tc) ** 0.5 * (self.m * (- (self.T / self.Tc) ** 0.5 + 1) + 1)
                ter_1 = 0.5 * self.m ** 2 * (self.T / self.Tc) ** 1.0 / self.T ** 2
                ter_2 = 0.5 * self.m * (self.T / self.Tc) ** 0.5 * (self.m * (- (self.T / self.Tc) ** 0.5 + 1) + 1) / self.T ** 2
                
                self.d2alfaT2 = ter_1 + ter_2
                self.a_ii = self.ac * self.alfa
                self.b_ii = self.bc
                
                self.da_iidT = self.ac * self.dalfadT
                d2adT2_puros = self.ac * self.d2alfaT2
    
            elif self.eq == "RKPR":
                # (RKPR)
                print ("No actualizada, intentalo de nuevo !!! ")            
    
            else:
                print ("Che boludo... Modelo no valido, intentalo de nuevo !!! ")
    
    
        def parametros(self):           
    
            if self.nC > 1:
                self.aij = np.ones((len(self.ni), len(self.ni)))
                self.bij = np.ones((len(self.ni), len(self.ni)))
                self.daijdT = np.ones((len(self.ni), len(self.ni)))
    
                for j in range(self.nC):
                    for i in range(self.nC):
                        self.aij[i, j] = (self.a_ii[i] * self.a_ii[j]) ** 0.5
                        self.bij[i, j] = (self.b_ii[i] + self.b_ii[j]) / 2
                        self.bij[i, j] = self.bij[i, j]
                        self.daijdT[i, j] = (self.da_iidT[i] * self.da_iidT[j]) ** 0.5
    
                for i in range(self.nC):
                    for  j in range(self.nC):
                        if i == j:
                            self.aij[i, j] = self.a_ii[i] * (1 - self.kij[i, j])
                            self.daijdT[i, j] = self.da_iidT[i] * (1 - self.kij[i, j])
                        elif i != j:
                            self.aij[i, j] = self.aij[i, j] * (1 - self.kij[i, j])
                            self.daijdT[i, j] = self.daijdT[i, j] * (1 - self.kij[i, j])
                       
            if self.nC == 1:
                return self.a_ii, self.b_ii, self.da_iidT
            else:
                return self.aij, self.bij, self.daijdT
    
        def parametro_D(self):
            if self.nC == 1:
                self.D = self.ni ** 2 * self.a_ii
                self.Di = 2 * self.ni * self.a_ii
            else:
                di = np.ones((len(self.ni), len(self.ni)))
                self.Di = np.ones((len(self.ni)))
                self.D = np.ones((len(self.ni)))
                for i in range(self.nC):
                    for j in range(self.nC):
                        di[i, j] = self.ni[j] * self.aij[i, j]
                        self.Di[i] = 2 * np.sum(di[i, :])
                self.D = 0.5 * np.sum(self.ni * self.Di)
    
            return self.D
        
        def parametro_delta_1(self):
            
            if self.nC == 1:
                self.D1m = np.zeros((len(self.ni)-1))
                self.dD1i = np.ones((len(self.ni)))
                self.dD1ij = np.ones((len(self.ni), len(self.ni)))
                
                for i in range(self.nC):
                    self.D1m = self.D1m + self.ni[i] * self.delta_1[i]
                
                self.D1m = self.D1m / self.nT
                
            else:
                self.D1m = np.zeros((len(self.ni)-1))
                self.dD1i = np.ones((len(self.ni)))
                self.dD1ij = np.ones((len(self.ni), len(self.ni)))
                
                for i in range(self.nC):
                    self.D1m = self.D1m + self.ni[i] * self.delta_1[i]
                
                self.D1m = self.D1m / self.nT
                
                for i in range(self.nC):
                    self.dD1i[i] = (self.delta_1[i] - self.D1m) / self.nT
                    for j in range(self.nC):
                        self.dD1ij[i,j] = (2.0 * self.D1m - self.delta_1[i] - self.delta_1[j]) / self.nT ** 2
                        
            return self.D1m, self.dD1i, self.dD1ij
                        
        def parametro_B(self):
            if self.nC == 1:
                self.B = self.ni * self.b_ii
            else:
                self.aux = np.zeros((len(self.ni)))
                for i in range(self.nC):
                    for j in range(self.nC):
                        self.aux[i] = self.aux[i] + self.ni[j] * self.bij[i, j]
    
                self.B = np.sum(self.ni * self.b_ii)
                #print("B = ", self.B)
    
            return self.B
    
        def presion(self):
            '''
            Con el metodo presion(), se calcula la Presión P(T, V, N) del sistema
            para una temperatura T, cantidad de moles N y un volumen V
            R = Constante universal de los gases
            nT = Número total de moles en el sistema
            Pcal = Peos = Presión calculada con la ecuación de estado
            Arv = Primera derivada parcial de la energía de Helmholz con respecto al
            volumen V, a T y N constantes
            '''
            self.gv = self.R * self.B / (self.V * (self.V - self.B))
            self.fv = - 1 / ((self.V + self.s1 * self.B) * (self.V + self.s2 * self.B))
            self.ArV = -self.nT * self.gv * self.T - self.D * self.fv
            self.Pcal = self.nT * self.R * self.T / self.V - self.ArV   
    
            return self.Pcal
    
        def dP_dV(self):
            self.dPdV = -self.ArV2 - self.R * self.T * self.nT / self.V ** 2
            return self.dPdV
    
        def Z_factor(self):
            self.Z = (self.P * self.V) / (self.nT * self.R * self.T)
            return self.Z
    
        def P_ideal(self):
            self.Pxi = (self.ni * self.P) / self.nT
            return self.Pxi
    
        def dF_dV(self):
            '''
            Primera derivada de F con respecto al volumen Ecu. (68)
            '''
            self.gv = self.R * self.B / (self.V * (self.V - self.B))
            self.fv = - 1 / ((self.V + self.s1 * self.B) * (self.V + self.s2 * self.B))
            self.ArV = -self.nT * self.gv * self.T - self.D * self.fv
            return self.ArV
    
        def dF_dVV(self):
            '''
            Segunda derivada de F con respecto al volumen Ecu. (74)
            '''
            self.gv2 = self.R * (1 / self.V ** 2 - 1 / (self.V - self.B) ** 2)
            self.fv2 = (- 1 / (self.V + self.s1 * self.B) ** 2 + 1 / (self.V + self.s2 * self.B) ** 2) / self.B / (self.s1 - self.s2)
            self.ArV2 = - self.nT * self.gv2 * self.T - self.D * self.fv2
            return self.ArV2
    
        def volumen_1(self):
            '''
            Calculo del volumen V(T,P,n) del fluido a una temperatura T, presión P
            y número de moles totales nT especificados.
            Se utiliza el método de Newton con derivada de la función analitica.
            Pendiente cambiar por una función de Scipy.
            '''
            self.V = 1.05 * self.B          
            lnP = np.log(self.P)
            Pite = self.presion()
            lnPcal = np.log(Pite)
            h = lnP - lnPcal
            errorEq = abs(h)
            i = 0
            s = 1.0
    
            while errorEq > self.ep:
                self.parametro_D()
                self.parametro_B()
                self.dF_dV()
                self.dF_dVV()
                dPite = self.dP_dV()
                Pite = self.presion()
                lnPcal = np.log(Pite)
                h = lnP - lnPcal
                dh = -dPite
                self.V = self.V - s * h / dh
                errorEq = abs(h)
                i += 1
                if i >= 900:
                    pass
                    #break
    
            return self.V
    
        def funcion_energia_F(self):
            self.g = self.R * np.log(1 - self.B / self.V)
            self.bv = self.B / self.V
            self.f = np.log((self.V + self.s1 * self.B) / (self.V + self.s2 * self.B)) / self.B / (self.s1 - self.s2)
            self.Ar = -self.nT * self.g * self.T - self.D * self.f
            return self.g, self.f, self.Ar, self.bv
        
        def tomar_B(self):
            print ("tomando B =", self.B)
            return self.B + 10
        
        def derivadas_delta_1(self):
            auxD2 = (1 + 2 / (1 + self.s1) ** 2)
            
            como_1 = (1 / (self.V + self.s1 * self.B) + 2 / (self.V + self.s2 * self.B) / (1 + self.s1) ** 2)
            como_2 = self.f * auxD2
            self.fD1 = como_1 - como_2
            self.fD1 = self.fD1/(self.s1 - self.s2)
            
            return self.fD1
    
        def primeras_derivadas1(self):
    
            if self.nC == 1:
                AUX = self.R * self.T / (self.V - self.B)
                self.fB = -(self.f + self.V * self.fv) / self.B
                self.FFB = self.nT * AUX - self.D * self.fB
                self.Di = 2 * self.nT * self.ac * self.alfa
                self.Bi = self.bc
                
                if self.eq != "RKPR":
                    self.Arn = -self.g * self.T + self.FFB * self.Bi - self.f * self.Di
                else:
                    self.Arn = -self.g * self.T + self.FFB * self.Bi - self.f * self.Di - self.D * self.fD1 * self.dD1i
            else:
                # Derivando la ecuación (64) se obtiene la ecuación eq (106)
                self.Bi = np.ones((len(self.ni)))
                for i in range(self.nC):
                    self.Bi[i] = (2 * self.aux[i] - self.B) / self.nT
    
                AUX = self.R * self.T / (self.V - self.B)
                self.fB = -(self.f + self.V * self.fv) / self.B
                self.FFB = self.nT * AUX - self.D * self.fB
                
                if self.eq != "RKPR":
                    self.Arn = -self.g * self.T + self.FFB * self.Bi - self.f * self.Di
                else:
                    auxD2 = (1 + 2 / (1 + self.s1) ** 2)
                    print("B delta1 = ", self.B)
                    co_1 = (1 / (self.V + self.s1 * self.B) + 2 / (self.V + self.s2 * self.B) / (1 + self.s1) ** 2)
                    co_2 = self.f * auxD2
                    self.fD1 = co_1 - co_2
                    self.fD1 = self.fD1/(self.s1 - self.s2)
                    self.Arn = -self.g * self.T + self.FFB * self.Bi - self.f * self.Di - self.D * self.fD1 * self.dD1i
    
            return self.Arn, self.Arn, self.Arn
        
    
        def coeficientes_fugacidad(self):
            self.Z = self.Z_factor()
            self.lnOi = self.Arn / (self.R * self.T) - np.log(self.Z)
            self.Oi = np.exp(self.lnOi)
            return self.Oi
    
        def fugacidad(self):
            self.Z = self.Z_factor()
            self.Pxi = self.P_ideal()
            self.lnFi = self.Arn / (self.R * self.T) - np.log(self.Z) + np.log(self.Pxi)
            self.Fi = np.exp(self.lnFi)
            self.PHILOG = self.Arn / (self.R * self.T) - np.log(self.Z)
            self.PHILOG_i = self.Arn - np.log(self.Z)
            self.FUGLOG = self.Arn / (self.R * self.T) + np.log(self.ni) + np.log((self.nT * self.R * self.T) / self.V)
            return self.Fi
    
        def exp_sol(self):
            '''
            Este método calcula el factor de corrección de la fugacidad del
            componente fluido para determinar la fugacidad del mismo componente
            en estado sólido.
            Fugacidad del sólido puro
            fi_s(T, P) = fi_l(T, P) * EXP(T, P)
            '''
            Tfus = 323.75
            # Temperatura de fusion de n-tetracosane
            # Unidad de Ti_f en Kelvin
            par_sol = np.array([[-176120.0, 8196.20, -55.911, 0.19357, -0.0002235],
                                [-1.66e6, 8.31e3, 0.0, 0.0, 0.0]])
            par_liq = np.array([[423160.0, 1091.9, 0.0, 0.0, 0.0],
                                [7.01e5, 1.47e3, 0.0, 0.0, 0.0]])
            #print ("par_sol", par_sol)
            #print ("par_liq", par_liq)
            # Las unidades de Cp están en J/Kmol.K
            Cp_solido = par_sol[:, 0] + par_sol[:, 1] * T + par_sol[:, 2] * T ** 2 + par_sol[:, 3] * T ** 3 + par_sol[:, 4] * T ** 4
            #print ("Cp_solido", Cp_solido)
            Cp_liquido= par_liq[:, 0] + par_liq[:, 1] * T + par_liq[:, 2] * T ** 2 + par_liq[:, 3] * T ** 3 + par_liq[:, 4] * T ** 4
            #print ("Cp_liquido", Cp_liquido)
            DeltaCp = (Cp_solido - Cp_liquido) * (1.0 / 1000)
            print ("Delta Cp", DeltaCp)
    
            #Unidades de Delta H de fusión en Kcal/mol
            DeltaH_f = np.array([13.12, 21.23]) * (1000 / 1.0) * (4.18 / 1.0)
            #print ("Delta H de fusion", DeltaH_f)
            T_f = np.array([323.75, 349.05])
            #print ("Temperaturas de fusion = ", T_f)
    
            Rp = 8.314
            A = (DeltaH_f / (Rp * Tfus)) * (1 - (Tfus / T))
            B = (DeltaCp / Rp) * (1 - (Tfus / T))
            C = (DeltaCp / Rp) * np.log(Tfus / T)
            self.EXP = np.exp(A - B - C)
    
            print ("A = ", A)
            print ("B = ", B)
            print ("C = ", C)
            print ("EXP = ", self.EXP)
    
            return self.EXP
        
        def exp_sol_1(self):
            '''
            Este método calcula el factor de corrección de la fugacidad del
            componente fluido para determinar la fugacidad del mismo componente
            en estado sólido.
            Fugacidad del sólido puro
            fi_s(T, P) = fi_l(T, P) * EXP(T, P)
            '''
            Tpt = 323.75
            Ppt = 1.38507E-8
            R = 8.314472
            AH = 54894000
            Av = -0.0376300841 #m3/kmol
            
            a = ((AH / (R * Tpt)) * (1 - (Tpt / self.T))) / 1000
            b = ((Av / (R * self.T)) * (self.P - Ppt)) * 100
            self.EXP_1 = a + b
            
            return self.EXP_1
        
        def exp_sol_3(self):
            '''
            Este método calcula el factor de corrección de la fugacidad del
            componente fluido para determinar la fugacidad del mismo componente
            en estado sólido.
            Fugacidad del sólido puro
            fi_s(T, P) = fi_l(T, P) * EXP(T, P)
            '''
            # [=] K
            # [=] bar
            # [m3 / Kmol]
            # Constante R [=] 0.08314472 bar.l/(mol.K)
            
            Tpt = 323.75
            Ppt = 3.2015002E-8
            #self.Avsl = -0.0565500835
            
            c1 = -14213.5004
            c2 = 605153.4382
            c3 = -591592.556
            
            R = 0.08314472
            
            A1 = c1 * (1 - Tpt / self.T)
            A2 = c2 * (-1 + Tpt / self.T + np.log(self.T / Tpt))
            A3 = c3 * (-1 + self.T / (2 * Tpt) + Tpt / (2 * self.T)) + (Tpt / self.T) * (self.P - Ppt)
            
            FE = (self.Avsl / (self.R * self.T)) * (A1 + A2 + A3)
            self.EXP_3 = np.exp(FE)    
            return self.EXP_3
            
    
        def fluido(self):
            ab = self.parametros()
            D = self.parametro_D()
            B = self.parametro_B()
            Vol_1 = self.volumen_1()
            F = self.funcion_energia_F()
            dF = self.primeras_derivadas1()
            Z = self.Z_factor()
            Zcal = (self.P * Vol_1) / (self.nT * self.R * self.T)
            Pq = self.presion()
            self.Fug = self.fugacidad()
            self.CoeFug = self.coeficientes_fugacidad()
            return self.Fug
    
        def solido(self):
            if self.nC == 1:
                Fug = self.fluido()
                #EXP = self.exp_sol()
                #EXP = self.exp_sol_1()
                EXP = self.exp_sol_3()
                
                FugS = Fug[0] * EXP
            else:
                print ("Aún no se qué hacer para una mezcla de sólidos !!!")
                FugS = 1
    
            return FugS
        
    #----------------
    def calculaFugacidad(x, Pe, nif, nCf, eq, TcDato, PcDato, wDAto, Avsl):
        #---------------------------------------------------------------------------
        # Temperatura en [=] K
        # Presión en [=] bar
        # Constante R [=] 0.08314472 bar.l/(mol.K)
        # x = variable que se cálcula, puede ser T ó P para el equilibrio sólido-fluido
        # Pe = Presión del sistema especificada
        # nif = número de moles del componente (i) en cada fase (f)
        # nCf = número de componentes en una fase (f)
        # eq = modelo de ecuación de estado, SRK, PR, RKPR
        # TcDato = Temperatura critica de la "base de datos"
        # PcDato = Presión critica de la "base de datos"
        # wDato = Factor acentrico de la "base de datos"
        # Avsl = Delta de volumen sólido-fluido
       
        # ep = Criterio de convergencia del método def volumen_1(self, P)
        
        T = x # 335.42 # x # 366.78 # 356.429 # 335.42 # 348.89 #327.0
        #print("Temperatura = ", T)
        P = Pe # 2575.0 # 2064.7 # 1524.4 #1164.2 # 865.0 
        # 560.3 # x #1054.6 #1560.3 # 2064.7 # 1524.4 # 560.3 # 1164.2 #865.0
        R = 0.08314472
        ep = 1e-5#1e-6
        #---------------------------------------------------------------------------    
        Tcm = TcDato
        Pcm = PcDato
        wm = wDato
            
        nC = nCf
        
        if nC == 1:
            #print ("...............................................................")
            
            #ni = nif
            ni = np.array([1.0])
            
            #print ("Número de moles = ", ni)
            # C24
            kij = 0.0
            lij = 0.0
            
            # Metano - Etano
            delta_1 = np.array([0.85])        
            k = np.array([1.50758])
            #C24
            Tc = Tcm[1]
            Pc = Pcm[1]
            w = wm[1]     
            print ("...............................................................")
        elif nC == 2:
            # metano - C24
            #ni = np.array([1-nif, nif])
            ni = nif #np.array([1-nif, nif])
            
            #ni = np.array([1 - 0.901, 0.901])
            #---------------------------------
            
            
            #ni = np.array([1 - 0.26, 0.26])
            
            #ni = np.array([1 - 0.104, 0.104])
            #print ("Número de moles = ", ni)
    
            kij = np.array([[0.000000, 0.083860],
                            [0.083860, 0.000000]])
            
            kij = np.array([[0.000000, 0.059600],
                            [0.059600, 0.000000]])
            
    
            lij = 0.0132
            
            
            #kij = np.array([[0.000000, 0.00],
            #                [0.00, 0.000000]])
            
            #lij = 0.0
            
            
            # Metano - C24
            delta_1 = np.array([0.85, 2.40])        
            k = np.array([1.50758, 4.90224])
            
            # metano sigma1 = 0.9253, sigma = 0.85, k = 1.49345, k = 1.50758
            # C24 sigma = 2.40 k = 4.90224
    
            Tc = Tcm
            Pc = Pcm
            w = wm
            print ("Temperatura Critica = ", Tc, "K")
            print ("Presión Critica = ", Pc, "bar")
            print ("Factor Acentrico = ", w)
            #print ("...............................................................")
    
        # Tempertura reducidad
        Tr = T / Tc
        # C24 puro
        V = 0.141604834257319
        nT = np.sum(ni)
    
        fugacidad = Fugacidad(eq, w, Tc, Pc, Tr, R, ep, ni, nT, nC, V, T, P, kij, lij, delta_1, k, Avsl)
        
        print(fugacidad.exp_sol_3())
                
        if nC == 1:
            SOL = fugacidad.solido()
            
            return SOL
        else:
            flu_1 = fugacidad.fluido()
            return flu_1
    
    
    #----------------

.. code:: python

    def equilibrioSF(x, Pe, nif, n1, n2, Avsl):
        
        # fugacidad del sólido puro
        FugS = calculaFugacidad(x, Pe, nif, n1, eq, TcDato, PcDato, wDato, Avsl)
        print(eq, TcDato, PcDato, wDato, Avsl)
        # fugacidad del fluido pesado en la mezcla fluida
        FugF = calculaFugacidad(x, Pe, nif, n2, eq, TcDato, PcDato, wDato, Avsl)
        
        # Función de igualdad de fugacidades del sólido y el fluido
        eqSF = np.abs(np.abs(np.log(FugS)) - np.abs(np.log(FugF[1])))
        print ("-"*80)
        print ("ln(Fugacidad Sólido) = ", np.log(FugS))
        print ("ln(Fugacidad Fluido) = ", np.log(FugF[1]))
        print ("ln(Fugacidad Sólido) - ln(Fugacidad Fluido) = ", eqSF)   
        
        return eqSF
    
    
    eq = 'PR'
    #Avsl = -0.0565500835
    #Avsl = -0.09605965500835
    
    #initial_temperature = [346.5] # T [=] K
    #initial_pressure = 136.9 # [=] bar
    
    #Tcal = fsolve(equilibrioSF,initial_temperature,args=(initial_pressure, 1, 2, Avsl), xtol=1e-4)
    #print(Tcal, "K")
    
    t_exp = [323.65, 326.04, 326.43, 328.12, 329.45, 329.89, 333.43, 335.12, 340.19, 344.58, 346.65, 352.53, 362.45, 362.76, 371.82, 379.74]
    temp = np.array(t_exp)
    
    p_exp = [1, 101.0, 136.9, 183.8, 266.2, 266.8, 426.9, 480.3, 718.9, 912.5, 1010.6, 1277.8, 1778.0, 1825.1, 2323.4, 2736.1]
    pres= np.array(p_exp)
    
    pos = np.arange(len(pres))
    Tcal = np.ones((len(pres)))
    Tcal
    
    Tres = np.array([ 322.65861561,  324.91946742,  325.73456905,  326.80151121,
            328.68045402,  328.69415114,  332.3526483 ,  333.57248076,
            338.99640222,  343.33723415,  345.50684642,  351.28742799,
            361.49784425,  362.4145721 ,  371.63445321,  378.63493779])
    
    Tcal - temp 
    
    Avsl = -0.32595074
    Avsl
    
    
    class Flash():
    
        def __init__(self, zi_F, temperature_f, pressure_f, TcDato_f, PcDato_f, wDato_f):
            self.zi = zi_F
            self.T = temperature_f
            self.P = pressure_f
            self.Tc = TcDato_f
            self.Pc = PcDato_f
            self.w = wDato_f        
            
        def wilson(self):
            # Ecuación wilson
            lnKi = np.log(self.Pc / self.P) + 5.373 * (1 + self.w) * (1 - self.Tc / self.T)
            self.Ki = np.exp(lnKi)
            return self.Ki
    
        def beta(self):
            # Estimación de la fracción de fase de vapor en el sistema
            self.Ki = self.wilson()
            #Bmin = np.divide((self.Ki * self.zi - 1), (self.Ki - 1))
            Bmin = (self.Ki * self.zi - 1) / (self.Ki - 1)
            
            #print (("Bmin_inter = ", Bmin))
            
            Bmax = (1 - self.zi) / (1 - self.Ki)
            #print (("Bmax_inter = ", Bmax))
            self.Bini = (np.max(Bmin) + np.min(Bmax)) / 2
            print("inib =", self.Bini)
            return self.Bini
    
        def rice(self):
            # Ecuación de Rachford-Rice para el equilibrio líqudo-vapor
            self.fg = np.sum(self.zi * (self.Ki - 1) / (1 - self.Bini + self.Bini * self.Ki))
            self.dfg = - np.sum(self.zi * (self.Ki - 1) ** 2 / (1 - self.Bini + self.Bini * self.Ki) ** 2)
            #print g, dg
            return self.fg, self.dfg
        
        def composicion_xy(self):
            # Ecuación de Rachford-Rice para calcular la composición del equilibrio líqudo-vapor
            self.xi = self.zi / (1 - self.Bini + self.Bini * self.Ki)
            self.yi = (self.zi * self.Ki) / (1 - self.Bini + self.Bini * self.Ki)
            self.li = (self.zi * (1 - self.Bini)) / (1 - self.Bini + self.Bini * self.Ki)
            self.vi = (self.zi * self.Bini * self.Ki) / (1 - self.Bini + self.Bini * self.Ki)
    
            return self.xi, self.yi, self.li, self.vi
    
        def flash_ideal(self):
            # Solución del flash (T,P,ni) isotermico para Ki_(T,P)
            self.Bini = self.beta()
            self.Ki = self.wilson()
            # print ("Ki_(P, T) = ", self.Ki)
            Eg = self.rice()
            errorEq = abs(Eg[0])
            # Especificaciones del método Newton precario, mientras se cambia por una librería Scipy
            i, s, ep = 0, 1, 1e-5
    
            while errorEq > ep:
                Eg = self.rice()
                self.Bini = self.Bini - s * Eg[0] / Eg[1]
                errorEq = abs(Eg[0])
                i += 1
                if i >= 50:
                    break
    
            xy = self.composicion_xy()        
            print ("-"*53, "\n", "-"*18, "Mole fraction", "-"*18, "\n","-"*53)
            print ("\n", "-"*13, "Zi phase composition", "-"*13, "\n")
            print ("{0} = {1} \n {2} = {3} \n {4}={5} \n {6}={7} \n".format(Componentes_f1.value, self.zi[0], Componentes_f2.value, self.zi[1], Componentes_f3.value, self.zi[2], Componentes_f4.value, self.zi[3]))
            print ("Sumatoria zi = {0}".format(np.sum(self.zi)))       
            print ("\n", "-"*13, "Liquid phase composition", "-"*13, "\n")
            print ("{0} = {1} \n {2} = {3} \n {4}={5} \n {6}={7} \n".format(Componentes_f1.value, self.xi[0], Componentes_f2.value, self.xi[1], Componentes_f3.value, self.xi[2], Componentes_f4.value, self.xi[3]))
            print ("Sumatoria xi = {0}".format(np.sum(self.xi)))
            print ("\n", "-"*14, "Vapor phase composition", "-"*13, "\n")
            print ("{0} = {1} \n {2} = {3} \n {4}={5} \n {6}={7} \n".format(Componentes_f1.value, self.yi[0], Componentes_f2.value, self.yi[1], Componentes_f3.value, self.yi[2], Componentes_f4.value, self.yi[3]))
            print ("Sumatoria yi = {0}".format(np.sum(self.yi)))
            print ("-"*53, "\n","Beta = {0}".format(self.Bini), "\n")
            print ("\n","Función R&R = {0}".format(Eg[0]), "\n")
            print ("\n","Derivada función R&R = {0}".format(Eg[1]), "\n", "-"*53)
    
    
            return #Eg[0], Eg[1], self.Bini
        
    class FlashHP(Fugacidad, Flash):
    
        def __init__(self, zF):
            Fugacidad.__init__(self, eq, w, Tc, Pc, Tr, R, ep, ni, nT, nC, V, T, P, kij, lij, delta_1, k, Avsl)
            self.zF = zF
            
          
        
        def flash_PT(self):
            # Solución del flash (T,P,ni) isotermico para Ki_(T,P,ni)
            flashID = self.flash_ideal()
            print ("flash (P, T, zi)")
            print ("g, dg, B = ", flashID)
            print ("-"*66)
    
            self.Bini = flashID[2]
            print ("Beta_r ini = ", self.Bini)
            moles = self.composicion_xy()
    
            self.xi, self.yi = moles[0], moles[1]
            nil, niv = moles[2], moles[3]
    
            fi_F = self.fugac()        
    
            self.Ki = fi_F[0] / fi_F[1]
    
            L = 1.0
    
            self.Ki = self.Ki * L
    
            Ki_1 = self.Ki
            print ("Ki_(P, T, ni) primera = ", self.Ki)
    
            print ("-"*66)
    
            #self.Ki = np.array([1.729, 0.832, 0.640])
    
            #self.Ki = self.wilson(self.Pc, self.Tc, self.w, self.T)
            #print "Ki_(P, T) = ", self.Ki
    
            while 1:
                i, s = 0, 0.1
    
                while 1:
                    Eg = self.rice()
                    print (Eg)
                    self.Bini = self.Bini - s * Eg[0] / Eg[1]
                    print (self.Bini)
                    errorEq = abs(Eg[0])
                    i += 1
                    #print i
    
                    #if self. Bini < 0 or self.Bini > 1:
                        #break
                    #    self.Bini = 0.5
                    if i >= 50:
                        pass
                        #break
                    if errorEq < 1e-5:
                        break
    
                print ("Resultado Real = ", Eg)
                print (" Beta r = ", self.Bini)
    
                moles = self.composicion_xy(zi, self.Ki, self.Bini)
                self.xi, self.yi = moles[0], moles[1]
    
                #xy = self.composicion_xy(zi, self.Ki, self.Bini)
    
                print ("C1 -i-C4 n-C4")
                print ("-"*13, "Composición de fase líquida", "-"*13)
                print ("xi = ", moles[0])
                print ("Sxi = ", np.sum(moles[0]))
                print ("-"*13, "Composición de fase vapor", "-"*13)
                print ("yi = ", moles[1])
                print ("Syi = ", np.sum(moles[1]))
    
                fi_F = self.fugac()
    
                self.Ki = fi_F[0] / fi_F[1]
                Ki_2 = self.Ki
                dKi = abs(Ki_1 - Ki_2)
                Ki_1 = Ki_2
                print ("Ki_(P, T, ni) = ", self.Ki)
    
                fun_Ki = np.sum(dKi)
                print ("fun_Ki = ", fun_Ki)
    
                if fun_Ki < 1e-5:
                    break
    
            return flashID
    
    url = 'Lectura Juan.xlsx'
    
    class DataGPEC():
        
        def __init__(self, url):
            self.url = url
            
        def leerGPEC_1(self):
            """
            El siguiente script python, se puede mejorar generalizando la lectura de etiquetas,
            mientras se pasa la transición GPEC librería
            """
            marcas = ['VAP', 'CRI', 'CEP']
            
            GPEC = pd.read_excel(url)
            
            """
            Revisar las etiquetas, nombre, roturlos de las figurar generadas con este script Python
            para que sean acordes a las variables que se desean gráficar, mientras se automatiza este
            proceso.
            """
            
            
            #------------------------------------------------------------------------------
            DatosGPEC = pd.DataFrame(GPEC)
            VAP = DatosGPEC.loc[(DatosGPEC['T(K)'] == marcas[0])]
            etiquetaVAP = VAP.index.get_values()
            inicioVAP = etiquetaVAP[0]+1
            finalVAP = etiquetaVAP[1]-2
            
            #------------------------------------------------------------------------------
            self.TemperaturaVAP = np.array([DatosGPEC.ix[inicioVAP:finalVAP,0]], dtype=np.float)
            self.PresionVAP = np.array([DatosGPEC.ix[inicioVAP:finalVAP,1]], dtype=np.float)
            self.VolumenLiqVAP = np.array([DatosGPEC.ix[inicioVAP:finalVAP,2]], dtype=np.float)
            self.VolumenVapVAP = np.array([DatosGPEC.ix[inicioVAP:finalVAP,3]], dtype=np.float)
            #------------------------------------------------------------------------------
            CRI = DatosGPEC.loc[(DatosGPEC['T(K)'] == marcas[1])]
            etiquetaCRI = CRI.index.get_values()
            inicioCRI = etiquetaCRI[0]+1
            finalCRI = etiquetaCRI[1]-2
            #------------------------------------------------------------------------------
            self.TemperaturaCRI = np.array([DatosGPEC.ix[inicioCRI:finalCRI,0]], dtype=np.float)
            self.PresionCRI = np.array([DatosGPEC.ix[inicioCRI:finalCRI,1]], dtype=np.float)
            self.VolumenLiqCRI = np.array([DatosGPEC.ix[inicioCRI:finalCRI,2]], dtype=np.float)
            self.VolumenVapCRI = np.array([DatosGPEC.ix[inicioCRI:finalCRI,3]], dtype=np.float)
            #------------------------------------------------------------------------------
            """
            En la segunda línea critica se tiene como referencia el final de la primera línea critica
            y la etiqueta CEP
            """
            
            CEP = DatosGPEC.loc[(DatosGPEC['T(K)'] == marcas[2])]
            etiquetaCEP = CEP.index.get_values()
            inicioCRI_2 = etiquetaCRI[1]+1
            finalCRI_2 = etiquetaCEP[0]-2
            
            self.TemperaturaCRI_2 = np.array([DatosGPEC.ix[inicioCRI_2:finalCRI_2,0]], dtype=np.float)
            self.PresionCRI_2 = np.array([DatosGPEC.ix[inicioCRI_2:finalCRI_2,1]], dtype=np.float)
            self.VolumenLiqCRI_2 = np.array([DatosGPEC.ix[inicioCRI_2:finalCRI_2,2]], dtype=np.float)
            self.VolumenVapCRI_2 = np.array([DatosGPEC.ix[inicioCRI_2:finalCRI_2,3]], dtype=np.float)
            
            return self.TemperaturaCRI_2
        
        def presionVapor(self):
            clear_output()
            pyplot.close("all")
            pyplot.scatter(self.TemperaturaVAP,self.PresionVAP, color = 'red', label = 'Presión de Vapor')
            pyplot.title('Temperatura-Presión')
            pyplot.legend(loc="upper left") 
            pyplot.xlabel('Temperatura [=] K')
            pyplot.ylabel('Presión [=] bar')
            
        def densidadPresion(self):
            clear_output()
            pyplot.close("all")
            pyplot.scatter(self.VolumenLiqVAP,self.PresionVAP, color = 'red', label = 'Líquido')
            pyplot.scatter(self.VolumenVapVAP,self.PresionVAP, color = 'blue', label = 'Vapor')
            pyplot.title('Diagrama Densidad-Presión')
            pyplot.legend(loc="upper right") 
            pyplot.xlabel('Densidad [=] -')  
            pyplot.ylabel('Presión [=] bar')
            
        def diagramaTPcritico(self):
            clear_output()
            pyplot.close("all")
            pyplot.scatter(self.TemperaturaCRI,self.PresionCRI, color = 'red', label = 'Presión Critica')
            pyplot.title('Diagrama Temperatura Cri-Presión Cri')
            pyplot.legend(loc="upper left") 
            pyplot.xlabel('Temperatura [=] K')  
            pyplot.ylabel('Presión [=] bar')
            
        def diagramaDensidadCri(self):
            clear_output()
            pyplot.close("all")
            pyplot.scatter(self.VolumenLiqCRI,self.PresionCRI, color = 'red', label = 'Líquido')
            pyplot.scatter(self.VolumenVapCRI,self.PresionCRI, color = 'blue', label = 'Vapor')
            pyplot.title('Diagrama Densidad Critica')
            pyplot.legend(loc="upper right") 
            pyplot.xlabel('Densidad [=] -')  
            pyplot.ylabel('Presión [=] bar')
            
        def diagramaCritico_2(self):
            clear_output()
            pyplot.close("all")
            fig_2= pyplot.scatter(self.TemperaturaCRI_2,self.PresionCRI_2)
            pyplot.scatter(self.TemperaturaCRI_2,self.PresionCRI_2, color = 'red', label = 'Presión de Critica 2')
            pyplot.title('Diagrama Critico 2')
            pyplot.legend(loc="upper left") 
            pyplot.xlabel('Temperatura [=] K')  
            pyplot.ylabel('Presión [=] bar')
    #------------------------------------------------------------------------------    

Interfaz "gráfica"
==================

.. code:: python

    Componentes_1 = widgets.SelectMultiple(
        description="Component 1",
        options=list(Etiquetas))
    
    Componentes_2 = widgets.SelectMultiple(
        description="Component 2",
        options=list(Etiquetas))
    
    button = widgets.Button(description="Upload Data")
    
    def cargarDatos(b):
        clear_output()
        print("Component 1: ", Componentes_1.value)
        Nombre = Componentes_1.value
        Propiedades = data2.loc[Nombre]
        Factor_Acentrico_1 = Propiedades[0]
        Temperatura_Critica_1 = Propiedades[1]
        Presion_Critica_1 = Propiedades[2]
        Z_Critico_1 = Propiedades[3]
    
        #print(Propiedades)
        print ("Acentric Factor = ", Factor_Acentrico_1)
        print ("Critical Temperature = ", Temperatura_Critica_1, "K")
        print ("Critical Pressure = ", Presion_Critica_1, "bar")
        print ("Z_Critical = ", Z_Critico_1, "\n")
        
        
        print("Component 2: ", Componentes_2.value)
        Nombre = Componentes_2.value
        Propiedades = data2.loc[Nombre]
        Factor_Acentrico_2 = Propiedades[0]
        Temperatura_Critica_2 = Propiedades[1]
        Presion_Critica_2 = Propiedades[2]
        Z_Critico_2 = Propiedades[3]
    
        #print(Propiedades)
        print ("Acentric Factor  = ", Factor_Acentrico_2)
        print ("Critical Temperature = ", Temperatura_Critica_2, "K")
        print ("Critical Pressure = ", Presion_Critica_2, "bar")
        print ("Z_Critical = ", Z_Critico_2)
        
        global TcDato, PcDato, wDato
        
        TcDato = np.array([Temperatura_Critica_1, Temperatura_Critica_2])
        PcDato = np.array([Presion_Critica_1, Presion_Critica_2])
        wDato = np.array([Factor_Acentrico_1, Factor_Acentrico_2])
    
    
    button.on_click(cargarDatos)
    #display(button)
    
    page1 = widgets.VBox(children=[Componentes_1, Componentes_2, button], padding=4)
    
    
    #VBox([VBox([Button(description='Press'), Dropdown(options=['a', 'b']), Button(description='Button')]), 
    #      VBox([Button(), Checkbox(), IntText()]), 
    #      VBox([Button(), IntSlider(), Button()])], background_color='#EEE')
    
    
    ecuacionEstado = widgets.Dropdown(description='Fluid :', padding=4, options=['SRK', 'PR', 'RKPR'])
    modeloSolido = widgets.Dropdown(description='Solid :', padding=4, options=['Model I', 'Model II', 'Model III'])
    
    button = widgets.Button(description="Upload Models")
    
    def cargarModelos(b):
        clear_output()
        global eq    
        eq = ecuacionEstado.value
        
        print("Component 1: ", Componentes_1.value)
        print("Component 2: ", Componentes_2.value)
        
        print("Fluid Model : ", ecuacionEstado.value)
        print("Solid Model : ", modeloSolido.value)
        
        
    
    
    button.on_click(cargarModelos)
    
    page2 = widgets.Box(children=[ecuacionEstado, modeloSolido, button], padding=4)
    
    
    Temp_ini = widgets.Text(description='Initial', padding=4, value="0.0")
    Temp_fin = widgets.Text(description='Final', padding=4, value="0.0")
    
    Pres_ini = widgets.Text(description='Initial', padding=4, value="0.0")
    Pres_fin = widgets.Text(description='Final', padding=4, value="0.0")
    
    n1 = widgets.Text(description='Mole light component', padding=4, value="0.0")
    n2 = widgets.Text(description='Mole heavy component', padding=4, value="0.0")
    
    #button = widgets.Button(description="Cargar Condiciones")
    
    titulo = widgets.HTML(value="<C><H1> System Conditions <H1>")
    tempe_info = widgets.HTML(value="<C><H3> Temperature <H3>")
    press_info = widgets.HTML(value="<C><H3> Pressure <H3>")
    fluid_info = widgets.HTML(value="<C><H3> Mole fracction in the fluid <H3>")
    
    button = widgets.Button(description="Upload Conditions")
    
    def cargarParametros(b):
        clear_output()
        
        global initial_temperature, initial_pressure, nif
        
        initial_temperature = float(Temp_ini.value)
        initial_pressure = float(Pres_ini.value)
        nif = np.array([float(n1.value), float(n2.value)])
            
        print("Component 1: ", Componentes_1.value)
        print("Component 2: ", Componentes_2.value)
        
        print("Fluid Model : ", ecuacionEstado.value)
        print("solid Model : ", modeloSolido.value)
        
        print("Initial_temperature = ", initial_temperature, type(initial_temperature))
        print("Final_temperature = ", Temp_fin.value)
        
        print("Initial_pressure =", initial_pressure, type(initial_pressure))
        print("Final_pressure =", Pres_fin.value)
        
        print("Mole fraccion light component n1 =", n1.value)
        print("Mole fraccion heavy component n2 =", n2.value)
        
        print("Mole fracction in the fluid ", nif) 
        
        print(initial_temperature, type(initial_temperature))
        
    
    button.on_click(cargarParametros)
    
    page3 = widgets.Box(children=[titulo, tempe_info, Temp_ini, Temp_fin, press_info, Pres_ini, Pres_fin, fluid_info, n1, n2, button], padding=4)
    
    
    
    
    
    
    button = widgets.Button(description="Solid-Fluid")
    #display(button)
    
    
    nnCC_1 = 1
    nnCC_2 = 2
    
    def calcularSolidoFluido(b):
        clear_output()
        #Tcal = fsolve(equilibrioSF,guess,args=(Pe, nnCC_1, nnCC_2), xtol=1e-4)
        
        #initial_temperature = [346.5] # T [=] K
        #initial_pressure = 137.9 # [=] bar
        
        Tcal = fsolve(equilibrioSF,initial_temperature,args=(initial_pressure, nif, 1, 2, Avsl), xtol=1e-4)
        print("Temperature ESF = ", Tcal, "K")
        
    
    button.on_click(calcularSolidoFluido)
    #display(button)
    
    
    page4 = widgets.Box(children=[button], padding=4)
    
    button = widgets.Button(description="Diagram Solid-Fluid")
    
    def DiagramaSolidoFluido(b):
        clear_output()
        #Tcal = fsolve(equilibrioSF,guess,args=(Pe, 1, 2), xtol=1e-4)
        #Tcal = fsolve(equilibrioSF,guess,args=(Pe, nnCC_1, nnCC_2), xtol=1e-4)
        initial_temperature =346.5 # T [=] K
        initial_pressure = 136.9 # [=] bar
        
        #346.5 136.9
        # n1, n2 = 1, 2 por defecto para el equilibrio sólido-fluido
        Tcal = fsolve(equilibrioSF,initial_temperature,args=(initial_pressure, nif, 1, 2, Avsl), xtol=1e-4)
    
        print(Tcal, "K")
        
        pyplot.scatter(Tres,pres, color = 'red', label = 'PR')
        pyplot.scatter(temp,pres, label = 'Data')
        pyplot.title('Temperature Equilibrium Solid Liquid')
        pyplot.legend(loc="upper left") 
        pyplot.xlabel('Temperature [=] K')  
        pyplot.ylabel('Pressure [=] bar')
        
    
    button.on_click(DiagramaSolidoFluido)
    
    
    page5 = widgets.Box(children=[button], padding=4)
    
    
    
    
    DatosTemperatura_Exp = np.array([323.65, 326.04, 326.43])
    #DatosTemperatura_Exp = np.array([323.65, 326.04, 326.43, 328.12])
    
    DatosPresionp_Exp = np.array([1.0, 101.0, 136.9])
    #DatosPresionp_Exp = np.array([1.0, 101.0, 136.9, 183.8])
    
    
    posicion = np.arange(len(DatosPresionp_Exp))
    TemperaturasModelo = np.ones((len(DatosPresionp_Exp)))
    TemperaturasModelo
    
    Avsl = -0.32595074
    
    
    button = widgets.Button(description="Regression of Parameters")
    
    def regresionParametros(b):
        clear_output()
        
        def minimizarVSL(Avsl):
            for T, P, i in zip(DatosTemperatura_Exp, DatosPresionp_Exp, posicion):
                print ("Initial Temperature = ", T, "K", "Pressure = ", P, "bar", "Experimental Data = ", i+1)
                initial_temperature = T # T [=] K
                initial_pressure = P # [=] bar
                # tol=
                TemperaturasModelo[i] = fsolve(equilibrioSF,initial_temperature,args=(initial_pressure, nif, 1, 2, Avsl), xtol=1e-4)
                
            funcionObjetivo = np.sum((DatosTemperatura_Exp - TemperaturasModelo) ** 2)
            print("modelTemperature = ", TemperaturasModelo)
            print("Objective Function = ", funcionObjetivo)
            
            return funcionObjetivo
        
        opt = sp.optimize.minimize(minimizarVSL, Avsl, method='L-BFGS-B')
        
        print("optimal parameter", opt)
    
    
    button.on_click(regresionParametros)
    
    
    page6 = widgets.Box(children=[button], padding=4)
    
    tabs = widgets.Tab(children=[page1, page2, page3, page4, page5, page6])
    #display(tabs)
    
    tabs.set_title(0, 'Components')
    tabs.set_title(1, 'Models')
    tabs.set_title(2, 'Conditions')
    tabs.set_title(3, 'Results')
    tabs.set_title(4, 'Experimental Data')
    tabs.set_title(5, 'Regression of Parameters')
    
    #--------------------- flash Isothermal------------------------------
    
    Componentes_f1 = widgets.SelectMultiple(
        description="Component 1",
        options=list(Etiquetas))
    
    Componentes_f2 = widgets.SelectMultiple(
        description="Component 2",
        options=list(Etiquetas))
    
    Componentes_f3 = widgets.SelectMultiple(
        description="Component 3",
        options=list(Etiquetas))
    
    Componentes_f4 = widgets.SelectMultiple(
        description="Component 4",
        options=list(Etiquetas))
    
    button = widgets.Button(description="Upload Data")
    
    def cargarDatos(b):
        clear_output()
        print("Component 1: ", Componentes_f1.value)
        Nombre = Componentes_f1.value
        Propiedades = data2.loc[Nombre]
        Factor_Acentrico_1 = Propiedades[0]
        Temperatura_Critica_1 = Propiedades[1]
        Presion_Critica_1 = Propiedades[2]
        Z_Critico_1 = Propiedades[3]
    
        #print(Propiedades)
        print ("Acentric Factor = ", Factor_Acentrico_1)
        print ("Critical Temperature = ", Temperatura_Critica_1, "K")
        print ("Critical Pressure = ", Presion_Critica_1, "bar")
        print ("Z_Critical = ", Z_Critico_1, "\n")
        
        
        print("Component 2: ", Componentes_f2.value)
        Nombre = Componentes_f2.value
        Propiedades = data2.loc[Nombre]
        Factor_Acentrico_2 = Propiedades[0]
        Temperatura_Critica_2 = Propiedades[1]
        Presion_Critica_2 = Propiedades[2]
        Z_Critico_2 = Propiedades[3]
    
        #print(Propiedades)
        print ("Acentric Factor  = ", Factor_Acentrico_2)
        print ("Critical Temperature = ", Temperatura_Critica_2, "K")
        print ("Critical Pressure = ", Presion_Critica_2, "bar")
        print ("Z_Critical = ", Z_Critico_2, "\n")
        
        print("Component 3: ", Componentes_f3.value)
        Nombre = Componentes_f3.value
        Propiedades = data2.loc[Nombre]
        Factor_Acentrico_3 = Propiedades[0]
        Temperatura_Critica_3 = Propiedades[1]
        Presion_Critica_3 = Propiedades[2]
        Z_Critico_3 = Propiedades[3]
    
        #print(Propiedades)
        print ("Acentric Factor = ", Factor_Acentrico_3)
        print ("Critical Temperature = ", Temperatura_Critica_3, "K")
        print ("Critical Pressure = ", Presion_Critica_3, "bar")
        print ("Z_Critical = ", Z_Critico_3, "\n")
        
        
        print("Component 4: ", Componentes_f4.value)
        Nombre = Componentes_f4.value
        Propiedades = data2.loc[Nombre]
        Factor_Acentrico_4 = Propiedades[0]
        Temperatura_Critica_4 = Propiedades[1]
        Presion_Critica_4 = Propiedades[2]
        Z_Critico_4 = Propiedades[3]
    
        #print(Propiedades)
        print ("Acentric Factor  = ", Factor_Acentrico_4)
        print ("Critical Temperature = ", Temperatura_Critica_4, "K")
        print ("Critical Pressure = ", Presion_Critica_4, "bar")
        print ("Z_Critical = ", Z_Critico_4, "\n")
        
        global TcDato_f, PcDato_f, wDato_f
        
        TcDato_f = np.array([Temperatura_Critica_1, Temperatura_Critica_2, Temperatura_Critica_3, Temperatura_Critica_4])
        PcDato_f = np.array([Presion_Critica_1, Presion_Critica_2, Presion_Critica_3, Presion_Critica_4])
        wDato_f = np.array([Factor_Acentrico_1, Factor_Acentrico_2, Factor_Acentrico_3, Factor_Acentrico_4])
    
    
    button.on_click(cargarDatos)
    #display(button)
    
    page_f1 = widgets.VBox(children=[Componentes_f1, Componentes_f2, Componentes_f3, Componentes_f4, button], padding=4)
    
    
    #------------------ page_f2
    ecuacionEstado_f = widgets.Dropdown(description='Fluid :', padding=4, options=['SRK', 'PR', 'RKPR'])
    
    button = widgets.Button(description="Upload Models")
    
    def cargarModelos(b):
        clear_output()
        global eq    
        eq = ecuacionEstado.value
        
        print("Component 1: ", Componentes_f1.value)
        print("Component 2: ", Componentes_f2.value)
        print("Component 3: ", Componentes_f3.value)
        print("Component 4: ", Componentes_f4.value)    
        print("Fluid Model : ", ecuacionEstado_f.value)
        
        
    button.on_click(cargarModelos)
    
    page_f2 = widgets.Box(children=[ecuacionEstado_f, button], padding=4)
    
    #------------------ page_f2
    
    #------------------ page_f3
    Temp_ini_f = widgets.Text(description='Initial', padding=4, value="0.0")
    
    Pres_ini_f = widgets.Text(description='Initial', padding=4, value="0.0")
    
    n1_f = widgets.Text(description='Component 1', padding=4, value="0.0")
    n2_f = widgets.Text(description='Component 2', padding=4, value="0.0")
    n3_f = widgets.Text(description='Component 3', padding=4, value="0.0")
    n4_f = widgets.Text(description='Component 4', padding=4, value="0.0")
    
    
    titulo = widgets.HTML(value="<C><H1> System Conditions <H1>")
    tempe_info = widgets.HTML(value="<C><H3> Temperature <H3>")
    press_info = widgets.HTML(value="<C><H3> Pressure <H3>")
    fluid_info = widgets.HTML(value="<C><H3> Mole fracction in the fluid <H3>")
    
    button = widgets.Button(description="Upload Conditions")
    
    def cargarParametros(b):
        clear_output()
        
        global zi_F, temperature_f, pressure_f, nif
        
        temperature_f = float(Temp_ini_f.value)
        pressure_f = float(Pres_ini_f.value)
        zi_F = np.array([float(n1_f.value), float(n2_f.value), float(n3_f.value), float(n4_f.value)])
        nif = np.array([float(n1_f.value), float(n2_f.value), float(n3_f.value), float(n4_f.value)])
            
        print("Component 1: ", Componentes_f1.value)
        print("Component 2: ", Componentes_f2.value)
        print("Component 3: ", Componentes_f3.value)
        print("Component 4: ", Componentes_f4.value)
        
        print("Fluid Model : ", ecuacionEstado_f.value)
        
        print("Temperature_f = ", temperature_f, type(temperature_f))
        
        print("Pressure_f = ", pressure_f, type(pressure_f))
        
        print("Mole fraccion component 1 = ", n1_f.value)
        print("Mole fraccion component 2 = ", n2_f.value)
        print("Mole fraccion component 3 = ", n3_f.value)
        print("Mole fraccion component 4 = ", n4_f.value)
        
        print("Mole fracction in the fluid = ", zi_F, type(zi_F)) 
        
        print(temperature_f, type(temperature_f))
        
    
    button.on_click(cargarParametros)
    
    page_f3 = widgets.Box(children=[titulo, tempe_info, Temp_ini_f, press_info, Pres_ini_f, fluid_info, n1_f, n2_f, n3_f, n4_f, button], padding=4)
    #------------------ page_f3
    
    #------------------ page_f4
    
    button = widgets.Button(description="Flash Calculation")
    
    def calcularFlashPT(b):
        clear_output()
        #Tcal = fsolve(equilibrioSF,guess,args=(Pe, nnCC_1, nnCC_2), xtol=1e-4)
        
        #initial_temperature = [346.5] # T [=] K
        #initial_pressure = 137.9 # [=] bar
        fhid = Flash(zi_F, temperature_f, pressure_f, TcDato_f, PcDato_f, wDato_f)
        fhid.flash_ideal()
        
    
    button.on_click(calcularFlashPT)
    #display(button)
    
    
    page_f4 = widgets.Box(children=[button], padding=4)
    
    
    
    #------------------ page_f4
    flash = widgets.Tab(children=[page_f1, page_f2, page_f3, page_f4])
    
    flash.set_title(0, 'Components')
    flash.set_title(1, 'Models')
    flash.set_title(2, 'Conditions')
    flash.set_title(3, 'Results')
    #tabs.set_title(4, 'Experimental Data')
    #tabs.set_title(5, 'Regression of Parameters')
    
    #--------------------- GPEC ------------------------------
    
    name_GPEC = widgets.Text(description='File name', padding=4, value=" ")
    url = name_GPEC.value
    
    titulo = widgets.HTML(value="<C><H1> Data GPEC <H1>")
    
    button_1 = widgets.Button(description="UpData GPEC")
    
    def upGPEC(b):
        clear_output()
        
        DGPEC = DataGPEC(url)
        DGPEC.leerGPEC_1()
        print ("Upload {0}".format(url))
            
    
    button_1.on_click(upGPEC)
    
    button_2 = widgets.Button(description="Vapor pressure")
    
    def diagram_1(b):
        clear_output()
        DGPEC = DataGPEC(url)
        DGPEC.leerGPEC_1()
        DGPEC.presionVapor()
        
    button_2.on_click(diagram_1)
    
    
    button_3 = widgets.Button(description="Diagram Density-Pressure")
    
    def diagram_2(b):
        clear_output()
        DGPEC = DataGPEC(url)
        DGPEC.leerGPEC_1()
        DGPEC.densidadPresion()
        
    
    button_3.on_click(diagram_2)
    
    page_G1 = widgets.Box(children=[titulo, name_GPEC, button_1, button_2, button_3], padding=4)
    
    gpec = widgets.Tab(children=[page_G1])
    
    gpec.set_title(0, 'Upload Data')
    
    
    
    accord = widgets.Accordion(children=[tabs, flash, gpec], width=400)
    display(accord)
    
    
    
    
    accord.set_title(0, 'Pure Solid-Binary Fluid')
    accord.set_title(1, 'Isothermal Flash Calculation')
    accord.set_title(2, 'Data GPEC')
    #accord.set_title(3, 'Regression of Parameters Solid-Fluid')
    #accord.set_title(4, 'Pure Fluid')
    
    
    
    #346.5 136.9
    # Lectura Juan.xlsx

.. code:: python

    url = 'Lectura Juan.xlsx'

