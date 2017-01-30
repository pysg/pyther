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
    


