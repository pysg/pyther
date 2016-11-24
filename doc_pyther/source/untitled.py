
class Componentes():
    """
    Las variables aux_ se utilizan para presentar de forma más clara y acotada
    las expresiones necesarias en los calculos. Estas, se numeran de acuerdo al orden de
    aparición dentro de una clase.
    
    """
    
    def __init__(self, propierties_FQ, sys_conditions):
    	self.propierties_FQ = propierties_FQ
    	self.sys_conditions = sys_conditions


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
        
        if self.sys_conditions[3] == "SRK":
            # Soave-Redlich-Kwong (SRK)
            self.component = self.cal_SRK_model()
        elif self.sys_conditions[3] == "PR":
            # Peng-Robinson (PR)
            self.component = self.cal_PR_model()
        elif self.sys_conditions[3] == "RKPR":
            # (RKPR)
            #self.component = self.cal_RKPR_model()
            print ("No actualizada, intentalo de nuevo !!! ")
        else:
            print ("Che boludo... Modelo no valido, intentalo de nuevo !!! ")

        return self.component
            
            