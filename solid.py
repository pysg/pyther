





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
