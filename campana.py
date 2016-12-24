import numpy as np
#import solution_matrix as sm


#c    Starting from critical point    c    c
#    WRITE(nout,15) TCmod(icomp),Pc(icomp),Dc(icomp),Dc(icomp)
#    T = 0.9999 * TCmod(icomp)
#    Vc = 1/DC(icomp)
#    Vv = 1.03*Vc
#    Zc = Pc(icomp)*Vc/RGAS/T
#    aaa = min(0.89 + (Zc - 0.2) / 2, 0.95)
#    Vl = aaa*Vc
#    NS = 3
#    delXS = 0.10

temperature = 0.9999 * Tc
critical_volume = 1 / DC
vapor_volume = 1.03 * critical_volume
Zc = Pc * critical_volume / RGAS / temperature
liquid_volume = min(0.89 + (Zc - 0.2) / 2, 0.95) * critical_volume
NS = 3
delXS = 0.10

#------------------------------------------------------
XVAR = np.log([temperature, liquid_volume, vapor_volume])
DFDS=0.0
DFDS[3]=1.0
RJAC=0.0
RJAC[3, NS]=1.0

NITER=0
T=exp(XVAR[1])
Vl=exp(XVAR[2])
Vv=exp(XVAR[3])

FMAXOLD=8.0
FMAX=7.0
DMAXOLD=8.0
DMAX=7.0
F[3]=0.0
delX=0.0
NV=0
#------------------------------------------------------








T = 123
TOLF = 1e-3
NS = 2
DFDS = np.zeros(3)
DFDS[2] = 1
F = np.ones(3)
F[2] = 0
RJAC = np.ones([3, 3])
RJAC[2, NS] = 0

delx_s = 0
NV = 0

Pl = 0.12
Pv = 3
FUGx = 2.3
FUGy = 3.7

condition_1 = Pl < 0.5 * Pv and (Pv - Pl) > 1e-8
condition_2 = Pl > 1.5 * Pv and (Pl - Pv) > 1e-8
print(Pl)
print(condition_1)
print(condition_2)

if condition_1:
	print('condition_1 = {0}'.format(condition_1))

if condition_2:
	print('condition_2 = {0}'.format(condition_2))

if condition_1 or condition_2:
	print(condition_1, condition_2)
	NV += 1
 
if Pl < 0:
	F[0] = TOLF
else:
	F[0] = np.log(Pl/Pv)

F[1] = FUGx - FUGy


#RJAC[0, 0] = T * (DPDTx / Pl - DPDTy / Pv)
#RJAC[0, 1] = Vl * DPDVx / Pl
#RJAC[0, 2] = -Vv * DPDVy / Pv

#RJAC[1, 0] = T * (FUGTx -FUGTy)
#RJAC[1, 1] = Vl * FUGVx
#RJAC[1, 2] = -Vv * FUGVy 

RJAC[0, 0] = 23
RJAC[0, 1] = 45
RJAC[0, 2] = -23

RJAC[1, 0] = 67
RJAC[1, 1] = 0.34
RJAC[1, 2] = -34.4 






print(F)
print(RJAC)
print(RJAC[2,2])

b = F
a = RJAC

print(a)

x = np.linalg.solve(a, b)
print(x)

print('sm.A = {0}'.format(sm.A))


class Component(object):		

	def function_Ar_cal(self):
		
		self.bv = self.B / self.V
		self.f = np.log((self.V + self.s1 * self.B) / (self.V + self.s2 * self.B)) / self.B / (self.s1 - self.s2)
		self.g = self.R * np.log(1 - self.B / self.V)

		self.AUX = self.R * self.T / (self.V - self.B)
		self.fB = -(self.f + self.V * self.fv) / self.B
		self.FFB = self.nT * AUX - self.D * self.fB
		self.Di = 2 * self.nT * self.ac * self.alfa
		self.Bi = self.bc

		self.Ar = -self.nT * self.g * self.T - self.D * self.f
		'''Primera derivada de F con respecto al volumen Ecu. (68)'''
		self.gv = self.R * self.B / (self.V * (self.V - self.B))
		self.fv = - 1 / ((self.V + self.s1 * self.B) * (self.V + self.s2 * self.B))
		self.ArV = -self.nT * self.gv * self.T - self.D * self.fv
		''' Segunda derivada de F con respecto al volumen Ecu. (74) '''
		self.gv2 = self.R * (1 / self.V ** 2 - 1 / (self.V - self.B) ** 2)
		self.fv2 = (- 1 / (self.V + self.s1 * self.B) ** 2 + 1 / (self.V + self.s2 * self.B) ** 2) / self.B / (self.s1 - self.s2)
		self.ArV2 = - self.nT * self.gv2 * self.T - self.D * self.fv2
		''' pressure '''
		self.Pcal = self.nT * self.R * self.T / self.V - self.ArV
		self.dPdV = -self.ArV2 - self.R * self.T * self.nT / self.V ** 2

		if self.eq != "RKPR":
			self.Arn = -self.g * self.T + self.FFB * self.Bi - self.f * self.Di
		else:
			self.Arn = -self.g * self.T + self.FFB * self.Bi - self.f * self.Di - self.D * self.fD1 * self.dD1i

		ArT = -nT * g - dDdT * f
		ArTV = -nT * gv - dDdT * fV
		ArTn = -g + (nT * AUX/T - dDdT * fB) * dBi - f * dDiT - dDdT * fD1 * dD1i
		ArVn = - gv * T + FFBV * dBi - fv * dDi - D * fVD1 * dD1i

		return self.g, self.f, self.Ar


component = Component()
component.function_Ar_cal()



class Component(object):		

	def function_Ar_cal(self):
		
		self.bv = self.B / self.V
		self.f = np.log((self.V + self.s1 * self.B) / (self.V + self.s2 * self.B)) / self.B / (self.s1 - self.s2)
		self.g = self.R * np.log(1 - self.B / self.V)

		self.AUX = self.R * self.T / (self.V - self.B)
		self.fB = -(self.f + self.V * self.fv) / self.B
		self.FFB = self.nT * AUX - self.D * self.fB
		self.Di = 2 * self.nT * self.ac * self.alfa
		self.Bi = self.bc

		self.Ar = -self.nT * self.g * self.T - self.D * self.f
		'''Primera derivada de F con respecto al volumen Ecu. (68)'''
		self.gv = self.R * self.B / (self.V * (self.V - self.B))
		self.fv = - 1 / ((self.V + self.s1 * self.B) * (self.V + self.s2 * self.B))
		self.ArV = -self.nT * self.gv * self.T - self.D * self.fv
		''' Segunda derivada de F con respecto al volumen Ecu. (74) '''
		self.gv2 = self.R * (1 / self.V ** 2 - 1 / (self.V - self.B) ** 2)
		self.fv2 = (- 1 / (self.V + self.s1 * self.B) ** 2 + 1 / (self.V + self.s2 * self.B) ** 2) / self.B / (self.s1 - self.s2)
		self.ArV2 = - self.nT * self.gv2 * self.T - self.D * self.fv2
		''' pressure '''
		self.Pcal = self.nT * self.R * self.T / self.V - self.ArV
		self.dPdV = -self.ArV2 - self.R * self.T * self.nT / self.V ** 2

		if self.eq != "RKPR":
			self.Arn = -self.g * self.T + self.FFB * self.Bi - self.f * self.Di
		else:
			self.Arn = -self.g * self.T + self.FFB * self.Bi - self.f * self.Di - self.D * self.fD1 * self.dD1i

		ArT = -nT * g - dDdT * f
		ArTV = -nT * gv - dDdT * fV
		ArTn = -g + (nT * AUX/T - dDdT * fB) * dBi - f * dDiT - dDdT * fD1 * dD1i
		ArVn = - gv * T + FFBV * dBi - fv * dDi - D * fVD1 * dD1i

		return self.g, self.f, self.Ar


component = Component()
component.function_Ar_cal()



# this is a example for the calculate for one composi
component.function_Ar_cal()




# there is a form to convert a componet object in a molec

