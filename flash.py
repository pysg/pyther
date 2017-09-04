def etiquetar():

    rotulo_Liquido = "Composición de fase líquido"
    rotulo_Vapor = "Composición de fase vapor"
    rotulo_Separador = "-" * 11

    etiqueta_liquido = "{0}{1}{0}".format(rotulo_Separador, rotulo_Liquido)
    etiqueta_vapor = "{0}{1}{0}".format(rotulo_Vapor, rotulo_Vapor)

    return etiqueta_liquido, etiqueta_vapor


print(etiquetar()[0])


def function():
    print("Metano, Butano, Hexano")
    etiqueta_liquido = "Composición de fase líquida"
    etiqueta_vapor = "Composición de fase vapor"

    print(" -*{} {etiqueta_liquido}".format(etiqueta_liquido))
    print("xi = ", xy[0])
    print("Sxi = ", np.sum(xy[0]))
    print("-" * 20)
    print("yi = ", xy[1])
    print("Syi = ", np.sum(xy[1]))
    pass



class ClassName(object):
    """docstring for ClassName"""
    def __init__(self, arg):
        super(ClassName, self).__init__()
        self.arg = arg
        
    def function_0():
        pass
    def function_1():
        pass


    def wilson(self, Pc, Tc, w, T):
        # Ecuación wilson
        lnKi = np.log(Pc / self.P) + 5.373 * (1 + w) * (1 - Tc / self.T)
        self.Ki = np.exp(lnKi)
        return self.Ki

    def beta(self, zi):
        self.zi = zi
        self.Ki = self.wilson(Pc, Tc, w, T)
        Bmin = np.divide((self.Ki * self.zi - 1), (self.Ki - 1))
        #print (("Bmin_inter = ", Bmin))
        Bmax = np.divide((1 - self.zi), (1 - self.Ki))
        #print (("Bmax_inter = ", Bmax))
        self.Bini = (np.max(Bmin) + np.min(Bmax)) / 2
        return self.Bini

    def rice(self, zi, Ki, Bini):
        self.zi = zi
        self.Bini = Bini
        self.Ki = Ki
        self.fg = np.sum(self.zi * (self.Ki - 1) / (1 - self.Bini + self.Bini * self.Ki))
        self.dfg = - np.sum(self.zi * (self.Ki - 1) ** 2 / (1 - self.Bini + self.Bini * self.Ki) ** 2)
        #print g, dg
        return self.fg, self.dfg

    def flash_ideal(self):
        self.Bini = self.beta(zi)
        self.Ki = self.wilson(self.Pc, self.Tc, self.w, self.T)
        print("Ki_(P, T) = ", self.Ki)
        Eg = self.rice(zi, self.Ki, self.Bini)
        errorEq = abs(Eg[0])
        i, s = 0, 1

        while errorEq > ep:
            Eg = self.rice(zi, self.Ki, self.Bini)
            self.Bini = self.Bini - s * Eg[0] / Eg[1]
            errorEq = abs(Eg[0])
            i += 1
            if i >= 50:
                break

        xy = self.composicion_xy(zi, self.Ki, self.Bini)
        print ("Metano, Butano, Hexano")
        print ("-------------Composición de fase líquida------------------------")
        print ("xi = ", xy[0])
        print ("Sxi = ", np.sum(xy[0]))
        print ("-------------Composición de fase vapor------------------------")
        print ("yi = ", xy[1])
        print ("Syi = ", np.sum(xy[1]))

        return Eg[0], Eg[1], self.Bini

    def flash_PT(self):
        flashID = self.flash_ideal()
        print ("flash (P, T, zi)")
        print ("g, dg, B = ", flashID)
        print ("---------------------------------------------------------------")

        self.Bini = flashID[2]
        print ("Beta_r ini = ", self.Bini)
        moles = self.composicion_xy(zi, self.Ki, self.Bini)

        self.xi, self.yi = moles[0], moles[1]
        nil, niv = moles[2], moles[3]

        fi_F = self.fugac()

        #print "nil = ", nil, np.sum(nil)
        #print "Snil = ", np.sum(nil)
        #print "niv", niv, np.sum(niv)
        #print "Sniv = ", np.sum(niv)

        #nF = 2
        #CoeFugi = np.ones((nF, nC))

        #for i in range(nF):
        #    if i == 1:
        #        self.ni = nil
        #        self.nT = np.sum(self.ni)
        #    elif i == 2:
        #        self.ni = niv
        #        self.nT = np.sum(self.ni)

        #    Flug_i = self.fluido(self.P)
        #    CoeFugi[i, :] = Flug_i[1]

        #print CoeFugi
        #self.Ki = CoeFugi[0, :] / CoeFugi[1, :]

        self.Ki = fi_F[0] / fi_F[1]

        L = 1.0

        self.Ki = self.Ki * L

        Ki_1 = self.Ki
        print ("Ki_(P, T, ni) primera = ", self.Ki)

        print ("----------------------------------------------------------------")

        #self.Ki = np.array([1.729, 0.832, 0.640])

        #self.Ki = self.wilson(self.Pc, self.Tc, self.w, self.T)
        #print "Ki_(P, T) = ", self.Ki

        while 1:
            i, s = 0, 0.1

            while 1:
                Eg = self.rice(zi, self.Ki, self.Bini)
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
            print ("-----------Composición de fase líquida----------------------")
            print ("xi = ", moles[0])
            print ("Sxi = ", np.sum(moles[0]))
            print ("-----------Composición de fase vapor------------------------")
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

    def fugac(self):
        '''
        Esta función que se llama fugacidad, calcula los coeficientes de fugacidad
        con una ecuación de estado para mezclas multicomponente
        T = temperatura en Kelvin
        Y = fracción molar del componente i en la fase de vapor
        X = fracción molar del componente i en la fase líquida
        '''
        self.T
        self.P
        #----------------------------------------------------------------------
        self.Fw = 0.48 + (1.574 * self.w) - (0.176 * self.w ** 2)
        a = ((0.42748 * (self.R * self.Tc) ** 2) / self.Pc) * ((1 + self.Fw * (1 - (self.Tr ** 0.5))) ** 2)
        b = (0.08664 * self.R * self.Tc) / self.Pc
        #----------------------------------------------------------------------
        #print Fw, "Parametro a:", a, b

        #Yf = np.array([Y4,Y5,Y6])
        Yf = self.yi
        Xf = self.xi
        #print Yf, Xf
        #----------------------- Vapor -----------------------------------------
        amv = np.sum(Yf * a ** 0.5) ** 2
        aml = np.sum(Xf * a ** 0.5) ** 2
        #-------------------------------
        #print "amv = ", amv
        #print "aml = ", aml

        bmv = np.sum(Yf * b)
        bml = np.sum(Xf * b)
        #print "bmv = ", bmv
        #print "bml = ", bml

        Av = (amv * self.P) / ((self.R * self.T) ** 2)
        Bv = (bmv * self.P) / (self.R * self.T)
        #-------------------- Liquido -------------------
        Al = (aml * self.P) / ((self.R * self.T) ** 2)
        Bl = (bml * self.P) / (self.R * self.T)
        #print "Av", Av
        #print "Bv", Bv
        #print "Al", Al
        #print "Av", Bl

        Zfv = [1, -1, (Av - Bv - Bv ** 2), (- Av * Bv)]
        ZfvR = np.roots(Zfv)
        Zv = np.max(ZfvR)
        #print "Zv = ", Zv

        Zfl = [1, -1, (Al - Bl - Bl ** 2), (- Al * Bl)]
        ZflR = np.roots(Zfl)
        Zl = np.min(ZflR)
        #print "Zl = ", Zl
        #----------------------------------------------------------------------
        lnfiv = (b / bmv) * (Zv - 1) - np.log(Zv - Bv) + (Av / Bv)  * ((b / bmv) - (2 * ((a / amv) ** 0.5))) * np.log((Zv + Bv) / Zv)
        self.fiv = np.exp(lnfiv)
        print ("fiv = ", self.fiv)
        lnfil = (b / bml) * (Zl - 1) - np.log(Zl - Bl) + (Al / Bl)  * ((b / bml) - (2 * ((a / aml) ** 0.5))) * np.log((Zl + Bl) / Zl)
        self.fil = np.exp(lnfil)
        print ("fil = ", self.fil)

        return self.fil, self.fiv


    def composicion_xy(self, zi, Ki, Bini):
        self.zi = zi
        self.Ki = Ki
        self.Bini = Bini
        self.xi = zi / (1 - self.Bini + self.Bini * self.Ki)
        self.yi = (zi * self.Ki) / (1 - self.Bini + self.Bini * self.Ki)
        self.li = (zi * (1 - self.Bini)) / (1 - self.Bini + self.Bini * self.Ki)
        self.vi = (zi * self.Bini * self.Ki) / (1 - self.Bini + self.Bini * self.Ki)

        return self.xi, self.yi, self.li, self.vi













def fugacity(self):

        self.Fw = 0.48 + (1.574 * self.w) - (0.176 * self.w ** 2)
        a = ((0.42748 * (self.R * self.Tc) ** 2) / self.Pc) * ((1 + self.Fw * (1 - (self.Tr ** 0.5))) ** 2)
        b = (0.08664 * self.R * self.Tc) / self.Pc

        Yf = self.yi
        Xf = self.xi

        # vapor
        amv = np.sum(Yf * a ** 0.5) ** 2
        bmv = np.sum(Yf * b)

        Av = (amv * self.P) / ((self.R * self.T) ** 2)
        Bv = (bmv * self.P) / (self.R * self.T)

        Zfv = [1, -1, (Av - Bv - Bv ** 2), (- Av * Bv)]
        ZfvR = np.roots(Zfv)
        Zv = np.max(ZfvR)

        # líquido
        aml = np.sum(Xf * a ** 0.5) ** 2
        bml = np.sum(Xf * b)
        Al = (aml * self.P) / ((self.R * self.T) ** 2)
        Bl = (bml * self.P) / (self.R * self.T)

        Zfl = [1, -1, (Al - Bl - Bl ** 2), (- Al * Bl)]
        ZflR = np.roots(Zfl)
        Zl = np.min(ZflR)

        # coeficiente de fugacidad
        lnfiv = (b / bmv) * (Zv - 1) - np.log(Zv - Bv) + (Av / Bv)  * ((b / bmv) - (2 * ((a / amv) ** 0.5))) * np.log((Zv + Bv) / Zv)
        self.fiv = np.exp(lnfiv)
        print("fiv = ", self.fiv)
        lnfil = (b / bml) * (Zl - 1) - np.log(Zl - Bl) + (Al / Bl)  * ((b / bml) - (2 * ((a / aml) ** 0.5))) * np.log((Zl + Bl) / Zl)
        self.fil = np.exp(lnfil)
        print("fil = ", self.fil)

        return self.fil, self.fiv



import numpy as np
import pyther as pt


class Flash(object):
    """
    Flash_TP is a Class for to calculate the flash with a temperature T and
    pressure P for a specified composition zi. In this case, the algorithm
    used, is a combination between the ideal flash with Ki(T,P) without
    dependence of the composition and the flash with Ki(T,P,zi) with dependence
    of the composition from the two phases in equilibrium.
    """

    def __init__(self, arg):
        self.Tc = arg[0]
        self.Pc = arg[1]
        self.w = arg[2]
        self.T = arg[3]
        self.P = arg[4]
        self.zi = arg[5]
        self.R = pt.RGAS

    def Ki_wilson(self):
        """Equation of wilson for to calculate the Ki(T,P)"""
        variable_0 = 5.373 * (1 + self.w) * (1 - self.Tc / self.T)
        lnKi = np.log(self.Pc / self.P) + variable_0
        self.Ki = np.exp(lnKi)
        return self.Ki

    def beta_initial(self):
        self.Ki = self.Ki_wilson()
        self.Bmin = np.divide((self.Ki * self.zi - 1), (self.Ki - 1))
        # print (("Bmin_inter = ", Bmin))
        self.Bmax = np.divide((1 - self.zi), (1 - self.Ki))
        # print (("Bmax_inter = ", Bmax))
        self.Binit = (np.max(self.Bmin) + np.min(self.Bmax)) / 2
        return self.Binit

    def rachford_rice(self):

        denominador = (1 - self.Binit + self.Binit * self.Ki)
        numerador = self.zi * (self.Ki - 1)
        self.function_rachford_rice = np.sum(numerador / denominador)
        # Derivate of function_rachford_rice with respect to Beta
        variable_1 = self.zi * (self.Ki - 1) ** 2 / denominador ** 2
        self.d_functions_rachford_rice = - np.sum(variable_1)
        return self.function_rachford_rice, self.d_functions_rachford_rice

    def composition_xy(self):
        denominador = (1 - self.Binit + self.Binit * self.Ki)
        self.xi = self.zi / denominador
        self.yi = (self.zi * self.Ki) / denominador
        return self.xi, self.yi

    def beta_newton(self):
        iteration, step, tolerance = 0, 1, 1e-5
        while True:
            self.Binit = self.Binit - step * self.rachford_rice()[0] / self.rachford_rice()[1]
            iteration += 1
            if abs(self.rachford_rice()[0]) <= tolerance or (iteration >= 50):
                break
        print(iteration)

        return self.Binit

    def flash_ideal_method_1(self):
        self.Bini = self.beta_initial()
        self.Ki = self.Ki_wilson()
        print("Ki_(P, T) = ", self.Ki)
        Eg = self.rachford_rice()
        errorEq = abs(Eg[0])
        i, s = 0, 1

        while errorEq > 1e-5:
            Eg = self.rachford_rice()
            self.Bini = self.Bini - s * Eg[0] / Eg[1]
            errorEq = abs(Eg[0])
            i += 1
            if i >= 50:
                break

        xy = self.composition_xy()
        print("C1, Ci4, C4")
        print("-" * 12, "Composición de fase líquida", "-" * 12)
        print("xi = ", xy[0])
        print("Sxi = ", np.sum(xy[0]))
        print("-" * 12, "Composición de fase vapor", "-" * 12)
        print("yi = ", xy[1])
        print("Syi = ", np.sum(xy[1]))

        return Eg[0], Eg[1], self.Bini

    def flash_ideal(self):
        self.Binit = self.beta_initial()
        self.Ki = self.Ki_wilson()
        self.Binit = self.beta_newton()
        self.xy = self.composition_xy()
        return self.rachford_rice()[0], self.rachford_rice()[1], self.Binit, self.xy, self.Ki

    def fugacity(self):

        self.m = 0.48 + (1.574 * self.w) - (0.176 * self.w ** 2)
        self.Tr = self.T / self.Tc
        alpha = (1 + self.m * (1 - (self.Tr ** 0.5))) ** 2
        ac = 0.42748 * (self.R * self.Tc) ** 2 / self.Pc
        a = ac * alpha
        b = 0.08664 * self.R * self.Tc / self.Pc

        Yf = self.yi
        Xf = self.xi

        # vapor
        amv = np.sum(Yf * a ** 0.5) ** 2
        bmv = np.sum(Yf * b)

        Av = (amv * self.P) / ((self.R * self.T) ** 2)
        Bv = (bmv * self.P) / (self.R * self.T)

        Zv = np.max(np.roots([1, -1, (Av - Bv - Bv ** 2), (- Av * Bv)]))

        # líquido
        aml = np.sum(Xf * a ** 0.5) ** 2
        bml = np.sum(Xf * b)
        Al = (aml * self.P) / ((self.R * self.T) ** 2)
        Bl = (bml * self.P) / (self.R * self.T)

        Zl = np.min(np.roots([1, -1, (Al - Bl - Bl ** 2), (- Al * Bl)]))

        # coeficiente de fugacidad
        aav = (a / amv)
        bbv = (b / bmv)

        aal = (a / aml)
        bbl = (b / bml)

        factor_1 = (bbv - (2 * (aav ** 0.5))) * np.log((Zv + Bv) / Zv)
        ln_phi_v = bbv * (Zv - 1) - np.log(Zv - Bv) + (Av / Bv) * factor_1
        self.phi_v = np.exp(ln_phi_v)

        factor_2 = (bbl - (2 * (aal ** 0.5))) * np.log((Zl + Bl) / Zl)
        ln_phi_l = bbl * (Zl - 1) - np.log(Zl - Bl) + (Al / Bl) * factor_2
        self.phi_l = np.exp(ln_phi_l)

        print("phi_v = ", self.phi_v)
        print("phi_l = ", self.phi_l)

        return self.phi_l, self.phi_v

    def flash_real(self):
        self.Binit = self.flash_ideal()[2]
        self.Ki = self.flash_ideal()[4]
        Ki_1 = self.Ki
        print("Ki_(P, T) inicial = ", self.Ki)
        tolerance = 1e-5

        while True:
            self.xi, self.yi = self.composition_xy()
            self.Ki = self.fugacity()[0] / self.fugacity()[1]
            self.Binit = self.beta_newton()

            Ki_2 = self.Ki
            dKi = abs(Ki_1 - Ki_2)
            Ki_1 = Ki_2

            if np.sum(dKi) <= tolerance:
                break

        return self.xi, self.yi, self.Binit


def main():

    print("-" * 79)

    # component = 'METHANE'
    # component = "ETHANE"
    # component = "3-METHYLHEPTANE"
    # component = "n-PENTACOSANE"

    # component = "ISOBUTANE"

    # component = ["METHANE", "n-TETRACOSANE", "n-PENTACOSANE", "ETHANE", "ISOBUTANE", "PROPANE", "3-METHYLHEPTANE"]

    # component = "METHANE"
    # component =  "ETHANE"
    # component = "n-TETRACOSANE"
    #component = "ISOBUTANE"
    #component = "n-BUTANE"
    component =  "PROPANE"
    #component = "HEXANE"

    properties_data = pt.Data_parse()
    properties_component = properties_data.selec_component(component)
    pt.print_properties_component(component, properties_component)

    dinputs = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
                        properties_component[1]['Omega'], properties_component[1]['Vc']])

    print(dinputs)
    print('-' * 79)


main()


c1 = np.array([1.90564000e+02, 4.53890000e+01, 1.15000000e-02, 9.86000000e-02])
c2 = np.array([3.05320000e+02, 4.80830000e+01, 9.95000000e-02, 1.45500000e-01])
c3 = np.array([3.69830000e+02, 4.19240000e+01, 1.52300000e-01, 2.00000000e-01])
c4 = np.array([4.25120000e+02, 3.74640000e+01, 2.00200000e-01, 2.55000000e-01])

ci4 = np.array([4.08140000e+02, 3.60030000e+01, 1.80800000e-01, 2.62700000e-01])
c24 = np.array([804.0, 9.672, 1.071, 1.41])


c3 = np.array([369.8, 42.49, 0.152])
ci4 = np.array([408.1, 36.48, 0.177])
c4 = np.array([425.2, 37.97, 0.193])



# Tc = np.array([c2[0], c3[0], ci4[0]])
# Pc = np.array([c2[1], c3[1], ci4[1]])
# w = np.array([c2[2], c3[2], ci4[2]])

# Tc = np.array([c1[0], c2[0], c3[0]])
# Pc = np.array([c1[1], c2[1], c3[1]])
# w = np.array([c1[2], c2[2], c3[2]])

#Tc = np.array([c2[0], c3[0], c4[0]])
#Pc = np.array([c2[1], c3[1], c4[1]])
#w = np.array([c2[2], c3[2], c4[2]])

Tc = np.array([c3[0], ci4[0], c4[0]])
Pc = np.array([c3[1], ci4[1], c4[1]])
w = np.array([c3[2], ci4[2], c4[2]])


T = 320.0
P = 8.0
zi = np.array([0.23, 0.67, 0.10])

argumentos = [Tc, Pc, w, T, P, zi]
flash = Flash(argumentos)
b = flash.flash_ideal()
d = flash.flash_real()

beta = b[2]

q = b[3]
xi = q[0]
yi = q[1]

print("*" * 70)
print("C3 -i-C4 n-C4")
print("---------- Composition of liquid phase ----------")
print("xi = {0} and Sxi ={1}". format(xi, np.sum(xi)))
print("---------- Composition dof vapor phase ----------")
print("yi = {0} and Syi ={1}". format(yi, np.sum(yi)))
print("-" * 70)
print("Beta(P, T) =", beta)
print("*" * 70)

print(d)


print(constans["Tc"])

Tc = np.float64(constans["Tc"])
Pc = np.float64(constans["Pc"])
w = np.float64(constans["Omega"])

print("Tc_1 = ", Tc)
print("Pc1 = ", Pc)
print("w_1 = ", w)

Tc_1 = [369.83, 408.14, 425.12]
Pc_1 = [41.924, 36.003, 37.464]
w_1 = [0.1523, 0.1808, 0.2002]

# -----------------------------------------------
c3 = np.array([369.8, 42.49, 0.152])
ci4 = np.array([408.1, 36.48, 0.177])
c4 = np.array([425.2, 37.97, 0.193])

#Tc = np.array([c3[0], ci4[0], c4[0]])
#Pc = np.array([c3[1], ci4[1], c4[1]])
#w = np.array([c3[2], ci4[2], c4[2]])
# -----------------------------------------------

#Tc = np.array([369.83, 408.14, 425.12])
#Pc = np.array([41.924, 36.003, 37.464])
#w = np.array([0.1523, 0.1808, 0.2002])

#(array([ 0.21574563,  0.68079758,  0.10345679]), array([ 0.37253378,  0.56203169,  0.06543454]), 0.090913934405840488)
#(array([ 0.21574563,  0.68079758,  0.10345679]), array([ 0.37253378,  0.56203169,  0.06543454]), 0.090913934405849051)

