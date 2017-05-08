import numpy as np


class Flash_TP(object):
    """
    Flash_TP is a Class for to calculate the flash with a temperature T and
    pressure P for a specified composition zi. In this case, the algorithm
    used, is a combination between the ideal flash with Ki(T,P) without
    dependence of the composition and the flash with Ki(T,P,zi) with dependence
    of the composition from the two phases in equilibrium.
    """

    def __init__(self, arg):
        # super(ClassName, self).__init__()
        self.Tc = arg[0]
        self.Pc = arg[1]
        self.w = arg[2]
        self.T = arg[3]
        self.zi = arg[4]

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
        # print g, dg
        return self.function_rachford_rice, self.d_functions_rachford_rice

    def flash_ideal(self):
        self.Binit = self.beta_initial()
        self.Ki = self.Ki_wilson()
        print("Ki_(P,T) = ", self.Ki)
        Eg = self.rachford_rice()
        errorEq = abs(Eg[0])
        i, s = 0, 1

        while errorEq > ep:
            Eg = self.rachford_rice(zi, self.Ki, self.Binit)
            self.Binit = self.Binit - s * Eg[0] / Eg[1]
            errorEq = abs(Eg[0])
            i += 1
            if i >= 50:
                break

        xy = self.composicion_xy(zi, self.Ki, self.Binit)


        return Eg[0], Eg[1], self.Binit

    def composition_xy(self, zi, Ki, Binit):
        self.zi = zi
        self.Ki = Ki
        self.Binit = Binit
        denominador = (1 - self.Binit + self.Binit * self.Ki)
        self.xi = zi / denominador
        self.yi = (zi * self.Ki) / denominador
        self.li = (zi * (1 - self.Binit)) / denominador
        self.vi = (zi * self.Binit * self.Ki) / denominador

        return self.xi, self.yi, self.li, self.vi

    def flash_PT(self):
        flashID = self.flash_ideal()
        print("flash (P, T, zi)")
        print("g, dg, B = ", flashID)
        print("-" * 20)

        self.Binit = flashID[2]
        print("beta_initial_r ini = ", self.Binit)
        moles = self.composicion_xy(zi, self.Ki, self.Binit)

        self.xi, self.yi = moles[0], moles[1]
        nil, niv = moles[2], moles[3]

        fi_F = self.fugac()
        self.Ki = fi_F[0] / fi_F[1]
        L = 1.0
        self.Ki = self.Ki * L

        Ki_1 = self.Ki
        print("Ki_(P, T, ni) primera = ", self.Ki)

        print("-" * 20)

        while 1:
            i, s = 0, 0.1

            while 1:
                Eg = self.rachford_rice(zi, self.Ki, self.Binit)
                print(Eg)
                self.Binit = self.Binit - s * Eg[0] / Eg[1]
                print(self.Binit)
                errorEq = abs(Eg[0])
                i += 1
                # print i

                #if self. Binit < 0 or self.Binit > 1:
                    #break
                #    self.Binit = 0.5
                if i >= 50:
                    pass
                    # break
                if errorEq < 1e-5:
                    break

            print("Resultado Real = ", Eg)
            print(" beta_initial r = ", self.Binit)

            moles = self.composicion_xy(zi, self.Ki, self.Binit)
            self.xi, self.yi = moles[0], moles[1]

            # xy = self.composicion_xy(zi, self.Ki, self.Binit)

            print("C1 -i-C4 n-C4")
            print("----------Composición de fase líquida----------")
            print("xi = ", moles[0])
            print("Sxi = ", np.sum(moles[0]))
            print("----------Composición de fase vapor----------")
            print("yi = ", moles[1])
            print("Syi = ", np.sum(moles[1]))

            fi_F = self.fugac()

            self.Ki = fi_F[0] / fi_F[1]
            Ki_2 = self.Ki
            dKi = abs(Ki_1 - Ki_2)
            Ki_1 = Ki_2
            print("Ki_(P, T, ni) = ", self.Ki)

            fun_Ki = np.sum(dKi)
            print("fun_Ki = ", fun_Ki)

            if fun_Ki < 1e-5:
                break

        return flashID


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










