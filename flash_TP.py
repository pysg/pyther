import numpy as np
import pandas as pd
import pyther as pt


class Flash(object):
    """
    Flash_TP is a Class for to calculate the flash with a temperature T and
    pressure P for a specified composition zi. In this case, the algorithm
    used, is a combination between the ideal flash with Ki(T,P) without
    dependence of the composition and the flash with Ki(T,P,zi) with dependence
    of the composition from the two phases in equilibrium.
    """

    def __init__(self, *args):

        self.components = args[0]
        self.constans = self.constans_to_flash()
        self.Tc = np.float64(self.constans["Tc"])
        self.Pc = np.float64(self.constans["Pc"])
        self.w = np.float64(self.constans["Omega"])
        self.T = args[1]
        self.P = args[2]
        self.zi = args[3]
        self.R = pt.RGAS

    def Ki_wilson(self):
        """Equation of wilson for to calculate the Ki(T,P)"""
        variable_0 = 5.373 * (1 + self.w) * (1 - self.Tc / self.T)
        lnKi = np.log(self.Pc / self.P) + variable_0
        self.Ki = np.exp(lnKi)
        return self.Ki

    def Ki_wilson_add(self):

        self.Ki = self.Ki_wilson()

        i = 0

        while True:

            g0 = np.sum(self.zi * self.Ki) - 1.0
            g1 = 1.0 - np.sum(self.zi / self.Ki)

            if (g0 < 0):
                self.Ki = 1.1 * self.Ki
            elif (g1 > 0):
                self.Ki = 0.9 * self.Ki

            i += 1
            print(i)
            print("go = {0} and g1 = {1}".format(g0, g1))

            if (g0 > 0 or g1 < 0):
                break

            if i >= 190:
                break

        return self.Ki

    def beta_initial(self):
        self.Ki = self.Ki_wilson_add()

        self.Bmin_values = np.divide((self.Ki * self.zi - 1), (self.Ki - 1))
        self.Bmin_values = np.array(list(filter(lambda x: 0 < x < 1, self.Bmin_values)))
        print(self.Bmin_values)
        if len(self.Bmin_values) == 0:
            self.Bmax_values = np.array([0])
        self.Bmin = np.max(self.Bmin_values)

        self.Bmax_values = np.divide((1 - self.zi), (1 - self.Ki))
        self.Bmax_values = np.array(list(filter(lambda x: 0 < x < 1, self.Bmax_values)))
        print(self.Bmax_values)
        if len(self.Bmax_values) == 0:
            self.Bmax_values = np.array([1])
        self.Bmax = np.min(self.Bmax_values)

        self.Binit = (self.Bmin + self.Bmax) / 2
        return self.Binit

    def rachford_rice(self):
        denominador = 1 + self.Binit * (self.Ki - 1)
        numerador = self.zi * (self.Ki - 1)
        self.function_rachford_rice = np.sum(numerador / denominador)
        # Derivate of function_rachford_rice with respect to Beta
        variable_1 = self.zi * (self.Ki - 1) ** 2 / denominador ** 2
        self.d_functions_rachford_rice = - np.sum(variable_1)
        return self.function_rachford_rice, self.d_functions_rachford_rice

    def composition_xy(self):
        denominador = 1 + self.Binit * (self.Ki - 1)
        self.yi = self.zi * self.Ki / denominador
        # self.xi = self.zi / denominador
        self.xi = self.yi / self.Ki
        return self.xi, self.yi

    def beta_newton(self):
        iteration, step, tolerance = 0, 0.1, 1e-6

        #while self.Beta < self.Bmin or self.Beta > self.Bmax:
        #    step = step / 2
        #    self.Beta = self.Beta - step

        while True:
            self.Binit = self.Binit - step * self.rachford_rice()[0] / self.rachford_rice()[1]

            self.advance = self.rachford_rice()[0] / self.rachford_rice()[1]

            i = 0
            while (self.Binit < self.Bmin) or (self.Binit > self.Bmax):
                self.Binit = self.Binit - self.advance / 2
                print(self.Bmin, self.Binit, self.Bmax)
                i += 1
                if i >= 50:
                    break

            iteration += 1
            if abs(self.rachford_rice()[0]) <= tolerance or (iteration >= 2000):
                break
        return self.Binit

    def isothermal_ideal(self):
        self.Binit = self.beta_initial()
        #self.Ki = self.Ki_wilson()
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

        # Vapor
        amv = np.sum(self.yi * a ** 0.5) ** 2
        bmv = np.sum(self.yi * b)
        Av = (amv * self.P) / ((self.R * self.T) ** 2)
        Bv = (bmv * self.P) / (self.R * self.T)
        print("a = ", a)
        print("yi =", self.yi)
        print("amv = ", amv)

        Zv = np.max(np.roots([1, -1, (Av - Bv - Bv ** 2), (- Av * Bv)]))
        #Zv = np.max(np.roots([1, -(1 - Bv), (Av - 3*Bv**2 - 2*Bv), (- Av * Bv - Bv**2-Bv**3)]))

        aav = (a / amv)
        bbv = (b / bmv)

        # Liquid
        aml = np.sum(self.xi * a ** 0.5) ** 2
        bml = np.sum(self.xi * b)

        print("xi =", self.xi)
        Al = (aml * self.P) / ((self.R * self.T) ** 2)
        Bl = (bml * self.P) / (self.R * self.T)
        Zl = np.min(np.roots([1, -1, (Al - Bl - Bl ** 2), (- Al * Bl)]))
        #Zl = np.max(np.roots([1, -(1 - Bl), (Al - 3*Bl**2 - 2*Bl), (- Al * Bl - Bl**2-Bl**3)]))

        aal = (a / aml)
        bbl = (b / bml)

        # Fugacity Coefficient

        factor_1 = (bbv - (2 * (aav ** 0.5))) * np.log((Zv + Bv) / Zv)
        ln_phi_v = bbv * (Zv - 1) - np.log(Zv - Bv) + (Av / Bv) * factor_1
        self.phi_v = np.exp(ln_phi_v)

        factor_2 = (bbl - (2 * (aal ** 0.5))) * np.log((Zl + Bl) / Zl)
        ln_phi_l = bbl * (Zl - 1) - np.log(Zl - Bl) + (Al / Bl) * factor_2
        self.phi_l = np.exp(ln_phi_l)

        return self.phi_l, self.phi_v

    def isothermal(self):
        self.Binit, self.Ki = self.isothermal_ideal()[2], self.isothermal_ideal()[4]
        Ki_1 = self.Ki
        tolerance = 1e-5

        while True:
            self.xi, self.yi = self.composition_xy()
            print("xi_iso = ", self.xi)
            print("yi_iso = ", self.yi)

            self.Ki = self.fugacity()[0] / self.fugacity()[1]
            print("Ki_iso = ", self.Ki)
            self.Binit = self.beta_newton()

            print("Beta_iso = ", self.Binit)
            print("avance_iso = ", self.advance)

            Ki_2 = self.Ki
            dKi = abs(Ki_1 - Ki_2)
            Ki_1 = Ki_2

            if np.max(dKi) <= tolerance:
                break

        return self.xi, self.yi, self.Binit

    def constans_to_flash(self):

        properties_data = pt.Data_parse()
        properties_component = properties_data.selec_component(self.components)

        constans = properties_component[1].loc[:, ["Omega", "Tc", "Pc"]]

        return constans


def main():

    # components = ["PROPANE", "ISOBUTANE", "n-BUTANE"]
    components = ["METHANE", "PROPANE", "n-PENTANE", "n-DECANE", "n-HEXADECANE"]

    T = 400.0
    #T = 383.15
    P = 30.0
    zi = np.array([0.822, 0.088, 0.050, 0.020, 0.020])

    #components = ["METHANE", "ETHANE", "PROPANE", "ISOBUTANE", "n-BUTANE"]
    #T = 320.0
    #P = 9.0
    #zi = np.array([0.05, 0.1, 0.23, 0.52, 0.10])

    flash_1 = Flash(components, T, P, zi)

    fk = flash_1.Ki_wilson()
    fkk = flash_1.Ki_wilson_add()

    print(fk, fkk)

    bk = flash_1.beta_initial()

    print(flash_1.Bmin, bk, flash_1.Bmax)
    print(flash_1.Bmin_values, flash_1.Bmax_values)

    b = flash_1.isothermal_ideal()

    #print(b)


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

    d = flash_1.isothermal()
    print("Beta_real = ", d[2])

    datos = np.array([zi, xi, yi]).T
    etiqueta_col = ["zi", "xi", "yi"]

    resultados_flash = pd.DataFrame(datos, components, etiqueta_col)

    print(resultados_flash)

    frames = [resultados_flash, resultados_flash]

    result = pd.concat(frames)

    print(result)


if __name__ == '__main__':
    main()
