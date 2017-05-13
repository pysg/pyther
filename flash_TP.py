import numpy as np
import pandas as pd

from .pure_data import Data_parse
from .constans import RGAS


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
        self.R = RGAS

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
        return self.Binit

    def isothermal_ideal(self):
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

        # Vapor
        amv = np.sum(self.yi * a ** 0.5) ** 2
        bmv = np.sum(self.yi * b)
        Av = (amv * self.P) / ((self.R * self.T) ** 2)
        Bv = (bmv * self.P) / (self.R * self.T)
        Zv = np.max(np.roots([1, -1, (Av - Bv - Bv ** 2), (- Av * Bv)]))

        aav = (a / amv)
        bbv = (b / bmv)

        # Liquid
        aml = np.sum(self.xi * a ** 0.5) ** 2
        bml = np.sum(self.xi * b)
        Al = (aml * self.P) / ((self.R * self.T) ** 2)
        Bl = (bml * self.P) / (self.R * self.T)
        Zl = np.min(np.roots([1, -1, (Al - Bl - Bl ** 2), (- Al * Bl)]))

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
        self.Binit = self.isothermal_ideal()[2]
        self.Ki = self.isothermal_ideal()[4]
        Ki_1 = self.Ki
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

    def table_flash(self):
        self.datos = np.array([self.zi, self.isothermal()[0], self.isothermal()[1]]).T
        self.etiqueta_colums = ["zi", "xi", "yi"]

        self.resultados_flash = pd.DataFrame(self.datos, self.components, self.etiqueta_colums)

        return self.resultados_flash


    def constans_to_flash(self):

        properties_data = Data_parse()
        properties_component = properties_data.selec_component(self.components)

        constans = properties_component[1].loc[:, ["Omega", "Tc", "Pc"]]

        return constans


def main():

    components = ["PROPANE", "ISOBUTANE", "n-BUTANE"]

    T = 320.0
    P = 8.0
    zi = np.array([0.23, 0.67, 0.10])

    flash_1 = Flash(components, T, P, zi)

    b = flash_1.isothermal_ideal()
    d = flash_1.isothermal()

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

    print(flash_1.table_flash())


if __name__ == '__main__':
    main()
