import numpy as np
import pyther as pt


class Flash_TP(object):
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
        self.li = (self.zi * (1 - self.Binit)) / denominador
        self.vi = (self.zi * self.Binit * self.Ki) / denominador

        return self.xi, self.yi, self.li, self.vi

    def beta_newton(self):
        iteration, step, tolerance = 0, 1, 1e-5
        while True:
            self.Binit = self.Binit - step * self.rachford_rice()[0] / self.rachford_rice()[1]
            iteration += 1
            if abs(self.rachford_rice()[0]) <= tolerance or (iteration >= 50):
                break
        print(iteration)

        return self.Binit

    def flash_ideal(self):
        self.Binit = self.beta_initial()
        self.Ki = self.Ki_wilson()
        self.Binit = self.beta_newton()
        self.xy = self.composition_xy()
        return self.rachford_rice()[0], self.rachford_rice()[1], self.Binit, self.xy, self.Ki

    def flash_PT(self):

        self.Binit = self.flash_ideal()[2]
        self.Ki = self.flash_ideal()[4]
        Ki_1 = self.Ki
        print("Ki_(P, T) inicial = ", self.Ki)
        tolerance = 1e-5

        while True:
            self.xi, self.yi, nil, niv = self.composicion_xy()
            self.Ki = self.fugac()[0] / self.fugac()[1]
            self.Binit = self.beta_newton()

            Ki_2 = self.Ki
            dKi = abs(Ki_1 - Ki_2)
            Ki_1 = Ki_2

            if np.sum(dKi) <= tolerance:
                break

        return self.xi, self.yi, self.Binit


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


def main():

    print("-" * 79)

    #component = 'METHANE'
    #component = "ETHANE"
    #component = "3-METHYLHEPTANE"
    #component = "n-PENTACOSANE"
    
    #component = "ISOBUTANE"

    #component = ["METHANE", "n-TETRACOSANE", "n-PENTACOSANE", "ETHANE", "ISOBUTANE", "PROPANE", "3-METHYLHEPTANE"]

    #component = "METHANE"
    #component =  "ETHANE"
    component = "n-TETRACOSANE"

    
    properties_data = pt.Data_parse()
    #properties_component = properties_data.selec_component(dppr_file, component)
    properties_component = properties_data.selec_component(component)

    pt.print_properties_component(component, properties_component)

    dinputs = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
                        properties_component[1]['Omega'], properties_component[1]['Vc']])

    #print(dinputs)
    print('-' * 79)


main()


c1 = np.array([1.90564000e+02, 4.53890000e+01, 1.15000000e-02, 9.86000000e-02])
c2 = np.array([3.05320000e+02, 4.80830000e+01, 9.95000000e-02, 1.45500000e-01])
c24 = np.array([804.0, 9.672, 1.071, 1.41])


Tc = np.array([c1[0], c2[0], c24[0]])
Pc = np.array([c1[1], c2[1], c24[1]])
w = np.array([c1[2], c2[2], c24[2]])
T = 300.0
P = 2.0
zi = np.array([0.4, 0.4, 0.4])

argumentos = [Tc, Pc, w, T, P, zi]

flash = Flash_TP(argumentos)

b = flash.flash_ideal()

print(b[0:3])

q = b[3]
x = q[0]
y = q[1]
print(x)
print(y)

# (-1.3325942660458168e-06, -5.1253037025537775, 0.6577960750582591)
# (-1.3325942660458168e-06, -5.1253037025537775, 0.65779633506124491)


print("C1 -i-C4 n-C4")
print("----------Composición de fase líquida----------")
# print("xi = {0} and Sxi ={1}". format(xi, np.sum(xi)))
print("----------Composición de fase vapor----------")
# print("yi = {0} and Syi ={1}". format(yi, np.sum(yi)))
# print("Ki_(P, T, ni) final =", Ki)
