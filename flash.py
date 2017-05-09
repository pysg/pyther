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
        print ("Ki_(P, T) = ", self.Ki)
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
