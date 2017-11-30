**************************************
14. Cálculo del flash Isotermico (T, P)
**************************************

Se presenta una implementación del calculo del flash isotermico bifasico
utilizando la ecuación de estado *Peng-Robinsong (PR)* [2] junto con las
reglas de mezclado de *Van Der Waalls* [2].

El cálculo del flash isotermico bifasico es un cálculo básico en la
introducción de los procesos de separación porque es el esqeuma
tecnologíco de separación más simple, en el que ingresa una corriente de
fluido a un "tanque" calentado por un flujo de calor en el que se
obtiene una corriente de salida por cada fase presente en el sistema. En
el caso bifasico, una corriente de líquido y otra de vapor, tal como se
muestra en la figura 1.

|pilares\_open\_science| Figura 1. Esquema del cálculo del flash
isotermico

.. image:: _static/flash_diagrama.jpg
  :width: 1200

14.1 Modelo flash líquido-vapor
----------------------------------------

El modelo del flash isotermico bifasico, corresponde al balance de
materia global y por componente en el tanque separador que se muestra en
la figura (1), junto con la condición de equilibrio de fases
líquido-vapor.

Coeficiente de distribución :math:`K_i`

.. math::  Ki = \frac {yi} {xi} 

Aproximación de wilson para el coeficiente de distribución :math:`K_i`

.. math::  lnK_i = ln \left(\frac {Pc_i} {P}\right ) + 5.373(1 + w_i)(1 - \frac {Tc_i} {T}) 

*Rachford-Rice* :math:`g(\beta)`

.. math::  g(\beta) = \sum \limits_{i=1}^{C} (y_i - x_i) 

.. math::  g(\beta) = \sum \limits_{i=1}^{C} \frac {K_i - 1} {1 - \beta + \beta K_i} 

Derivada de la función *Rachford-Rice* :math:`g(\beta)`

.. math::  \frac {dg} {d \beta} = \sum \limits_{i=1}^{C} z_i \frac {(K_i - 1)^2} {(1 - \beta + \beta K_i)^2} < 0 

Valores límites de la función *Rachford-Rice* :math:`g(\beta)`

.. math::  g(0) = \sum \limits_{i=1}^{C} (z_i K_i - 1) > 0 

.. math::  g(1) = \sum \limits_{i=1}^{C} (1 - \frac {z_i} {K_i}) < 0 

Ecuaciones para calcular las fracciones molares de cada fase

.. math::  y_i \frac{K_i z_i} {1 - \beta + \beta K_i} 

.. math::  x_i = \frac{z_i} {1 - \beta + \beta K_i} 

Relaciones que determinan los valores mínimos y máximos para
:math:`\beta`

.. math::  1 - \beta + \beta K_i >= K_i z_i 

.. math::  \beta \geq \frac {K-i z_i - 1} {K_i - 1} 

.. math::  1 - \beta + \beta K_i >= z_i 

.. math::  \beta \leq \frac {z_i - 1} {1 - K_i} 

Valores extremos de la fracción de vapor en el sistema :math:`\beta`

.. math::  \beta_{min} = 0 

.. math::  \beta_{max} = 1 

14.2 Algoritmo
------------

-  Especificar la Presión :math:`P`, Temperatura :math:`T` y número de
   moles :math:`N` de cada componente del sistema
-  Calcular el coeficiente de distribución :math:`K_i^{wilson}` a partir
   de la relación de Wilson
-  Calcular el valor de :math:`\beta_{min}`
-  Calcular el valor de :math:`\beta_{max}`
-  Calcular el promedio de beta, usando Beta minimo y Beta máximo
-  Resolver la ecuación de *Rachford-Rice* :math:`g(\beta)`, para
   calcular :math:`\beta` con una tolerancia de :math:`1x10^{-6}`
-  Calcular las fracciones molares del líquido :math:`x_i` y del vapor
   :math:`y_i`
-  Calcular los coeficientes de fugacidad :math:`\hat{\phi_i}` para las
   fracciones molares del líquido :math:`x_i` y del vapor :math:`y_i`
-  Calcular el coeficiente de distribución :math:`K_i` a partir de los
   coeficientes de fugacidad del componente i :math:`\hat{\phi_i}`
-  Volver a resolver la ecuación de *Rachford-Rice* :math:`g(\beta)`,
   para calcular :math:`\beta` con una tolerancia de :math:`1x10^{-6}`
-  Verificar la convergencia del sistema con una tolerancia de
   :math:`1x10^{-6}` para :math:`\Delta K_i =  \left | K_{i}^{j+1} - K_{i}^{j} \right|`,
   siendo está situación la convergencia del procedimiento.

14.2.1 Implementación
------------------

En la implementación del cálculo del flash isotermico, se tiene 3 partes
importantes:

-  Cálculo de los coeficientes de distribución por medio de la ecuación
   de Wilson
-  Cálculo de los valores mínimos y máximos para la fracción
   :math:`\beta`
-  Cálculo del *step* para calcular la fracción :math:`\beta`

Ecuación de Wilson
~~~~~~~~~~~~~~~~~~

.. code:: ipython3

        def Ki_wilson(self):
            """Equation of wilson for to calculate the Ki(T,P)"""
            variable_0 = 5.373 * (1 + self.w) * (1 - self.Tc / self.T)
            lnKi = np.log(self.Pc / self.P) + variable_0
            self.Ki = np.exp(lnKi)
            return self.Ki

Cálculo de los valores mínimos y máximos para la fracción :math:`\beta`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

        def beta_initial(self):
            self.Ki = self.Ki_wilson()
            self.Bmin = (self.Ki * self.zi - 1) / (self.Ki - 1)
            self.Bmax = (1 - self.zi) / (1 - self.Ki)
            self.Binit = (np.max(self.Bmin) + np.min(self.Bmax)) / 2
            return self.Binit

Cálculo del *step* para calcular la fracción :math:`\beta`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def beta_newton(self):
        iteration, step, tolerance = 0, 1, 1e-5
        while True:
            self.Binit = self.Binit - step * self.rachford_rice()[0] / self.rachford_rice()[1]
            iteration += 1
            while self.Binit < self.Bmin or self.Binit > self.Bmax:
                step = step / 2
            if abs(self.rachford_rice()[0]) <= tolerance or (iteration >= 50):
                break
        return self.Binit

14.3. Resultados
------------------

A continuación se muestran los resultados numéricos del calculo del
flash isotermico bifasico para una mezcla de los componentes
(C3-Ci4-C4), que corresponde al cálculo del flash isotermico propuesto
por (Elliott & Lira, 2012) el ejemplo 10.7 de su libro Introductory
Chemical engineering thermodynamics. En la tabla 1, se presentan las
especificaciones de la presión P, temperatura T y flujo F junto con las
fracciones molares del líquido, del vapor y la fracción de fase
resultanten usando como modelo termodinámico la ecuación de estado
*Peng-robinson (PR)* y las reglas de mezclado de *Van Der Waalls*.

En la tabla 1., se presenta el resultado del cálculo del flash
isotermico utilizando solo el :math:`K_i^{wilson}`

Tabla.1 flash isotermico :math:`K_i(T, P)` Mezcla ideal

+---------------+-----------------+-----------------+
| Presión Bar   | Temperatura K   | Flujo F mol/h   |
+===============+=================+=================+
| 8             | 320             | 1               |
+---------------+-----------------+-----------------+

| Componente \| :math:`z_i` \| líquido :math:`x_i` \| Vapor :math:`y_i`
\|
| :---------:\| ---------- ------------ ------------\|
|     C3 \| 0.23 \|0.18357118 \|0.37209837 \|
|     Ci4 \| 0.67 \|0.70479988 \|0.56349276 \|
|     C4 \| 0.10 \|0.11162895 \|0.06440887 \|

+--------------------------+------------------------------------------------+-----------------------+
| función g                | derivada función :math:`\frac{dg}{d \beta }`   | :math:`\beta`         |
+==========================+================================================+=======================+
| 6.1017797856749434e-07   | -0.20663315922997191                           | 0.24627123315157093   |
+--------------------------+------------------------------------------------+-----------------------+

mientras que en la tabla 2, se muestra el resultado del cálculo del
flash isotermico utilizando el resultado de :math:`K_i^{wilson}` como
valor inicial para el procedimiento del cálculo del flash isotermico
incluyento el cálculo de los coeficientes de fugacidad
:math:`\hat{\phi_i}` con la ecuación de estado PR.

Tabla.2 Flash isotermico :math:`K_i(T, P, x_i, y_i)` **(PR)**

| Componente \| :math:`z_i` \| líquido :math:`x_i` \| Vapor :math:`y_i`
\|
| :---------:\| ---------- ------------ ------------\|
|     C3 \| 0.23 \|0.20070242 \|0.35071046 \|
|     Ci4 \| 0.67 \|0.69183981 \|0.5800167 \|
|     C4 \| 0.10 \|0.10745949 \|0.06926579 \|

+---------------------------+------------------------------------------------+-----------------+
| función g                 | derivada función :math:`\frac{dg}{d \beta }`   | :math:`\beta`   |
+===========================+================================================+=================+
| -9.7482523918959729e-06   | -0.13108663002971882                           | 0.19530673657   |
+---------------------------+------------------------------------------------+-----------------+

De esta forma, se observa que el algoritmo empleando la ecuación de
estado **Peng-Robinson (PR)** converge en a una solución *cercana* de la
solución que utiliza la aproximación de wilson para el coeficiente de
distribución **Ki**, mostrando ser efieciente para casos simples como el
presente en este capítulo.

14.3.1 Efecto de la temperatura y presión sobre :math:`\beta`
----------------------------------------------------------

Para el mismo sistema que se presentó en las tabla 1 y 2, en la figura 2
se muestra la solución del cálculo del flash isotermico para un rango de
presión y temperatura en el cual la fracción vaporizada :math:`\beta`
varia entre 0 y 1. En este caso, al aumentar la presión :math:`\beta`
disminuye mientras que el efecto de la temperatura es el contrario.

|pilares\_open\_science| Figura 2. Efecto de la temperatura y presión
sobre :math:`\beta`

.. |pilares\_open\_science| image:: img/index.png

14.4 Conclusiones
---------------

-  Se implemento el cálculo del flash isotermico bifasico utilizando la
   ecuación de estado *Peng-Robinsong (PR)* tomando las recomendaciones
   planteadas en el curso de termodinámica de fluidos para mejorar la
   convergencia del cálculo.

-  Se encontró que se utilizan en promedio 3 iteraciones para calcular
   el valor :math:`\beta` en cada paso que se mantienen constantes los
   valores :math:`K_i`.

14.5 Referencias
--------------

1. Curso de especialización en Termodinámica de fluidos. Ph.D Martín
   Cismondí. Marzo-Junio (2017)

2. Introductory Chemical engineering thermodynamics. J. Richard Elliott
   , Carl T. Lira. Prentice Hall (2012)
