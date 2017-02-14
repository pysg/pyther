6. Thermodynamics correlations for pure components
==================================================

En esta sección se muestra la class *Thermodynamic\_correlations()* la
cual permite realizar el cálculo de propiedades termodinámicas de
sustancias puras como una función de la temperatura. En este caso se
pueden tener 6 situaciones para cada una de las 13 propiedades
termofísicas soportadas:

1. Especificar una sustancia pura sin especificar una temperatura. En
   este caso por defecto la propiedad termodinámica se calcula entre el
   intervalo mínimo y máximo de validez experimental para cada
   correlación.

2. Especificar una sustancia pura y especificar una temperatura.

3. Especificar una sustancia pura y especificar varias temperaturas.

4. Especificar varias sustancias puras sin especificar una temperatura.

5. Especificar varias sustancias puras y especificar una temperatura.

6. Especificar varias sustancias puras y especificar varias temperaturas

la clase *Thermodynamics\_correlations* es usada para calcular 13
propiedades termodinámicas de sustancias puras en función de la
temperatura y se sigue la siguente convención para presentar identificar
las propiedades termodinámicas

property thermodynamics = name property, units, correlation, equation

The thermodynamic correlations are:

-**Solid\_Density** = "Solid Density", "[kmol/m^3]",
"A+B*T+C*\ T\ :sup:`2+D\ *T^3+E*\ T`\ 4", 0

-**Liquid\_Density** = "Liquid Density", "[kmol/m^3]",
"A/B:sup:`(1+(1-T/C)`\ D)", 1

-**Vapour\_Pressure** = "Vapour Pressure", "[Bar]",
"exp(A+B/T+C*ln(T)+D*\ T^E) \* 1e-5", 2

-**Heat\_of\_Vaporization** = "Heat of Vaporization", "[J/kmol]",
"A\*(1-Tr):sup:`(B+C*Tr+D*\ Tr`\ 2)", 3

-**Solid\_Heat\_Capacity** = "Solid Heat Capacity", "[J/(kmol\*K)]",
"A+B*T+C*\ T\ :sup:`2+D\ *T^3+E*\ T`\ 4", 4

-**Liquid\_Heat\_Capacity** = "Liquid Heat Capacity", "[J/(kmol\*K)]",
"A:sup:`2/(1-Tr)+B-2\ *A*\ C\ *(1-Tr)-A*\ D\ *(1-Tr):sup:`2-C`\ 2*\ (1-Tr)`\ 3/3-C\ *D*\ (1-Tr):sup:`4/2-D`\ 2\*(1-Tr)^5/5",
5

-**Ideal\_Gas\_Heat\_Capacity** = "Ideal Gas Heat Capacity"
"[J/(kmol\*K)]", "A+B*(C/T/sinh(C/T))^2+D*\ (E/T/cosh(E/T))^2", 6

-**Second\_Virial\_Coefficient** = "Second Virial Coefficient",
"[m^3/kmol]", "A+B/T+C/T:sup:`3+D/T`\ 8+E/T^9", 7

-**Liquid\_Viscosity** = "Liquid Viscosity", "[kg/(m\*s)]",
"exp(A+B/T+C*ln(T)+D*\ T^E)", 8

-**Vapour\_Viscosity** = "Vapour Viscosity", "[kg/(m\*s)]",
"A\*T:sup:`B/(1+C/T+D/T`\ 2)", 9

-**Liquid\_Thermal\_Conductivity** = "Liquid Thermal Conductivity",
"[J/(m*s*\ K)]", "A+B*T+C*\ T\ :sup:`2+D\ *T^3+E*\ T`\ 4", 10

-**Vapour\_Thermal\_Conductivity** = "Vapour Thermal Conductivity",
"[J/(m*s*\ K)]", "A\*T:sup:`B/(1+C/T+D/T`\ 2)", 11

-**Surface\_Tension** = "Surface Tension", "[kg/s^2]",
"A\*(1-Tr):sup:`(B+C*Tr+D*\ Tr`\ 2)", 12

Para empezar se importan las librerías que se van a utilizar, que en
este caso son numpy, pandas, pyther y especificar que las figuras
generadas se muesten dentro del jupyter notebook



.. code:: python

    import numpy as np
    import pandas as pd
    import pyther as pt
    import matplotlib.pyplot as plt
    %matplotlib inline


1. Especificar una sustancia pura sin especificar una temperatura.
==================================================================

Luego se carga el archivo que contine las constantes de las
correlaciones de las propiedades termodinamicas, que se llama en este
caso *"PureFull\_mod\_properties.xls"* y se asigna a la variable
*dppr\_file*.

Creamos un objeto llamado **thermodynamic\_correlations** y se pasan
como parametros las variables **component** y
**property\_thermodynamics** que en el ejemplo se especifica para el
componente METHANE la Vapour\_Pressure

.. code:: python

    dppr_file = "PureFull_mod_properties.xls"
    
    thermodynamic_correlations = pt.Thermodynamic_correlations(dppr_file)
    
    component = ['METHANE']
    property_thermodynamics = "Vapour_Pressure"
    
    Vapour_Pressure = thermodynamic_correlations.property_cal(component, property_thermodynamics)
    print("Vapour Pressure = {0}". format(Vapour_Pressure))



.. parsed-literal::

    ----------------------------------------------------------------------
    Pure substance without temperature especific: ['METHANE']
    ----------------------------------------------------------------------
    Vapour Pressure = [  0.11687017   0.13272851   0.15029231   0.1696935    0.19106965
       0.21456383   0.24032459   0.26850587   0.29926689   0.33277204
       0.36919081   0.40869762   0.45147173   0.49769708   0.54756216
       0.60125987   0.65898737   0.72094595   0.78734085   0.85838113
       0.93427952   1.01525227   1.101519     1.19330257   1.29082892
       1.39432695   1.50402838   1.62016765   1.74298174   1.87271013
       2.00959463   2.1538793    2.30581038   2.46563616   2.63360692
       2.80997486   2.99499402   3.18892023   3.39201109   3.60452587
       3.82672552   4.05887261   4.30123136   4.55406758   4.8176487
       5.09224376   5.37812343   5.67556002   5.98482753   6.30620166
       6.63995987   6.98638141   7.34574742   7.71834093   8.10444699
       8.5043527    8.91834734   9.34672242   9.78977179  10.24779173
      10.7210811   11.20994139  11.71467689  12.23559478  12.7730053
      13.32722183  13.89856107  14.48734319  15.09389194  15.71853484
      16.36160334  17.02343294  17.70436342  18.40473898  19.1249084
      19.86522527  20.62604814  21.40774072  22.21067207  23.03521683
      23.88175537  24.75067404  25.64236538  26.55722832  27.4956684
      28.45809802  29.44493665  30.45661106  31.49355559  32.55621234
      33.64503148  34.76047146  35.90299928  37.07309076  38.2712308
      39.49791367  40.75364324  42.03893333  43.35430794  44.7003016 ]


para realizar un gráfico simple de la propiedad termodinámica se utiliza
el método **graphical(temperature, property\_thermodynamics,
label\_property\_thermodynamics, units)**.

En donde se pasan como argumentos la temperatura a la cual se claculó la
propiedad termodinamica, los valores calculados de la propiedad
termodinamica, el label de la propiedad termodinámica y las unidades
correspondientes de temperatura y la propiedad termodinámica en cada
caso.

.. code:: python

    temperature_vapour = thermodynamic_correlations.temperature
    units = thermodynamic_correlations.units
    print(units)
    
    thermodynamic_correlations.graphical(temperature_vapour, Vapour_Pressure, property_thermodynamics, units)


.. parsed-literal::

    ('K', '[Pa]')



.. image:: output_9_1.png


2. Especificar una sustancia pura y una temperatura.
====================================================

Siguiendo con la sustacia pura *METHANE* se tiene el segundo caso en el
cual ademas de especificiar el componente se especifica también solo un
valor de temperatura, tal como se muestra en la variable *temperature*.

Dado que cada correlación de propiedad termodinámica tiene un rango
mínimo y máximo de temperatura en la cual es valida, al especificar un
valor de temperatura se hace una verificación para determinar si la
temperatura ingresada se encuentra entre el intervalo aceptado para cada
componente y cada propiedad termodinámica. En caso contrario la
temperatura se clasifica como invalida y no se obtiene valor para la
propiedad termodinámica seleccionada.

.. code:: python

    component = ['METHANE']
    property_thermodynamics = "Vapour_Pressure"
    temperature = [180.4]
    
    Vapour_Pressure = thermodynamic_correlations.property_cal(component, property_thermodynamics, temperature)
    print("Vapour Pressure = {0} {1}". format(Vapour_Pressure, units[1]))



.. parsed-literal::

    ----------------------------------------------------------------------
    Pure substance with a temperature especific: ['METHANE']
    ----------------------------------------------------------------------
    Temperature_enter = [180.4]
    Temperature_invalid = []
    Temperature_valid = [180.4]
    ----------------------------------------------------------------------
    Vapour Pressure = [ 33.32655377] [Pa]


3. Especificar una sustancia pura y especificar varias temperaturas.
====================================================================

Ahora se tiene la situación de contar con un solo componente "METHANE"
sin embargo, esta vez se especifica varios valores para la temperatura
en las cuales se quiere determinar el correspondiente valor de una
proiedad termodinámica, que como en los casos anteriores es la
*Vapour\_Pressure*.

.. code:: python

    component = ['METHANE']
    property_thermodynamics = "Vapour_Pressure"
    temperature = [180.4, 181.4, 185.3, 210, 85]
    
    Vapour_Pressure = thermodynamic_correlations.property_cal(component, "Vapour_Pressure", temperature)
    print("Vapour Pressure = {0} {1}". format(Vapour_Pressure, units[1]))


.. parsed-literal::

    ----------------------------------------------------------------------
    Pure substance with a temperature especific: ['METHANE']
    ----------------------------------------------------------------------
    Temperature_enter = [180.4, 181.4, 185.3, '210 K is a temperature not valid', '85 K is a temperature not valid']
    Temperature_invalid = ['210 K is a temperature not valid', '85 K is a temperature not valid']
    Temperature_valid = [180.4, 181.4, 185.3]
    ----------------------------------------------------------------------
    Vapour Pressure = [ 33.32655377  34.43422601  39.01608023] [Pa]


Se debe notar que al ingresar una serie de valores de temperatura, en
este caso 5 valores, se obtienen solo 3 valores de la propiedad
termodinámica. Esto se debe a que para este caso 2 valores de
temperatura no se encuentran en el valor mínimo y máximo en donde es
valida la correlación termodinámica. Por tanto, esto se avisa por medio
del mensaje: *Temperature\_invalid = ['210 K is a temperature not
valid', '85 K is a temperature not valid']*

4. Especificar varias sustancias puras sin especificar una temperatura.
=======================================================================

Otra de las posibilidades que se puede tener es la opción de especificar
varios componentes para una misma propiedad termodinámica sin que se
especifique una o más valores de temperatura. En esta opción se pueden
ingresar multiples componentes sin un limite, siempre y cuando estén en
la base de datos con la que se trabaja o en dado caso sean agregados a
la base de datos nuevas correlaciones para sustancias puras *Ver sección
base de datos*. Para este ejemplo se utiliza una *list components* con 3
sustancias puras por cuestiones de visibilidad de las gráficas de
*Vapour\_Pressure*.

.. code:: python

    components = ["METHANE", "n-TETRACOSANE", "ISOBUTANE"]
    property_thermodynamics = "Vapour_Pressure"
    
    Vapour_Pressure = thermodynamic_correlations.property_cal(components, property_thermodynamics)
    temperature_vapour = thermodynamic_correlations.temperature

por medio del método *multi\_graphical(components, temperature,
property\_thermodynamics)* al cual se pasan los parámetros
correspondiente a las sustancias puras, la temperatura a la cual se
realiza el calculo de la propiedad termodinámica y los valores de la
propiedad termodinámica de cada sustancia pura, para obtener la
siguiente figura.

.. code:: python

    
    thermodynamic_correlations.multi_graphical(components, temperature_vapour, Vapour_Pressure)



.. image:: output_21_0.png


sin embargo como se menciono anteriormente, es posible calcular una
propiedad termodinámica para un gran número de sustancias puras y luego
realizar las gráficas correspondientes dependiendo de las necesidades de
visualización entre otros criterios. Para ejemplificar esto, ahora se
tienen 7 sustancias puras y se quiere gŕaficar la propiedad
termodinámica de solo: *n-PENTACOSANE, ETHANE y el ISOBUTANE*.

.. code:: python

    components = ["METHANE", "n-TETRACOSANE", "n-PENTACOSANE", "ETHANE", "ISOBUTANE", "PROPANE", "3-METHYLHEPTANE"]
    property_thermodynamics = "Vapour_Pressure"
    
    Vapour_Pressure = thermodynamic_correlations.property_cal(components, property_thermodynamics)
    temperature_vapour = thermodynamic_correlations.temperature

.. code:: python

    thermodynamic_correlations.multi_graphical(components[2:5], temperature_vapour[2:5], Vapour_Pressure[2:5])



.. image:: output_24_0.png


5. Especificar varias sustancias puras y una temperatura.
=========================================================

Como en el caso anterios, en este ejemplo se espcifican 3 sustancias
puras pero con la especificación de un solo valor de temperatura. Esta
temperatura será común para las sustancias puras con las que se trabaje
por tanto puede darse el caso de que sea una temperatura valida para
algunas sustancias puras mientras que para otras no dependiendo del
intervalo de valides de cada correlación termodinámica.

.. code:: python

    dppr_file = "PureFull_mod_properties.xls"
    
    thermodynamic_correlations = pt.Thermodynamic_correlations(dppr_file)
    
    components = ["METHANE", "n-TETRACOSANE", "ISOBUTANE"]
    property_thermodynamics = "Vapour_Pressure"
    temperature = [180.4]
    
    Vapour_Pressure = thermodynamic_correlations.property_cal(components, property_thermodynamics, temperature)
    print("Vapour Pressure = {0} {1}". format(Vapour_Pressure, units[1]))
    



.. parsed-literal::

    ----------------------------------------------------------------------
    Pure substances with a temperature especific: ['METHANE', 'n-TETRACOSANE', 'ISOBUTANE']
    ----------------------------------------------------------------------
    [180.4]
    Temperature_enter = [[180.4], ['180.4 K is a temperature not valid'], [180.4]]
    Temperature_invalid = [[], ['180.4 K is a temperature not valid'], []]
    Temperature_valid = [array([ 180.4]), array([], dtype=float64), array([ 180.4])]
    vapour_Pressure =  [array([ 33.32655377]) array([], dtype=float64) array([ 0.0074373])] (3,)
    3
    Vapour Pressure = [array([ 33.32655377]) array([], dtype=float64) array([ 0.0074373])] [Pa]


en este caso se tiene como resultado un con 2 valores de presión de
vapor, uno para METHANE y otro para ISOBUTANE, mientras que se obtiene
un array vacio en el caso "de n-TETRACOSANE, puesto que la temperatura
de 180 K especificada no se encuentra como valida.

para verificar tanto los valores de las constantes como los valores
mínimos y máximos de cada correlación termodinámica para cada una de las
sustancias puras que se especifique se utiliza el atributo
*component\_constans* tal como se muestra a continuación

.. code:: python

    thermodynamic_correlations.component_constans





.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>A</th>
          <th>B</th>
          <th>C</th>
          <th>D</th>
          <th>E</th>
          <th>T Min [K]</th>
          <th>T Max [K]</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>METHANE</th>
          <td>39.205</td>
          <td>-1324.4</td>
          <td>-3.4366</td>
          <td>3.1019e-05</td>
          <td>2</td>
          <td>90.69</td>
          <td>190.56</td>
        </tr>
        <tr>
          <th>n-TETRACOSANE</th>
          <td>211.42</td>
          <td>-21711</td>
          <td>-26.255</td>
          <td>7.7485e-06</td>
          <td>2</td>
          <td>323.75</td>
          <td>804</td>
        </tr>
        <tr>
          <th>ISOBUTANE</th>
          <td>100.18</td>
          <td>-4841.9</td>
          <td>-13.541</td>
          <td>0.020063</td>
          <td>1</td>
          <td>113.54</td>
          <td>408.14</td>
        </tr>
      </tbody>
    </table>
    </div>



6. Especificar varias sustancias puras y especificar varias temperaturas
========================================================================

En esta opción se puede manipular varias sustancias puras de forma
simultanea con la especificación de varios valores de temperaturas, en
donde cada valor de temperatura especificado será común para cada
sustancia pura, de tal forma que se obtendra valores adecuados para
aquellos valores de temperatura que sean validos para cada caso
considerado.

.. code:: python

    import numpy as np
    import pandas as pd
    import pyther as pt
    import matplotlib.pyplot as plt
    %matplotlib inline

.. code:: python

    dppr_file = "PureFull_mod_properties.xls"
    
    thermodynamic_correlations = pt.Thermodynamic_correlations(dppr_file)
    
    #components = ["METHANE", "n-TETRACOSANE", "ISOBUTANE"]
    components = ["METHANE", "n-TETRACOSANE", "n-PENTACOSANE", "ETHANE", "ISOBUTANE", "PROPANE", "3-METHYLHEPTANE"]
    property_thermodynamics = "Vapour_Pressure"
    temperature = [180.4, 181.4, 185.3, 210, 800]
    
    Vapour_Pressure = thermodynamic_correlations.property_cal(components, property_thermodynamics, temperature)
    print("Vapour Pressure = {0}". format(Vapour_Pressure))


.. parsed-literal::

    ----------------------------------------------------------------------
    Pure substances with a temperature especific: ['METHANE', 'n-TETRACOSANE', 'n-PENTACOSANE', 'ETHANE', 'ISOBUTANE', 'PROPANE', '3-METHYLHEPTANE']
    ----------------------------------------------------------------------
    [180.4, 181.4, 185.3, 210, 800]
    Temperature_enter = [[180.4, 181.4, 185.3, '210 K is a temperature not valid', '800 K is a temperature not valid'], ['180.4 K is a temperature not valid', '181.4 K is a temperature not valid', '185.3 K is a temperature not valid', '210 K is a temperature not valid', 800], ['180.4 K is a temperature not valid', '181.4 K is a temperature not valid', '185.3 K is a temperature not valid', '210 K is a temperature not valid', 800], [180.4, 181.4, 185.3, 210, '800 K is a temperature not valid'], [180.4, 181.4, 185.3, 210, '800 K is a temperature not valid'], [180.4, 181.4, 185.3, 210, '800 K is a temperature not valid'], [180.4, 181.4, 185.3, 210, '800 K is a temperature not valid']]
    Temperature_invalid = [['210 K is a temperature not valid', '800 K is a temperature not valid'], ['180.4 K is a temperature not valid', '181.4 K is a temperature not valid', '185.3 K is a temperature not valid', '210 K is a temperature not valid'], ['180.4 K is a temperature not valid', '181.4 K is a temperature not valid', '185.3 K is a temperature not valid', '210 K is a temperature not valid'], ['800 K is a temperature not valid'], ['800 K is a temperature not valid'], ['800 K is a temperature not valid'], ['800 K is a temperature not valid']]
    Temperature_valid = [array([ 180.4,  181.4,  185.3]), array([800]), array([800]), array([ 180.4,  181.4,  185.3,  210. ]), array([ 180.4,  181.4,  185.3,  210. ]), array([ 180.4,  181.4,  185.3,  210. ]), array([ 180.4,  181.4,  185.3,  210. ])]
    7
    Vapour Pressure = [array([ 33.32655377,  34.43422601,  39.01608023]) array([ 9.23391967])
     array([ 7.9130031])
     array([ 0.80394112,  0.85063572,  1.05335836,  3.33810867])
     array([ 0.0074373 ,  0.00816353,  0.01160766,  0.07565701])
     array([ 0.05189654,  0.05605831,  0.07505225,  0.35872729])
     array([  2.09878094e-07,   2.50494222e-07,   4.89039104e-07,
             1.75089920e-05])]


como se muestra en los resultados anteriores, se comienza a complicar la
manipulación de los datos conforme incrementa el número de sustancias
puras y temperaturas involucradas en el analisis, por tal motivo
conviene utilizar las bondades de librerías especializadas para el
procesamiento de datos como *Pandas* para obtener resultados más
eficientes.

El método *data\_temperature(components, temperature, Vapour\_Pressure,
temp\_enter)* presenta un DataFrame con los resultados obtenidos luego
de calcular la propiedad termodinámica indicada, señalan que para las
temperaturas invalidas en el intervalo de aplicación de la correlación
termodinámica, el resultado será *NaN*, tal como se muestra con el
ejemplo a continuación.

.. code:: python

    temp_enter = thermodynamic_correlations.temperature_enter
    thermodynamic_correlations.data_temperature(components, temperature, Vapour_Pressure, temp_enter)




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>180.4 K</th>
          <th>181.4 K</th>
          <th>185.3 K</th>
          <th>210 K</th>
          <th>800 K</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>METHANE</th>
          <td>3.332655e+01</td>
          <td>3.443423e+01</td>
          <td>3.901608e+01</td>
          <td>NaN</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>n-TETRACOSANE</th>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>9.233920</td>
        </tr>
        <tr>
          <th>n-PENTACOSANE</th>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>7.913003</td>
        </tr>
        <tr>
          <th>ETHANE</th>
          <td>8.039411e-01</td>
          <td>8.506357e-01</td>
          <td>1.053358e+00</td>
          <td>3.338109</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>ISOBUTANE</th>
          <td>7.437302e-03</td>
          <td>8.163530e-03</td>
          <td>1.160766e-02</td>
          <td>0.075657</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>PROPANE</th>
          <td>5.189654e-02</td>
          <td>5.605831e-02</td>
          <td>7.505225e-02</td>
          <td>0.358727</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>3-METHYLHEPTANE</th>
          <td>2.098781e-07</td>
          <td>2.504942e-07</td>
          <td>4.890391e-07</td>
          <td>0.000018</td>
          <td>NaN</td>
        </tr>
      </tbody>
    </table>
    </div>



7. Future work
==============

-  Actualmente PyTher se encuentra implementando la opción de multiples
   propiedades termodinámicas de forma simultanea para el caso de
   multiples sustancias puras con multiples opciones de temepratura.

-  Dar soporte a la manipulación de bases de datos por parte de usuarios
   para agregar, modificar, eliminar, renombrar sustancias puras y/o
   correlaciones termodinámicas.

8. References
=============

Numpy

