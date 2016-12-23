5a. PyTher: Model, parameters and properties component
******************************************************
******************************************************

Temperatura critica
Presión critica
factor acentrico Omega
Volumen critico

 





importar las linrerías requeridas, en este caso se trata de las
librerías numpy, pandas junto con pyther

.. code:: python

    import numpy as np
    import pandas as pd
    import pyther as pt

En los ejemplos siguientes se utilizan los datos termodísicos de la base
de datos DPPR. Para el caso se tiene como especificación la ecuación de
estado RKPR y las constantes criticas para el componente
3-METHYLHEPTANE a continuación.

.. code:: python

    properties_data = pt.Data_parse()
    
    dppr_file = "PureFull.xls"
    component = "3-METHYLHEPTANE"
    
    NMODEL = "RKPR"
    ICALC = "constants_eps"
    
    properties_component = properties_data.selec_component(dppr_file, component)
    pt.print_properties_component(component, properties_component)
    dinputs = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
                        properties_component[1]['Omega'], properties_component[1]['Vc']])
    
    component_eos = pt.models_eos_cal(NMODEL, ICALC, dinputs)
    
    #ac = component_eos[0]
    print(component_eos)
    
    



.. parsed-literal::

    Component = 3-METHYLHEPTANE
    Acentric_factor = 0.3718
    Critical_Temperature = 563.67 K
    Critical_Pressure = 25.127 Bar
    Critical_Volume = 0.464 cm3/mol
    Compressibility_factor_Z = 0.252
    del1ini = 6.038268203938681
    Zc = 0.24877058378575795
    The NMODEL is eos_RKPR and method ICALC is constants_eps
    params = [ac, b, rm, del1]
    [46.430578671675555, 0.10935115096084358, 2.5860921512475117, 6.0431253541984447]


De esta forma se observa el calculo simple de los parámetros para la
sustancia pura 3-METHYLHEPTANE\_RKPR

A continuación se realiza el mismo tipo de calculo pero tomando una
serie de 9 sustancias puras, que se pueden extender facilmente a n
sustancias, para obtener sus parámetros de nuevo con la ecuación de
estado RKPR.

.. code:: python

    properties_data = pt.Data_parse()
    
    dppr_file = "PureFull.xls"
    components = ["ISOBUTANE", "CARBON DIOXIDE", 'METHANE', "ETHANE", "3-METHYLHEPTANE", "n-PENTACOSANE",
                  "NAPHTHALENE", "m-ETHYLTOLUENE", "2-METHYL-1-HEXENE"]
    
    NMODEL = "RKPR"
    ICALC = "constants_eps"
    component_eos_list = np.zeros( (len(components),4) )
    
    for index, component in enumerate(components):
        
        properties_component = properties_data.selec_component(dppr_file, component)
        pt.print_properties_component(component, properties_component)
        dinputs = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
                            properties_component[1]['Omega'], properties_component[1]['Vc']])
        
        component_eos = pt.models_eos_cal(NMODEL, ICALC, dinputs)
        component_eos_list[index] = component_eos 
    
        
    components_table = pd.DataFrame(component_eos_list, index=components, columns=['ac', 'b', 'rm', 'del1'])
    
    print(components_table)
    



.. parsed-literal::

    Component = ISOBUTANE
    Acentric_factor = 0.18080000000000002
    Critical_Temperature = 408.14 K
    Critical_Pressure = 36.003 Bar
    Critical_Volume = 0.2627 cm3/mol
    Compressibility_factor_Z = 0.28200000000000003
    del1ini = 3.9722378008963446
    Zc = 0.27871152548257544
    The NMODEL is eos_RKPR and method ICALC is constants_eps
    params = [ac, b, rm, del1]
    Component = CARBON DIOXIDE
    Acentric_factor = 0.22360000000000002
    Critical_Temperature = 304.21 K
    Critical_Pressure = 72.865 Bar
    Critical_Volume = 0.094 cm3/mol
    Compressibility_factor_Z = 0.274
    del1ini = 4.462908059336361
    Zc = 0.2707937660977233
    The NMODEL is eos_RKPR and method ICALC is constants_eps
    params = [ac, b, rm, del1]
    Component = METHANE
    Acentric_factor = 0.0115
    Critical_Temperature = 190.564 K
    Critical_Pressure = 45.389 Bar
    Critical_Volume = 0.09860000000000001 cm3/mol
    Compressibility_factor_Z = 0.28600000000000003
    del1ini = 3.7519407434981633
    Zc = 0.2824567739174239
    The NMODEL is eos_RKPR and method ICALC is constants_eps
    params = [ac, b, rm, del1]
    Component = ETHANE
    Acentric_factor = 0.0995
    Critical_Temperature = 305.32 K
    Critical_Pressure = 48.083 Bar
    Critical_Volume = 0.14550000000000002 cm3/mol
    Compressibility_factor_Z = 0.279
    del1ini = 4.161423913263858
    Zc = 0.2755907402334964
    The NMODEL is eos_RKPR and method ICALC is constants_eps
    params = [ac, b, rm, del1]
    Component = 3-METHYLHEPTANE
    Acentric_factor = 0.3718
    Critical_Temperature = 563.67 K
    Critical_Pressure = 25.127 Bar
    Critical_Volume = 0.464 cm3/mol
    Compressibility_factor_Z = 0.252
    del1ini = 6.038268203938681
    Zc = 0.24877058378575795
    The NMODEL is eos_RKPR and method ICALC is constants_eps
    params = [ac, b, rm, del1]
    Component = n-PENTACOSANE
    Acentric_factor = 1.1053
    Critical_Temperature = 812 K
    Critical_Pressure = 9.376 Bar
    Critical_Volume = 1.46 cm3/mol
    Compressibility_factor_Z = 0.20500000000000002
    del1ini = 10.600246415857843
    Zc = 0.20275882073834256
    The NMODEL is eos_RKPR and method ICALC is constants_eps
    params = [ac, b, rm, del1]
    Component = NAPHTHALENE
    Acentric_factor = 0.3022
    Critical_Temperature = 748.35 K
    Critical_Pressure = 39.98 Bar
    Critical_Volume = 0.41300000000000003 cm3/mol
    Compressibility_factor_Z = 0.269
    del1ini = 4.8204311891035925
    Zc = 0.2653709654843225
    The NMODEL is eos_RKPR and method ICALC is constants_eps
    params = [ac, b, rm, del1]
    Component = m-ETHYLTOLUENE
    Acentric_factor = 0.3226
    Critical_Temperature = 637.15 K
    Critical_Pressure = 28.029 Bar
    Critical_Volume = 0.49 cm3/mol
    Compressibility_factor_Z = 0.263
    del1ini = 5.246526144851435
    Zc = 0.2592551086535563
    The NMODEL is eos_RKPR and method ICALC is constants_eps
    params = [ac, b, rm, del1]
    Component = 2-METHYL-1-HEXENE
    Acentric_factor = 0.3094
    Critical_Temperature = 538 K
    Critical_Pressure = 28.325 Bar
    Critical_Volume = 0.398 cm3/mol
    Compressibility_factor_Z = 0.255
    del1ini = 5.784189965441039
    Zc = 0.2520206003977051
    The NMODEL is eos_RKPR and method ICALC is constants_eps
    params = [ac, b, rm, del1]
                               ac         b        rm       del1
    ISOBUTANE           15.743219  0.064343  2.205509   4.000470
    CARBON DIOXIDE       4.409808  0.022801  2.280728   4.492210
    METHANE              2.696405  0.024259  1.282178   3.777713
    ETHANE               6.649597  0.035503  1.673541   4.190762
    3-METHYLHEPTANE     46.430579  0.109351  2.586092   6.043125
    n-PENTACOSANE      289.947431  0.320522  4.581358  10.628260
    NAPHTHALENE         49.312554  0.099495  2.591582   4.847168
    m-ETHYLTOLUENE      51.786960  0.117115  2.565531   5.267361
    2-METHYL-1-HEXENE   37.214555  0.094214  2.338038   5.794610


Como se observa, los resultados obtenidos son organizados en un
DataFrame permitiendo agilizar la manipulación de los datos de una serie
de sustancias puras.

.. code:: python

    components_table




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>ac</th>
          <th>b</th>
          <th>rm</th>
          <th>del1</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>ISOBUTANE</th>
          <td>15.743219</td>
          <td>0.064343</td>
          <td>2.205509</td>
          <td>4.000470</td>
        </tr>
        <tr>
          <th>CARBON DIOXIDE</th>
          <td>4.409808</td>
          <td>0.022801</td>
          <td>2.280728</td>
          <td>4.492210</td>
        </tr>
        <tr>
          <th>METHANE</th>
          <td>2.696405</td>
          <td>0.024259</td>
          <td>1.282178</td>
          <td>3.777713</td>
        </tr>
        <tr>
          <th>ETHANE</th>
          <td>6.649597</td>
          <td>0.035503</td>
          <td>1.673541</td>
          <td>4.190762</td>
        </tr>
        <tr>
          <th>3-METHYLHEPTANE</th>
          <td>46.430579</td>
          <td>0.109351</td>
          <td>2.586092</td>
          <td>6.043125</td>
        </tr>
        <tr>
          <th>n-PENTACOSANE</th>
          <td>289.947431</td>
          <td>0.320522</td>
          <td>4.581358</td>
          <td>10.628260</td>
        </tr>
        <tr>
          <th>NAPHTHALENE</th>
          <td>49.312554</td>
          <td>0.099495</td>
          <td>2.591582</td>
          <td>4.847168</td>
        </tr>
        <tr>
          <th>m-ETHYLTOLUENE</th>
          <td>51.786960</td>
          <td>0.117115</td>
          <td>2.565531</td>
          <td>5.267361</td>
        </tr>
        <tr>
          <th>2-METHYL-1-HEXENE</th>
          <td>37.214555</td>
          <td>0.094214</td>
          <td>2.338038</td>
          <td>5.794610</td>
        </tr>
      </tbody>
    </table>
    </div>



En el siguiente ejemplo se utiliza la ecuación RKPR pero esta vez con la
especificación de la temperatura y densidad de líquido saturado para el
CARBON DIOXIDE y de esta forma encontrar el valor del parámetro *delta*
que verifica la especificación realizada para la densidad de líquido
saturado.

.. code:: python

    properties_data = pt.Data_parse()
    
    dppr_file = "PureFull.xls"
    component = "CARBON DIOXIDE"
    
    NMODEL = "RKPR"
    ICALC = "density"
    
    properties_component = properties_data.selec_component(dppr_file, component)
    pt.print_properties_component(component, properties_component)
    #dinputs = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
    #                    properties_component[1]['Omega'], properties_component[1]['Vc']])
    
    T_especific = 270.0
    RHOLSat_esp = 21.4626
    # valor initial of delta_1
    delta_1 = 1.5
    
    dinputs = np.array([properties_component[1]['Tc'], properties_component[1]['Pc'],
                        properties_component[1]['Omega'], delta_1, T_especific, RHOLSat_esp])
    
    
    component_eos = pt.models_eos_cal(NMODEL, ICALC, dinputs)
    
    print(component_eos)


.. parsed-literal::

    Component = CARBON DIOXIDE
    Acentric_factor = 0.22360000000000002
    Critical_Temperature = 304.21 K
    Critical_Pressure = 72.865 Bar
    Critical_Volume = 0.094 cm3/mol
    Compressibility_factor_Z = 0.274
    The NMODEL is eos_RKPR and method ICALC is density
    The parameter delta1(rho,T) = [ 2.65756708]
    [ 2.65756708]


