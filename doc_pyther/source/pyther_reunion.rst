
importar las linrerías requeridas, en este caso se trata de la numpy
junto con pyther

.. code:: python

    import numpy as np
    import pyther as pt

.. code:: python

    properties_data = pt.Data_parse()
    
    dppr_file = "PureFull.xls"
    #component = 'METHANE'
    #component = "ETHANE"
    #component = "3-METHYLHEPTANE"
    #component = "n-PENTACOSANE"
    component = "ISOBUTANE"
    component = "CARBON DIOXIDE"
    
    NMODEL = "RKPR"
    ICALC = "constants_eps"
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
    
    #ac = component_eos[0]
    print(component_eos)
    print('-' * 79)
    



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
    -------------------------------------------------------------------------------


