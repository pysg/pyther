10. Termodinámica de fluidos 
****************************
****************************


En la termodinámica clásica los efectos de la capilaridad y la la segragación por gravedad no son tenidos encuenta en el análisis PVT del equilibrio de fases


Phase  equilibrium  of  fluids  confined  in  porous  media  from an  extended  Peng–Robinson  equation  of  state

Para abordar la termodinámica de fluidos confinados en shale media, se han destacado principalmente 3 enfoques: la simulación molecular, los modelos de adsorción convencionales y el creciente campo de la extensión de ecuaciones de estado [Ref_x].

1. En el primer método se cuenta con una buena aproximación del cálculo de las propiedades termodinámicas, sin embargo, es costoso computacionalmente realizar este tipo de cálculos para sistemas complejos y/o que requieran realizarse como parte de un simulador que involucre cálculos de otra índole más allá de un perfil de concentración.

2. En el segundo método, se enfoca en determinar las propiedades termodinámicas del sistema subdividiendo el problema, por un lado la fases que se encuentra confinada en el interior del medio poroso se modela utilizando el enfoque convencional de las isotermas de adsorción, mientras que las propiedades de la fases presenten en el bulk, se modela con enfoques, por ejemplo las tradicionales como las ecuaciones, dando lugar a  problemas de consistencia en el modelo puesto que no se puede modelar el sistema completo con un solo set de parámetros, lo cual es una inconsistencia.

3. El tercer método se concentra en extender los modelos de las ecuaciones de estado, que han sido ampliamente probadas para calcular las propiedades termodinámicas para fases no confinadas, para ampliar su rango de aplicación con fases confinadas en nanoporos, en donde el equilibrio de fases involucra una mayor complejidad que el tratamiento convencional puede representar.



Modelo de tensión superficial:

.. math:: P^v - \frac {2 \sigma Cos(\theta)} {R} - P^l = 0

.. math:: ln f_i^l(T,P^l,x_i) - lnf_i^v(T, P^v, y_i) = 0

.. math:: \sum\limits_{i=1}^{c} {x_i-1} = 0

.. math:: \sum\limits_{i=1}^{c} {y_i-1} = 0

.. math:: a = \sum\limits_{i=1}^{c} {\sum\limits_{j=1}^{c} {x_ix_j \sqrt{a_i a_j} (1-k_{ij}) }}

.. math:: a = \sum\limits_{j=1}^{c} {x_ib_i}

.. math:: \sigma^ \frac{1}{4} = P (\rho^l - \rho^v)

.. math:: \sigma = kT_c \left( \frac{N_A}{V_c}\right) (4.35 + 4.14w) t^{1.26}

.. math:: t = 1 - \frac{T}{T_c} 

.. math:: \sigma^ \frac{1}{4} = \sum\limits_{i=1}^{c} {P_i (x_i\rho^l - y_i\rho^v)}

.. math:: \sigma^ \frac{1}{E} = \sum\limits_{i=1}^{c} {P_i (x_i\rho^l - y_i\rho^v)}

.. math:: E = 3.583+0.16(\rho^l - \rho^v)

.. math:: \frac {\sigma_R}{\sigma_plane} =  \frac {1}{1+\frac{2\sigma}{R}}

.. math:: \frac{P^C}{P^D} = \frac{-V_mix(P^D,y)}{V^l(y)} Z_av(y)ln\chi



El presente trabajo doctoral se oriente sobre este tercer enfoque, el del análisis de la extensión de ecuaciones de estado para representar el comportamiento termodinámico de las fases en el bulk y en el interior de los nanoporos. En este sentido, se destaca el trabajo que realizado por [], en cual  muestra una extensión de

.. math:: \prod_{i=1}^{3}x_{i}

.. math:: Q(T,V,N_1,N_2,...,N_nc) = \prod_{i=1}^{c} \left(\frac{q_i^{n_i}} {\lambda_i^{3N_i} N_i!} \right) V_f^N exp\left( \int_{\infty}^{T}\frac{E_{conf}}{kT^2} \right)



























































