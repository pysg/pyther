11. Estabilidad Material de las Fases
************************************
************************************

A temperatura y presión constantes una mezcla de composición z será estable termodinámicamente sólo si se encuentra en su estado o configuración de menor energía de Gibbs posible para tales condiciones. Por lo tanto, ninguna separación posible de una cantidad infinitesimal de una fase de composición w podrá disminuir aun mas la energía libre del sistema. Esto puede expresarse matemáticamente en lo que se conoce como la condición del plano tangente de Gibbs, que se muestra más adelante.

La condición de igualdad de fugacidades de los componentes de una mezcla de diferentes fases solo es una de las condiciones necesarias para establecer el equilibrio termodinámico. De esta forma, una mezcla es estable a :math:`(T, P)` si y solo si la energía libre de Gibbs total del sistema se encuentra en un **minimo global**. El cambio de la energía libre de Gibbs por la tranferencia de :math:` \delta n` moles de un componente i de una fase líquida a una fase de vapor es:

.. math:: \delta G = (\mu_i^v - \mu_i^l) \delta n_i

Por tanto, en el minimo global :math:` \delta G` debe ser cero para cualquier tranferencia de materia, mantienedo el cumplimiento de la igualdad del potencial químico o fugacidad como una condición necesaria del equilibrio de fases termodiámico.

Al considerar una fase de composición z, con potencial químico :math:` \delta (z)` y asumiendo que una cantidad infinitisimal :math:` \delta e` de una nueva fase de composición molar w es formada, entonces elcambio en la energía libre de Gibbs del sistema está expresada como sigue:

.. math:: \delta G = dn \sum\limits_{j} {w_i(\mu_i(w) - \mu_i(z))} \geq 0

donde una cantidad :math:` w_i \delta e` del componentes i es transferido, resultando en que una condición necesaria para la estabilidad termodinámica de la fase de composición z es que :math: ` \delta G` no sea negativo para cualquier :math:` \delta e` positivo, expresado en la siguiente desigualdad:

.. math:: \sum\limits_{i=1}^{c} {w_i(\mu_i(w) - \mu_i(z))} \geq 0 

para cualquier composición w, este resultado es denominada como la condición del plano tangente de Gibbs para evaluar la estabilidad termodinámica.


11.1 Resolución de la condición de estabilidad
---------------------------------------------

Para iniciar se considera una mezcla de C componentes de composición z a una temperatura T y presión P especificada, para escribir la condición suficiente de estabilidad de la mezcla como la función de la distancia del plano tangente `TPD(w)` por sus siglas en inglés:

.. math:: TPD(w) = \sum\limits_{i=1}^{c} {w_i(\mu_i(w) - \mu_i(z))} \geq 0 

que como se mencionó anteriormente, se exige que este sea no negativo para una composición w de una nueva fase en formación. Normalmente, conviene reescribir la TPD en terminos de la fugacidad al utilizar modelos de ecuaciones de estado (puede ser en terminos de otras variables termodinámicas) como sigue:

expresión del potencial químico a T, P y composición w

.. math:: \mu_i(T,P,w) = \mu_i^*(T,P_o) + RTln \left( \frac{f_i(T,P,w)}{P_o}\right)

y reemplazando el termino de la fugacidad del componente i en la mezcla a T, P y composición w

.. math:: \mu_i(T,P,w) = \mu_i^*(T,P_o) + RT(ln (w_i) + ln \frac{P)}{P_o} ln \phi_i(T,P,w)))

con lo cual se obtiene la expresión de la función de la distancia del plano tangente reducido tpd como se muestra a continuación:

.. math:: tpd(w) = \frac{TPD(w)}{RT} = \sum\limits_{i=1}^{C} {w_i(ln(w_i) + ln(\phi_i(w)) - ln(z_i) - ln (\phi_i(z)) )}

agrupando terminos 

.. math:: tpd(w) = \frac{TPD(w)}{RT} = \sum\limits_{i=1}^{C} {w_i \left( ln(w_i) + ln(\phi_i(w)) - d_i\right)}

donde 

.. math:: d_i = ln(z_i) + ln(\phi_i(z))


Por tanto, una apromaximación computacional puede ser basada en el hecho de que la condición del plano tangete es no negativa siempre, si y solo si es no negativa para todos los minimos, en ese sentido la recomendación de Michelsen & Mollurup [1]_, para implementar la evaluación de la condición del plano tangente son:

1. Localizar todos los minimos locales de la distancia del plano tangente.

2. Verificar que el valor de **tpd** es no negavita en todos los minimos. En caso de encontrar un valor negativo de **tpd** durante el procedimiento en alguno de los minimos locales de la función, la mezcla se evaluara como inestable.

11.1.1 Formas de resolver la función tpd
---------------------------------------

En primera instancia se puede mencionar los métodos de optimización para encontrar los minimos de la función tpd, sin embargo, en está sección se presentara brevemente la estrategía de expresar este problema como un problema de un sistema de ecuaciónes algebraicas no lineales.


.. math:: tm(W) = 1 + \sum\limits_{i}^{C} {W_i(ln(W_i + ln \phi_i(W) - d_i -1)}


.. math:: \frac{\partial tm}{ \partial W_i} = lnW_i + ln \phi_i(W) - d_i


.. math:: tm(W)^{SP} = 1 - W_T

.. math:: W_T = \sum\limits_{i}^{C} {W_i}


.. math:: tm(W) = 1 + W_T \sum\limits_{i}^{C} {w_i(ln(W_T + ln w_i + ln \phi_i - d_i - 1)}

.. math:: tm(W) = (1 - W_T + W_TlnW_T) + W_Ttpd(w)


Método de solución

.. math:: lnW_{i}^{k+1} = d_i - ln\phi_i(W^{k})

.. math:: \phi_i(W) = \phi_i(W_i)

.. math:: w_i = \frac{W_i}{W_T}


.. math:: g_i = ln W_i + ln \phi_i - d_i

Matriz Hessiana

.. math:: H_{ij} = \frac{\partial g_i} {\partial W_j} = \frac{1}{W_i} \sigma_{ij} + \frac{\partial ln \phi_i}{\partial W_i}


corrector de Newton

.. math:: H \Delta W + g = 0

.. math:: W^{k+1} = w^{k} \Delta W






