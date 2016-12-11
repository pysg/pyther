33. Articulo: Data science to phase equilibrium
*************************
*************************

3.1 Introduccion
----------------

En esta sección se presenta una breve explicación de la primera parte del diagrama de clases que corresponde al cálculo de la fugacidad de un componente puro. Por tanto, los nombres de clases, atributos y métodos están sujetos a modificaciones, actualizaciones y correciones, junto con la estructura de las relaciones entre clases presentadas aquí, así como también las que se puedan adicionar como resultado de la continua revisión del diagrama de clases.

En este trabajo se presenta una librería open source pyther que busca facilitar la manipulacion, procesamiento y visualización de datos relevantes para el equilibrio termodinámico de fases de sustancias puras.


Recientemente se vienen dando grandes avances en materia del desarrollo de software de interes cientifico por parte de investigadores de diferentes áreas con poca formación en desarrollo de software y buenas practicas de programación. en ese sentido se vienen dando una series de herramientas que facilitan esta incorporación del desarroolo de software para outsiders que cuenten facilmente con más simplicidad en materia de codigo para enfocarse en la implmentación de servicios y microservicios sin detenerse tanto en la parte del desarrollo de software pero al mismo tiempo que esto no implique ni resulte en codigo mal implementado, codigo espagueti, no escalable ni facil de mantener por carencia de buenas practicas de implementación y ausencia de documentación.


El enfoque que se viene dando es la automatización de diferentes tareas repetitivas que consumen tiempo que los investigadores pueden destinar al analisis de los resultados en lugar de hacer estas tareas monotonas cada vez que requiere hacer un ensayo.


En este trabajo se toma un enfoque de desarrollo de una librería que facilite la manipulación, procesamiento y visualización de datos relevantes al equilibrio de fases y fluidos supercriticos.

También se enfoca en desarrollar una documentación que permita llegar esta librería a interesados que no estén familiarizados con este tipo de herramientas, de forma que lo puedan poner a prueba y aportar mejoras y correcciónes.

Otro de los puntos importantes que se tratan en este trabajo es la capcidad de desarrollo de harramientas computaciones que predigan los puntos de interes termodinámico 


este desarrollo parte de la manipulación de las propiedades termofisicas de sustancias puras para luego implementar la capcidad de mezclas y sistemas binarios


en este contexto se propician herramientas que buscan hacer el desarrollo, la implementación y sobre todo la difusión de información teneindo un fuerte enfasis en lo interactivo



Otra de las psoibilidades que hace fuerte el ecosistema cientifico de Python es la posiblidad de utilizar librerías especializadas para diferentes campos de la cienca, como por ejemplo la librería Pandas que se enfoca en la manipulación de datos, la librería mathplotlib enfica en la visualización de datos.

En el caso de la termodinámica se encuentran diferentes tipos de sotfware esfocados en el equilibrio de fases principalmente desarrolado en lenguajes como FORTRAN manteniendo librerías de calculo numérico que no son libres.

En caso de las sustancias puras, en este trabajo se presenta la ecuación de estado RKPR la cual es una ecuación de estado de 3 parámetros 




Figura 1. Diagrama de Clases Pyther

.. image:: _static/diagrama_clase_8.jpeg
	:width: 1200

La primera parte del diagrama de clases corresponde a:

1. **DatosComponentesPuros**
2. **CondicionesSistema**
3. **Componente**
4. **ParametrosBD**
5. **PropiedadesVolumetricas**
6. **ModulosMM**
7. **PropiedadesTermodinamicas**

La segunda parte del diagrama de clases que será comentado en el siguiente avance corresponde a:

8. **SolidoPuro**
9. **Solido Fluido**
10. **RegresionParametros**
11. **Flash_i**
12. **Flash_Fi**
13. **Estabilidad_Material**
14. **Interfaz Gráfica**

3.2 Clase DatosComponentesPuros
-------------------------------

En la primera clase **DatosComponentesPuros** se tiene:

- Atributos

DIPPR = Este atributo es una variable tipo string que corresponde al nombre que tiene el archivo que actualmente hace de "base de datos" provisional y se verifica que el nombre del archivo coincida con el preestablecido **DPPR** para mostrar por pantalla si se ha cargado o no los datos correctamente. Cuando se adicione la posibilidad de otras "bases de datos", en esta clase se deberá contar con más atributos para manipularlas adecuadamente.

- Métodos

LeerBaseDatos() = Carga los datos del archivo "DPPR" en una variable del sistema para poder manipular dichos datos a conveniencia.

AgregarBaseDatos() = Carga los datos de un archivo con nombre diferente al archivo por defecto "DPPR". Nota: Falta generalizar el formato en el que se pretatarian los diferentes archivos con datos supuestos para que se puedan manipular dentro del sistema.

ModificarBaseDatos() = Crea una copia del archivo "DPPR" en el cual se modifica uno o más valores de los registros del archivo o adiciona un campo nuevo cuyo nombre es especificado por el usuario. Falta generalizar la opción dehacer una agrupación de componentes de acuerdo a un criterio para crear dichos "nuevos" pseudocomponentes.

CrearBaseDatos() = Crea un archivo con datos obtenidos durante la realización de cálculos, por ejemplo la regresión de parametros o puntos importantes de diagrama de fases por mencionar algunas posibilidades para que dicha información se almacene de forma estructurada para su uso en calculos posteriores sin requerir realizar de nuevo el calculo. Actualmente en pruebas. 

3.3 Clase CondicionesSistema
----------------------------

En la segunda clase **CondicionesSistema**

Esta clase tiene como objetivo capturar del usuario las condiciones del sistema al cual se realizará los cálculos, como lo son temperatura, presión, fracción molar, volumén (según sea el caso), el modelo (ecuación de estado/modelo sólido puro) y el componente. 

- Atributos

Se tienen los siguientes atributos

1 Temperatura 
2 Presión 
3 Fracción Molar 
4 Volumen 
5 Modelo
6 Componentes


3.4 Clase Componente
--------------------

Esta clase tiene como objetivo la definición del o los componentes que se manejaran para realizar un cálculo con base a los registros (que se identifican con el nombre de una sustancia química) seleccionados de la clase **DatosComponentesPuros** a las condiciones establecidas en la clase **CondicionesSistema**. Luego se crea cada componente de acuerdo al modelo especificado en la clase **CondicionesSistema**), por ejemplo METHANE-SRK.

- Atributos

propiedadesFQ = Corresponde a un array que contiene las propiedades (temperatura critica, presión critica y factor acentrico) que se definió en la selección del nombre de la sustancia química que se quiere utilizar.

CondicionesSistema = Corresponde a un array que contiene la definición de la temperatura, presión fracción molar, modelo y nombre de la sustancia química que se quiere utilizar.

- Métodos

ModeloSRK
ModeloPR
ModeloRKPR

Los métodos (ModeloSRK, ModeloPR, ModeloRKPR) corresponden al cálculo de los parametros requeridos para los modelos SKR, PR, RKPR según sea el caso que se especifique en la clase **CondicionesSistema**. 

3.5 Clase ParametrosBD
----------------------
 
Esta clase obtiene la información del o los **componentes**, por ejemplo ("METHANE SRK"), para calcular los parámetros B y D correspodientes.

- Atributos

componente = es un array que contiene los parámetros necesarios para cálcular las variables B y D 

- Métodos

Parametro B = Calcula el parametro B con la información provista en **componente** 
Parametro D = Calcula el parametro D con la información provista en **componente**

3.6 Clase PropiedadesVolumetricas
---------------------------------

Esta clase tiene como objetivo la manipulación de la ecuación de estado cúbica para determinar la presión, temperatura o volúmen según sea el caso de las especificaciones dadas en la clase **CondicionesSistema**. Por ejemplo, al especificar P, T y ni, encontrar el V en dichas condiciones y un modelo y parametros determinados. Esta clase se separa de de la clase **Modulos MM** (se muestra a continuación) para aprovechar el enfoque modular y acceder al calculo de propiedades volumetricas de forma independiente del calculo de propiedades termodinámicas y sus correspondientes modulos (funcion de helmholtz, primeras derivas y segundas derivadas), según sean requeridas (las propiedades volumetricas). 

- Atributos

Parametro B = parametro B determinado en la clase **ParametrosBD** 
Parametro D = parametro D determinado en la clase **ParametrosBD**
Optimizador = corresponde a la selección y especificación de los parámetos requeridos para acceder y ejecutar un método ńumérico de resolución de ecuaciones no lineales de la librería Scipy. 

- Métodos

Volumen = cálcua el volumén con una ecuación de estado para una P, T y ni especificados
Temperatura = cálcua la temperatura con una ecuación de estado para una P, V y ni especificados (Falta por implementar). 
Presión = cálcua la presión con una ecuación de estado para una T, V y ni especificados

En caso de especificiar el volumen V, se calcula la presión P para la temperatura T y ni especificada. Para el caso contrario de especificar la presión P, se determina el volumen V para la temperatura T y ni especificada.

3.7 Clase ModulosMM
-------------------

Esta clase se tiene como objetivo calcular la función de energía de Helmholtz siguiendo el enfoque modular de Michelsen & Mollerup, partiendo de los parametros B y D obtenidos en la clase **ParametrosBD** y la propiedad volumetrica "volumen" o "presión" según sea el caso especificado (Esta clase tiene la capacidad de navegar y acceder a los otros atributos como lo son la temperatura, fracción molar). En esta clase se tienen tres métodos, que calculan la función de energía de Helmholtz ya mencionada, las primeras derivadas de esta función con respecto a las variables como son: emperatTura, Presión, Volumen y Número de moles (para el caso de la fracción molar hay relaciones que permiten obtener las derivadas en función de las fracciones molares a partir de las derivadas del númerod de moles), así mismo para el caso de las segundas derivadas de la función de energía de Helmholtz.

- Atributos

Parametro B = parametro B determinado en la clase **ParametrosBD** 
Parametro D = parametro D determinado en la clase **ParametrosBD**
Volumen = corresponde al volumén calculado con una ecuación de estado para una P, T y ni especificados
Presión = corresponde a la presión con una ecuación de estado para una T, V y ni especificados

En esta clase los atributos de presión P, volumen V se acceden desde la clase **PropiedadesVolumetricas**y como ya se ha mencionado estos pueden ser una especificación o calculados según sea el caso.

- Métodos

funciónHelmholtz = este método calcula la función de energía de Helmholtz con los parametros indicados para la especificación del modelo (por ejemplo METHANE SKR) y las condiciones del sistema.

primerasDerivadas = este método calcula las primeras derivadas de la función de energía de Helmholtz con respecto a las variables como son: Temperatura, Presión, Volumen y Número de moles (para el caso de la fracción molar hay relaciones que permiten obtener las derivadas en función de las fracciones molares a partir de las derivadas del númerod de moles), a las vcon los parametros indicados para la especificación del modelo (por ejemplo METHANE SKR) y las condiciones del sistema.

segundasDerivadas = este método calcula las segundas derivadas de la función de energía de Helmholtz con respecto a las variables como son: Temperatura, Presión, Volumen y Número de moles (para el caso de la fracción molar hay relaciones que permiten obtener las derivadas en función de las fracciones molares a partir de las derivadas del númerod de moles), , a las vcon los parametros indicados para la especificación del modelo (por ejemplo METHANE SKR) y las condiciones del sistema.

3.8 Clase PropiedadesTermodinamicas
-----------------------------------

En esta clase se tiene los métodos para calcular las propiedades termodinámicas siguiendo el enfoque modular de Michelsen & Mollerup. Esta clase no tiene atributos y sus métodos corresponden a las propiedades termodinámicas como: Fugacidad, Entalpía y Entropía. (Se está implementando para el método de la energía libre de Gibbs)

- Atributos

No tiene atributos.

- Métodos

Fugacidad = este método calcula la fungacidad de un componente puro o mezcla multicomponente, según sea la especificación (puro o multicomponente) siguiendo el enfoque modular de Michelsen & Mollerup partiendo de los métodos de la clase **ModulosMM**, que ya contienen toda la información pertinente para realizar el calculo de la propiedad termodinámica.

Entalpía = este método calcula la entalpía de un componente puro o mezcla multicomponente, según sea la especificación (puro o multicomponente) siguiendo el enfoque modular de Michelsen & Mollerup partiendo de los métodos de la clase **ModulosMM** para el calculo de las primeras y segundas derivadas de la función de energía de Helmholtz, que ya contienen toda la información pertinente para realizar el calculo de la propiedad termodinámica.

Entropía = este método calcula la entropía de un componente puro o mezcla multicomponente, según sea la especificación (puro o multicomponente) siguiendo el enfoque modular de Michelsen & Mollerup partiendo de los métodos de la clase **ModulosMM** para el calculo de las primeras y segundas derivadas de la función de energía de Helmholtz, que ya contienen toda la información pertinente para realizar el calculo de la propiedad termodinámica.

.. Note:: para el caso de las propiedades termodinámica aún no se han terminado de realizar las pruebas que corroboren que los calculos implementados tienen resultados correctos. 

3.9 Clase Estabilidad_Material
------------------------------

En esta clase falta por empezar a documentarla.





