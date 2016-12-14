5b. PyTher: Análisis computacíonal de consistencia termodinámica ELV
********************************************************************
********************************************************************

A. Salazar * , M. Cismondí


**Resumen**. En este trabajo se presenta la herramienta IPyTherm, la cual se
enfoca en cálculos termodinámicos del comportamiento de fases a través de
la plataforma Jupyter para realizar el análisis computacional de la
consistencia termodinámica de datos experimentales entre fases líquido-
vapor, permitiendo una manipulación eficiente de los datos experimentales
para determinar su calidad de forma programática e interactiva.

**Palabras clave**: PyTher, Termodinámica computacional, consistencia termodinámica, Python, Análisis de datos.

1. **Introducción**

Desde hace bastante tiempo se viene trabajando en la generación, recopilación y
procesamiento de los datos en el ámbito científico en distintas áreas de forma
programática, sin embargo, desde hace 20 años se presenta un crecimiento dramático de
datos que incluyen datos experimentales científicos reportados en literatura
especializada de dominio público (Frenkel, 2013), que sirve como punto de partida para
otras investigaciones, por ejemplo las involucradas en moldeamiento, simulación y
optimización. Por tanto, en el campo de la termodinámica los datos referidos a las áreas
de termofísica y termoquímica que son una fuente importante de datos tanto para la
investigación científica en áreas fundamentales como el desarrollo de tecnología y
aplicaciones en nuevos diseños de procesos y productos. Recientemente el
Thermodynamics Research Center (TRC) del US National Institute of Standards and
Technology (NIST), publico estadísticas referentes al crecimiento de los datos de
propiedades termofísicas y termoquímicas (Frenkel , 2015), reportadas en las 5
principales revistas especializadas en esta temática (Journal of Chemical and
Engineering Data, The Journal of Chemical Thermodynamics, Fluid Phase Equilibria,
Thermochimica Acta, and the International Journal of Thermophysics), mostrando que
la cantidad de datos se ha duplicado en los últimos 10 años y viene presentando un
crecimiento anual del 7% en el volumen de datos reportados. Esto se debe entre varias
cosas, por el aumento en la eficiencia y capacidad tecnológica de la medición de datos
experimentales de propiedades termofísicas y termoquímicas junto con la
automatización de sistemas de control y adquisición de datos para la medición de
presión, temperatura, concentración entre otras variables, lo que resulta en un aumento
en la productividad en la adquisición de datos, sin embargo, este aumento de
productividad no ha venido acompañada con el aumento de la capacidad de evaluar la
“calidad” de los datos medidos y reportados en la literatura especializada, debido a que
los equipos comerciales que tradicionalmente son empleados para realizar las
mediciones han sido desarrollados sin involucrar suficientemente personal altamente
calificado en cada temática específica, además del uso de software que utiliza metadatos
para completar de forma “engañosa” la información de propiedades termodinámicas
(Frenkel , 2015), que además tiene un factor agravante que es la dificultad de la
adecuada verificación por parte de los pares evaluadores de la gran cantidad de artículos
presentados para su publicación con un tiempo insuficiente para corroborar la calidad de
los datos experimentales reportados (Chirico et al, 2013; Frenkel et al, 2006).

En este trabajo se presenta la herramienta IPyTherm para el procesamiento y
visualización de datos experimentales del equilibrio de fases líquido-vapor, la cual se
basa en la tecnología de la plataforma IPython que en su tercera versión recibe el
nombre de Jupyter. Esta plataforma se desarrolla bajo el concepto del “peper
ejecutable” (Pérez and Granger, 2007; Pérez, 2013), puesto que frecuentemente en el
desarrollo de una investigación científica actual se requiere de la computación,
procesamiento, visualización y presentación de una gran cantidad de información y
datos que habitualmente se realiza con diferentes herramientas computacionales que no
siempre están adaptadas para funcionar juntas lo que implica un esfuerzo considerable,
tener que enfocarse en llevar datos de un formato a otro para poder avanzar en el
procesamiento, que principio no hace parte del objetivo de la investigación científica
que se está realizando, resultando en un proceso improductivo por el costo de tiempo
que involucra la manipulación de herramientas de cálculo científico tradicionalmente
implementado en lenguajes como FORTRAN, el cual es limitado para el procesamiento
y visualización de grandes cantidades de datos (M. Gaitan et al. 2012).


