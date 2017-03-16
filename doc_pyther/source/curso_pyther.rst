13. Curso de postgrado: Termodinámica de fluidos
************************************************
************************************************

-  **Parte 1: Python básico**
-  **Parte 2: Termodinámica del equilibrio de fases**
-  **Parte 3: Termodinámica Propiedades de exces**
-  **Parte 4: Thermodynamics and Data**
-  **Parte 5: Proyecto final**

13.1 Contenido del curso
------------------------

-  **Parte 1: Python básico**
*****************************

1. ¿Qué es Python?
   
   - ¿Por qué Python?
   - ¿Por qué Python para ciencia e ingeniería?
   - Jupyter: Literate Computing environment
   - IPython notebook

2. Ptyhon básico

   -  Instalación y primeros pasos
   -  Tipos de datos en Python
   -  Control de flujo

      -  La sentencia if
      -  La sentencia for
      -  La función range()
      -  Las sentencias break, continue, y else en lazos
      -  La sentencia pass
      -  Definiendo funciones
      -  Más sobre definición de funciones
      -  Argumentos con valores por omisión
      -  Palabras claves como argumentos
      -  Listas de argumentos arbitrarios
      -  Desempaquetando una lista de argumentos
      -  Expresiones lambda
      -  Cadenas de texto de documentación
      -  Anotación de funciones

   -  Estructuras de datos

      -  Más sobre listas
      -  Listas por comprensión anidadas
      -  La instrucción del
      -  Tuplas y secuencias
      -  Conjuntos
      -  Diccionarios

   -  Módulos

      -  Ejecutando módulos como scripts
      -  Búsqueda de los módulos
      -  Archivos "compilados" de Python
      -  Módulos estándar
      -  La función dir()
      -  Paquetes

         -  Importando \* desde un paquete
         -  Referencias internas en paquetes
         -  Paquetes en múltiples directorios

   -  Entrada y salida

      -  Formateo elegante de la salida
      -  Leyendo y escribiendo archivos
      -  Métodos de los objetos Archivo

   -  Errores y excepciones

      -  Errores de sintaxis
      -  Excepciones
      -  Manejando excepciones
      -  Levantando excepciones
      -  Excepciones definidas por el usuario

   -  Entornos Virtuales y Paquetes

      -  Introducción
      -  Creando Entornos Virtuales
      -  Manejando paquetes con pip y anaconda

   -  Programación orientada a objetos básica

      -  Un primer vistazo a las clases
      -  Sintaxis de definición de clases
      -  Objetos clase
      -  Objetos instancia
      -  Objetos método
      -  Variables de clase y de instancia
      -  Herencia
      -  Herencia múltiple
      -  Variables privadas

3. Computación cientifica en el ecosistema Python

   -  NumPy y SciPy para computación cientifica e ingenieria
   -  matplotlib y Bokeh para visualización de datos
   -  Pandas para análisis de datos

-  **Parte 2: Termodinámica del equilibrio de fases**
******************************************************

1. Calculation of Thermodynamic Properties

   -  The Helmholtz function
   -  Thermodynamic properties of the Helmholtz function .
   -  Test of fugacity coefficients and partial derivatives.
   -  Calculation of the partial derivatives of F

      -  First order derivatives .
      -  Seond order derivatives

   -  Sympy for symbolic mathematics of the Helmholtz function

2. Thermodynamic Properties from a Cubic Equation of State

   -  The cubic equation of state .
   -  The pure component parameters

      -  The critical point
      -  Subcritical temperatures
      -  Supercritical temperatures
      -  The temperature dependence of b

   -  Mixtures.

      -  The Helmholtz function
      -  The Lorentz-Berthelot combining rules.

   -  Derivatives of the Helmholtz function .

      -  The derivatives of g(() and f(()

   -  Calculation of the volume.

      -  Calculation of the volume in a two-parameter equation of state.

   -  Elimination of the gas constant.
   -  Numerical example with Jupyter

3. The Isothermal Two-Phase Flash

   -  Successive substitution and the Rachford-Rice equation
   -  Convergence analysis
   -  Initial estimates .
   -  Accelerated direct substitution
   -  Gibbs energy minimization by second order methods
   -  Strategy for a flash algorithm
   -  Tangent plane analysis
   -  Locating the minima of tm
   -  Procedures for minimizing tm
   -  Hybrid models .
   -  Liquid-liquid equilibrium
   -  When speed counts.
   -  Numerical example with Jupyter

4. The Multiphase Isothermal Flash

   -  Successive substitution
   -  Pure phases and solids
   -  Acceleration of successive substitution.
   -  Gibbs energy minimization by second order methods
   -  Stability analysis
   -  Near-critical phases
   -  Numerical example with Jupyter

5. Saturation Points and Phase Envelopes

   -  Ideal solution based methods
   -  Constructing the phase envelope
   -  Step selection and stepsize control .
   -  Unusual phase envelopes .
   -  Phase diagrams for binary mixtures

      -  General conditions .
      -  Properties of the solutions
      -  Thre&-phase equilibrium
      -  Binary equilibrium lines
      -  Isolated regions .
      -  Low temperature LLE .

   -  Numerical example with Jupyter

6. Phase diagrams for ternary an cuaternary mixtures

   -  Ternary example.
   -  Cuaternary example

7. Equilibrium Solid-Fluid

   -  Modeling of pure solid
   -  Fugacity of Pure Solid

      -  Model I
      -  Model II
      -  Model III

   -  Equilibrium Solid-liquid

      -  Binary systems

8. Equilibrio Químico

   -  Constante de equilibrio químico
   -  Le Chatelier's principle

      -  Temperature
      -  Composition
      -  Pressure
      -  Inert molecules

   -  Reacción química simple

      -  Compotamiento ideal
      -  Modelos termodinámico de Actividad
      -  Modelos de Fugacidad: Ecuaciones de estado

   -  Reacciones químicas multiples
   -  Principles of Separation & Reaction

      -  System controlled by the kinetic and system ontrolled by the
         equilibrium
      -  Reactive PT-Flash calculation
      -  Thermodynamic Equilibrium for the Esterification reactions of
         Carboxylic Acid

-  **Parte 3: Termodinámica Propiedades de exceso**
***************************************************

1. Propiedades de exceso

   -  Funciones de exceso y coeficientes de actividad.
   -  Expresiones empíricas para funciones exceso (Porter-Margules,
      Redlich-Kister, etc.).
   -  Modelo de van Laar y ecuación de van der Waals.
   -  Teoría de soluciones regulares.
   -  Determinación de equilibrio v-l por el método de Barker.
   -  Introducción a la teoría quasiquímica.
   -  Composiciones locales. Método de Scott.
   -  Modelos de Wilson, NRTL y UNIQUAC.
   -  Concepto de solución de grupos.
   -  Modelo UNIFAC.
   -  Cálculo de equilibrio v-l con funciones de exceso.

-  **Parte 4: Thermodynamics and Data**
***************************************

1. Base de datos de propiedades termodinámicas

   -  DIPPR data thermodynamics
   -  Properties thermodynamics
   -  Pandas for Data thermodynamics

2. Analisis de consistencia termodinamica

   -  Sistemas a baja presión
   -  Sistemas con alta presión

3. Regression of parameters

   - Modeling objective function
   - Regression of parameters with Scipy
   - Casos de de aplicación

-  **Parte 5: Proyecto final**
******************************

1. Evaluación del curso

   - Definición del caso a trabajar
   - Formulación de la solución
   - Implementación con Jupyter
   - Presentación de resultados y conclusiones


- **Información**
-----------------

**Horario de Clases:** Viernes de 2 pm a 4 pm en el anfiteatro C de la FCEFyN.

**Modalidad:** Guias de Problemas - Prácticas con computación.

**Evaluación:** Los alumnos rendirán un examen final y una sustentación de un proyecto.

**Profesores:** Martín Cismondi. 