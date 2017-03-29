1. Instalación de PyTher
**************************
**************************

# Requisitos

Para realizar la instalación de *PyTher* se requiere tener pre-instalado
*Jupyter Notebook* y *Python*.

Instalación de Jupyter utilizando Anaconda
==========================================

Se recomienda instalar
`*Anaconda* <https://www.continuum.io/downloads>`__ porque de forma
simple no solo instala `*Python* <https://www.python.org/>`__ y
`*Jupyter Notebook* <http://jupyter.org/>`__ sino que también un gran
número de librerías para computación científica.

Pasos de instalación:

1. Descargar Anaconda. Se recomienda la descarga de Anaconda superior a
   Python 3.X.

2. Instalar la versión de Anaconda que descargo, siguiendo las
   instrucciones según sea el
   `caso <https://www.continuum.io/downloads#windows>`__

3. Muy bien, ya se instalo Jupyter Notebook. Ahora vamos a probarlo en
   una `línea de
   comandos <https://es.wikipedia.org/wiki/S%C3%ADmbolo_del_sistema>`__
   y se ejecuta:

jupyter notebook

Luego de tener abierto *Jupyter Notebook* se puede realizar la
instalación de *PyTher* desde una celda del mismo *Jupyter Notebook*
utilizando *PyPi* con la instrucción:

.. code:: python

    !pip install pyther


.. parsed-literal::

    Requirement already satisfied: pyther in ./anaconda3/lib/python3.5/site-packages


*NO olvidar el símbolo **!** inicial*

Luego de instalar PyTher, se puede probar con una importación simple de
la librería con el sigiente ejemplo:

.. code:: python

    import pyther as pt

.. code:: python

    print("PyTher version: ", pt.__version__)


.. parsed-literal::

    PyTher version:  0.6


En este caso se accedió al atributo ***version*** de PyTher para
verificar su correcta instalación.

De esta forma, ya se encuentra disponible la librería PyTher para ser
utilizada con los ejemplos que vienen más adelante.
