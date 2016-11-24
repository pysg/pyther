

http://thomas-cokelaer.info/tutorials/sphinx/docstring_python.html
http://rest-sphinx-memo.readthedocs.org/en/latest/ReST.html

http://pythonbooks.revolunet.com/

http://interactivepython.org/runestone/static/pythonds/index.html

http://pybonacci.org/2015/01/14/introduccion-a-machine-learning-con-python-parte-1/
https://www.cs.uct.ac.za/mit_notes/Python/Object_Oriented_Programming.html
https://www.toptal.com/python/python-design-patterns


Ecuación ( 5.3) página V-3 

.. math:: f_{puro,2}^s(T,P) = \hat f_2 \left(T,P,y_2 \right)

con la función objetivo 

.. math:: F = f_{puro,2}^s(T,P) - \hat f_2 \left(T,P,y_2 \right)

En este caso se tiene una ecuación con una incognita, que puede ser la :math:`T` cuando se especifica la :math:`P` y viceversa.

Se importa el método **fsolve** de la librería **scipy** como sigue:

**from scipy.optimize import fsolve**

el cual es utilizado para resolver la ecuación de igualdad de fugacidades del componente pesado puro en el sólido y el mismo componente en el fluido calculada en la función **equilibrioSF**, que utiliza 





Ecuación (5.7) página V-6 

.. math:: f_{puro,2}^l(T,P) = \hat f_2 \left(T,y_2=1, v_{2, puro} \right)


.. math:: f(T,P)_{puro,i}^s = f(T,P)_{puro,i}^l exp \left( \frac{\Delta v_i^{s-l}}{RT_{i,f}} \left( \lambda_1 + \lambda_2 + \lambda_3 \right)   \right)

.. math:: \lambda_1 = C_1\left( 1- \frac{T_{i,f}}{T}\right)

.. math:: \lambda_2 = C_2\left( -1- \frac{T_{i,f}}{T}  + ln \left(\frac{T}{T_{i,f}}\right) \right)

.. math:: \lambda_3 = C_3\left( -1 + \frac{T}{2T_{i,f}} + \frac{T_{i,f}}{2T} \right) + \frac{T_{i,f}}{T} \left(P-P_{i,f} \right)

.. math:: P = P_{i,f} + C_1\left( 1- \frac{T}{T_{i,f}}\right) + C_2\left(\frac{T}{T_{i,f}} -1 + \frac{T}{T_{i,f}}ln \left(\frac{T_{i,f}}{T}\right) \right) + C_3\left(\frac{T}{T_{i,f}} - \frac{T^2}{2T_{i,f}^2} - \frac{1}{2}  \right)

.. math:: C_1 = \frac{\Delta h_{i,f}}{\Delta v^{s-l}}

.. math:: C_2 = \frac{AT_{i,f}}{\Delta v^{s-l}}

.. math:: C_3 = \frac{BT_{i,f}^2}{\Delta v^{s-l}}

en este caso el equilibrio sólido-fluido se calcula utilizando una ecuación de estado para describir el fluido y para el caso del sólido se emplea un factor de corrección que multiplica a la fugacidad considerando el fluido puro a las condiciones de T y P especificadas.





Por ahora cumplo con lo que estaba pendiente... Respecto al cálculo de parámetros para RKPR, esto es lo que llegué a encontrar y adapté un poco, partiendo de lo que alguna vez le pasé a Gaitan para implementar en GPEC. 
Se alude al programa CubicParam que yo desarrollé en fortran, y los archivos conparin y conparout.dat. Después puedo pasarte el código, y la idea es que vos desarrolles las mismas posibilidades en Python para IPyTherm.
Alguno de estos días, quizás el miércoles, lo charlamos por si tenés alguna duda.
Saludos

You can choose what to specify, besides Tc, Pc and the acentric factor:

- Zrat = campo con el default 1.168
- RhoLsat at T = campo con el default 0.7*Tc
- RhoLsat at Tr = campo con el default 0.70 (y esta será la opción marcada por default)
- Delta1 parameter
i
Si alguien marca la primera opción, escribirás el CONPARIN.DAT igual que  hasta ahora. Ejemplo:
0  3                                 ICALC,NMODEL
617.7  21.1  0.49233  0.822  1.37    Tc, Pc, omega, Vceos(L/mol)  C10

Si queda marcada la tercera (que de ahora en más será lo "normal") o si el usuario marcó la segunda, el CONPARIN deberá tener esta organización:

3  3                             ICALC,NMODEL
304.21  73.83  0.2236      Tc, Pc, omega
270.0  21.4626                 T(K), RhoLsat (L/mol)

Y el CONPAROUT.DAT que  recibirás será este:

  304.2100   73.8300  0.104682  0.22360     Tc,Pc,Vceos,OM
    3.9809  0.026440  2.509688   2.04173     ac,b,del1,rk
 
Si se elige la cuarta opción, el CONPARIN.DAT será por ejemplo:

2  3  0                             ICALC,NMODEL,IVAP
C64    948.8  3.91     2.723    1.350    ac,b,k,del1

El ejecutable que lee conparin y escribe conparout es "CubicParam.exe".













