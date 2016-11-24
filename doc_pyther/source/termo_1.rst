**************************************
15. Termo (T, P)
**************************************



.. code:: python

    import scipy as sp
    from scipy import optimize
    from scipy.optimize import fsolve
    import numpy as np
    from matplotlib import pyplot
    %matplotlib inline
    import pandas as pd
    from numpy import linalg as LA
    from IPython.html import widgets
    from IPython.display import display
    from IPython.display import clear_output
    
    # encoding: utf-8
    
    from pandas import read_csv



.. parsed-literal::

    /home/andres.python/anaconda3/lib/python3.4/site-packages/IPython/html.py:14: ShimWarning: The `IPython.html` package has been deprecated. You should import from `notebook` instead. `IPython.html.widgets` has moved to `ipywidgets`.
      "`IPython.html.widgets` has moved to `ipywidgets`.", ShimWarning)


.. code:: python

    
    f = pd.read_excel("PureFull.xls")
    f.head()
    data2 = pd.DataFrame(f)
    data2 = data2.set_index('Name')
    data2 = data2.ix[:, 1:12]
    Etiquetas = data2.index.get_values()
    #Etiquetas.shape
    
    print(type(Etiquetas))


.. parsed-literal::

    <class 'numpy.ndarray'>


.. code:: python

    Componentes_1 = widgets.SelectMultiple(description="Component 1", options=list(Etiquetas))
    display(Componentes_1)  

.. code:: python

    titulo = widgets.HTML(value="<C><H1> Welcome to PyTher <H1>")
    titulo

.. code:: python

    # have to comment out the magic IPython functions because of a bug
    #%pylab inline
    import numpy as np
    import pandas as pd

.. code:: python

    print (np.random.randn(20))


.. parsed-literal::

    [-1.84827553  0.19748324 -0.84097815 -0.76943984 -0.59399768 -1.16957456
     -0.20641787  0.26595243 -0.69885132  1.40674578  0.82962834  1.35402756
     -3.0501609  -0.59944992  0.44388827 -0.03891357 -1.99067461 -0.3815753
      0.41178857  0.62843074]


.. code:: python

    pd.Series(np.random.randn(100)).plot()




.. parsed-literal::

    <matplotlib.axes._subplots.AxesSubplot at 0x7f781bb606d8>

.. image:: programando.jpg


