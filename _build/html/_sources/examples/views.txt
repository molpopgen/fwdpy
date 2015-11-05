
Example of taking 'views' from simulated populations
====================================================

.. code:: python

    from __future__ import print_function
    import fwdpy as fp
    import pandas as pd
    from background_selection_setup import *

Get the mutations that are segregating in each population:

.. code:: python

    mutations = [fp.view_mutations(i) for i in pops]

Look at the raw data in the first element of each list:

.. code:: python

    for i in mutations:
        print(i[0])


.. parsed-literal::

    {'g': 9928, 'h': 1.0, 'neutral': False, 'pos': -0.9881998631171882, 's': -0.05000000074505806, 'n': 1}
    {'g': 9990, 'h': 1.0, 'neutral': False, 'pos': -0.9669022976886481, 's': -0.05000000074505806, 'n': 1}
    {'g': 9997, 'h': 1.0, 'neutral': False, 'pos': -0.981396772665903, 's': -0.05000000074505806, 'n': 1}
    {'g': 9976, 'h': 1.0, 'neutral': False, 'pos': -0.9961579141672701, 's': -0.05000000074505806, 'n': 29}


Let's make that nicer, and convert each list of dictionaries to a Pandas
DataFrame object:

.. code:: python

    mutations2 = [pd.DataFrame(i) for i in mutations]

.. code:: python

    for i in mutations2:
        print(i.head())


.. parsed-literal::

          g  h  n neutral       pos     s
    0  9928  1  1   False -0.988200 -0.05
    1  9996  1  2   False -0.945347 -0.05
    2  9982  1  4   False -0.936773 -0.05
    3  9986  1  5   False -0.926391 -0.05
    4  9997  1  3   False -0.912301 -0.05
          g  h   n neutral       pos     s
    0  9990  1   1   False -0.966902 -0.05
    1  9981  1   5   False -0.928520 -0.05
    2  9993  1   2   False -0.905435 -0.05
    3  9978  1   2   False -0.876393 -0.05
    4  9994  1  11   False -0.861297 -0.05
          g  h  n neutral       pos     s
    0  9997  1  1   False -0.981397 -0.05
    1  9999  1  1   False -0.961600 -0.05
    2  9998  1  2   False -0.958168 -0.05
    3  9985  1  1   False -0.860568 -0.05
    4  9992  1  3   False -0.852956 -0.05
          g  h   n neutral       pos     s
    0  9976  1  29   False -0.996158 -0.05
    1  9992  1   1   False -0.985585 -0.05
    2  9985  1   3   False -0.958816 -0.05
    3  9999  1   1   False -0.924767 -0.05
    4  9992  1   8   False -0.878102 -0.05


The columns are:

-  g = the generation when the mutation first arose
-  h = the dominance
-  n = the number of copies of the mutation in the population. You can
   use this to get its frequency.
-  neutral = a boolean
-  pos = the position of the mutation
-  s = selection coefficient/effect size

We can do all the usual subsetting, etc., using regular pandas tricks.
For example, let's get the neutral mutations for each population:

.. code:: python

    nmuts = [i[i.neutral == True] for i in mutations2]
    for i in nmuts:
        print(i.head())


.. parsed-literal::

           g  h    n neutral       pos  s
    42  9996  0    2    True  0.001781  0
    43  9636  0   37    True  0.007266  0
    44  9763  0  539    True  0.011170  0
    45  9937  0   23    True  0.011725  0
    46  9999  0    1    True  0.031578  0
           g  h    n neutral       pos  s
    56  9425  0  132    True  0.004553  0
    57  8659  0  663    True  0.011580  0
    58  9985  0    3    True  0.013945  0
    59  9979  0    6    True  0.025695  0
    60  9954  0   53    True  0.038241  0
           g  h     n neutral       pos  s
    51  6991  0  1872    True  0.004533  0
    52  5616  0   128    True  0.010608  0
    53  9994  0     1    True  0.012792  0
    54  9991  0    13    True  0.021528  0
    55  9999  0     1    True  0.021921  0
           g  h     n neutral       pos  s
    52  9973  0    28    True  0.000875  0
    53  9561  0    66    True  0.002232  0
    54  8814  0  1658    True  0.006304  0
    55  8590  0  1658    True  0.017952  0
    56  8317  0   342    True  0.018824  0


We can also take views of gametes:

.. code:: python

    gametes = [fp.view_gametes(i) for i in pops]

The format is really ugly. v Each gamete is a dict with two elements:

-  'neutral' is a list of mutations *not* affecting fitness. The format
   is the same as for the mutation views above.
-  'selected' is a list of mutations that *do* affect fitness. The
   format is the same as for the mutation views above.

.. code:: python

    for i in gametes:
        print(i[0])


.. parsed-literal::

    {'neutral': [{'g': 9763, 'h': 0.0, 'neutral': True, 'pos': 0.01117014978080988, 's': 0.0, 'n': 539}, {'g': 7495, 'h': 0.0, 'neutral': True, 'pos': 0.04206915246322751, 's': 0.0, 'n': 1491}, {'g': 8137, 'h': 0.0, 'neutral': True, 'pos': 0.05866316333413124, 's': 0.0, 'n': 1491}, {'g': 7999, 'h': 0.0, 'neutral': True, 'pos': 0.07763492036610842, 's': 0.0, 'n': 1444}, {'g': 8444, 'h': 0.0, 'neutral': True, 'pos': 0.11814258061349392, 's': 0.0, 'n': 1407}, {'g': 9838, 'h': 0.0, 'neutral': True, 'pos': 0.12312695779837668, 's': 0.0, 'n': 225}, {'g': 8830, 'h': 0.0, 'neutral': True, 'pos': 0.21917402558028698, 's': 0.0, 'n': 1420}, {'g': 7304, 'h': 0.0, 'neutral': True, 'pos': 0.255468935938552, 's': 0.0, 'n': 1474}, {'g': 7907, 'h': 0.0, 'neutral': True, 'pos': 0.2975613162852824, 's': 0.0, 'n': 1474}, {'g': 8598, 'h': 0.0, 'neutral': True, 'pos': 0.3026043844874948, 's': 0.0, 'n': 1474}, {'g': 8955, 'h': 0.0, 'neutral': True, 'pos': 0.303952107205987, 's': 0.0, 'n': 1437}, {'g': 9796, 'h': 0.0, 'neutral': True, 'pos': 0.5242417093832046, 's': 0.0, 'n': 223}, {'g': 9361, 'h': 0.0, 'neutral': True, 'pos': 0.6134322723373771, 's': 0.0, 'n': 1351}, {'g': 7541, 'h': 0.0, 'neutral': True, 'pos': 0.6526365035679191, 's': 0.0, 'n': 1605}, {'g': 7170, 'h': 0.0, 'neutral': True, 'pos': 0.6547259239014238, 's': 0.0, 'n': 1605}, {'g': 8355, 'h': 0.0, 'neutral': True, 'pos': 0.6865988140925765, 's': 0.0, 'n': 1472}, {'g': 8185, 'h': 0.0, 'neutral': True, 'pos': 0.704608783358708, 's': 0.0, 'n': 1472}, {'g': 9747, 'h': 0.0, 'neutral': True, 'pos': 0.7162117944099009, 's': 0.0, 'n': 585}, {'g': 9848, 'h': 0.0, 'neutral': True, 'pos': 0.7177490566391498, 's': 0.0, 'n': 222}, {'g': 8584, 'h': 0.0, 'neutral': True, 'pos': 0.8198388759046793, 's': 0.0, 'n': 1490}, {'g': 7656, 'h': 0.0, 'neutral': True, 'pos': 0.828043204266578, 's': 0.0, 'n': 1554}, {'g': 9702, 'h': 0.0, 'neutral': True, 'pos': 0.8485063153784722, 's': 0.0, 'n': 610}, {'g': 6614, 'h': 0.0, 'neutral': True, 'pos': 0.8667287053540349, 's': 0.0, 'n': 1445}, {'g': 2480, 'h': 0.0, 'neutral': True, 'pos': 0.9092434337362647, 's': 0.0, 'n': 1885}, {'g': 2970, 'h': 0.0, 'neutral': True, 'pos': 0.9112986992113292, 's': 0.0, 'n': 1885}, {'g': 3901, 'h': 0.0, 'neutral': True, 'pos': 0.9496243973262608, 's': 0.0, 'n': 1885}], 'selected': [], 'n': 120}
    {'neutral': [{'g': 8659, 'h': 0.0, 'neutral': True, 'pos': 0.011580322636291385, 's': 0.0, 'n': 663}, {'g': 9423, 'h': 0.0, 'neutral': True, 'pos': 0.17388591263443232, 's': 0.0, 'n': 253}, {'g': 8293, 'h': 0.0, 'neutral': True, 'pos': 0.26215685391798615, 's': 0.0, 'n': 573}, {'g': 9157, 'h': 0.0, 'neutral': True, 'pos': 0.28463278082199395, 's': 0.0, 'n': 253}, {'g': 9181, 'h': 0.0, 'neutral': True, 'pos': 0.3756080197636038, 's': 0.0, 'n': 264}, {'g': 9935, 'h': 0.0, 'neutral': True, 'pos': 0.40037358715198934, 's': 0.0, 'n': 177}, {'g': 8526, 'h': 0.0, 'neutral': True, 'pos': 0.4443175867199898, 's': 0.0, 'n': 574}, {'g': 9494, 'h': 0.0, 'neutral': True, 'pos': 0.508453646209091, 's': 0.0, 'n': 253}, {'g': 8460, 'h': 0.0, 'neutral': True, 'pos': 0.5480429537128657, 's': 0.0, 'n': 574}, {'g': 8230, 'h': 0.0, 'neutral': True, 'pos': 0.5577912081498653, 's': 0.0, 'n': 574}, {'g': 9412, 'h': 0.0, 'neutral': True, 'pos': 0.6276867836713791, 's': 0.0, 'n': 253}, {'g': 8506, 'h': 0.0, 'neutral': True, 'pos': 0.6890833463985473, 's': 0.0, 'n': 562}, {'g': 7736, 'h': 0.0, 'neutral': True, 'pos': 0.878798944177106, 's': 0.0, 'n': 946}, {'g': 7836, 'h': 0.0, 'neutral': True, 'pos': 0.9040587430354208, 's': 0.0, 'n': 946}, {'g': 9360, 'h': 0.0, 'neutral': True, 'pos': 0.9324489531572908, 's': 0.0, 'n': 255}, {'g': 7694, 'h': 0.0, 'neutral': True, 'pos': 0.9534462296869606, 's': 0.0, 'n': 940}], 'selected': [], 'n': 138}
    {'neutral': [{'g': 6991, 'h': 0.0, 'neutral': True, 'pos': 0.004532672232016921, 's': 0.0, 'n': 1872}, {'g': 7371, 'h': 0.0, 'neutral': True, 'pos': 0.09877051902003586, 's': 0.0, 'n': 1853}, {'g': 9507, 'h': 0.0, 'neutral': True, 'pos': 0.13153062527999282, 's': 0.0, 'n': 517}, {'g': 8609, 'h': 0.0, 'neutral': True, 'pos': 0.1337847807444632, 's': 0.0, 'n': 1786}, {'g': 9521, 'h': 0.0, 'neutral': True, 'pos': 0.16495748492889106, 's': 0.0, 'n': 517}, {'g': 9299, 'h': 0.0, 'neutral': True, 'pos': 0.2192109227180481, 's': 0.0, 'n': 713}, {'g': 7001, 'h': 0.0, 'neutral': True, 'pos': 0.2334199717734009, 's': 0.0, 'n': 1786}, {'g': 7767, 'h': 0.0, 'neutral': True, 'pos': 0.26152393594384193, 's': 0.0, 'n': 1786}, {'g': 7895, 'h': 0.0, 'neutral': True, 'pos': 0.28069930942729115, 's': 0.0, 'n': 1786}, {'g': 7034, 'h': 0.0, 'neutral': True, 'pos': 0.28539203526452184, 's': 0.0, 'n': 1786}, {'g': 8290, 'h': 0.0, 'neutral': True, 'pos': 0.3596646406222135, 's': 0.0, 'n': 1786}, {'g': 9580, 'h': 0.0, 'neutral': True, 'pos': 0.4241732864174992, 's': 0.0, 'n': 514}, {'g': 9362, 'h': 0.0, 'neutral': True, 'pos': 0.4602559239137918, 's': 0.0, 'n': 710}, {'g': 5261, 'h': 0.0, 'neutral': True, 'pos': 0.4855107609182596, 's': 0.0, 'n': 1581}, {'g': 9587, 'h': 0.0, 'neutral': True, 'pos': 0.5322959397453815, 's': 0.0, 'n': 514}, {'g': 8881, 'h': 0.0, 'neutral': True, 'pos': 0.5442884352523834, 's': 0.0, 'n': 1581}, {'g': 8226, 'h': 0.0, 'neutral': True, 'pos': 0.5709147118031979, 's': 0.0, 'n': 1581}, {'g': 9288, 'h': 0.0, 'neutral': True, 'pos': 0.6156564524862915, 's': 0.0, 'n': 710}, {'g': 8683, 'h': 0.0, 'neutral': True, 'pos': 0.6614894329104573, 's': 0.0, 'n': 1570}, {'g': 6045, 'h': 0.0, 'neutral': True, 'pos': 0.7605923339724541, 's': 0.0, 'n': 1342}, {'g': 8117, 'h': 0.0, 'neutral': True, 'pos': 0.8810987235046923, 's': 0.0, 'n': 1342}], 'selected': [], 'n': 185}
    {'neutral': [{'g': 8814, 'h': 0.0, 'neutral': True, 'pos': 0.006303782807663083, 's': 0.0, 'n': 1658}, {'g': 8590, 'h': 0.0, 'neutral': True, 'pos': 0.017951839603483677, 's': 0.0, 'n': 1658}, {'g': 8824, 'h': 0.0, 'neutral': True, 'pos': 0.04905842989683151, 's': 0.0, 'n': 1658}, {'g': 9660, 'h': 0.0, 'neutral': True, 'pos': 0.14721597940661013, 's': 0.0, 'n': 505}, {'g': 7888, 'h': 0.0, 'neutral': True, 'pos': 0.22597101191058755, 's': 0.0, 'n': 1700}, {'g': 9193, 'h': 0.0, 'neutral': True, 'pos': 0.2749505708925426, 's': 0.0, 'n': 557}, {'g': 7714, 'h': 0.0, 'neutral': True, 'pos': 0.2756908859591931, 's': 0.0, 'n': 1689}, {'g': 7682, 'h': 0.0, 'neutral': True, 'pos': 0.5999884365592152, 's': 0.0, 'n': 1683}, {'g': 6908, 'h': 0.0, 'neutral': True, 'pos': 0.6519452347420156, 's': 0.0, 'n': 1683}, {'g': 7858, 'h': 0.0, 'neutral': True, 'pos': 0.8093185543548316, 's': 0.0, 'n': 1777}, {'g': 8205, 'h': 0.0, 'neutral': True, 'pos': 0.8573803363833576, 's': 0.0, 'n': 1683}], 'selected': [], 'n': 187}


OK, let's clean that up. We'll focus on the selected mutations for each
individual, and turn everything into a pd.DataFrame.

We're only going to do this for the first simulated population.

.. code:: python

    smuts = [i['selected'] for i in gametes[0]]

We now have a list of lists stored in 'smuts'.

.. code:: python

    smutsdf = pd.DataFrame()
    ind=0
    ##Add the non-empty individuals to the df
    for i in smuts:
        if len(i)>0:
            smutsdf = pd.concat([smutsdf,pd.DataFrame(i,index=[ind]*len(i))])
        ind += 1

.. code:: python

    smutsdf.head()




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>g</th>
          <th>h</th>
          <th>n</th>
          <th>neutral</th>
          <th>pos</th>
          <th>s</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>25</th>
          <td>9914</td>
          <td>1</td>
          <td>22</td>
          <td>False</td>
          <td>1.984707</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>31</th>
          <td>9974</td>
          <td>1</td>
          <td>18</td>
          <td>False</td>
          <td>1.939535</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>34</th>
          <td>9933</td>
          <td>1</td>
          <td>12</td>
          <td>False</td>
          <td>1.877543</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>37</th>
          <td>9989</td>
          <td>1</td>
          <td>10</td>
          <td>False</td>
          <td>-0.699703</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>37</th>
          <td>9972</td>
          <td>1</td>
          <td>11</td>
          <td>False</td>
          <td>1.890960</td>
          <td>-0.05</td>
        </tr>
      </tbody>
    </table>
    </div>



That's much better. We can use the index to figure out which individual
has which mutations, and their effect sizes, etc.

Finally, we can also take views of diploids. Let's get the first two
diploids in each population:

.. code:: python

    dips = [fp.view_diploids(i,[0,1]) for i in pops]

Again, the format here is ugly. Each diploid view is a dictionary:

.. code:: python

    for key in dips[0][0]:
        print(key)


.. parsed-literal::

    chrom1
    e
    w
    g
    chrom0


The values are:

-  chrom0, chrom1 are gamete views, just like what we dealt with above
-  g = genetic component of phenotype
-  e = random component of phenotype
-  w = fitness

Please note that g, e, and w, may or may not be set by a particular
simulation. Their use is optional.
