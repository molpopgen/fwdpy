
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

    {'g': 9999, 'h': 1.0, 'neutral': False, 'pos': 1.2961536727380008, 's': -0.05000000074505806, 'label': 3, 'n': 1}
    {'g': 9974, 'h': 0.0, 'neutral': True, 'pos': 0.3750512143597007, 's': 0.0, 'label': 1, 'n': 17}
    {'g': 9999, 'h': 0.0, 'neutral': True, 'pos': 0.802856310736388, 's': 0.0, 'label': 1, 'n': 1}
    {'g': 9832, 'h': 0.0, 'neutral': True, 'pos': 0.3733968627639115, 's': 0.0, 'label': 1, 'n': 68}


Let's make that nicer, and convert each list of dictionaries to a Pandas
DataFrame object:

.. code:: python

    mutations2 = [pd.DataFrame(i) for i in mutations]

.. code:: python

    for i in mutations2:
        print(i.head())


.. parsed-literal::

          g    h  label   n neutral       pos     s
    0  9999  1.0      3   1   False  1.296154 -0.05
    1  9999  0.0      1   1    True  0.808702  0.00
    2  9985  0.0      1  20    True  0.851213  0.00
    3  9997  1.0      2   1   False -0.743512 -0.05
    4  9997  1.0      3   3   False  1.256208 -0.05
          g    h  label    n neutral       pos    s
    0  9974  0.0      1   17    True  0.375051  0.0
    1  9954  0.0      1   53    True  0.038241  0.0
    2  9999  0.0      1    1    True  0.922048  0.0
    3  9986  0.0      1   16    True  0.663966  0.0
    4  9637  0.0      1  148    True  0.516812  0.0
          g    h  label   n neutral       pos     s
    0  9999  0.0      1   1    True  0.802856  0.00
    1  9987  1.0      3   2   False  1.719877 -0.05
    2  9993  0.0      1  11    True  0.044257  0.00
    3  9999  0.0      1   1    True  0.319559  0.00
    4  9961  1.0      2  19   False -0.092556 -0.05
          g    h  label   n neutral       pos     s
    0  9832  0.0      1  68    True  0.373397  0.00
    1  9995  1.0      2   2   False -0.317365 -0.05
    2  9998  1.0      3   1   False  1.072465 -0.05
    3  9868  0.0      1  12    True  0.045664  0.00
    4  9749  0.0      1  42    True  0.644652  0.00


The columns are:

-  g = the generation when the mutation first arose
-  h = the dominance
-  n = the number of copies of the mutation in the population. You can
   use this to get its frequency.
-  neutral = a boolean
-  pos = the position of the mutation
-  s = selection coefficient/effect size
-  label = The label assigned to a mutation. These labels can be
   associated with Regions and Sregions. Here, 1 is a mutation from the
   neutral region, 2 a selected mutation from the 'left' region and 3 a
   selected mutation from the 'right' regin.

We can do all the usual subsetting, etc., using regular pandas tricks.
For example, let's get the neutral mutations for each population:

.. code:: python

    nmuts = [i[i.neutral == True] for i in mutations2]
    for i in nmuts:
        print(i.head())


.. parsed-literal::

          g    h  label    n neutral       pos    s
    1  9999  0.0      1    1    True  0.808702  0.0
    2  9985  0.0      1   20    True  0.851213  0.0
    5  9916  0.0      1   14    True  0.048614  0.0
    6  7673  0.0      1  478    True  0.442744  0.0
    8  9925  0.0      1    6    True  0.238877  0.0
          g    h  label    n neutral       pos    s
    0  9974  0.0      1   17    True  0.375051  0.0
    1  9954  0.0      1   53    True  0.038241  0.0
    2  9999  0.0      1    1    True  0.922048  0.0
    3  9986  0.0      1   16    True  0.663966  0.0
    4  9637  0.0      1  148    True  0.516812  0.0
          g    h  label     n neutral       pos    s
    0  9999  0.0      1     1    True  0.802856  0.0
    2  9993  0.0      1    11    True  0.044257  0.0
    3  9999  0.0      1     1    True  0.319559  0.0
    5  8486  0.0      1   147    True  0.087042  0.0
    7  8290  0.0      1  1786    True  0.359665  0.0
          g    h  label   n neutral       pos    s
    0  9832  0.0      1  68    True  0.373397  0.0
    3  9868  0.0      1  12    True  0.045664  0.0
    4  9749  0.0      1  42    True  0.644652  0.0
    5  9968  0.0      1  20    True  0.316820  0.0
    6  9970  0.0      1   7    True  0.605128  0.0


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

    {'neutral': [{'g': 9999, 'h': 1.0, 'neutral': False, 'pos': 1.2961536727380008, 's': -0.05000000074505806, 'label': 3, 'n': 1}, {'g': 9999, 'h': 0.0, 'neutral': True, 'pos': 0.8087020928505808, 's': 0.0, 'label': 1, 'n': 1}, {'g': 9985, 'h': 0.0, 'neutral': True, 'pos': 0.8512128263246268, 's': 0.0, 'label': 1, 'n': 20}, {'g': 9997, 'h': 1.0, 'neutral': False, 'pos': -0.7435115049593151, 's': -0.05000000074505806, 'label': 2, 'n': 1}, {'g': 9997, 'h': 1.0, 'neutral': False, 'pos': 1.2562080402858555, 's': -0.05000000074505806, 'label': 3, 'n': 3}, {'g': 9916, 'h': 0.0, 'neutral': True, 'pos': 0.04861398716457188, 's': 0.0, 'label': 1, 'n': 14}, {'g': 7673, 'h': 0.0, 'neutral': True, 'pos': 0.4427443742752075, 's': 0.0, 'label': 1, 'n': 478}, {'g': 9942, 'h': 1.0, 'neutral': False, 'pos': -0.6328390480484813, 's': -0.05000000074505806, 'label': 2, 'n': 9}, {'g': 9925, 'h': 0.0, 'neutral': True, 'pos': 0.23887674184516072, 's': 0.0, 'label': 1, 'n': 6}, {'g': 9992, 'h': 1.0, 'neutral': False, 'pos': 1.7460428576450795, 's': -0.05000000074505806, 'label': 3, 'n': 3}, {'g': 9951, 'h': 0.0, 'neutral': True, 'pos': 0.9244193523190916, 's': 0.0, 'label': 1, 'n': 9}, {'g': 9998, 'h': 0.0, 'neutral': True, 'pos': 0.40817153407260776, 's': 0.0, 'label': 1, 'n': 0}, {'g': 3901, 'h': 0.0, 'neutral': True, 'pos': 0.9496243973262608, 's': 0.0, 'label': 1, 'n': 1885}, {'g': 9987, 'h': 1.0, 'neutral': False, 'pos': -0.8368993068579584, 's': -0.05000000074505806, 'label': 2, 'n': 5}, {'g': 7838, 'h': 0.0, 'neutral': True, 'pos': 0.4856802138965577, 's': 0.0, 'label': 1, 'n': 479}, {'g': 9999, 'h': 1.0, 'neutral': False, 'pos': -0.08335124771110713, 's': -0.05000000074505806, 'label': 2, 'n': 1}, {'g': 9993, 'h': 1.0, 'neutral': False, 'pos': 1.1813256666064262, 's': -0.05000000074505806, 'label': 3, 'n': 8}, {'g': 9999, 'h': 0.0, 'neutral': True, 'pos': 0.9340733217541128, 's': 0.0, 'label': 1, 'n': 1}, {'g': 8617, 'h': 0.0, 'neutral': True, 'pos': 0.655554321128875, 's': 0.0, 'label': 1, 'n': 121}, {'g': 9997, 'h': 1.0, 'neutral': False, 'pos': -0.04492441425099969, 's': -0.05000000074505806, 'label': 2, 'n': 0}, {'g': 9994, 'h': 1.0, 'neutral': False, 'pos': 1.5325499204918742, 's': -0.05000000074505806, 'label': 3, 'n': 2}], 'selected': [{'g': 9999, 'h': 1.0, 'neutral': False, 'pos': 1.2961536727380008, 's': -0.05000000074505806, 'label': 3, 'n': 1}], 'n': 1}
    {'neutral': [{'g': 9974, 'h': 0.0, 'neutral': True, 'pos': 0.3750512143597007, 's': 0.0, 'label': 1, 'n': 17}, {'g': 9954, 'h': 0.0, 'neutral': True, 'pos': 0.03824054426513612, 's': 0.0, 'label': 1, 'n': 53}, {'g': 9999, 'h': 0.0, 'neutral': True, 'pos': 0.9220477596390992, 's': 0.0, 'label': 1, 'n': 1}, {'g': 9986, 'h': 0.0, 'neutral': True, 'pos': 0.663966491818428, 's': 0.0, 'label': 1, 'n': 16}, {'g': 9637, 'h': 0.0, 'neutral': True, 'pos': 0.5168115310370922, 's': 0.0, 'label': 1, 'n': 148}, {'g': 9983, 'h': 1.0, 'neutral': False, 'pos': -0.843893475830555, 's': -0.05000000074505806, 'label': 2, 'n': 24}, {'g': 9989, 'h': 0.0, 'neutral': True, 'pos': 0.290936284000054, 's': 0.0, 'label': 1, 'n': 8}, {'g': 9996, 'h': 0.0, 'neutral': True, 'pos': 0.32976379804313183, 's': 0.0, 'label': 1, 'n': 0}, {'g': 9976, 'h': 1.0, 'neutral': False, 'pos': -0.5398803392890841, 's': -0.05000000074505806, 'label': 2, 'n': 0}, {'g': 9992, 'h': 0.0, 'neutral': True, 'pos': 0.6648217977490276, 's': 0.0, 'label': 1, 'n': 6}, {'g': 9998, 'h': 1.0, 'neutral': False, 'pos': 1.7505667991936207, 's': -0.05000000074505806, 'label': 3, 'n': 0}, {'g': 9999, 'h': 0.0, 'neutral': True, 'pos': 0.27948754257522523, 's': 0.0, 'label': 1, 'n': 1}, {'g': 9981, 'h': 1.0, 'neutral': False, 'pos': -0.9285198170691729, 's': -0.05000000074505806, 'label': 2, 'n': 5}, {'g': 8751, 'h': 0.0, 'neutral': True, 'pos': 0.5460440656170249, 's': 0.0, 'label': 1, 'n': 216}, {'g': 9986, 'h': 1.0, 'neutral': False, 'pos': 1.4597016391344368, 's': -0.05000000074505806, 'label': 3, 'n': 8}], 'selected': [], 'n': 15}
    {'neutral': [{'g': 9999, 'h': 0.0, 'neutral': True, 'pos': 0.802856310736388, 's': 0.0, 'label': 1, 'n': 1}, {'g': 9987, 'h': 1.0, 'neutral': False, 'pos': 1.7198766963556409, 's': -0.05000000074505806, 'label': 3, 'n': 2}, {'g': 9993, 'h': 0.0, 'neutral': True, 'pos': 0.044257442466914654, 's': 0.0, 'label': 1, 'n': 11}, {'g': 9999, 'h': 0.0, 'neutral': True, 'pos': 0.31955947587266564, 's': 0.0, 'label': 1, 'n': 1}, {'g': 9961, 'h': 1.0, 'neutral': False, 'pos': -0.09255639626644552, 's': -0.05000000074505806, 'label': 2, 'n': 19}, {'g': 8486, 'h': 0.0, 'neutral': True, 'pos': 0.08704224601387978, 's': 0.0, 'label': 1, 'n': 147}, {'g': 9988, 'h': 1.0, 'neutral': False, 'pos': -0.6743454595562071, 's': -0.05000000074505806, 'label': 2, 'n': 1}, {'g': 8290, 'h': 0.0, 'neutral': True, 'pos': 0.3596646406222135, 's': 0.0, 'label': 1, 'n': 1786}, {'g': 9996, 'h': 1.0, 'neutral': False, 'pos': 1.1176744250115007, 's': -0.05000000074505806, 'label': 3, 'n': 0}, {'g': 9999, 'h': 1.0, 'neutral': False, 'pos': -0.6342870420776308, 's': -0.05000000074505806, 'label': 2, 'n': 1}, {'g': 8881, 'h': 0.0, 'neutral': True, 'pos': 0.5442884352523834, 's': 0.0, 'label': 1, 'n': 1581}, {'g': 9993, 'h': 0.0, 'neutral': True, 'pos': 0.9674568464979529, 's': 0.0, 'label': 1, 'n': 4}, {'g': 9999, 'h': 1.0, 'neutral': False, 'pos': -0.8462745684664696, 's': -0.05000000074505806, 'label': 2, 'n': 1}, {'g': 9999, 'h': 1.0, 'neutral': False, 'pos': -0.11604839516803622, 's': -0.05000000074505806, 'label': 2, 'n': 1}, {'g': 9997, 'h': 1.0, 'neutral': False, 'pos': -0.45378030952997506, 's': -0.05000000074505806, 'label': 2, 'n': 0}, {'g': 9992, 'h': 0.0, 'neutral': True, 'pos': 0.16772447689436376, 's': 0.0, 'label': 1, 'n': 5}, {'g': 9884, 'h': 0.0, 'neutral': True, 'pos': 0.2799046675208956, 's': 0.0, 'label': 1, 'n': 17}, {'g': 9864, 'h': 0.0, 'neutral': True, 'pos': 0.9488759697414935, 's': 0.0, 'label': 1, 'n': 36}, {'g': 9998, 'h': 0.0, 'neutral': True, 'pos': 0.48674530582502484, 's': 0.0, 'label': 1, 'n': 1}, {'g': 9995, 'h': 1.0, 'neutral': False, 'pos': 1.7994547698181123, 's': -0.05000000074505806, 'label': 3, 'n': 4}], 'selected': [], 'n': 19}
    {'neutral': [{'g': 9832, 'h': 0.0, 'neutral': True, 'pos': 0.3733968627639115, 's': 0.0, 'label': 1, 'n': 68}, {'g': 9995, 'h': 1.0, 'neutral': False, 'pos': -0.31736462796106935, 's': -0.05000000074505806, 'label': 2, 'n': 2}, {'g': 9998, 'h': 1.0, 'neutral': False, 'pos': 1.072465442121029, 's': -0.05000000074505806, 'label': 3, 'n': 1}, {'g': 9868, 'h': 0.0, 'neutral': True, 'pos': 0.04566416423767805, 's': 0.0, 'label': 1, 'n': 12}, {'g': 9749, 'h': 0.0, 'neutral': True, 'pos': 0.6446520905010402, 's': 0.0, 'label': 1, 'n': 42}, {'g': 9968, 'h': 0.0, 'neutral': True, 'pos': 0.31681952835060656, 's': 0.0, 'label': 1, 'n': 20}, {'g': 9970, 'h': 0.0, 'neutral': True, 'pos': 0.6051280729006976, 's': 0.0, 'label': 1, 'n': 7}, {'g': 9992, 'h': 1.0, 'neutral': False, 'pos': 1.4546588633675128, 's': -0.05000000074505806, 'label': 3, 'n': 2}, {'g': 9998, 'h': 0.0, 'neutral': True, 'pos': 0.2997995924670249, 's': 0.0, 'label': 1, 'n': 2}, {'g': 9997, 'h': 1.0, 'neutral': False, 'pos': 1.260333125013858, 's': -0.05000000074505806, 'label': 3, 'n': 5}, {'g': 9989, 'h': 1.0, 'neutral': False, 'pos': 1.2523613322991878, 's': -0.05000000074505806, 'label': 3, 'n': 9}, {'g': 9992, 'h': 1.0, 'neutral': False, 'pos': 1.2467866179067641, 's': -0.05000000074505806, 'label': 3, 'n': 1}, {'g': 8760, 'h': 0.0, 'neutral': True, 'pos': 0.5242203767411411, 's': 0.0, 'label': 1, 'n': 63}], 'selected': [], 'n': 11}


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
          <th>label</th>
          <th>n</th>
          <th>neutral</th>
          <th>pos</th>
          <th>s</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>9999</td>
          <td>1.0</td>
          <td>3</td>
          <td>1</td>
          <td>False</td>
          <td>1.296154</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>3</th>
          <td>9999</td>
          <td>1.0</td>
          <td>3</td>
          <td>1</td>
          <td>False</td>
          <td>1.296154</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>6</th>
          <td>9999</td>
          <td>1.0</td>
          <td>3</td>
          <td>1</td>
          <td>False</td>
          <td>1.296154</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>8</th>
          <td>9999</td>
          <td>1.0</td>
          <td>3</td>
          <td>1</td>
          <td>False</td>
          <td>1.296154</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>10</th>
          <td>9999</td>
          <td>1.0</td>
          <td>3</td>
          <td>1</td>
          <td>False</td>
          <td>1.296154</td>
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

    sh0
    sh1
    e
    g
    w
    n0
    n1
    chrom1
    chrom0


The values are:

-  chrom0, chrom1 are gamete views, just like what we dealt with above
-  g = genetic component of phenotype
-  e = random component of phenotype
-  w = fitness
-  n0 and n1 are the number of selected variants on chrom0 and chrom1,
   respectively.
-  sh0 and sh1 are the sum of :math:`s \times h` for all selected
   mutations on chrom0 and chrom1, respectively

Please note that g, e, and w, may or may not be set by a particular
simulation. Their use is optional.
