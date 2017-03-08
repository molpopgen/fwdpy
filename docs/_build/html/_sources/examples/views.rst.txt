
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

    {'g': 9974, 'h': 1.0, 'pos': 1.2837485196068883, 'n': 8, 's': -0.05000000074505806, 'ftime': 4294967295, 'label': 0, 'neutral': False}
    {'g': 9997, 'h': 1.0, 'pos': -0.8252747415099293, 'n': 4, 's': -0.05000000074505806, 'ftime': 4294967295, 'label': 0, 'neutral': False}
    {'g': 9866, 'h': 0.0, 'pos': 0.12148048472590744, 'n': 3, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}
    {'g': 9998, 'h': 1.0, 'pos': -0.6509648731444031, 'n': 1, 's': -0.05000000074505806, 'ftime': 4294967295, 'label': 0, 'neutral': False}


Let's make that nicer, and convert each list of dictionaries to a Pandas
DataFrame object:

.. code:: python

    mutations2 = [pd.DataFrame(i) for i in mutations]

.. code:: python

    for i in mutations2:
        print(i.head())


.. parsed-literal::

            ftime     g    h  label   n neutral       pos     s
    0  4294967295  9974  1.0      0   8   False  1.283749 -0.05
    1  4294967295  9970  0.0      0  13    True  0.361552  0.00
    2  4294967295  9996  1.0      0   1   False -0.741675 -0.05
    3  4294967295  9996  0.0      0   2    True  0.812076  0.00
    4  4294967295  9909  0.0      0  78    True  0.566514  0.00
            ftime     g    h  label    n neutral       pos     s
    0  4294967295  9997  1.0      0    4   False -0.825275 -0.05
    1  4294967295  9928  0.0      0   28    True  0.841211  0.00
    2  4294967295  8471  0.0      0  542    True  0.394552  0.00
    3  4294967295  9993  0.0      0   10    True  0.930931  0.00
    4  4294967295  9998  0.0      0    4    True  0.129942  0.00
            ftime     g    h  label   n neutral       pos     s
    0  4294967295  9866  0.0      0   3    True  0.121480  0.00
    1  4294967295  9999  1.0      0   1   False -0.287123 -0.05
    2  4294967295  9999  0.0      0   1    True  0.101539  0.00
    3  4294967295  9974  0.0      0  48    True  0.443054  0.00
    4  4294967295  9998  1.0      0   1   False -0.722185 -0.05
            ftime     g    h  label     n neutral       pos     s
    0  4294967295  9998  1.0      0     1   False -0.650965 -0.05
    1  4294967295  9986  1.0      0     5   False -0.427709 -0.05
    2  4294967295  9999  0.0      0     1    True  0.709887  0.00
    3  4294967295  8598  0.0      0  1546    True  0.227702  0.00
    4  4294967295  9938  0.0      0   111    True  0.898534  0.00


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

            ftime     g    h  label   n neutral       pos    s
    1  4294967295  9970  0.0      0  13    True  0.361552  0.0
    3  4294967295  9996  0.0      0   2    True  0.812076  0.0
    4  4294967295  9909  0.0      0  78    True  0.566514  0.0
    5  4294967295  9883  0.0      0  79    True  0.167945  0.0
    6  4294967295  9952  0.0      0   6    True  0.296263  0.0
            ftime     g    h  label    n neutral       pos    s
    1  4294967295  9928  0.0      0   28    True  0.841211  0.0
    2  4294967295  8471  0.0      0  542    True  0.394552  0.0
    3  4294967295  9993  0.0      0   10    True  0.930931  0.0
    4  4294967295  9998  0.0      0    4    True  0.129942  0.0
    6  4294967295  9938  0.0      0   37    True  0.526806  0.0
            ftime     g    h  label     n neutral       pos    s
    0  4294967295  9866  0.0      0     3    True  0.121480  0.0
    2  4294967295  9999  0.0      0     1    True  0.101539  0.0
    3  4294967295  9974  0.0      0    48    True  0.443054  0.0
    6  4294967295  3177  0.0      0  1416    True  0.715795  0.0
    7  4294967295  9087  0.0      0   432    True  0.266955  0.0
            ftime     g    h  label     n neutral       pos    s
    2  4294967295  9999  0.0      0     1    True  0.709887  0.0
    3  4294967295  8598  0.0      0  1546    True  0.227702  0.0
    4  4294967295  9938  0.0      0   111    True  0.898534  0.0
    5  4294967295  9993  0.0      0     1    True  0.516694  0.0
    6  4294967295  9727  0.0      0   234    True  0.311537  0.0


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

    {'neutral': [{'g': 9258, 'h': 0.0, 'pos': 0.01718478580005467, 'n': 316, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8257, 'h': 0.0, 'pos': 0.2615648997016251, 'n': 1711, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8488, 'h': 0.0, 'pos': 0.30679267970845103, 'n': 1711, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9268, 'h': 0.0, 'pos': 0.37840621965005994, 'n': 108, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9174, 'h': 0.0, 'pos': 0.385991113493219, 'n': 108, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9840, 'h': 0.0, 'pos': 0.5386483389884233, 'n': 101, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 7088, 'h': 0.0, 'pos': 0.5690092630684376, 'n': 1710, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 7159, 'h': 0.0, 'pos': 0.6020533256232738, 'n': 1710, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 7814, 'h': 0.0, 'pos': 0.6173335334751755, 'n': 1710, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9272, 'h': 0.0, 'pos': 0.6502099095378071, 'n': 108, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8276, 'h': 0.0, 'pos': 0.6987745801452547, 'n': 1710, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9518, 'h': 0.0, 'pos': 0.9861047994345427, 'n': 1573, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}], 'selected': [{'g': 9974, 'h': 1.0, 'pos': 1.2837485196068883, 'n': 8, 's': -0.05000000074505806, 'ftime': 4294967295, 'label': 0, 'neutral': False}], 'n': 1}
    {'neutral': [{'g': 8764, 'h': 0.0, 'pos': 0.0003967492375522852, 'n': 1653, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9278, 'h': 0.0, 'pos': 0.09578495984897017, 'n': 524, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8726, 'h': 0.0, 'pos': 0.15159097942523658, 'n': 1447, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9016, 'h': 0.0, 'pos': 0.22119286027736962, 'n': 532, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9735, 'h': 0.0, 'pos': 0.36381777515634894, 'n': 335, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9791, 'h': 0.0, 'pos': 0.5702131441794336, 'n': 338, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9309, 'h': 0.0, 'pos': 0.5829871296882629, 'n': 485, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9079, 'h': 0.0, 'pos': 0.6875497191213071, 'n': 499, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9223, 'h': 0.0, 'pos': 0.6895107077434659, 'n': 499, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8448, 'h': 0.0, 'pos': 0.7205626938957721, 'n': 1445, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9259, 'h': 0.0, 'pos': 0.98554897448048, 'n': 499, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}], 'selected': [{'g': 9997, 'h': 1.0, 'pos': -0.8252747415099293, 'n': 4, 's': -0.05000000074505806, 'ftime': 4294967295, 'label': 0, 'neutral': False}], 'n': 4}
    {'neutral': [{'g': 8752, 'h': 0.0, 'pos': 0.09111091843806207, 'n': 1312, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8029, 'h': 0.0, 'pos': 0.3158982314635068, 'n': 841, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9410, 'h': 0.0, 'pos': 0.3186325046699494, 'n': 507, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8726, 'h': 0.0, 'pos': 0.32298249402083457, 'n': 841, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9141, 'h': 0.0, 'pos': 0.33764092018827796, 'n': 841, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 7835, 'h': 0.0, 'pos': 0.3412583307363093, 'n': 841, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 7017, 'h': 0.0, 'pos': 0.3460872450377792, 'n': 1273, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8387, 'h': 0.0, 'pos': 0.3504885055590421, 'n': 841, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9460, 'h': 0.0, 'pos': 0.36058812984265387, 'n': 507, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 5030, 'h': 0.0, 'pos': 0.4330524855758995, 'n': 1426, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 6864, 'h': 0.0, 'pos': 0.47038055770099163, 'n': 1426, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 2161, 'h': 0.0, 'pos': 0.4902536778245121, 'n': 1426, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 2569, 'h': 0.0, 'pos': 0.49638046440668404, 'n': 1426, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8604, 'h': 0.0, 'pos': 0.50460350420326, 'n': 755, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 5706, 'h': 0.0, 'pos': 0.5091155949048698, 'n': 1426, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 5841, 'h': 0.0, 'pos': 0.5238643542397767, 'n': 1442, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 4601, 'h': 0.0, 'pos': 0.6289758402854204, 'n': 1416, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 3177, 'h': 0.0, 'pos': 0.7157951188273728, 'n': 1416, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 4704, 'h': 0.0, 'pos': 0.7261011672671884, 'n': 1352, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8372, 'h': 0.0, 'pos': 0.7275875445920974, 'n': 1352, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 5598, 'h': 0.0, 'pos': 0.7377090277150273, 'n': 1355, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 2753, 'h': 0.0, 'pos': 0.755001342156902, 'n': 1355, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 4733, 'h': 0.0, 'pos': 0.7580284273717552, 'n': 1355, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 3797, 'h': 0.0, 'pos': 0.7767451286781579, 'n': 1355, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9630, 'h': 0.0, 'pos': 0.8485506123397499, 'n': 518, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8715, 'h': 0.0, 'pos': 0.8772767938207835, 'n': 1355, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}], 'selected': [], 'n': 16}
    {'neutral': [{'g': 7698, 'h': 0.0, 'pos': 0.14117315271869302, 'n': 1943, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8440, 'h': 0.0, 'pos': 0.16768103023059666, 'n': 1943, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8712, 'h': 0.0, 'pos': 0.25365707697346807, 'n': 454, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 6884, 'h': 0.0, 'pos': 0.2973291710950434, 'n': 455, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9346, 'h': 0.0, 'pos': 0.353145023342222, 'n': 130, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8134, 'h': 0.0, 'pos': 0.3580396131146699, 'n': 455, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 6057, 'h': 0.0, 'pos': 0.4042419760953635, 'n': 455, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 7426, 'h': 0.0, 'pos': 0.406794270966202, 'n': 455, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 5825, 'h': 0.0, 'pos': 0.4344785118009895, 'n': 455, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 5333, 'h': 0.0, 'pos': 0.46032888046465814, 'n': 455, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9064, 'h': 0.0, 'pos': 0.49631250905804336, 'n': 455, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9860, 'h': 0.0, 'pos': 0.5083352492656559, 'n': 114, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9112, 'h': 0.0, 'pos': 0.7280409946106374, 'n': 1120, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 8680, 'h': 0.0, 'pos': 0.7307664987165481, 'n': 1352, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}, {'g': 9358, 'h': 0.0, 'pos': 0.9505829166155308, 'n': 1027, 's': 0.0, 'ftime': 4294967295, 'label': 0, 'neutral': True}], 'selected': [], 'n': 1}


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
          <th>ftime</th>
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
          <td>4294967295</td>
          <td>9974</td>
          <td>1.0</td>
          <td>0</td>
          <td>8</td>
          <td>False</td>
          <td>1.283749</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>2</th>
          <td>4294967295</td>
          <td>9995</td>
          <td>1.0</td>
          <td>0</td>
          <td>2</td>
          <td>False</td>
          <td>1.290498</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>4</th>
          <td>4294967295</td>
          <td>9996</td>
          <td>1.0</td>
          <td>0</td>
          <td>1</td>
          <td>False</td>
          <td>-0.741675</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>10</th>
          <td>4294967295</td>
          <td>9992</td>
          <td>1.0</td>
          <td>0</td>
          <td>3</td>
          <td>False</td>
          <td>-0.626625</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>11</th>
          <td>4294967295</td>
          <td>9998</td>
          <td>1.0</td>
          <td>0</td>
          <td>1</td>
          <td>False</td>
          <td>-0.073575</td>
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
