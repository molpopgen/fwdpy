A note for R users
==============================

Many users may prefer to do downstram analysis, plotting, etc., using R.  The basic data object in R is a "data.table".  Python's Pandas_ library provides pandas.DataFrame, which is an analagous data structure.

These pandas DataFrame objects may be written to a text file as follow:

.. code-block:: python

		#Assume x is a pandas.DataFrame.
		#The index=False means to not print row names,
		#which you typically don't want in a file that
		#will be process in R
		x.to_csv("file.txt",sep="\t",index=False)

Then, in R:

.. code-block:: r

		#You can now dplyr and ggplot2 this thing:
		x=read.table("file.txt",header=T)

R users should read the tutorials at the Pandas_ site for how to manipulate DataFrame objects.

.. _Pandas: http://pandas.pydata.org/
