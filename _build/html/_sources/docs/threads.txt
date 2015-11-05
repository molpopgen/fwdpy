How many threads to use?
===========================

Many of the simulation functions in this package allow you to simulate more than one replicate at a time.   Each replicate will be simulated using a different thread of execution.

In order to maximize performance, you want to use the "right" number of threads.  However, that number is hard to know, as it depends on many things, including:

* The population size(s) being simulated.
* The number of cores on your machine.

For example, our HPC system at UCI is mostly made up of 64-core AMD Opteron nodes.  If I simulate single demes of size :math:`N=10^3`, I can simulate 64 replicates at a time (using 64 cores), and use one hundred percent of a node.  However, changing to :math:`N=10^4` results in less efficient execution when using 64 cores, effectively slowing each replicate down by about 30 percent.  However, I can run two 32-core jobs on a node and get back to 100 percent CPU use.  I would suggest some experimentation to figure out what is most effective on your system.  Utilities like top and htop will be invaluable here.  
