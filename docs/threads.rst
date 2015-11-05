How many threads to use?
===========================

Many of the simulation functions in this package allow you to simulate more than one replicate at a time.   Each replicate will be simulated using a different thread of execution.

In order to maximize performance, you want to use the "right" number of threads.  However, that number is hard to know, as it depends on many things, including:

* The population size(s) being simulated.
* The number of cores on your machine.

If you have a 64-core machine, you may or may not see a 100% load on that machine with 64 cores.  On the UCI HPC system, we may see loads more like 95-97%.  If you want to squeeze out the remaining few percent, run 32 core jobs and merge the output.
