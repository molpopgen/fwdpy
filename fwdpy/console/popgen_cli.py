#
# Copyright (C) 2015-2016 Kevin Thornton <krthornt@uci.edu>
#
# This file is part of fwdpy.
#
# fwdpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fwdpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy.  If not, see <http://www.gnu.org/licenses/>.
#
from __future__ import print_function

import argparse
import sys
import numpy 
import fwdpy
import fwdpy.fitness
import fwdpy.demography

fwdpy_citation_text= """
If you use fwdpy in your work, please cite the following article(s):

Thornton, K (2014) A C++ Template Library for Efficient Forward-Time
Population Genetic Simulation of Large Populations, 
Genetics 198: 157-166; doi: 10.1534/genetics.114.165019
"""

def sregion_type(s):
    return s

def get_parser():
    """
    Returns an argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description="Command line interface to population-genetic simulations using fwdpy.",
        epilog=fwdpy_citation_text)


    parser.add_argument(
    "-V", "--version", action='version',
    version='%(prog)s {}'.format(fwdpy.__version__))
    # This is required to get uniform behaviour in Python2 and Python3
    #subparsers = parser.add_subparsers(dest="subcommand")
    #subparsers.required = True

    #Arguments for "regions"
    group = parser.add_argument_group("Mutation rates and effect sizes")
    group.add_argument('--nregion','-nr',type=float,nargs=3,default=None,
            action="append",
            metavar=('start','stop','4Nu/site'),help="Parameters are start, stop, theta per site")
    group.add_argument('--recregion','-rr',type=float,nargs=3,default=None,
            action="append",
            metavar=('start','stop','4Nr/site'),help="Parameters are start, stop, rho per site")
    group.add_argument('--sregion','-sr',type=sregion_type,nargs='*',default=None,
            action="append",
            metavar=('type','[params]'),help="TBD")

    #Population size changes
    group = parser.add_argument_group("Population size history")
    group.add_argument('--burnin','-b',type=int,default=None,
            help="Initial number of generations to simulation. Default of None will be converted to 10*popsize")
    group.add_argument('--epoch','-e',type=sregion_type,nargs='*',default=None,
            action='append',
            metavar=('TYPE POPSIZE NGENS'),help="Change population size.  TYPE must be either 'growth' or 'constant'. In the case of growth, the population size is adjusted to POPSIZE over NGENS generations.  In the case of 'constant', the population size is changed immediately to POPSIZE and evolved for NGENS generations.")

    #How many threads to use
    group = parser.add_argument_group("Number of threads, etc.")
    group.add_argument('--nthreads','-T',type=int,default=1,action='store',
            help="How many threads to use.")
    group.add_argument('--nreps','-R',type=int,default=1,
            help="Number of replicates to simulate per thread.")

    group = parser.add_argument_group("Random number seeds")
    group.add_argument('--seed','-S',type=int,default=None,
            help="Random number seed.  If nothing is provided, a random seed will be generated")

    #Positional arguments
    parser.add_argument('popsize',type=int,default=None)

    return parser

def validate_parser(parser):
    error=False
    if parser.popsize is None:
        print("error: popsize cannot be None")
        error=True
    if parser.nregion is None and parser.sregion is None:
        print("error: must have at least one nregion or sregion")
        error=True

    if error is True:
        sys.exit(0)

class SimRunner(object):
    """
    Does the actual work
    """
    def __init__(self,args):
        self.nregions=[]
        self.sregions=[]
        self.recregions=[]
        self.popsize = int(args.popsize)
        self.nthreads=int(args.nthreads)
        self.nreps=int(args.nreps)
        if args.seed is not None:
            self.seed=int(args.seed)
        else:
            self.seed = numpy.random.randint(42000000)

        if args.burnin is None:
            self.burnin=10*self.popsize
        else:
            self.burnin = int(args.burnin)
        if args.nregion is not None:
            for i in args.nregion:
                self.nregions.append(fwdpy.Region(i[0],i[1],i[2]/(4.*float(self.popsize))))
        if args.recregion is not None:
            for i in args.recregion:
                self.recregions.append(fwdpy.Region(i[0],i[1],i[2]/(4.*float(self.popsize))))
        if args.sregion is not None:
            for i in args.sregion:
                sregion_type = i[0]
                if sregion_type == b'exp':
                    self.sregions.append(fwdpy.ExpS(beg=float(i[1]),end=float(i[2]),weight=float(i[3]),mean=float(i[4]/(4.*float(self.popsize))),h=float(i[5])))
                else:
                    print("sregion type "+sregion_type+" not recognized")
                    sys.exit(0)

        self.neutral_mut_rate=0.
        self.selected_mut_rate=0.
        self.recrate=0.
        for i in self.nregions:
            self.neutral_mut_rate += i.w
        for i in self.recregions:
            self.recrate += i.w
        for i in self.sregions:
            self.selected_mut_rate += i.w

        self.popvec=fwdpy.SpopVec(self.nthreads,self.popsize)
        self.nlist=numpy.array([self.popsize]*self.burnin,dtype=numpy.uint32)
        last_size = self.popsize
        if args.epoch is not None:
            for e in args.epoch:
                if e[0] == 'constant':
                    numpy.append(self.nlist,[int(e[1])]*e[2])
                    last_size=int(e[1])
                elif e[1] == 'growth':
                    numpy.append(self.nlist,fwdpy.demography.exponential_size_change(last_size,int(e[1]),int(e[2])))
                    last_size=int(e[1])

    def run(self):
        rng=fwdpy.GSLrng(self.seed)
        fwdpy.evolve_regions_sampler(rng,self.popvec,
                fwdpy.NothingSampler(len(self.popvec)),
                self.nlist,
                self.neutral_mut_rate,
                self.selected_mut_rate,
                self.recrate,
                self.nregions,
                self.sregions,
                self.recregions,
                0)

def popgen_cli_main(arg_list=None):
    parser=get_parser()
    args = parser.parse_args(arg_list)
    validate_parser(args)
    runner=SimRunner(args)
    runner.run()
