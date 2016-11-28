#
# Copyright (C) 2015-2016 Kevin Thornton <krthornt@uci.edu>
#
# This file is part of fwdpy.
#
# fwdpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
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

import fwdpy.console.dfe as dfe
import fwdpy.console.validate as validate
import fwdpy.console.demography as demography
import fwdpy.console.common as common

fwdpy_citation_text= """
If you use fwdpy in your work, please cite the following article(s):

Thornton, K (2014) A C++ Template Library for Efficient Forward-Time
Population Genetic Simulation of Large Populations, 
Genetics 198: 157-166; doi: 10.1534/genetics.114.165019
"""

sregions_help="""
The following DFE models are allowed:
(constant start stop 4Nv/site 4Ns h),
(exp start stop 4Nv/site mean_4Ns h),
(uniform start stop 4Nv/site 4Ns_lo 4Ns_hi h),
(gamma start stop 4Nv/site mean_4Ns shape h), 
where 4Nv/site is the scaled mutation rate.
"""

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
    parser.add_argument(
    "-v", "--verbose", action='store_true',help="Verbose output.  Extra info printed to standard error stream.")
    # This is required to get uniform behaviour in Python2 and Python3
    #subparsers = parser.add_subparsers(dest="subcommand")
    #subparsers.required = True

    dfe.add_popgen_dfe_parser(parser)
    demography.add_epoch_parser(parser)
    common.add_common_option_parser(parser)

    #Positional arguments
    parser.add_argument('popsize',type=int,default=None)
    return parser


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
        self.verbose=args.verbose
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
                self.nregions.append(fwdpy.Region(*i))
        if args.recregion is not None:
            for i in args.recregion:
                self.recregions.append(fwdpy.Region(*i))
        for i in args.exp:
            self.sregions.append(fwdpy.ExpS(*i))
        for i in args.constant:
            self.sregions.append(fwdpy.ConstantS(*i))
        for i in args.uniform:
            self.sregions.append(fwdpy.UniformS(*i))
        for i in args.gamma:
            self.sregions.append(fwdpy.GammaS(*i))

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
        self.demography_string = b'\tBurn in at N = ' + str(self.popsize)
        self.demography_string += ' for ' + str(self.burnin) + b' generations\n';
        if args.epoch is not None:
            for e in args.epoch:
                if e[0] == 'constant':
                    numpy.append(self.nlist,[int(e[1])]*int(e[2]))
                    self.demography_string += b'\tChange N from ' + format(last_size)
                    self.demography_string += b' to ' + format(e[1])
                    self.demography_string += b' for ' + format(e[2]) + b' generations\n';
                    last_size=int(e[1])
                elif e[0] == 'growth':
                    numpy.append(self.nlist,fwdpy.demography.exponential_size_change(last_size,int(e[1]),int(e[2])))
                    self.demography_string += b'\tExponential size change from N = ' + format(last_size)
                    self.demography_string += b' to N = ' +format(e[1]) + b' over ' + format(e[2]) + b' generations.\n'
                    last_size=int(e[1])
                else:
                    print("epoch type must be either 'constant' or 'growth'")
                    sys.exit(1)
    def __str__(self):
        rep=b'Simulation details:\n'
        rep+=b'Population size: ' + str(self.popsize) + b'\n'
        rep+=b'Seed: ' + str(self.seed) + b'\n'
        rep+=b'Number of threads: ' + str(self.nthreads) + b'\n'
        rep+=b'Number of replicates per thread: ' + str(self.nreps) + b'\n'
        rep+=b'Population size history:\n' + self.demography_string
        rep+=b'Mutation and recombination rates:\n'
        rep+=b'\tTotal neutral mutation rate (per gamete per generation): ' + str(self.neutral_mut_rate) + b'\n'
        rep+=b'\tTotal selected mutation rate (per gamete per generation): ' + str(self.selected_mut_rate) + b'\n'
        rep+=b'\tTotal recombination rate (per diploid per generation): ' + str(self.recrate) + b'\n'
        rep+=b'Neutral mutations:\n'
        for i in self.nregions:
            rep += b'\t' + str(i) + b'\n'
        rep += b'Distributions of fitness effects:\n'
        for i in self.sregions:
            rep += b'\t' + str(i) + b'\n'
        rep += b'Recombination rate variation:\n'
        for i in self.recregions:
            rep += b'\t' + str(i) + b'\n'
        return rep
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
    dfe.rescale_neutral_mut_rec_params(args)
    dfe.rescale_popgen_DFE(args)
    validate.validate_parser(args)
    runner=SimRunner(args)
    if args.verbose is True:
        print(runner,file=sys.stderr)
    runner.run()
