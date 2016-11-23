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

import fwdpy
import fitness
import fwdpyio

fwdpy_citation_text= """
If you use fwdpy in your work, please cite the following article(s):

Thornton, K (2014) A C++ Template Library for Efficient Forward-Time
Population Genetic Simulation of Large Populations, 
Genetics 198: 157-166; doi: 10.1534/genetics.114.165019
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
    # This is required to get uniform behaviour in Python2 and Python3
    subparsers = parser.add_subparsers(dest="subcommand")
    subparsers.required = True
    return parser

def popgen_cli_main(arg_list=None):
    print("I am here")
    parser=get_parser()
    args = parser.parse_args(arg_list)
