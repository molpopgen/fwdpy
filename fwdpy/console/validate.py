from __future__ import print_function
import sys

def validate_parser(parser):
    error=False
    if parser.popsize is None:
        print("error: popsize cannot be None")
        error=True
    if parser.nregion is None and (parser.exp == [] and parser.gamma==[] and parser.constant ==[] and parser.uniform == []):
        print("error: must specify at least one neutral region or distribution of fitness effects.")
        error=True

    if error is True:
        sys.exit(0)
