import argparse
class Epoch(argparse._AppendAction):
    def __call__(self,parser,namespace,values,option_string=None):
        epoch=[v for v in values]
        try:
            str(epoch[0])
        except:
            parser.error("first argument to epoch must be convertible to string")
        if (epoch[0] != b'growth' and epoch[0] != b'exp'):
            parser.error("first argument to epoch must be 'growth' or 'exp'")
        try:
            int(epoch[1])
        except:
            parser.error("second argument to epoch must be convertible to int")
        if int(epoch[1])<=0:
            parser.error("second argument to epoct must be non-negative integer")
        try:
            int(epoch[2])
        except:
            parser.error("third argument to epoch must be convertible to int")
        if int(epoch[2])<=0:
            parser.error("third argument to epoct must be non-negative integer")
        super(Epoch,self).__call__(parser,namespace,epoch,option_string)

def add_epoch_parser(parser):
    #Population size changes
    group = parser.add_argument_group("Population size history")
    group.add_argument('--burnin','-b',type=int,default=None,
            help="Initial number of generations to simulation. Default of None will be converted to 10*popsize")
    group.add_argument('--epoch','-e',nargs=3,default=None,
            action=Epoch,
            metavar=('TYPE','POPSIZE','NGENS'),help="Change population size.  TYPE must be either 'growth' or 'constant'. In the case of growth, the population size is adjusted to POPSIZE over NGENS generations.  In the case of 'constant', the population size is changed immediately to POPSIZE and evolved for NGENS generations.")
