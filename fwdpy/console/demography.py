def epoch_type(e):
    return e

def add_epoch_parser(parser):
    #Population size changes
    group = parser.add_argument_group("Population size history")
    group.add_argument('--burnin','-b',type=int,default=None,
            help="Initial number of generations to simulation. Default of None will be converted to 10*popsize")
    group.add_argument('--epoch','-e',type=epoch_type,nargs=3,default=None,
            action='append',
            metavar=('TYPE','POPSIZE','NGENS'),help="Change population size.  TYPE must be either 'growth' or 'constant'. In the case of growth, the population size is adjusted to POPSIZE over NGENS generations.  In the case of 'constant', the population size is changed immediately to POPSIZE and evolved for NGENS generations.")
