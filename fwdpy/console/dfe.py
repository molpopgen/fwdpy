import argparse

def rescale_neutral_mut_rec_params(args):
    """
    Convert scaled parameters into per-generation parameters for neutral
    mutation rates and recombination rates.
    """
    for i in args.nregion:
        i[2]/=float(4*args.popsize)
    for i in args.recregion:
        i[2]/=float(4*args.popsize)

def rescale_popgen_DFE(args):
    """
    Convert scaled params for DFE into per-generation parameters.
    """
    for i in args.exp:
        i[2]/=float(4*args.popsize)
        i[3]/=float(4*args.popsize)
    for i in args.constant:
        i[2]/=float(4*args.popsize)
        i[3]/=float(4*args.popsize)
    for i in args.uniform:
        i[2]/=float(4*args.popsize)
        i[3]/=float(4*args.popsize)
        i[4]/=float(4*args.popsize)
    for i in args.gamma:
        i[2]/=float(4*args.popsize)
        i[3]/=float(4*args.popsize)

def add_popgen_dfe_parser(parser):
    #Arguments for "regions"
    group = parser.add_argument_group("Mutation rates and distributiosn of effect sizes")
    group.add_argument('--neutral','-n',type=float,nargs=3,default=None,
            action="append",dest='nregion',
            metavar=('START','STOP','4Nu/SITE'),help="Parameters are start, stop, theta per site")
    group.add_argument('--recombination','-r',type=float,nargs=3,default=None,
            action="append",dest='recregion',
            metavar=('START','STOP','4Nr/SITE'),help="Parameters are start, stop, rho per site")
    group.add_argument('--exp','-ex',type=float,nargs=5,action='append',default=[],
            metavar=('START', 'STOP', '4Nv/SITE', 'MEAN_4Ns', 'DOMINANCE'))
    group.add_argument('--gamma','-ga',type=float,nargs=6,action='append',default=[],
            metavar=('START', 'STOP', '4Nv/SITE', '4Ns', 'SHAPE','DOMINANCE'))
    group.add_argument('--constant','-co',type=float,nargs=5,action='append',default=[],
            metavar=('START', 'STOP', '4Nv/SITE', '4Ns', 'DOMINANCE'))
    group.add_argument('--uniform','-un',type=float,nargs=6,action='append',default=[],
            metavar=('START', 'STOP', '4Nv/SITE', '4Ns_LO', '4Ns_HI','DOMINANCE'))
