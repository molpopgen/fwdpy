def add_common_option_parser(parser):
    #How many threads to use
    group = parser.add_argument_group("Number of replicates. The total number of simulations run will be nthreads*nreps.")
    group.add_argument('--nthreads','-T',type=int,default=1,action='store',
            help="How many threads to use.")
    group.add_argument('--nreps','-R',type=int,default=1,
            help="Number of replicates to simulate per thread.")

    group = parser.add_argument_group("Random number seeds")
    group.add_argument('--seed','-S',type=int,default=None,
            help="Random number seed.  If nothing is provided, a random seed will be generated")
