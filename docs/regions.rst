# Introduction

Several of the simulation routines allow the details of the mutation and recombination models to vary along a "sequence" or "region".  A user is able to specify the details of such variation by passing _lists_ to package functions.  For example, you are able to:

* Vary the neutral mutation rate along a sequence.
* Vary the distribution of selection coefficients (and the dominance associated with selected mutations) along a sequence.
* Vary the recombination rate along a sequence.

The implementation of such variation along a region is _discrete_.  A region is specified by a beginning, and end, and a weight, plus any additional data required to specify selection coefficients, dominance, etc.

## Background

The models are parameterized through R's S4 class system, which is a method of object-oriented programming in R.  If you are not familiar with it, Hadley Wickham's [tutorial](http://adv-r.had.co.nz/S4.html) is a good place to start.

## Mutation rates, recombination rates, and a weighting system

A simulation will typically have a mutation rate, $\mu$, which represents the mean of a Poisson number of mutations per gamete per generation), and a recombination rate, $r$, which again is the mean of Poisson number of crossover events (per diploid, per generation).  These parameters are the _total_ rates across an entire simulated region.  Variation in these parameters along the region are affected by a set of positions coupled with "weights", which the user specifies using S4 classes.

### The base S4 class: 'region'

A 'region' is an S4 class with the following slots:

* $b$, which is the beginning/start of the region. The type is "numeric" and the default value is 0.
* $e$, which is the end/stop of the region. The type is "numeric" and the default value is 1.
* $w$, which is a weighting factor associated with the region. The type is "numeric" and the default value is 1.
* c$$, which is a logical value which we will explain in more detail below.  The default value is TRUE.

These slots are also documented in the help page:

~~~{r}
help("region-class")
~~~

The four slots are used to inform the C++ code about the relative abundance of new mutations or recombination events will occur in what region.  Briefly, the number of events that occur in region $i$ are proportional to $w_i/\sum_i w$, _i.e._, the weight assigned to region $i$ divided by the sum of weights assigned to all regions.  The weights for mutation events and for recombination events are considered separately.  Thus, in order to model a correlation between mutational processes and recombination, it is up to the user to generate regions whose weights are correlated.

foRward allows the $w$ slot to be interpreted in one of two ways:

* It is affected by the length of a region $(e - b)$.  This is the default, which is specified by $c$ defaulting to TRUE.  Internally, the weight assigned to a region will be converted to $w(b-e)$.
* It is _not_ affected by the length of region.  Interally, the weight assigned is simply $w$.  This approach is specified by setting $c$ to FALSE.

These two approaches allow for considerable modeling flexibility.  For example, the former approach allows _w_ to be interpreted as a "per base-pair" rate.  Imagine that you wanted to simulate variation in recombination along discrete 100 kilobase chunks, and the rate of crossing-over _per base pair_ increases in each chunk, and includes an initial chunk with no recombination:

| Start | Stop | $r_{bp}$ |
|:----:|:----:|:----:|
| 1 | $10^5$ | 0 |
| $10^5$ | $2 \times 10^5$ | $10^{-8}$ |
| $2 \times 10^5$ | $3 \times 10^5$ | $10^{-7}$ |

This model boils down to the relative number of crossing overs per region occuring in the ratio $0 : 10^{-8} : 10^{-7}$.  This is easily represented using foRward's S4 classes:

```{r}
library(foRward)
rec.regions1 = list( region1 = new('region',b=1,e=1e5,w=0),
	    region2 = new('region',b=1e5, e=2e5, w=1e-8),
	    region3 = new('region',b=2e5, e=3e5, w=1e-7) )
#Let's calculate the total recombination rate (per diploid, per region, per generation)
littler = 0
for( i in 1:length(rec.regions1) )
{
 littler = littler + rec.regions1[[i]]@w*(rec.regions1[[i]]@e-rec.regions1[[i]]@b)
}
print(paste("rec rate is ",littler))
```

For this hypothetical example, the region lengths are all even, and thus an equivalent specification would be this:

```{r}
rec.regions2 = list( region1 = new('region',b=0,e=1,w=0,c=FALSE),
	    region2 = new('region',b=1, e=2, w=1e-8,c=FALSE),
	    region3 = new('region',b=2, e=3, w=1e-7,c=FALSE) )
```

Or this:

```{r}
rec.regions3 = list( region1 = new('region',b=0,e=1,w=0,c=FALSE),
	    region2 = new('region',b=1, e=2, w=1,c=FALSE),
	    region3 = new('region',b=2, e=3, w=10,c=FALSE) )
```

The latter two (rec.regions2 and rec.regions3) are equivalent to rec.regions1 assuming a total recombination rate of `r littler`.  The only difference will be the positions of those crossover events, which should tell you that your mutational "regions" and your recombination "regions" must overlap each other in sensible ways!

## Specific examples

### Mutations not affecting fitness ("neutral" mutations)

You specify regions where neutral mutations arise via the S4 class 'region'.  A region has a beginning, end, and a weight, specified by 'slots' _b_, _e_, and _w_, respectively.  The default values are 0, 1, and 1, respectively.  Thus, the following list would specify that 100% of neutral mutations occur on the continuous interval [0,1):

```{r}
require("foRward")
neutral.regions = list( new('region') )
##or, equivalently:
neutral.regions.explicit = list( new('region',b=0,e=1,w=1) )
print(neutral.regions)
```

The beginning and end positions can be whatever you like:

```{r}
##defaulting to weight of 1...
neutral.regions = list( new('region',b=0,e=100) )
```

To specify variation in the netural mutation process along a sequence, combine multiple regions in your list:

```{r}
neutral.regions = list( new('region',b=0,e=1,w=1),
		##This region is 10x "bigger", but we want the same mutational weight per hypothetical "base pair"
		new('region',b=2,e=12,w=1) )
```

Internally, the total "mutational weight" of the first region will be a function of its length, which is 1(1-0)=1.  The second region's total weight will be 1*(12-2)=10, and it will have 10xas many new mutations arising as the first region.

### Mutations affecting fitness

Type types of mutations affecting fitness that we consider will have two parameters associated with them:

* $s$, the selection coefficient
* $h$, the effect of the mutation in a heterozygote (a.k.a. the "dominance" of the mutation).

In a simulation, we may place a distribution on either $s$ itself or on the scaled selection parameter $\alpha = 2Ns$.  These two methods are represented by the S4 classes 'sregion' and 'twoNsregion', respectively.  These classes contain/extend the 'region' class described above, and thus inherit their _b_, _e_, and _w_ slots.  These new classes contain the following additional slots:

* 'sregion' and 'twoNsregion' contain a slot called _dominance_, which is itself a type inheriting the S4 class 'dominance', which we will cover in more detail below.
* 'twoNsregion' contains the slots _sign_ and _twoN_.  The former is a constant by which a selection coefficient will be multiplied.  A sign of -1 corresponds to a deleterious mutation, and 1 would be a beneficial mutation.  The slot _twoN_ is the value of 2N in the $\alpha = 2Ns$.

The following types extend 'sregion':

* 's.constant' represents a mutation model where selected mutations always have the same effect on fitness.
* 's.beta' represents a model where $f(s) = x\beta(\alpha,\beta)$, where $x$ is a scaling parameter that the user may specify.
* 's.gaussian' is a model where $f(s) = N(0,\sigma)$, where $\sigma$ is the standard deviation

The following types extend 'twoNsregion':

* 'twoNs.exp' is a model where $f(\alpha) = \lambda e^{-\lambda \alpha}$, parameterized by a mean equal to $1/\lambda$.
* 'twoNs.gamma' is a model where $f(\alpha ; \overline{\alpha},\beta) = \frac{(\beta/\overline{\alpha})\alpha^{\overline{\alpha}-1}e^{-\beta \alpha/\overline{\alpha}}}{\Gamma (\beta)}$, which is a Gamma distribution with mean $\alpha$ and scale parameter $\beta$.  See the following paper for the rationale for including this distribution: Eyre-Walker, A. (2010). Evolution in health and medicine Sackler colloquium: Genetic architecture of a complex trait and its implications for fitness and genome-wide association studies. Proceedings of the National Academy of Sciences, 107 Suppl 1(suppl 1), 1752-1756. http://doi.org/10.1073/pnas.0906182107

#### Dominance of mutations affecting fitness

All of the type describe above contain a slot called _dominance_, which represents how the dominance of mutations affecting fitness is calculated.  This slot is an S4 class extending the base type 'dominance', which has no slots.  The following S4 types representing dominance are currently supported:

* 'h.constant' is a model where the dominance of mutations is fixed.  It contains a single slot called _h_, which defaults to 1, which represents a mutation with an additive effect.  _All of the sregions and twoNsregions described above default to this type._
* 'h.beta' is a model where $f(h) = x\beta (\alpha,\beta)$, where $x$ is a scaling parameter that the user may specify.  The two shape parameters are specified by the slots _alpha_ and _beta_, respectively, and the slot _scaling_ specifies $x$.


### Crossover rate variation

Just like neutral mutations, intervals with different crossover rates are specified by different 'region' objects.  Let's set up the following concrete example:

* A region where crossovers occur between positions [0,1)
* Positions [0,0.45) and [0.55,1) have uniform recombintion rates
* Positions [0.45,0.55) are a recombination hotspot with 100x the background intensity (per "base pair").

The above model can be represented as:

```{r}
library(foRward)
recrates = list( new('region',b=0,e=0.45,w=1),
new('region',b=0.55,e=1,w=1),
##This is the hotspot:
new('region',b=0.45,e=0.55,w=100) )

print(recrates)
```

Internally, this is what will happen to the above input:

* The total weight on the first region will be $w = w \times (e-b)$ = 1*(0.45-0) = 0.45
* The weight on the second region will be 1*(1-0.55) = 0.45
* The weight on the hotspot will be 100*(0.55-0.45) = 10

This gives us what we want: the hotspot is 100x hotter "per base", and is 10% of the total region in length.  We therefore expect 10x as many crossovers in that region as in the flanking regions.

## Where to get more help/details/examples

For any of the types above, help files exist.  For class X, the help is called X-class.  For example:

~~~{r}
help(h.constant-class)
~~~

## How to set up a model

When setting up a model, it is important that you think in terms of conditional probabilities.  In other words, if the total rate to neutral variants is $\mu_n$, then the weights passed along to a function have the interpretations "Given that a neutral mutation occurs, the probability that it occurs in a certain interval is $x$", where $x$ is determined by the relative weight assigned to an interval.

The 'weights' that you assign are _relative_ and need not sum to 1.  Each weight must be $> 0$, though.  They are used to generate discrete probability distributions for sampling using code like the following:

```{r}
suppressWarnings(library(Rcpp))

Sys.setenv("PKG_LIBS"="-lgsl -lgslcblas","PKG_CXXFLAGS"="-std=c++11")

sourceCpp(code="#include <Rcpp.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Use Rcpp to create an R function
// [[Rcpp::export()]]
Rcpp::IntegerVector GSLdiscrete() {
    //GSL rng setup		    		 
    gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,0);
    
    //The weights
    std::vector<double> weights={10.,2.,0.001};
    //Generate a lookup table for a discrete distribution with the above weights
    gsl_ran_discrete_t * lookup = gsl_ran_discrete_preproc(weights.size(),&weights[0]);

    //Sample from it
    Rcpp::IntegerVector rv(1000000);
    for( int i = 0 ; i < rv.size() ; ++i )
    {
	rv[i] = gsl_ran_discrete(r,lookup);
    }
    //Cleanup and return
    gsl_ran_discrete_free(lookup);
    gsl_rng_free(r);
    return rv;
}")

x = GSLdiscrete()
table(x)
```

The output should be in proportions similar to the weights set in the above C++ function.

See [here](https://www.gnu.org/software/gsl/manual/html_node/General-Discrete-Distributions.html) for more details on the GSL discrete sampler.

## Example

```{r}
library(foRward)

n=commandArgs(trailing=TRUE)
SEED = 202

## Initialize the random number system
rng = init.gsl.rng(SEED)

##Some basic parameters
N=1e3
theta=100
rho=100

##All neutral muts are [0,1)
nregions = list( new("region",b=0,e=1,w=1) )

##Selected mutations.  All are additive, to keep this example simple.
smodels = list(
    ##Strongly deleterious mutations to the "left"
    new("s.constant",s=-0.1,b=-1,e=0,w=0.99/2),
    ##Weaker mutations (2Ns = 10 on average) to the "right"
    new("s.exp",mean=-10,b=1,e=2,w=0.99/2),
    ## 1% of selected mutations will be positively selected
    ## and uniform throughout the region.  The distribution
    ## of s will be beta(2,10)
    new("s.beta",alpha=2,beta=10,scaling=1,b=-1,e=2,w=0.01)
)

##Recombination models
rmodels = list(
    ##uniform throughtout the region 
    new("region",b=-1,e=1,w=1),
    ## 10x hotspot in middle
    new("region",b=0.45,e=0.55,w=10)
)

pop = evolve.regions(rng,
    ##pop size over time (short sim...)
    rep(N,1+N),
    ##neutral mut rate
    theta/(4*N),
    ##selected mut rate is 1/10th the neutral mutation rate
    0.1*theta/(4*N),
    rho/(4*N),
    nregions,
    smodels,
    rmodels)

#Take a sample from the population.
pop.sample = sample.single.deme(pop[[1]],rng,10)

print(pop.sample)
```

## A simple example of "exons"

Let's consider the following model of a "two-exon gene":

* Exon 1 spans positions 0.2 to 0.4.
* Exon 2 spans positions 0.6 to 1.
* There will be selection against amino acid replacements.

For the sake of simplicity, we will ignore:

* mutations in the introns and UTR regions

We can define the following region lengths:

* $l_1 = 0.2$ is the length of Exon 1.
* $l_2 = 0.4$ is the length of Exon 2.

Thus the relative weights assigned to Exon 1 and 2 must satisfy $w_2/w_1 = l_2/l_1 = 2.$.

Further, within an exon, $\approx 3/4$ of new mutations will be amino acid replacements.  Thus, 3/4 of the total mutational weight in a region will be on selected mutations.

Our model looks like this:

~~~{r}
nregions = list( new('region',b=0.2,e=0.4,w=0.25),
	 new('region',b=0.6,e=1,w=0.5) )
sregions = list( new('s.constant',b=0.2,e=0.4,w=0.75,s=-0.01),
	 new('s.constant',b=0.6,e=0.1,w=1.5,s=-0.01) )
~~~

If we pass _the same neutral and selected mutation rates to evolve.regions_, then the above model satisfies:

* The total number of mutations occurring in Exon 2 is 2x the number occuring in Exon 1.
* Within an expon, 3/4 of all new mutations are deleterious.
