#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
int
main(int argc, char **argv)
{
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);

    gsl_rng_free(r);
}
