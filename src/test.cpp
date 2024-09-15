#include <stdio.h>
//#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>

using namespace std;

int main (void)
{
const gsl_rng_type * T;
gsl_rng * r;
int i, n = 1000000;
gsl_rng_env_setup();
T = gsl_rng_default;
r = gsl_rng_alloc (T);

double F=2.;
double SD=2;
for (i = 0; i < n; i++)
{
double u =(gsl_ran_gaussian_tail (r, -F, SD)+F);
printf ("%1.5f\n", u);
}
gsl_rng_free (r);
return 0;
}
