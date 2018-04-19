#include <math.h>
#include "mt19937.h"

//Use this along with the provided mt19937.h 

// Random no. generator wihin Normal distribution using Marsaglia Polar method 
// Refer https://en.wikipedia.org/wiki/Marsaglia_polar_method
// Variance=1.0 and mean=0 for this case
double gaussian_rand()
{
    static unsigned int hasSpare = 0;
    static double spare;
    static double u, v, s;
    if(hasSpare)
    {
        hasSpare = 0;
        return spare;
    }
    hasSpare = 1;
    do
    {
        u = 2*dsfmt_genrand() - 1;
        v = 2*dsfmt_genrand() - 1;
        s = u*u + v*v;
    }
    while(s >= 1 || s == 0);
    s = sqrt(-2.0 * log(s) / s);
    spare = v*s;
    return u*s;
}
