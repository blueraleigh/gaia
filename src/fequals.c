#include <R.h>
#include "fequals.h"

#define SQRT_DBL_EPSILON 1.490116119384765696e-8

// test for a == b
// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
int
fequals(double a, double b)
{
    if (!R_FINITE(a) || !R_FINITE(b)) return 0;
    double diff = fabs(a - b);
    if (diff <= SQRT_DBL_EPSILON)
    {
        return 1;
    }
    a = fabs(a);
    b = fabs(b);
    double M = (b > a) ? b : a;
    return (diff <= M * SQRT_DBL_EPSILON) ? 1 : 0;
}
