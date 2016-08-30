#ningtegrate

##A set of functions for an adaptive guass quadrature integration.

1) Include the nintegrate header.

    #include "nintegrate.h"

2) Define your function to integrate with this prototype.

    double f(double, void*);

3) The nintegrate function will pass a void* as parameters to the function.

4) Call the function

    nintegrate(f, params, LOWER, UPPER);

