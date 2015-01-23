#include "nintegrate.h"

double change_interval (double (*func)(double, void*), void *params, double z, double a, double b) {
    return (b-a)/2*func(z*(b-a)/2+(b+a)/2, params);
}

double adaptive_gauss_quadrature (double (*func)(double, void*), void *params, double a, double b, int recursions, quadrature_weights *weights) {
    int i;
    double gauss_total;
    double kronrod_total;
    double func_eval[weights->n_kronrod];

    /* First we make all the function evaluations */
    for (i = 0; i < weights->n_kronrod; i++) {
        func_eval[i] = change_interval (func, params, weights->node[i], a, b);
    }

    /* Next calculate the gauss and kronrod sums */
    gauss_total = 0;
    for (i = 0; i < weights->n_gauss; i++) {
        gauss_total += func_eval[i] * weights->gauss[i];
    }
    kronrod_total = 0;
    for (i = 0; i < weights->n_kronrod; i++) {
        kronrod_total += func_eval[i] * weights->kronrod[i];
    }

    /* Finally, check for error and recursion cap */
    if (recursions < 0) // If recursions negative then there is no limit
        recursions = -2;

    if (200 * pow(fabs(kronrod_total - gauss_total), 1.5) < DBL_EPSILON || recursions == 0 )
        return kronrod_total;
    else {
        double S1 = adaptive_gauss_quadrature (func, params, a, (a+b)/2, recursions/2, weights);
        double S2 = adaptive_gauss_quadrature (func, params, (a+b)/2, b, recursions/2, weights);
        return S1 + S2;
    }
}


double nintegrate(double (*func)(double, void*), void *params, double a, double b) {
    return nintegrate_r(func, params, a, b, 1024);
}

double nintegrate_r(double (*func)(double, void*), void *params, double a, double b, int recursions) {
    /* We have to define the nodes and weights for the adaptive gauss quadrature */
    quadrature_weights weights;
    double gauss_weights[]= {
        0.4179591836734693877551020,
        0.3818300505051189449503698,
        0.3818300505051189449503698,
        0.2797053914892766679014678,
        0.2797053914892766679014678,
        0.1294849661688696932706114,
        0.1294849661688696932706114
    };
    double kronrod_weights[]= {
        0.2094821410847278280129992,
        0.1903505780647854099132564,
        0.1903505780647854099132564,
        0.1406532597155259187451896,
        0.1406532597155259187451896,
        0.0630920926299785532907007,
        0.0630920926299785532907007,
        0.2044329400752988924141620,
        0.2044329400752988924141620,
        0.1690047266392679028265834,
        0.1690047266392679028265834,
        0.1047900103222501838398763,
        0.1047900103222501838398763,
        0.0229353220105292249637320,
        0.0229353220105292249637320
    };
    double nodes[]= {
         0.0000000000000000000000000,
         0.4058451513773971669066064,
        -0.4058451513773971669066064,
         0.7415311855993944398638648,
        -0.7415311855993944398638648,
         0.9491079123427585245261897,
        -0.9491079123427585245261897,
         0.2077849550078984676006894,
        -0.2077849550078984676006894,
         0.5860872354676911302941448,
        -0.5860872354676911302941448,
         0.8648644233597690727897128,
        -0.8648644233597690727897128,
         0.9914553711208126392068547,
        -0.9914553711208126392068547
    };

    /* We put these in a structure to be passed to the recursive function */
    weights.n_gauss = 7;
    weights.n_kronrod = 15;
    weights.node = nodes;
    weights.gauss = gauss_weights;
    weights.kronrod = kronrod_weights;

    /* Now return the numerical integral */
    return adaptive_gauss_quadrature(func, params, a, b, recursions, &weights);
}
