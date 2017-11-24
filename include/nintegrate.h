/* nintegrate.h */

#ifndef __NINTEGRATE__
#define __NINTEGRATE__

#include <math.h>   // needed for pow() and fabs()
#include <float.h>  // needed for DBL_EPSILON

/** quadrature_weights
 * This is a structure to store the weights and nodes of the adaptive gauss quadrature.
 */
typedef struct __quadrature_weights quadrature_weights;
struct __quadrature_weights {
    int n_gauss;      ///< Number of gauss weights.
    int n_kronrod;    ///< Number of kronrod weights.
    const double *node;     ///< Array of nodes to evaluate at.
    const double *gauss;    ///< Array of gauss weights.
    const double *kronrod;  ///< Array of kronrod weights.
};

/** nintegrate_r()
 * Returns the numerical integral of f from a to b. This function calls the recursive function 
 * adaptive_gauss_quadrature which will recurse no deeper than log2(recursions). If recursions is
 * negative, then adaptive_gauss_quadrature will recurse until estimated error is within DBL_EPSILON.
 */
double nintegrate_r(
    double (*f)(double, void*), ///< Function to be integrated.
    void *params,               ///< Parameters for the function.
    double a,                   ///< Lower limit.
    double b,                   ///< Upper limit.
    int recursions              ///< Maximum number of recursions.
);

/** nintegrate()
 * Returns the numerical integral of f from a to b. Identical to nintegrate_r, except this always sets
 * the recursion limit to 1024.
 */
double nintegrate(
    double (*f)(double, void*), ///< Function to be integrated.
    void *params,               ///< Parameters for the function.
    double a,                   ///< Lower limit.
    double b                    ///< Upper limit.
);

/** change_interval()
 * This function scales f so that the interval of integration is -1 to 1,
 * then evaluates at z. This is only called by adaptive_gauss_quadrature. 
 */
double change_interval (
    double (*f)(double, void*), ///< Function to be scaled/evaluated.
    void *params,               ///< Parameters for the function.
    double z,                   ///< Node to evaluate at.
    double a,                   ///< Lower limit.
    double b                    ///< Upper limit.
);

/** adaptive_gauss_quadrature()
 * Returns the numerical integral of f using kronrod nodes and weights. This function estimates error
 * by comparing the kronrod sum and the gauss sum. If the error bigger than DBL_EPSILON, then this
 * function calls itself on both halves of the interval and returns the sum of the numerical integral
 * of both halves. This function will only recurse if recursions != 0.
 */
double adaptive_gauss_quadrature (
    double (*f)(double, void*), ///< Function to be integrated.
    void *params,               ///< Parameters for the function.
    double a,                   ///< Lower limit.
    double b,                   ///< Upper limit.
    int recursions,             ///< Number of recursions left.
    quadrature_weights *weights ///< Nodes and weights to use.
);

#endif
