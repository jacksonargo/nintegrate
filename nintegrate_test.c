// Tests the nintegrate function

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "nintegrate.h"

#define PI 3.14159265358979323846264338327950
#define E 2.718281828459045235360287471352662497757247093699959574966

#define SHOW_RESULTS(a, b) (printf("Actual:     %.*e\nCalculated: %.*e\n%%Error:     %.*e\n\n", \
      DBL_DIG, a, DBL_DIG, b, DBL_DIG, percent_error(a, b)))

double percent_error(double actual, double calculated) {
  return fabs(100*(actual-calculated)/actual);
}

double always1(double _x, void* params) {
  return 1;
}

double quadratic(double x, void* params) {
  return x*x;
}

double sin2(double x, void* params){
  return sin(x)*sin(x);
}

double ln_sin_xpx(double x, void *params) {
  return log(sin(x)+x);
}

double normdistribution(double x, void* param) {
  int u = 1;
  int s = 1;
  double coef = 1/(sqrt(2*PI*s*s));
  double exp = -1*(x-u)*(x-u)/(2*s*s);
  return coef*pow(E, exp);
}

double test_function(
  double (*func)(double x, void* params),
  void* params,
  double x1,
  double x2,
  char* description,
  double actual
) {
  double calculated;

  printf("Testing %s\n\tfor x in [%e, %e]\n", description, x1, x2);
  calculated = nintegrate(func, params, x1, x2);
  SHOW_RESULTS(actual, calculated);
  return percent_error(actual, calculated);
}

int main(int argc, char **argv) {
  // Test always1
  test_function(&always1, NULL, 0, 10, "f(x) = 1", 10);

  // Test quadratic
  test_function(&quadratic, NULL, 0, 10, "f(x) = x^2", 1.0/3*1000);

  // Test sin2
  test_function(&sin2, NULL, 0, PI, "f(x) = sin^2(x)", PI/2);

  // Test log(sin(x) + x)
  test_function(&ln_sin_xpx, NULL, 1, 2, "f(x) = log(x+sin(x))", 0.890373044577735);

  exit(0);
}
