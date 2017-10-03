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

int main(int argc, char **argv) {
  double calculated, actual;

  // Test always1
  printf("Testing f(x) = 1 for x in [0, 10]...\n");
  actual = 10;
  calculated = nintegrate(&always1, NULL, 0, 10);
  SHOW_RESULTS(actual, calculated);

  // Test quadratic
  printf("Testing f(x) = x^2 for x in [0, 10]...\n");
  actual = 1.0/3*1000;
  calculated = nintegrate(&quadratic, NULL, 0, 10);
  SHOW_RESULTS(actual, calculated);

  // Test sin2
  printf("Testing f(x) = sin^2(x) for x in [0, pi]...\n");
  actual = PI/2;
  calculated = nintegrate(&sin2, NULL, 0, PI);
  SHOW_RESULTS(actual, calculated);

  // Test log(sin(x) + x)
  printf("Testing f(x) = log(x+sin(x)) for x in [1, 2]\n");
  actual = 0.890373044577735;
  calculated = nintegrate(&ln_sin_xpx, NULL, 1, 2);
  SHOW_RESULTS(actual, calculated);
  exit(0);
}
