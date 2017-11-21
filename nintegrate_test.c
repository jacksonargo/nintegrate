// Tests the nintegrate function

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>  // needed for DBL_EPSILON
#include "nintegrate.h"

#define COLOR_RED     "\x1b[31m"
#define COLOR_GREEN   "\x1b[32m"
#define COLOR_YELLOW  "\x1b[33m"
#define COLOR_BLUE    "\x1b[34m"
#define COLOR_MAGENTA "\x1b[35m"
#define COLOR_CYAN    "\x1b[36m"
#define END_COLOR     "\x1b[0m"

#define PI 3.14159265358979323846264338327950
#define E 2.718281828459045235360287471352662497757247093699959574966

#define BAD_ERROR_THRESHOLD 1e-10

double always1(double _x, void* params) {
  return 1;
}

double quadratic(double x, void* params) {
  return x*x;
}

double sin2x(double x, void* params) {
  return sin(x)*sin(x);
}

double sinx2(double x, void* params) {
  return sin(x*x);
}

double sinsinx(double x, void *params) {
  return sin(sin(x));
}

double ln_sin_xpx(double x, void *params) {
  return log(sin(x)+x);
}

double ln_cos_x(double x, void *params) {
  return log(cos(x));
}

double ln_poly(double x, void *params) {
  return log(1 + 3*x + x*x)/x;
}

double poly1(double x, void *params) {
  return 1 + 1/sqrt(x) + 1/x;
}

double normdistribution(double x, void* param) {
  int u = 1;
  int s = 1;
  double coef = 1/(sqrt(2*PI*s*s));
  double exp = -1*(x-u)*(x-u)/(2*s*s);
  return coef*pow(E, exp);
}

double percent_error(double actual, double calculated) {
  return fabs(100*(actual-calculated)/actual);
}

void show_results(double actual, double calculated) {
  double pe = percent_error(actual, calculated);
  printf("Actual:     %.*e\n", DBL_DIG, actual);
  printf("Calculated: %.*e\n", DBL_DIG, calculated);
  if(pe > BAD_ERROR_THRESHOLD) printf(COLOR_RED);
  printf("%%Error:     %.*e\n\n", DBL_DIG, pe);
  if(pe > BAD_ERROR_THRESHOLD) printf(END_COLOR);
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
  show_results(actual, calculated);
  return percent_error(actual, calculated);
}

int main(int argc, char **argv) {
  double total_error = 0;
  printf("Running tests\n");
  printf("Machine epsilon: %.*e\n\n", DBL_DIG, DBL_EPSILON);
  // Test always1
  total_error += test_function(&always1, NULL, 0, 10, "f(x) = 1", 10);

  // Test quadratic
  total_error += test_function(&quadratic, NULL, 0, 10, "f(x) = x^2", 1.0/3*1000);

  // Test sin2x
  total_error += test_function(&sin2x, NULL, 0, PI, "f(x) = sin^2(x)", PI/2);

  // Test sinx2
  total_error += test_function(&sinx2, NULL, 0, PI, "f(x) = sin(x^2)", 0.77265171269006565);

  total_error += test_function(&sinsinx, NULL, 0, PI, "f(x) = sin(sin(x))", 1.786487481950061);

  // Test log(sin(x) + x)
  total_error += test_function(&ln_sin_xpx, NULL, 1, 2, "f(x) = log(x+sin(x))", 0.890373044577735);

  // Test log(cos(x))
  total_error += test_function(&ln_cos_x, NULL, 0, 0.5*PI, "f(x) = log(cos(x))", -1.088793045151791);

  total_error += test_function(&ln_poly, NULL, 0, 1, "f(x) = log(1+3x+x^2)/x", 2.108063708002621);

  total_error += test_function(&poly1, NULL, 1, 2, "f(x) = 1+1/sqrt(x) + 1/x", 2.52157430530614);

  printf("Tests complete!\n");
  if(total_error > BAD_ERROR_THRESHOLD) printf(COLOR_RED);
  printf("Total %%Error: %.*e\n", DBL_DIG, total_error);
  if(total_error > BAD_ERROR_THRESHOLD) printf(END_COLOR);
  exit(0);
}
