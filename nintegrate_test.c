// Tests the nintegrate function

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <signal.h>
#include <sys/resource.h>
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

struct stats_t {
  int n_calls;
  double pe;
  struct rusage usage_pre;
  struct rusage usage_post;
};

struct stats_t* init_stats_t() {
  struct stats_t *stats = malloc(sizeof(struct stats_t));
  stats->n_calls = 0;
  stats->pe = 0;
  getrusage(RUSAGE_SELF, &(stats->usage_pre));
  return stats;
}

void track_stats(void *params) {
  struct stats_t *f_stats = (struct stats_t*)(params);
  f_stats->n_calls++;
}

double always1(double _x, void* params) {
  track_stats(params);
  return 1;
}

double quadratic(double x, void* params) {
  track_stats(params);
  return x*x;
}

double sin2x(double x, void* params) {
  track_stats(params);
  return sin(x)*sin(x);
}

double sinx2(double x, void* params) {
  track_stats(params);
  return sin(x*x);
}

double sinsinx(double x, void *params) {
  track_stats(params);
  return sin(sin(x));
}

double ln_sin_xpx(double x, void *params) {
  track_stats(params);
  return log(sin(x)+x);
}

double ln_cos_x(double x, void *params) {
  track_stats(params);
  return log(cos(x));
}

double ln_poly(double x, void *params) {
  track_stats(params);
  return log(1 + 3*x + x*x)/x;
}

double poly1(double x, void *params) {
  track_stats(params);
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

double track_percent_error(double error) {
  static double total_error = 0;
  total_error += error;
  return total_error;
}

void show_results(struct stats_t *stats, double actual, double calculated) {
  double pe = percent_error(actual, calculated);
  struct rusage usage_pre = stats->usage_pre;
  struct rusage usage_post = stats->usage_post;
  int utime = usage_post.ru_utime.tv_usec - usage_pre.ru_utime.tv_usec;

  stats->pe = pe;

  printf("Actual: %.*e\n", DBL_DIG, actual);
  printf("Calced: %.*e\n", DBL_DIG, calculated);
  if(pe > BAD_ERROR_THRESHOLD) printf(COLOR_RED);
  printf("%%Error: %.*e\n\n", DBL_DIG, pe);
  if(pe > BAD_ERROR_THRESHOLD) printf(END_COLOR);

  printf("function calls: %i\n", stats->n_calls);

  printf("usr time: %i\n\n", utime);
}

void test_function(
  double (*func)(double x, void* params),
  double x1,
  double x2,
  char* description,
  double actual
) {
  int x, recursions[] = {0, 1024, -1};
  int n_r = sizeof(recursions)/sizeof(int);
  double calculated;
  struct stats_t *stats;

  printf("Testing %s\n\tfor x in [%e, %e]\n\n", description, x1, x2);

  for(x = 0; x < n_r; x++) {
    if(recursions[x] == -1)
      printf("Recurse until complete...\n\n");
    else
      printf("Limit %d recursions...\n\n", recursions[x]);
    stats = init_stats_t();
    calculated = nintegrate_r(func, stats, x1, x2, recursions[x]);
    getrusage(RUSAGE_SELF, &(stats->usage_post));
    show_results(stats, actual, calculated);
    if(recursions[x] == -1)
      track_percent_error(stats->pe);
    free(stats);
  }
}

void print_total_error(void) {
  double total_error = track_percent_error(0.0);
  if(total_error > BAD_ERROR_THRESHOLD) printf(COLOR_RED);
  printf("Total %%Error: %.*e\n\n", DBL_DIG, total_error);
  if(total_error > BAD_ERROR_THRESHOLD) printf(END_COLOR);
}

void print_ru(void) {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  printf("usr time: %i\n", usage.ru_utime.tv_usec);
  printf("sys time: %i\n", usage.ru_stime.tv_usec);
}

void sigint_handler(int sig) {
  printf("Tests interrupted!\n\n");
  print_total_error();
  print_ru();
  exit(0);
}

int main(int argc, char **argv) {
  clock_t start = clock(); // Make this the first operation.
  double total_error;

  signal(SIGINT, sigint_handler);

  printf("Running tests...\n");

  printf("System constants:\n");
  printf("Machine epsilon: %.*e\n\n", DBL_DIG, DBL_EPSILON);
  // Test always1
  test_function(&always1, 0, 10, "f(x) = 1", 10);

  // Test quadratic
  test_function(&quadratic, 0, 10, "f(x) = x^2", 1.0/3*1000);

  // Test sin2x
  test_function(&sin2x, 0, PI, "f(x) = sin^2(x)", PI/2);

  // Test sinx2
  test_function(&sinx2, 0, PI, "f(x) = sin(x^2)", 0.77265171269006565);

  test_function(&sinsinx, 0, PI, "f(x) = sin(sin(x))", 1.786487481950061);

  // Test log(sin(x) + x)
  test_function(&ln_sin_xpx, 1, 2, "f(x) = log(x+sin(x))", 0.890373044577735);

  // Test log(cos(x))
  test_function(&ln_cos_x, 0, 0.5*PI, "f(x) = log(cos(x))", -1.088793045151791);

  test_function(&ln_poly, 0, 1, "f(x) = log(1+3x+x^2)/x", 2.108063708002621);

  test_function(&poly1, 1, 2, "f(x) = 1+1/sqrt(x) + 1/x", 2.52157430530614);

  printf("Tests complete!\n\n");
  print_total_error();
  print_ru();
  exit(0);
}
