#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "acc_vecop.h"

/*
gcc -O2 -ansi -pedantic -Wall -fno-strict-aliasing main.c acc_vecop.c
*/

#define N (1 << 21)
#define M (1 <<  5)
#define W 1.0/7

#define _COMPILE_ASSERT_SYMBOL_INNER(line, msg) \
        __COMPILE_ASSERT_ ## line ## _ ## msg
#define _COMPILE_ASSERT_SYMBOL(line, msg) \
        _COMPILE_ASSERT_SYMBOL_INNER(line, msg)
#define COMPILE_ASSERT(exp, msg) \
        typedef char _COMPILE_ASSERT_SYMBOL(__LINE__, msg) [(exp) ? 1 : -1]

static double rand_dbl()
{
  unsigned int r[2] = { 0 };
  double R;

  r[1] |= (          0x000003ffU  << 20);
  r[1] |= ((rand() & 0x000fffffU) <<  0);
  r[0] |= ((rand() & 0xffffffffU) <<  0);

  memcpy(&R, &r, 8);
  
  return 2*R-3;
}

int float_to_int(float A)
{
  int iA = *(int*)&A;

  return iA < 0 ? 0x80000000 - iA : iA;
}

long double_to_long(double A)
{
  long iA = *(long*)&A;

  return iA < 0 ? 0x8000000000000000 - iA : iA;
}

float int_to_float(int iA)
{
  if(iA < 0)
    iA = 0x80000000 - iA;
  return *(float*)&iA;
}

double long_to_double(long iA)
{
  if(iA < 0)
    iA = 0x8000000000000000 - iA;
  return *(double*)&iA;
}

unsigned int ulp_dist(float A, float B)
{
  int aInt = *(int*)&A;
  int bInt = *(int*)&B;

  if (aInt < 0)
    aInt = 0x80000000 - aInt;

  if (bInt < 0)
    bInt = 0x80000000 - bInt;

  return aInt > bInt ? aInt - bInt : bInt - aInt;
}

unsigned long ulp_dist_d(double A, double B)
{
  long aInt = *(long*)&A;
  long bInt = *(long*)&B;

  if (aInt < 0)
    aInt = 0x8000000000000000 - aInt;


  if (bInt < 0)
    bInt = 0x8000000000000000 - bInt;

  return aInt > bInt ? aInt - bInt : bInt - aInt;
}

int msb(unsigned long x)
{
  int r = 0;
  while(x != 0)
  {
    x >>= 1;
    r += 1;
  }
  return r;
}

static int test_flt(void)
{
  float sum_a1, sum_a2, sum_a3;
  float sum_b1, sum_b2, sum_b3;
  float r, s, t;
  int k, l;

  float *v = malloc(N * sizeof(float));
  float *w = malloc(N * sizeof(float));
  float *z = malloc(N * sizeof(float));

  for (k = 0; k < N; ++k)
  {
    v[k] = rand_dbl();
    w[k] = rand_dbl();
    z[k] = rand_dbl();
  }

  for(l = 0; l < 0; ++l)
  {
    r = sum_flt_acc(z, N);
    s = dot_flt_acc(v, w, N);
    t = dot_flt_acc(w, w, N);

    for (k = 0; k < N; ++k)
    {
      v[k] -= s/t * w[k];
      z[k] -= r/N;
    }
  }
  
  for (k = 0; k < 4*M; ++k)
  {
    sum_a1 = sum_flt(z, N);
    sum_a2 = sum_flt_acc(z, N);
    sum_a3 = sum_flt_ac2(z, N);
  }

  for (k = 0; k < M; ++k)
  {
    sum_b1 = dot_flt(v, w, N);
    sum_b2 = dot_flt_acc(v, w, N);
    sum_b3 = dot_flt_ac2(v, w, N);
  }

  printf("%22s %22s %22s\n", "vec_sum", "acc_sum", "vec_sum_d");
  printf("%+.15e %+.15e %+.15e %d %d\n", sum_a1, sum_a2, sum_a3,
         msb(ulp_dist(sum_a1, sum_a3)), msb(ulp_dist(sum_a2, sum_a3)));
  printf("%+.15e %+.15e %+.15e %d %d\n", sum_b1, sum_b2, sum_b3,
         msb(ulp_dist(sum_b1, sum_b3)), msb(ulp_dist(sum_b2, sum_b3)));

  free(z);
  free(w);
  free(v);

  return 0;
}

static int test_dbl(void)
{
  double sum_a1, sum_a2, sum_a3;
  double sum_b1, sum_b2, sum_b3;
  double r, s, t;
  int k, l;

  double *v = malloc(N * sizeof(double));
  double *w = malloc(N * sizeof(double));
  double *z = malloc(N * sizeof(double));

  for (k = 0; k < N; ++k)
  {
    v[k] = rand_dbl();
    w[k] = rand_dbl();
    z[k] = rand_dbl();
  }

  for(l = 0; l < 0; ++l)
  {
    r = sum_dbl_acc(z, N);
    s = dot_dbl_acc(v, w, N);
    t = dot_dbl_acc(w, w, N);

    for (k = 0; k < N; ++k)
    {
      v[k] -= s/t * w[k];
      z[k] -= r/N;
    }
  }

  for (k = 0; k < 4*M; ++k)
  {
    sum_a1 = sum_dbl(z, N);
    sum_a2 = sum_dbl_acc(z, N);
    sum_a3 = sum_dbl_ac2(z, N);
  }

  for (k = 0; k < M; ++k)
  {
    sum_b1 = dot_dbl(v, w, N);
    sum_b2 = dot_dbl_acc(v, w, N);
    sum_b3 = dot_dbl_ac2(v, w, N);
  }

  printf("%22s %22s %22s\n", "vec_sum", "acc_sum", "vec_sum_d");
  printf("%+.15e %+.15e %+.15e %d %d\n", sum_a1, sum_a2, sum_a3,
         msb(ulp_dist_d(sum_a1, sum_a3)), msb(ulp_dist_d(sum_a2, sum_a3)));
  printf("%+.15e %+.15e %+.15e %d %d\n", sum_b1, sum_b2, sum_b3,
         msb(ulp_dist_d(sum_b1, sum_b3)), msb(ulp_dist_d(sum_b2, sum_b3)));

  free(z);
  free(w);
  free(v);

  return 0;
}

int main(void)
{
  test_flt();
  test_dbl();
  return 0;
}