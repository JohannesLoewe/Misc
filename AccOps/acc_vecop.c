#include <stdlib.h>
#include <stdio.h>

#include "acc_vecop.h"

static void TwoSum_dbl(double a, double b, double *x, double *y)
{
  double A = a;

  *x = A + b;
  *y = *x - A;
  *y = (A - (*x - *y) + (b - *y));
}

static void TwoSum_flt(float a, float b, float *x, float *y)
{
  float A = a;

  *x = A + b;
  *y = *x - A;
  *y = (A - (*x - *y) + (b - *y));
}

#define SPLIT_DBL(a, a1, a2, c) (c ) = 134217729.0 * (a); \
                                (a1) = (c) - ((c) - (a)); \
                                (a2) = (a) - (a1);

#define SPLIT_FLT(a, a1, a2, c) (c ) = 4097.0f * (a); \
                                (a1) = (c) - ((c) - (a)); \
                                (a2) = (a) - (a1);

static void TwoProduct_flt(float a, float b, float *x, float *y)
{
  float work, h[2], g[2];

  SPLIT_FLT(a, h[0], g[0], work);
  SPLIT_FLT(b, h[1], g[1], work);

  *x = a * b;
  *y = g[0]*g[1] - (((*x - h[0]*h[1]) - h[1]*g[0]) - h[0]*g[1]);
}

static void TwoProduct_dbl(double a, double b, double *x, double *y)
{
  double work, h[2], g[2];

  SPLIT_DBL(a, h[0], g[0], work);
  SPLIT_DBL(b, h[1], g[1], work);

  *x = a * b;
  *y = g[0]*g[1] - (((*x - h[0]*h[1]) - h[1]*g[0]) - h[0]*g[1]);
}

float sum_flt(float *v, int n)
{
  float s = 0.0f;
  int k;

  for (k = 0; k < n; ++k)
    s += v[k];

  return s;
}

float sum_flt_acc(float *v, int n)
{
  float s1 = 0.0, s2 = 0.0, w;
  int k;

  for (k = 0; k < n; ++k)
  {
    TwoSum_flt(s1, v[k], &s1, &w);
    s2 += w;
  }

  return s1 + s2;
}

float sum_flt_ac2(float *v, int n)
{
  double s = 0.0;
  int k;

  for (k = 0; k < n; ++k)
    s += (double)v[k];

  return s;
}

double sum_dbl(double *v, int n)
{
  double s = 0.0;
  int k;

  for (k = 0; k < n; ++k)
    s += v[k];

  return s;
}

double sum_dbl_acc(double *v, int n)
{
  double s1 = 0.0, s2 = 0.0, w;
  int k;

  for (k = 0; k < n; ++k)
  {
    TwoSum_dbl(s1, v[k], &s1, &w);
    s2 += w;
  }

  return s1 + s2;
}

double sum_dbl_ac2(double *v, int n)
{
  long double s = 0.0;
  int k;

  for (k = 0; k < n; ++k)
    s += (long double)v[k];

  return s;
}

double dot_dbl(double *x, double *y, int n)
{
  double s = 0.0;
  int k;

  for (k = 0; k < n; ++k)
    s += x[k]*y[k];

  return s;
}

double dot_dbl_acc(double *x, double *y, int n)
{
  double s = 0.0, p = 0.0;
  int k;

  for (k = 0; k < n; ++k)
  {
    double hi, lo;

    TwoProduct_dbl(x[k], y[k], &hi, &lo);
    TwoSum_dbl(p, hi, &p, &hi);

    s += hi + lo;
  }

  return s + p;
}

double dot_dbl_ac2(double *x, double *y, int n)
{
  long double s = 0.0;
  int k;

  for (k = 0; k < n; ++k)
    s += (long double)x[k] * (long double)y[k];

  return s;
}

float dot_flt(float *x, float *y, int n)
{
  float s = 0.0;
  int k;

  for (k = 0; k < n; ++k)
    s += x[k]*y[k];

  return s;
}

float dot_flt_acc(float *x, float *y, int n)
{
  float s = 0.0, p = 0.0;
  int k;

  for (k = 0; k < n; ++k)
  {
    float hi, lo;

    TwoProduct_flt(x[k], y[k], &hi, &lo);
    TwoSum_flt(p, hi, &p, &hi);

    s += hi + lo;
  }

  return s + p;
}

float dot_flt_ac2(float *x, float *y, int n)
{
  double s = 0.0;
  int k;

  for (k = 0; k < n; ++k)
    s += (double)x[k] * (double)y[k];

  return s;
}