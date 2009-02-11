/*----------------------------------------------------------------------
 *   discrepancies.c - C definition of the three different types of
 *                     discrepancies
 *
 *         code (C) Christine Choirat and Raffaello Seri
 *
 *    This file:
 *      created 2004-09-08
 *      last modified 2009-02-11
 *
 *--------------------------------------------------------------------*/

#include <math.h>

void fGDStar(int* piD, double* pvXk, double* pvXl, double* pdResult) {
  int j;
  double dResult;
  dResult = 1.0;
  for (j = 0;j < *piD; j++)
    dResult = dResult * (2 - (pvXk[j] < pvXl[j] ? pvXl[j] : pvXk[j]));
  pdResult[0] = dResult;
}

void fGDSymm(int* piD, double* pvXk, double* pvXl, double* pdResult) {
  int j;
  double dResult;
  dResult = 1.0;
  for (j = 0;j < *piD; j++)
    dResult = dResult * 2 * (1 - fabs(pvXk[j] - pvXl[j]));
  pdResult[0] = dResult;
}

void fGDCent(int* piD, double* pvXk, double* pvXl, double* pdResult) {
  int j;
  double dResult;
  dResult = 1.0;
  for (j = 0;j < *piD; j++)
    dResult = dResult * (1 + .5 * fabs(pvXk[j] - .5)
			 + .5 * fabs(pvXl[j] - .5)
			 - .5 * fabs(pvXk[j] - pvXl[j]));
  pdResult[0] = dResult;
}

void gGDStar(int* piD, double* pvXk, double* pdResult) {
  int j;
  double dResult;
  dResult = 1.0;
  for (j = 0;j < *piD; j++)
    dResult = dResult * (1.5 - .5 * (pvXk[j]) * (pvXk[j])) ;
  pdResult[0] = dResult;
}

void gGDSymm(int* piD, double* pvXk, double* pdResult) {
  int j;
  double dResult;
  dResult = 1.0;
  for (j = 0;j < *piD; j++)
    dResult = dResult * (1 + 2 * pvXk[j] - 2 * (pvXk[j]) * (pvXk[j])) ;
  pdResult[0] = dResult;
}

void gGDCent(int* piD, double* pvXk, double* pdResult) {
  int j;
  double dResult;
  dResult = 1.0;
  for (j = 0;j < *piD; j++)
    dResult = dResult * (1 + .5 * fabs(pvXk[j] - .5) - .5 * (pvXk[j] - .5) * (pvXk[j] - .5));
  pdResult[0] = dResult;
}

