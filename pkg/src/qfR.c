/*----------------------------------------------------------------------
 *   gfR.c - R wrapper to Robert Davies' qf function to compute
 *           linear combination of chi-squared random variables
 *           (from http://www.robertnz.net/ftp/qf.tar.gz).
 *
 *         code (C)  Christine Choirat and Raffaello Seri
 *
 *    This file:
 *      created 2004-09-08
 *      last modified 2009-02-11
 *
 *--------------------------------------------------------------------*/

/* Include Robert Davies' C file (renamed from "qfc.c" to "qfc.h"
   so that we don't need a specific Makefile). */

#include "qfc.h"

/* Interface to R */
void qfR(real* lb1, real* nc1, int* n1, int* r1, real* sigma, real* c1,
         int* lim1, real* acc, real* trace, int* ifault, real* qfval)
{
  qfval[0] = qf(lb1, nc1, n1, r1[0], sigma[0], c1[0], lim1[0], acc[0], trace, ifault);
}

