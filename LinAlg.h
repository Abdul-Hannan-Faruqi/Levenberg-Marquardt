#ifndef LinAlg_H_
#define LinAlg_H_

#include <math.h>
#include <stdbool.h>


void vec_op(int n, double k1, double k2, double *a, double *b, double *c, bool op);
void GE(int n, double A[][n], double *B, double *X);
double norm(int n, double *v);
double dot(int n, double *a, double *b);
void mult(int m, int p, double A[m][p], int q, double B[p][q], double C[m][q]);

#endif
