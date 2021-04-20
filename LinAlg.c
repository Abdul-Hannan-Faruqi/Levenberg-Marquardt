#include "LinAlg.h"

void vec_op(int n, double k1, double k2, double *a, double *b, double *c, bool op){
int i;

for(i=0;i<n;i++)
c[i] = k1*a[i] + ((!op)*k2*b[i]) - ((op)*k2*b[i]);
}

void GE(int n, double A[][n], double *B, double *X){

int row, col, i;
double factor;

for(col=0; col<n; col++){
   for(row=col+1; row<n; row++){
      factor = A[row][col]/A[col][col];
      for(i=col;i<n;i++)
         A[row][i] = A[row][i] - factor*A[col][i];
      B[row] = B[row] - factor*B[col];
   }
}

for(row = n-1; row>=0; row--){
double sum = 0;

   if(row!=n-1){
      for(col = row+1; col<n; col++)
         sum = sum + A[row][col]*X[col];
   }
   X[row] = (B[row] - sum)/A[row][row];
}

}

double norm(int n, double *v){
int i;
double l2n = 0;

for(i = 0; i<n; i++){
    l2n = l2n + pow(v[i],2);
}
l2n = sqrt(l2n);
return l2n;
}

double dot(int n, double *a, double *b){
int i;
double dp = 0;

for(i=0;i<n;i++)
dp = dp + a[i]*b[i];

return dp;
}

void mult(int m, int p, double A[m][p], int q, double B[p][q], double C[m][q]){
int i, j, k;
double sum = 0;

for(i=0;i<m;i++){
   for(j=0;j<q;j++){
   double sum = 0;
      for(k=0; k <p;k++)
         sum = sum + A[i][k]*B[k][j];
   C[i][j] = sum;
   }
}

}
