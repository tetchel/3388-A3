
/*            PURPOSE : Matrix and vector dynamic memory management.

        PREREQUISITES : NONE

*/

#include <stdio.h>
#include <stdlib.h>

void nrerror(char error_text[])

{ void exit() ;

  fprintf(stderr,"Run-time error...\n") ;
  fprintf(stderr,"%s\n",error_text) ;
  exit(1) ;
}


float *vector(int nl, int nh)

{ float *v ;

  v = (float *)malloc((unsigned) (nh - nl + 1)*sizeof(float)) ;
  if (!v) {
    nrerror("VECTOR: allocation failure") ;
  }
  return (v - nl) ;
}


int *ivector(int nl, int nh)

{ int *v ;

  v = (int *)malloc((unsigned) (nh - nl + 1)*sizeof(int)) ;
  if (!v) {
    nrerror("IVECTOR: allocation failure") ;
  }
  return (v - nl) ;
}


double *dvector(int nl,int nh)

{ double *v ;

  v = (double *)malloc((unsigned) (nh - nl + 1)*sizeof(double)) ;
  if (!v) {
    nrerror("DVECTOR: allocation failure") ;
  }
  return (v - nl) ;
}


float **matrix(int nrl, int nrh, int ncl, int nch)

{ int i ;
  float **m ;

  m = (float **)malloc((unsigned) (nrh - nrl +1)*sizeof(float)) ;
  if (!m) {
    nrerror("MATRIX: allocation failure") ;
  }
  m -= nrl ;

  for (i = nrl ; i <= nrh ; i++) {
    m[i] = (float *)malloc((unsigned) (nch - ncl + 1)*sizeof(float)) ;
    if (!m[i]) {
      nrerror("MATRIX: allocation failure") ;
    }
    m[i] -= ncl ;
  }
  return (m) ;
}


double **dmatrix(int nrl, int nrh, int ncl, int nch)

{ int i ;
  double **m ;

  m = (double **)malloc((unsigned) (nrh - nrl +1)*sizeof(double)) ;
  if (!m) {
    nrerror("DMATRIX: allocation failure") ;
  }
  m -= nrl ;

  for (i = nrl ; i <= nrh ; i++) {
    m[i] = (double *)malloc((unsigned) (nch - ncl + 1)*sizeof(double)) ;
    if (!m[i]) {
      nrerror("DMATRIX: allocation failure") ;
    }
    m[i] -= ncl ;
  }
  return (m) ;
}


int **imatrix(int nrl, int nrh, int ncl, int nch)

{ int i ;
  int **m ;

  m = (int **)malloc((unsigned) (nrh - nrl +1)*sizeof(int)) ;
  if (!m) {
    nrerror("IMATRIX: allocation failure") ;
  }
  m -= nrl ;

  for (i = nrl ; i <= nrh ; i++) {
    m[i] = (int *)malloc((unsigned) (nch - ncl + 1)*sizeof(int)) ;
    if (!m[i]) {
      nrerror("IMATRIX: allocation failure") ;
    }
    m[i] -= ncl ;
  }
  return (m) ;
}


float **submatrix(float **a, int oldrl,int oldrh, int oldcl, int oldch, int newrl, int newcl)

{ int i, j ;

  float **m ;

  m = (float **)malloc((unsigned) (oldrh - oldrl + 1)*sizeof(float*)) ;
  if (!m) {
    nrerror("SUBMATRIX: allocation failure") ;
  }
  m -= newrl ;

  for (i = oldrl, j = newrl ; i <= oldrh ; i++, j++) {
    m[j] = a[i] + oldcl - newcl ;
  }
  return (m) ;
}


void free_vector(float *v, int nl, int nh)

{ free((char *) (v + nl)) ; }


void free_ivector(int *v, int nl, int nh)

{ free((char *) (v + nl)) ; }


void free_dvector(double *v, int nl, int nh)

{ free((char *) (v + nl)) ; }


void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)

{ int i ;

  for (i = nrh ; i >= nrl ; i--) {
    free((char *) (m[i] + ncl)) ;
  }
  free((char *) (m + nrl)) ;
}


void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)

{ int i ;

  for (i = nrh ; i >= nrl ; i--) {
    free((char *) (m[i] + ncl)) ;
  }
  free((char *) (m + nrl)) ;
}


void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)

{ int i ;

  for (i = nrh ; i >= nrl ; i--) {
    free((char *) (m[i] + ncl)) ;
  }
  free((char *) (m + nrl)) ;
}


void free_submatrix(float **b, int nrl, int nrh,int ncl, int nch)

{ free((char *) (b + nrl)) ; }


float **convert_matrix(float *a, int nrl, int nrh, int ncl, int nch)

{ int i, j, nrow, ncol ;
  float **m ;

  nrow = nrh - nrl + 1 ;
  ncol = nch - ncl + 1 ;

  m = (float **)malloc((unsigned) (nrow)*sizeof(float *)) ;
  if (!m) {
    nrerror("CONVERT_MATRIX: allocation failure") ;
  }
  m -= nrl ;
  for (i = 0, j = nrl ; i <= nrow - 1 ; i++, j++) {
    m[j] = a + ncol*i - ncl ;
  }
  return(m) ;
}

void free_convert_matrix(float **b, int nrl, int nrh, int ncl, int nch)

{ free((char *) (b + nrl)) ; }
