#include <stdlib.h>
#include <math.h>

void cox_stuff(double *x, double *y, int *ic, double *offset, int *xlen,
	       int *ylen, double *fail_times, double *s, int *d, int *dd,
	       int *nf, double *nn, double *nno) {
  short is_new;
  int i, j, k;

  *nf=0;
  for (i=0; i<*ylen; i++) {
    if (ic[i]==1) {
      is_new=1;
      for (j=0; j<*nf; j++)
	if (y[i]==fail_times[j])
	  is_new=0;
      if (is_new==1) {
	fail_times[*nf] = y[i];
	*nf += 1;
      }
    }
  }
  for (i=0; i<*nf; i++) {
    nn[i]=0;
    nno[i]=0;
    for (j=0; j<*ylen; j++) {
      if (y[j] >= fail_times[i]) {
	nn[i]++;
	nno[i] += exp(offset[i]);
      }
    }
  }
  for (i=0; i<*nf; i++) {
    d[i]=0;
    for (k=0; k<*xlen; k++)
      s[(*xlen * i)+k]=0;
    for (j=0; j<*ylen; j++) {
      if ((ic[j]==1) && (y[j]==fail_times[i])) {
	d[i]++;
	for (k=0; k<*xlen; k++) {
	  s[(*xlen * i)+k] += x[(*xlen * j)+k];
	}
      }
    }
  }
  for (i=0; i<*nf; i++) {
    dd[i]=0;
    for (j=0; j<*ylen; j++) {
      if ((ic[j]==1) && (y[j]==fail_times[i]))
	dd[j]=d[i];
    }
  }
}

void cox_scor(double *x, double *y, int *ic, double *offset, int *xlen,
	       int *ylen, double *fail_times, double *s, int *d, int *dd,
	       int *nf, double *nn, double *nno, double *scor) {
  int i, j, k;

  for (k=0; k<*xlen; k++)
    scor[k]=0;
  for (i=0; i<*nf; i++) {
    for (k=0; k<*xlen; k++)
      scor[k]+=s[(*xlen * i)+k];
    for (j=0; j<*ylen; j++)
      if (y[j] >= fail_times[i])
	for (k=0; k<*xlen; k++)
	  scor[k]-= (d[i] * x[(*xlen * j)+k] * exp(offset[j]))/nno[i];
  }
}

void cox_var(double *x, double *y, int *ic, double *offset, int *xlen,
	       int *ylen, double *fail_times, double *ss, int *d, int *dd,
	       int *nf, double *nn, double *nno, double *sd) {
  int i, j, k;
  double *s, *sx;

  s = calloc(*xlen, sizeof(double));
  sx = calloc(*xlen, sizeof(double));
  for (k=0; k<*xlen; k++)
    sd[k]=0;
  for (i=0; i<*nf; i++) {
    for (k=0; k<*xlen; k++)
      s[k] = sx[k] = 0;
    for (j=0; j<*ylen; j++)
      if (y[j] >= fail_times[i])
	for(k=0; k<*xlen; k++) {
	  sx[k]+= (x[(*xlen * j)+k] * exp(offset[j]))/nno[i];
	  s[k]+= (x[(*xlen * j)+k] * x[(*xlen * j)+k] *
		  exp(offset[j]))/nno[i];
	}
    for (k=0; k<*xlen; k++)
      sd[k]+= d[i] * (s[k] - sx[k]*sx[k]);
  }
}

void cox_func(double *x, double *y, int *icens, int *xlen, int *ylen,
	      int *nf, double *scor, double *sd) {
  double *offset, *fail_times, *s, *nn, *nno;
  int i, *d, *dd;

  offset = calloc(*ylen, sizeof(double));
  for(i=0; i<*ylen; i++)
    offset[i]=0;
  fail_times = calloc(*nf, sizeof(double));
  s = calloc(((*nf)*(*xlen)), sizeof(double));
  d = calloc(*nf, sizeof(int));
  dd = calloc(*ylen, sizeof(int));
  nn = calloc(*nf, sizeof(double));
  nno = calloc(*nf, sizeof(double));
  cox_stuff(x, y, icens, offset, xlen, ylen, fail_times, s, d, dd, nf,
	    nn, nno);
  cox_scor(x, y, icens, offset, xlen, ylen, fail_times, s, d, dd, nf,
	   nn, nno, scor);
  cox_var(x, y, icens, offset, xlen, ylen, fail_times, s, d, dd, nf,
	  nn, nno, sd);
}
