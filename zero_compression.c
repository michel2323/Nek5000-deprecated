#include <stdio.h>

void zero_compress_(double *in, int *n_, char *out, int *m_, int *ierr) {
  char *ctmp=0;
  long int n=(long int) *n_;
  long int li=0;
  short sj=0;
  short sk=0;
  for(li=0;li<n*8;li=li+1) out[li]=0;
  long int c=0;
  short b=7;
  for(li=0;li<n;li=li+1) {
    if(in[li]==0) {
      b=b-1;
      if(b<0) {
        b=7;
        c=c+1;
      }
    }
    else {
      out[c] |= 1 << b; 
      b=b-1;
      if(b<0) {
        b=7;
        c=c+1;
      }

      ctmp=(char*) &in[li]; 
      for(sj=7;sj>=0;sj=sj-1) {
        for(sk=7;sk>=0;sk=sk-1) {
          out[c] |= ((*(ctmp+sj) >> sk) & 1) << b; 
          b=b-1;
          if(b<0) {
            b=7;
            c=c+1;
          }
        }
      }
    }
  }
  *m_=(int) (c+1);
  /*printf("Compressed to %ld of %ld.\n",c+1,n*8);*/
}

void zero_decompress_(char *in, double *out, long int *n_, int *ierr) {
  int c=0;
  int b=7;
  int bit=0;
  long int n=*n_;
  long int li=0;
  short sj=0;
  short sk=0;
  double tmp=0;
  char *ctmp=0;
  for(li=0;li<n;li=li+1) {
    bit = (in[c] >> b) & 1;
    b=b-1;
    if(b<0) {
      b=7;
      c=c+1;
    }
    if(bit==0) {
      tmp=0;
    }
    else {
      tmp=0;
      ctmp=(char*)&tmp;
      for(sk=7;sk>=0;sk=sk-1) {
        for(sj=7;sj>=0;sj=sj-1) {
          *(ctmp+sk) |= ((in[c] >> b) & 1) << sj; 
          b=b-1;
          if(b<0) {
            b=7;
            c=c+1;
          }
        }
      }
    }
    out[li]=tmp;
  }
}

void zero_compress(double *in, int *n_, char *out, int *m_, int *ierr) {
zero_compress_(in, n_, out, m_, ierr);
}
void zero_decompress(char *in, double *out, long int *n_, int *ierr) {
zero_decompress_(in, out, n_, ierr);
}

