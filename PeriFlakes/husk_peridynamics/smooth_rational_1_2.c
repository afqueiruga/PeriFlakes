#include "smooth_rational_1_2.h"

#include "math.h"
#include <gsl/gsl_linalg.h>

void smooth_rational_1_2_eval(int l_edge,
/*Inputs:*/
const real_t * restrict x,
const real_t * restrict delta,
const real_t * restrict y,
const real_t * restrict alpha,
/*Outputs:*/
real_t * restrict ys)
{
real_t m[1];{int i; for(i= 0;i<1;i++) m[i]=0.0;}

/* Evaluation of m */
m[0]= 10;

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of m */
m[0]+= alpha[i - 1]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   (1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0])/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
}
}

/* Evaluation of ys */
ys[0]= 10*y[0]/m[0];
ys[1]= 10*y[1]/m[0];

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of ys */
ys[0]+= alpha[i - 1]*y[2*i]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   (1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0])/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/m[0];
ys[1]+= alpha[i - 1]*y[2*i + 1]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   (1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0])/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/m[0];
}
}
}

void smooth_rational_1_2_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
   smooth_rational_1_2_eval(l_edge,in, in+(int)(l_edge + 1), in+(int)(l_edge + 1)+(int)(1), in+(int)(l_edge + 1)+(int)(1)+(int)(l_edge + 1), 
out);
}

void kmap_smooth_rational_1_20(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)((int)((1.0L/2.0L)*l_edge + 1.0L/2.0L));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_smooth_rational_1_21(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)(1));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_smooth_rational_1_22(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)((1.0L/2.0L)*l_edge - 1.0L/2.0L));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L) + 1*i];
    }
  }
}
void kmap_smooth_rational_1_23(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)((int)(1));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
kernel_t kernel_smooth_rational_1_2 = {
.nmap = 4,
.maps = {
kmap_smooth_rational_1_20,
kmap_smooth_rational_1_21,
kmap_smooth_rational_1_22,
kmap_smooth_rational_1_23
},

.ninp = 4,
.inp = {
{ .field_number=0, .map_num=0, .name="x" },
{ .field_number=1, .map_num=1, .name="delta" },
{ .field_number=2, .map_num=0, .name="y" },
{ .field_number=3, .map_num=2, .name="alpha" }
},

.noutp = 1,
.outp = {
{ .rank = 1, .nmap = 1, .map_nums = { 3 }, .name="ys" }
},
.eval=smooth_rational_1_2_eval_wr,
.name="smooth_rational_1_2"
};

