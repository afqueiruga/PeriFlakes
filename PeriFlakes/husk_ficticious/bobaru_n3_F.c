#include "bobaru_n3_F.h"

#include "math.h"
#include <gsl/gsl_linalg.h>

void bobaru_n3_F_eval(int l_edge,
/*Inputs:*/
const real_t * restrict load,
const real_t * restrict p_E,
const real_t * restrict p_nu,
const real_t * restrict y,
const real_t * restrict x,
/*Outputs:*/
real_t * restrict K,
real_t * restrict R)
{
/* Evaluation of R */
R[0]= -load[0]*((p_nu[0] + 1)*(x[0] - x[2])/(p_E[0]*sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + (-pow(p_nu[0], 2) + 1)*(x[1] - x[3])/(p_E[0]*sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2)))) - p_nu[0]*(x[3] - x[5] - y[3] + y[5])/sqrt(pow(-x[2] + x[4], 2) + pow(-x[3] + x[5], 2)) + (x[0] - x[2] - y[0] + y[2])/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
R[1]= -load[1]*((p_nu[0] + 1)*(x[1] - x[3])/(p_E[0]*sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + (-pow(p_nu[0], 2) + 1)*(x[0] - x[2])/(p_E[0]*sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2)))) - p_nu[0]*(x[2] - x[4] - y[2] + y[4])/sqrt(pow(-x[2] + x[4], 2) + pow(-x[3] + x[5], 2)) + (x[1] - x[3] - y[1] + y[3])/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
R[2]= 0;
R[3]= 0;
R[4]= 0;
R[5]= 0;

/* Evaluation of K */
K[0]= -1/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
K[1]= 0;
K[2]= pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), -1.0L/2.0L);
K[3]= p_nu[0]/sqrt(pow(-x[2] + x[4], 2) + pow(-x[3] + x[5], 2));
K[4]= 0;
K[5]= -p_nu[0]/sqrt(pow(-x[2] + x[4], 2) + pow(-x[3] + x[5], 2));
K[6]= 0;
K[7]= -1/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
K[8]= p_nu[0]/sqrt(pow(-x[2] + x[4], 2) + pow(-x[3] + x[5], 2));
K[9]= pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), -1.0L/2.0L);
K[10]= -p_nu[0]/sqrt(pow(-x[2] + x[4], 2) + pow(-x[3] + x[5], 2));
K[11]= 0;
K[12]= 0;
K[13]= 0;
K[14]= 0;
K[15]= 0;
K[16]= 0;
K[17]= 0;
K[18]= 0;
K[19]= 0;
K[20]= 0;
K[21]= 0;
K[22]= 0;
K[23]= 0;
K[24]= 0;
K[25]= 0;
K[26]= 0;
K[27]= 0;
K[28]= 0;
K[29]= 0;
K[30]= 0;
K[31]= 0;
K[32]= 0;
K[33]= 0;
K[34]= 0;
K[35]= 0;
}

void bobaru_n3_F_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
   bobaru_n3_F_eval(l_edge,in, in+(int)(2), in+(int)(2)+(int)(1), in+(int)(2)+(int)(1)+(int)(1), in+(int)(2)+(int)(1)+(int)(1)+(int)(6), 
out, out+(int)(36));
}

void kmap_bobaru_n3_F0(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)((int)(1));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_bobaru_n3_F1(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)(1));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_bobaru_n3_F2(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)((int)(3));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
kernel_t kernel_bobaru_n3_F = {
.nmap = 3,
.maps = {
kmap_bobaru_n3_F0,
kmap_bobaru_n3_F1,
kmap_bobaru_n3_F2
},

.ninp = 5,
.inp = {
{ .field_number=0, .map_num=0, .name="load" },
{ .field_number=1, .map_num=1, .name="p_E" },
{ .field_number=2, .map_num=1, .name="p_nu" },
{ .field_number=3, .map_num=2, .name="y" },
{ .field_number=4, .map_num=2, .name="x" }
},

.noutp = 2,
.outp = {
{ .rank = 2, .nmap = 1, .map_nums = { 2 }, .name="K" },
{ .rank = 1, .nmap = 1, .map_nums = { 2 }, .name="R" }
},
.eval=bobaru_n3_F_eval_wr,
.name="bobaru_n3_F"
};

