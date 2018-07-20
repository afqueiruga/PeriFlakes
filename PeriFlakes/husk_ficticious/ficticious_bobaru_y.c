#include "ficticious_bobaru_y.h"

#include "math.h"
#include <gsl/gsl_linalg.h>

void ficticious_bobaru_y_eval(int l_edge,
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
R[0]= -load[0]*(-pow(p_nu[0], 2) + 1)/p_E[0] + p_nu[0]*(-x[5] + x[7] + y[5] - y[7])/fabs(x[4] - x[6]) + (-x[0] + x[2] + y[0] - y[2])/fabs(x[1] - x[3]);
R[1]= -load[1]*(p_nu[0] + 1)/p_E[0] + p_nu[0]*(-x[4] + x[6] + y[4] - y[6])/fabs(x[4] - x[6]) + (-x[1] + x[3] + y[1] - y[3])/fabs(x[1] - x[3]);
R[2]= 0;
R[3]= 0;
R[4]= 0;
R[5]= 0;
R[6]= 0;
R[7]= 0;

/* Evaluation of K */
K[0]= 1.0/fabs(x[1] - x[3]);
K[1]= 0;
K[2]= -1/fabs(x[1] - x[3]);
K[3]= 0;
K[4]= 0;
K[5]= p_nu[0]/fabs(x[4] - x[6]);
K[6]= 0;
K[7]= -p_nu[0]/fabs(x[4] - x[6]);
K[8]= 0;
K[9]= 1.0/fabs(x[1] - x[3]);
K[10]= 0;
K[11]= -1/fabs(x[1] - x[3]);
K[12]= p_nu[0]/fabs(x[4] - x[6]);
K[13]= 0;
K[14]= -p_nu[0]/fabs(x[4] - x[6]);
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
K[36]= 0;
K[37]= 0;
K[38]= 0;
K[39]= 0;
K[40]= 0;
K[41]= 0;
K[42]= 0;
K[43]= 0;
K[44]= 0;
K[45]= 0;
K[46]= 0;
K[47]= 0;
K[48]= 0;
K[49]= 0;
K[50]= 0;
K[51]= 0;
K[52]= 0;
K[53]= 0;
K[54]= 0;
K[55]= 0;
K[56]= 0;
K[57]= 0;
K[58]= 0;
K[59]= 0;
K[60]= 0;
K[61]= 0;
K[62]= 0;
K[63]= 0;
}

void ficticious_bobaru_y_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
   ficticious_bobaru_y_eval(l_edge,in, in+(int)(2), in+(int)(2)+(int)(1), in+(int)(2)+(int)(1)+(int)(1), in+(int)(2)+(int)(1)+(int)(1)+(int)(8), 
out, out+(int)(64));
}

void kmap_ficticious_bobaru_y0(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)(1);
  if(edge) verts[0]=0;
}
void kmap_ficticious_bobaru_y1(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)(1));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_ficticious_bobaru_y2(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)((int)(4));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
kernel_t kernel_ficticious_bobaru_y = {
.nmap = 3,
.maps = {
kmap_ficticious_bobaru_y0,
kmap_ficticious_bobaru_y1,
kmap_ficticious_bobaru_y2
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
.eval=ficticious_bobaru_y_eval_wr,
.name="ficticious_bobaru_y"
};

