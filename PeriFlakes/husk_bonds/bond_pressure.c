#include "bond_pressure.h"

#include "math.h"
#include <gsl/gsl_linalg.h>

void bond_pressure_eval(int l_edge,
/*Inputs:*/
const real_t * restrict x,
const real_t * restrict p,
/*Outputs:*/
real_t * restrict R)
{
/* Evaluation of R */
R[0]= p[0]*(-x[0] + x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
R[1]= p[0]*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
R[2]= -p[0]*(-x[0] + x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
R[3]= -p[0]*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
}

void bond_pressure_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
   bond_pressure_eval(l_edge,in, in+(int)(4), 
out);
}

void kmap_bond_pressure0(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)((int)(2));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_bond_pressure1(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)(1));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(2) + 1*i];
    }
  }
}
kernel_t kernel_bond_pressure = {
.nmap = 2,
.maps = {
kmap_bond_pressure0,
kmap_bond_pressure1
},

.ninp = 2,
.inp = {
{ .field_number=0, .map_num=0, .name="x" },
{ .field_number=1, .map_num=1, .name="p" }
},

.noutp = 1,
.outp = {
{ .rank = 1, .nmap = 1, .map_nums = { 0 }, .name="R" }
},
.eval=bond_pressure_eval_wr,
.name="bond_pressure"
};

