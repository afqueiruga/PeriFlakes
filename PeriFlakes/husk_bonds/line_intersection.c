#include "line_intersection.h"

#include "math.h"
#include <gsl/gsl_linalg.h>

void line_intersection_eval(int l_edge,
/*Inputs:*/
const real_t * restrict x,
const real_t * restrict pts,
/*Outputs:*/
real_t * restrict test)
{

if( (-(pts[0] - pts[2])*(pts[1] - x[1]) + (pts[0] - x[0])*(pts[1] - pts[3]))/((pts[0] - pts[2])*(x[1] - x[3]) - (pts[1] - pts[3])*(x[0] - x[2])) >= 0 && ((pts[0] - x[0])*(x[1] - x[3]) - (pts[1] - x[1])*(x[0] - x[2]))/((pts[0] - pts[2])*(x[1] - x[3]) - (pts[1] - pts[3])*(x[0] - x[2])) >= 0 && (-(pts[0] - pts[2])*(pts[1] - x[1]) + (pts[0] - x[0])*(pts[1] - pts[3]))/((pts[0] - pts[2])*(x[1] - x[3]) - (pts[1] - pts[3])*(x[0] - x[2])) <= 1 && ((pts[0] - x[0])*(x[1] - x[3]) - (pts[1] - x[1])*(x[0] - x[2]))/((pts[0] - pts[2])*(x[1] - x[3]) - (pts[1] - pts[3])*(x[0] - x[2])) <= 1 ) {
/* Evaluation of test */
test[0]= 1;
} else {
/* Evaluation of test */
test[0]= 0;
}

}

void line_intersection_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
   line_intersection_eval(l_edge,in, in+(int)(4), 
out);
}

void kmap_line_intersection0(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)((int)(2));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_line_intersection1(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 4;
  *n = (int)(1);
  if(edge) verts[0]=0;
}
void kmap_line_intersection2(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)(1));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(2) + 1*i];
    }
  }
}
kernel_t kernel_line_intersection = {
.nmap = 3,
.maps = {
kmap_line_intersection0,
kmap_line_intersection1,
kmap_line_intersection2
},

.ninp = 2,
.inp = {
{ .field_number=0, .map_num=0, .name="x" },
{ .field_number=1, .map_num=1, .name="pts" }
},

.noutp = 1,
.outp = {
{ .rank = 1, .nmap = 1, .map_nums = { 2 }, .name="test" }
},
.eval=line_intersection_eval_wr,
.name="line_intersection"
};

