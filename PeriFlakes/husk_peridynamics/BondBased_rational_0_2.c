#include "BondBased_rational_0_2.h"

#include "math.h"
#include <gsl/gsl_linalg.h>

void BondBased_rational_0_2_eval(int l_edge,
/*Inputs:*/
const real_t * restrict x,
const real_t * restrict p_Vol,
const real_t * restrict p_E,
const real_t * restrict delta,
const real_t * restrict y,
const real_t * restrict alpha,
/*Outputs:*/
real_t * restrict K,
real_t * restrict R)
{
real_t m[1];{int i; for(i= 0;i<1;i++) m[i]=0.0;}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of m */
m[0]+= alpha[i - 1]*(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
}
}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of R */
R[0]+= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
R[1]+= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[1] + y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
/* Evaluation of R */
R[2*i]-= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
R[2*i + 1]-= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[1] + y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
}
}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of K */
K[0]+= alpha[i - 1]*p_E[0]*p_Vol[0]*pow(-y[0] + y[2*i], 2)*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(y[0] - y[2*i])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) - alpha[i - 1]*p_E[0]*p_Vol[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[1]+= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(y[1] - y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[l_edge + 1]+= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(y[0] - y[2*i])*(-y[1] + y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[l_edge + 2]+= alpha[i - 1]*p_E[0]*p_Vol[0]*pow(-y[1] + y[2*i + 1], 2)*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[1] + y[2*i + 1])*(y[1] - y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) - alpha[i - 1]*p_E[0]*p_Vol[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
/* Evaluation of K */
K[2*i]+= alpha[i - 1]*p_E[0]*p_Vol[0]*pow(-y[0] + y[2*i], 2)*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(y[0] - y[2*i])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*i + 1]+= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(y[1] - y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L));
K[2*i + l_edge + 1]+= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + alpha[i - 1]*p_E[0]*p_Vol[0]*(y[0] - y[2*i])*(-y[1] + y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L));
K[2*i + l_edge + 2]+= alpha[i - 1]*p_E[0]*p_Vol[0]*pow(-y[1] + y[2*i + 1], 2)*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[1] + y[2*i + 1])*(y[1] - y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
/* Evaluation of K */
K[2*i*(l_edge + 1)]-= alpha[i - 1]*p_E[0]*p_Vol[0]*pow(-y[0] + y[2*i], 2)*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(y[0] - y[2*i])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) - alpha[i - 1]*p_E[0]*p_Vol[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*i*l_edge + 2*i + 1]-= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(y[1] - y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*i*l_edge + 2*i + l_edge + 1]-= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(y[0] - y[2*i])*(-y[1] + y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*i*l_edge + 2*i + l_edge + 2]-= alpha[i - 1]*p_E[0]*p_Vol[0]*pow(-y[1] + y[2*i + 1], 2)*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[1] + y[2*i + 1])*(y[1] - y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) - alpha[i - 1]*p_E[0]*p_Vol[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
/* Evaluation of K */
K[2*i*(l_edge + 1) + 2*i]-= alpha[i - 1]*p_E[0]*p_Vol[0]*pow(-y[0] + y[2*i], 2)*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(y[0] - y[2*i])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*i*l_edge + 4*i + 1]-= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(y[1] - y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L));
K[2*i*l_edge + 4*i + l_edge + 1]-= alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + alpha[i - 1]*p_E[0]*p_Vol[0]*(y[0] - y[2*i])*(-y[1] + y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L));
K[2*i*l_edge + 4*i + l_edge + 2]-= alpha[i - 1]*p_E[0]*p_Vol[0]*pow(-y[1] + y[2*i + 1], 2)*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-y[1] + y[2*i + 1])*(y[1] - y[2*i + 1])*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + alpha[i - 1]*p_E[0]*p_Vol[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1.0/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(m[0]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
}
}
}

void BondBased_rational_0_2_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
   BondBased_rational_0_2_eval(l_edge,in, in+(int)(l_edge + 1), in+(int)(l_edge + 1)+(int)(1), in+(int)(l_edge + 1)+(int)(1)+(int)(1), in+(int)(l_edge + 1)+(int)(1)+(int)(1)+(int)(1), in+(int)(l_edge + 1)+(int)(1)+(int)(1)+(int)(1)+(int)(l_edge + 1), 
out, out+(int)(pow(l_edge + 1, 2)));
}

void kmap_BondBased_rational_0_20(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)((int)((1.0L/2.0L)*l_edge + 1.0L/2.0L));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_BondBased_rational_0_21(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)(1));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_BondBased_rational_0_22(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)((1.0L/2.0L)*l_edge - 1.0L/2.0L));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L) + 1*i];
    }
  }
}
kernel_t kernel_BondBased_rational_0_2 = {
.nmap = 3,
.maps = {
kmap_BondBased_rational_0_20,
kmap_BondBased_rational_0_21,
kmap_BondBased_rational_0_22
},

.ninp = 6,
.inp = {
{ .field_number=0, .map_num=0, .name="x" },
{ .field_number=1, .map_num=1, .name="p_Vol" },
{ .field_number=2, .map_num=1, .name="p_E" },
{ .field_number=3, .map_num=1, .name="delta" },
{ .field_number=4, .map_num=0, .name="y" },
{ .field_number=5, .map_num=2, .name="alpha" }
},

.noutp = 2,
.outp = {
{ .rank = 2, .nmap = 1, .map_nums = { 0 }, .name="K" },
{ .rank = 1, .nmap = 1, .map_nums = { 0 }, .name="R" }
},
.eval=BondBased_rational_0_2_eval_wr,
.name="BondBased_rational_0_2"
};

