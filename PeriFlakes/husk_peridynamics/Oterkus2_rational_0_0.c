#include "Oterkus2_rational_0_0.h"

#include "math.h"
#include <gsl/gsl_linalg.h>

void Oterkus2_rational_0_0_eval(int l_edge,
/*Inputs:*/
const real_t * restrict x,
const real_t * restrict p_E,
const real_t * restrict p_nu,
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
m[0]+= alpha[i - 1]*sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
));
}
}

real_t a_integ[1];{int i; for(i= 0;i<1;i++) a_integ[i]=0.0;}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of a_integ */
a_integ[0]+= alpha[i - 1]*pow(-x[0] + x[2*i], 2)*pow(-x[1] + x[2*i + 1], 2)*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2));
}
}

real_t b_denom[1];{int i; for(i= 0;i<1;i++) b_denom[i]=0.0;}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of b_denom */
b_denom[0]+= alpha[i - 1]*pow(-x[0] + x[2*i], 2)*pow(-x[1] + x[2*i + 1], 2)*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2));
}
}

real_t theta[1];{int i; for(i= 0;i<1;i++) theta[i]=0.0;}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of theta */
theta[0]+= 2*alpha[i - 1]*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/m[0];
}
}

real_t theta_Dy0[2];{int i; for(i= 0;i<2;i++) theta_Dy0[i]=0.0;}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of theta_Dy0 */
theta_Dy0[0]+= 2*alpha[i - 1]*(y[0] - y[2*i])*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + 2*alpha[i - 1]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((-x[0] + x[2*i])*pow(-y[0] + y[2*i], 2)/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) - (-x[0] + x[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/m[0];
theta_Dy0[1]+= 2*alpha[i - 1]*(y[1] - y[2*i + 1])*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + 2*alpha[i - 1]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])*pow(-y[1] + y[2*i + 1], 2)/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) - (-x[1] + x[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/m[0];
}
}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of R */
R[0]+= 2*alpha[i - 1]*(-y[0] + y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
R[1]+= 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
/* Evaluation of R */
R[2*i]-= 2*alpha[i - 1]*(-y[0] + y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
R[2*i + 1]-= 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
}
}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of K */
K[0]+= 2*alpha[i - 1]*pow(-y[0] + y[2*i], 2)*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[0] + y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*pow(-y[0] + y[2*i], 2)/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) - (-x[0] + x[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)))/m[0] + (1.0L/2.0L)*p_E[0]*(y[0] - y[2*i])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) - 2*alpha[i - 1]*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*theta_Dy0[0]*(-y[0] + y[2*i])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[1]+= 2*alpha[i - 1]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[0] + y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])*pow(-y[1] + y[2*i + 1], 2)/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) - (-x[1] + x[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(y[1] - y[2*i + 1])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*theta_Dy0[1]*(-y[0] + y[2*i])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[l_edge + 1]+= 2*alpha[i - 1]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*pow(-y[0] + y[2*i], 2)/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) - (-x[0] + x[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)))/m[0] + (1.0L/2.0L)*p_E[0]*(y[0] - y[2*i])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*theta_Dy0[0]*(-y[1] + y[2*i + 1])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[l_edge + 2]+= 2*alpha[i - 1]*pow(-y[1] + y[2*i + 1], 2)*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])*pow(-y[1] + y[2*i + 1], 2)/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) - (-x[1] + x[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(y[1] - y[2*i + 1])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) - 2*alpha[i - 1]*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*theta_Dy0[1]*(-y[1] + y[2*i + 1])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
/* Evaluation of K */
K[2*i]+= 2*alpha[i - 1]*(-y[0] + y[2*i])*(y[0] - y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[0] + y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(y[0] - y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[0] + x[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(y[0] - y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)))/m[0] + (1.0L/2.0L)*p_E[0]*(-y[0] + y[2*i])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
K[2*i + 1]+= 2*alpha[i - 1]*(-y[0] + y[2*i])*(y[1] - y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[0] + y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(y[1] - y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])*(y[1] - y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-y[1] + y[2*i + 1])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
K[2*i + l_edge + 1]+= 2*alpha[i - 1]*(y[0] - y[2*i])*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(y[0] - y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[0] + x[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(y[0] - y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)))/m[0] + (1.0L/2.0L)*p_E[0]*(-y[0] + y[2*i])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
K[2*i + l_edge + 2]+= 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(y[1] - y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(y[1] - y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])*(y[1] - y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-y[1] + y[2*i + 1])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
/* Evaluation of K */
K[2*i*(l_edge + 1)]-= 2*alpha[i - 1]*pow(-y[0] + y[2*i], 2)*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[0] + y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*pow(-y[0] + y[2*i], 2)/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) - (-x[0] + x[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)))/m[0] + (1.0L/2.0L)*p_E[0]*(y[0] - y[2*i])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) - 2*alpha[i - 1]*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*theta_Dy0[0]*(-y[0] + y[2*i])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*i*l_edge + 2*i + 1]-= 2*alpha[i - 1]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[0] + y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])*pow(-y[1] + y[2*i + 1], 2)/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) - (-x[1] + x[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(y[1] - y[2*i + 1])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*theta_Dy0[1]*(-y[0] + y[2*i])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*i*l_edge + 2*i + l_edge + 1]-= 2*alpha[i - 1]*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*pow(-y[0] + y[2*i], 2)/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) - (-x[0] + x[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)))/m[0] + (1.0L/2.0L)*p_E[0]*(y[0] - y[2*i])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*theta_Dy0[0]*(-y[1] + y[2*i + 1])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*i*l_edge + 2*i + l_edge + 2]-= 2*alpha[i - 1]*pow(-y[1] + y[2*i + 1], 2)*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])*pow(-y[1] + y[2*i + 1], 2)/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) - (-x[1] + x[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(y[1] - y[2*i + 1])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) - 2*alpha[i - 1]*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*theta_Dy0[1]*(-y[1] + y[2*i + 1])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
/* Evaluation of K */
K[2*i*(l_edge + 1) + 2*i]-= 2*alpha[i - 1]*(-y[0] + y[2*i])*(y[0] - y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[0] + y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(y[0] - y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[0] + x[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(y[0] - y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)))/m[0] + (1.0L/2.0L)*p_E[0]*(-y[0] + y[2*i])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
K[2*i*l_edge + 4*i + 1]-= 2*alpha[i - 1]*(-y[0] + y[2*i])*(y[1] - y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[0] + y[2*i])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(y[1] - y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])*(y[1] - y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-y[1] + y[2*i + 1])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
K[2*i*l_edge + 4*i + l_edge + 1]-= 2*alpha[i - 1]*(y[0] - y[2*i])*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(y[0] - y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[0] + x[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(y[0] - y[2*i])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)))/m[0] + (1.0L/2.0L)*p_E[0]*(-y[0] + y[2*i])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
K[2*i*l_edge + 4*i + l_edge + 2]-= 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(y[1] - y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L) + 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])*(y[1] - y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])*(y[1] - y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*pow(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-y[1] + y[2*i + 1])/(b_denom[0]*(2.0*p_nu[0] + 2.0)*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)) + 2*alpha[i - 1]*(theta[0]*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))/m[0] + (1.0L/2.0L)*p_E[0]*(-sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) + sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)))/(b_denom[0]*(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2));
{
int j;
for(j=(int)(1);j<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);j++) {
/* Evaluation of K */
K[2*j]+= 2*alpha[i - 1]*(-y[0] + y[2*i])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*(2*alpha[j - 1]*(-y[0] + y[2*j])*((-x[0] + x[2*j])*(-y[0] + y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + 2*alpha[j - 1]*(-sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) + sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2)))*((-x[0] + x[2*j])*(-y[0] + y[2*j])*(y[0] - y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[0] + x[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(y[0] - y[2*j])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/m[0])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*j + 1]+= 2*alpha[i - 1]*(-y[0] + y[2*i])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*(2*alpha[j - 1]*(-y[1] + y[2*j + 1])*((-x[0] + x[2*j])*(-y[0] + y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + 2*alpha[j - 1]*(-sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) + sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2)))*((-x[0] + x[2*j])*(-y[0] + y[2*j])*(y[1] - y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])*(y[1] - y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/m[0])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*j + l_edge + 1]+= 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*(2*alpha[j - 1]*(-y[0] + y[2*j])*((-x[0] + x[2*j])*(-y[0] + y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + 2*alpha[j - 1]*(-sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) + sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2)))*((-x[0] + x[2*j])*(-y[0] + y[2*j])*(y[0] - y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[0] + x[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(y[0] - y[2*j])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/m[0])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*j + l_edge + 2]+= 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*(2*alpha[j - 1]*(-y[1] + y[2*j + 1])*((-x[0] + x[2*j])*(-y[0] + y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + 2*alpha[j - 1]*(-sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) + sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2)))*((-x[0] + x[2*j])*(-y[0] + y[2*j])*(y[1] - y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])*(y[1] - y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/m[0])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
/* Evaluation of K */
K[2*i*(l_edge + 1) + 2*j]-= 2*alpha[i - 1]*(-y[0] + y[2*i])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*(2*alpha[j - 1]*(-y[0] + y[2*j])*((-x[0] + x[2*j])*(-y[0] + y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + 2*alpha[j - 1]*(-sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) + sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2)))*((-x[0] + x[2*j])*(-y[0] + y[2*j])*(y[0] - y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[0] + x[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(y[0] - y[2*j])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/m[0])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*i*l_edge + 2*i + 2*j + 1]-= 2*alpha[i - 1]*(-y[0] + y[2*i])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*(2*alpha[j - 1]*(-y[1] + y[2*j + 1])*((-x[0] + x[2*j])*(-y[0] + y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + 2*alpha[j - 1]*(-sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) + sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2)))*((-x[0] + x[2*j])*(-y[0] + y[2*j])*(y[1] - y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])*(y[1] - y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/m[0])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*i*l_edge + 2*i + 2*j + l_edge + 1]-= 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*(2*alpha[j - 1]*(-y[0] + y[2*j])*((-x[0] + x[2*j])*(-y[0] + y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + 2*alpha[j - 1]*(-sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) + sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2)))*((-x[0] + x[2*j])*(-y[0] + y[2*j])*(y[0] - y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[0] + x[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(y[0] - y[2*j])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/m[0])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
K[2*i*l_edge + 2*i + 2*j + l_edge + 2]-= 2*alpha[i - 1]*(-y[1] + y[2*i + 1])*(-2*a_integ[0]*p_E[0]/(b_denom[0]*(2.0*p_nu[0] + 2.0)) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((-x[0] + x[2*i])*(-y[0] + y[2*i])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))) + (-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])/(sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2))))*(2*alpha[j - 1]*(-y[1] + y[2*j + 1])*((-x[0] + x[2*j])*(-y[0] + y[2*j])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))) + 2*alpha[j - 1]*(-sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) + sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2)))*((-x[0] + x[2*j])*(-y[0] + y[2*j])*(y[1] - y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*j + 1])*(-y[1] + y[2*j + 1])*(y[1] - y[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*pow(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2), 3.0L/2.0L)) + (-x[1] + x[2*j + 1])/(sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))*sqrt(pow(-y[0] + y[2*j], 2) + pow(-y[1] + y[2*j + 1], 2))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/m[0])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   1
)
: (
   0.0
))/(m[0]*sqrt(pow(-y[0] + y[2*i], 2) + pow(-y[1] + y[2*i + 1], 2)));
}
}
}
}
}

void Oterkus2_rational_0_0_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
   Oterkus2_rational_0_0_eval(l_edge,in, in+(int)(l_edge + 1), in+(int)(l_edge + 1)+(int)(1), in+(int)(l_edge + 1)+(int)(1)+(int)(1), in+(int)(l_edge + 1)+(int)(1)+(int)(1)+(int)(1), in+(int)(l_edge + 1)+(int)(1)+(int)(1)+(int)(1)+(int)(l_edge + 1), 
out, out+(int)(pow(l_edge + 1, 2)));
}

void kmap_Oterkus2_rational_0_00(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)((int)((1.0L/2.0L)*l_edge + 1.0L/2.0L));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_Oterkus2_rational_0_01(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)(1));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_Oterkus2_rational_0_02(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)((1.0L/2.0L)*l_edge - 1.0L/2.0L));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L) + 1*i];
    }
  }
}
kernel_t kernel_Oterkus2_rational_0_0 = {
.nmap = 3,
.maps = {
kmap_Oterkus2_rational_0_00,
kmap_Oterkus2_rational_0_01,
kmap_Oterkus2_rational_0_02
},

.ninp = 6,
.inp = {
{ .field_number=0, .map_num=0, .name="x" },
{ .field_number=1, .map_num=1, .name="p_E" },
{ .field_number=2, .map_num=1, .name="p_nu" },
{ .field_number=3, .map_num=1, .name="delta" },
{ .field_number=4, .map_num=0, .name="y" },
{ .field_number=5, .map_num=2, .name="alpha" }
},

.noutp = 2,
.outp = {
{ .rank = 2, .nmap = 1, .map_nums = { 0 }, .name="K" },
{ .rank = 1, .nmap = 1, .map_nums = { 0 }, .name="R" }
},
.eval=Oterkus2_rational_0_0_eval_wr,
.name="Oterkus2_rational_0_0"
};

