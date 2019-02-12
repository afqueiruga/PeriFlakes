#include "Fbased_rational_2_1.h"

#include "math.h"
#include <gsl/gsl_linalg.h>

void Fbased_rational_2_1_eval(int l_edge,
/*Inputs:*/
const real_t * restrict x,
const real_t * restrict p_Vol,
const real_t * restrict p_E,
const real_t * restrict p_nu,
const real_t * restrict delta,
const real_t * restrict y,
const real_t * restrict alpha,
/*Outputs:*/
real_t * restrict K,
real_t * restrict R)
{
real_t mK[4];{int i; for(i= 0;i<4;i++) mK[i]=0.0;}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of mK */
mK[0]+= alpha[i - 1]*p_Vol[0]*pow(-x[0] + x[2*i], 2)*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
mK[1]+= alpha[i - 1]*p_Vol[0]*(-x[0] + x[2*i])*(-x[1] + x[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
mK[2]+= alpha[i - 1]*p_Vol[0]*(-x[0] + x[2*i])*(-x[1] + x[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
mK[3]+= alpha[i - 1]*p_Vol[0]*pow(-x[1] + x[2*i + 1], 2)*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
}
}

real_t mKi[4];{int i; for(i= 0;i<4;i++) mKi[i]=0.0;}

/* Evaluation of mKi */
mKi[0]= 1.0/mK[0] + mK[1]*mK[2]/(pow(mK[0], 2)*(mK[3] - mK[1]*mK[2]/mK[0]));
mKi[1]= -mK[1]/(mK[0]*(mK[3] - mK[1]*mK[2]/mK[0]));
mKi[2]= -mK[2]/(mK[0]*(mK[3] - mK[1]*mK[2]/mK[0]));
mKi[3]= 1.0/(mK[3] - mK[1]*mK[2]/mK[0]);

real_t mN[4];{int i; for(i= 0;i<4;i++) mN[i]=0.0;}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of mN */
mN[0]+= alpha[i - 1]*p_Vol[0]*(-x[0] + x[2*i])*(-y[0] + y[2*i])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
mN[1]+= alpha[i - 1]*p_Vol[0]*(-x[1] + x[2*i + 1])*(-y[0] + y[2*i])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
mN[2]+= alpha[i - 1]*p_Vol[0]*(-x[0] + x[2*i])*(-y[1] + y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
mN[3]+= alpha[i - 1]*p_Vol[0]*(-x[1] + x[2*i + 1])*(-y[1] + y[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
}
}

real_t N00_dy0[2];{int i; for(i= 0;i<2;i++) N00_dy0[i]=0.0;}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of N00_dy0 */
N00_dy0[0]+= -alpha[i - 1]*p_Vol[0]*(-x[0] + x[2*i])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
N00_dy0[1]+= 0;
}
}

real_t N01_dy0[2];{int i; for(i= 0;i<2;i++) N01_dy0[i]=0.0;}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of N01_dy0 */
N01_dy0[0]+= 0;
N01_dy0[1]+= -alpha[i - 1]*p_Vol[0]*(-x[0] + x[2*i])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
}
}

real_t N10_dy0[2];{int i; for(i= 0;i<2;i++) N10_dy0[i]=0.0;}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of N10_dy0 */
N10_dy0[0]+= -alpha[i - 1]*p_Vol[0]*(-x[1] + x[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
));
N10_dy0[1]+= 0;
}
}

real_t N11_dy0[2];{int i; for(i= 0;i<2;i++) N11_dy0[i]=0.0;}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of N11_dy0 */
N11_dy0[0]+= 0;
N11_dy0[1]+= -alpha[i - 1]*p_Vol[0]*(-x[1] + x[2*i + 1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
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
R[0]+= p_Vol[0]*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*p_E[0]*(1.0*mKi[0]*mN[0] + 1.0*mKi[2]*mN[1] - 1)/(2.0*p_nu[0] + 2.0) + (-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*(1.0*mKi[0]*mN[0] + 1.0*mKi[1]*mN[2] + 1.0*mKi[2]*mN[1] + 1.0*mKi[3]*mN[3] - 2))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 2.0*alpha[i - 1]*mKi[1]*p_E[0]*(0.5*mKi[0]*mN[2] + 0.5*mKi[1]*mN[0] + 0.5*mKi[2]*mN[3] + 0.5*mKi[3]*mN[1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(alpha[i - 1]*mKi[2]*(2.0*p_E[0]*(1.0*mKi[0]*mN[0] + 1.0*mKi[2]*mN[1] - 1)/(2.0*p_nu[0] + 2.0) + (-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*(1.0*mKi[0]*mN[0] + 1.0*mKi[1]*mN[2] + 1.0*mKi[2]*mN[1] + 1.0*mKi[3]*mN[3] - 2))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 2.0*alpha[i - 1]*mKi[3]*p_E[0]*(0.5*mKi[0]*mN[2] + 0.5*mKi[1]*mN[0] + 0.5*mKi[2]*mN[3] + 0.5*mKi[3]*mN[1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)));
R[1]+= p_Vol[0]*((-x[0] + x[2*i])*(2.0*alpha[i - 1]*mKi[0]*p_E[0]*(0.5*mKi[0]*mN[2] + 0.5*mKi[1]*mN[0] + 0.5*mKi[2]*mN[3] + 0.5*mKi[3]*mN[1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*p_E[0]*(1.0*mKi[1]*mN[2] + 1.0*mKi[3]*mN[3] - 1)/(2.0*p_nu[0] + 2.0) + (-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*(1.0*mKi[0]*mN[0] + 1.0*mKi[1]*mN[2] + 1.0*mKi[2]*mN[1] + 1.0*mKi[3]*mN[3] - 2))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(2.0*alpha[i - 1]*mKi[2]*p_E[0]*(0.5*mKi[0]*mN[2] + 0.5*mKi[1]*mN[0] + 0.5*mKi[2]*mN[3] + 0.5*mKi[3]*mN[1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*p_E[0]*(1.0*mKi[1]*mN[2] + 1.0*mKi[3]*mN[3] - 1)/(2.0*p_nu[0] + 2.0) + (-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*(1.0*mKi[0]*mN[0] + 1.0*mKi[1]*mN[2] + 1.0*mKi[2]*mN[1] + 1.0*mKi[3]*mN[3] - 2))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))));
/* Evaluation of R */
R[2*i]-= p_Vol[0]*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*p_E[0]*(1.0*mKi[0]*mN[0] + 1.0*mKi[2]*mN[1] - 1)/(2.0*p_nu[0] + 2.0) + (-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*(1.0*mKi[0]*mN[0] + 1.0*mKi[1]*mN[2] + 1.0*mKi[2]*mN[1] + 1.0*mKi[3]*mN[3] - 2))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 2.0*alpha[i - 1]*mKi[1]*p_E[0]*(0.5*mKi[0]*mN[2] + 0.5*mKi[1]*mN[0] + 0.5*mKi[2]*mN[3] + 0.5*mKi[3]*mN[1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(alpha[i - 1]*mKi[2]*(2.0*p_E[0]*(1.0*mKi[0]*mN[0] + 1.0*mKi[2]*mN[1] - 1)/(2.0*p_nu[0] + 2.0) + (-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*(1.0*mKi[0]*mN[0] + 1.0*mKi[1]*mN[2] + 1.0*mKi[2]*mN[1] + 1.0*mKi[3]*mN[3] - 2))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 2.0*alpha[i - 1]*mKi[3]*p_E[0]*(0.5*mKi[0]*mN[2] + 0.5*mKi[1]*mN[0] + 0.5*mKi[2]*mN[3] + 0.5*mKi[3]*mN[1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)));
R[2*i + 1]-= p_Vol[0]*((-x[0] + x[2*i])*(2.0*alpha[i - 1]*mKi[0]*p_E[0]*(0.5*mKi[0]*mN[2] + 0.5*mKi[1]*mN[0] + 0.5*mKi[2]*mN[3] + 0.5*mKi[3]*mN[1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*p_E[0]*(1.0*mKi[1]*mN[2] + 1.0*mKi[3]*mN[3] - 1)/(2.0*p_nu[0] + 2.0) + (-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*(1.0*mKi[0]*mN[0] + 1.0*mKi[1]*mN[2] + 1.0*mKi[2]*mN[1] + 1.0*mKi[3]*mN[3] - 2))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(2.0*alpha[i - 1]*mKi[2]*p_E[0]*(0.5*mKi[0]*mN[2] + 0.5*mKi[1]*mN[0] + 0.5*mKi[2]*mN[3] + 0.5*mKi[3]*mN[1])*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*p_E[0]*(1.0*mKi[1]*mN[2] + 1.0*mKi[3]*mN[3] - 1)/(2.0*p_nu[0] + 2.0) + (-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*(1.0*mKi[0]*mN[0] + 1.0*mKi[1]*mN[2] + 1.0*mKi[2]*mN[1] + 1.0*mKi[3]*mN[3] - 2))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))));
}
}

{
int i;
for(i=(int)(1);i<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);i++) {
/* Evaluation of K */
K[0]+= N00_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[1], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[2]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N01_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(alpha[i - 1]*mKi[2]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[3], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0))) + N10_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N11_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))));
K[1]+= N00_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[1], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[2]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N01_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(alpha[i - 1]*mKi[2]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[3], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0))) + N10_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N11_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))));
K[l_edge + 1]+= N00_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0))) + N01_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N10_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*pow(mKi[0], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N11_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*pow(mKi[2], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))));
K[l_edge + 2]+= N00_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0))) + N01_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N10_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*pow(mKi[0], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N11_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*pow(mKi[2], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))));
/* Evaluation of K */
K[2*i]+= 0;
K[2*i + 1]+= 0;
K[2*i + l_edge + 1]+= 0;
K[2*i + l_edge + 2]+= 0;
/* Evaluation of K */
K[2*i*(l_edge + 1)]-= N00_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[1], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[2]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N01_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(alpha[i - 1]*mKi[2]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[3], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0))) + N10_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N11_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))));
K[2*i*l_edge + 2*i + 1]-= N00_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[1], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[2]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N01_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(alpha[i - 1]*mKi[2]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[3], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0))) + N10_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N11_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))));
K[2*i*l_edge + 2*i + l_edge + 1]-= N00_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0))) + N01_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N10_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*pow(mKi[0], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N11_dy0[0]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*pow(mKi[2], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))));
K[2*i*l_edge + 2*i + l_edge + 2]-= N00_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0))) + N01_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N10_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*pow(mKi[0], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)))) + N11_dy0[1]*p_Vol[0]*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*pow(mKi[2], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))));
/* Evaluation of K */
K[2*i*(l_edge + 1) + 2*i]-= 0;
K[2*i*l_edge + 4*i + 1]-= 0;
K[2*i*l_edge + 4*i + l_edge + 1]-= 0;
K[2*i*l_edge + 4*i + l_edge + 2]-= 0;
{
int j;
for(j=(int)(1);j<(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L);j++) {
/* Evaluation of K */
K[2*j]+= alpha[j - 1]*pow(p_Vol[0], 2)*(-x[0] + x[2*j])*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[1], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[2]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
)) + alpha[j - 1]*pow(p_Vol[0], 2)*(-x[1] + x[2*j + 1])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
));
K[2*j + 1]+= alpha[j - 1]*pow(p_Vol[0], 2)*(-x[0] + x[2*j])*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(alpha[i - 1]*mKi[2]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[3], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
)) + alpha[j - 1]*pow(p_Vol[0], 2)*(-x[1] + x[2*j + 1])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
));
K[2*j + l_edge + 1]+= alpha[j - 1]*pow(p_Vol[0], 2)*(-x[0] + x[2*j])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
)) + alpha[j - 1]*pow(p_Vol[0], 2)*(-x[1] + x[2*j + 1])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*pow(mKi[0], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
));
K[2*j + l_edge + 2]+= alpha[j - 1]*pow(p_Vol[0], 2)*(-x[0] + x[2*j])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
)) + alpha[j - 1]*pow(p_Vol[0], 2)*(-x[1] + x[2*j + 1])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*pow(mKi[2], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
));
/* Evaluation of K */
K[2*i*(l_edge + 1) + 2*j]-= alpha[j - 1]*pow(p_Vol[0], 2)*(-x[0] + x[2*j])*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[1], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[2]*(2.0*mKi[0]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[0]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
)) + alpha[j - 1]*pow(p_Vol[0], 2)*(-x[1] + x[2*j + 1])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
));
K[2*i*l_edge + 2*i + 2*j + 1]-= alpha[j - 1]*pow(p_Vol[0], 2)*(-x[0] + x[2*j])*((-x[0] + x[2*i])*(alpha[i - 1]*mKi[0]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(alpha[i - 1]*mKi[2]*(2.0*mKi[2]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*pow(mKi[3], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
)) + alpha[j - 1]*pow(p_Vol[0], 2)*(-x[1] + x[2*j + 1])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
));
K[2*i*l_edge + 2*i + 2*j + l_edge + 1]-= alpha[j - 1]*pow(p_Vol[0], 2)*(-x[0] + x[2*j])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[1]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[0]*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
)) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0)))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
)) + alpha[j - 1]*pow(p_Vol[0], 2)*(-x[1] + x[2*j + 1])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*pow(mKi[0], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[1]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[1]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
));
K[2*i*l_edge + 2*i + 2*j + l_edge + 2]-= alpha[j - 1]*pow(p_Vol[0], 2)*(-x[0] + x[2*j])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[1]*mKi[2]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*mKi[2]*mKi[3]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + 1.0*alpha[i - 1]*mKi[2]*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
)) + alpha[j - 1]*pow(p_Vol[0], 2)*(-x[1] + x[2*j + 1])*((-x[0] + x[2*i])*(1.0*alpha[i - 1]*mKi[0]*mKi[2]*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[1]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))) + (-x[1] + x[2*i + 1])*(1.0*alpha[i - 1]*pow(mKi[2], 2)*p_E[0]*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))/(2.0*p_nu[0] + 2.0) + alpha[i - 1]*mKi[3]*(2.0*mKi[3]*p_E[0]/(2.0*p_nu[0] + 2.0) + 1.0*mKi[3]*(-0.666666666666667*p_E[0]/(2.0*p_nu[0] + 2.0) + p_E[0]/(-6.0*p_nu[0] + 3.0)))*((sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*i], 2) + pow(-x[1] + x[2*i + 1], 2))
)
: (
   0.0
))))*((sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2)) <= delta[0]) ? (
   pow(1.0 - sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))/delta[0], 2)/sqrt(pow(-x[0] + x[2*j], 2) + pow(-x[1] + x[2*j + 1], 2))
)
: (
   0.0
));
}
}
}
}
}

void Fbased_rational_2_1_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
   Fbased_rational_2_1_eval(l_edge,in, in+(int)(l_edge + 1), in+(int)(l_edge + 1)+(int)(1), in+(int)(l_edge + 1)+(int)(1)+(int)(1), in+(int)(l_edge + 1)+(int)(1)+(int)(1)+(int)(1), in+(int)(l_edge + 1)+(int)(1)+(int)(1)+(int)(1)+(int)(1), in+(int)(l_edge + 1)+(int)(1)+(int)(1)+(int)(1)+(int)(1)+(int)(l_edge + 1), 
out, out+(int)(pow(l_edge + 1, 2)));
}

void kmap_Fbased_rational_2_10(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)((int)((1.0L/2.0L)*l_edge + 1.0L/2.0L));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_Fbased_rational_2_11(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)(1));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_Fbased_rational_2_12(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)((int)((1.0L/2.0L)*l_edge - 1.0L/2.0L));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)((1.0L/2.0L)*l_edge + 1.0L/2.0L) + 1*i];
    }
  }
}
kernel_t kernel_Fbased_rational_2_1 = {
.nmap = 3,
.maps = {
kmap_Fbased_rational_2_10,
kmap_Fbased_rational_2_11,
kmap_Fbased_rational_2_12
},

.ninp = 7,
.inp = {
{ .field_number=0, .map_num=0, .name="x" },
{ .field_number=1, .map_num=1, .name="p_Vol" },
{ .field_number=2, .map_num=1, .name="p_E" },
{ .field_number=3, .map_num=1, .name="p_nu" },
{ .field_number=4, .map_num=1, .name="delta" },
{ .field_number=5, .map_num=0, .name="y" },
{ .field_number=6, .map_num=2, .name="alpha" }
},

.noutp = 2,
.outp = {
{ .rank = 2, .nmap = 1, .map_nums = { 0 }, .name="K" },
{ .rank = 1, .nmap = 1, .map_nums = { 0 }, .name="R" }
},
.eval=Fbased_rational_2_1_eval_wr,
.name="Fbased_rational_2_1"
};

