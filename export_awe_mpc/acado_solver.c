/*
 *    This file was auto-generated using the ACADO Toolkit.
 *    
 *    While ACADO Toolkit is free software released under the terms of
 *    the GNU Lesser General Public License (LGPL), the generated code
 *    as such remains the property of the user who used ACADO Toolkit
 *    to generate this code. In particular, user dependent data of the code
 *    do not inherit the GNU LGPL license. On the other hand, parts of the
 *    generated code that are a direct copy of source code from the
 *    ACADO Toolkit or the software tools it is based on, remain, as derived
 *    work, automatically covered by the LGPL license.
 *    
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *    
 */


#include "acado_common.h"




/******************************************************************************/
/*                                                                            */
/* ACADO code generation                                                      */
/*                                                                            */
/******************************************************************************/


int acado_modelSimulation(  )
{
int ret;

int lRun1;
ret = 0;
for (lRun1 = 0; lRun1 < 60; ++lRun1)
{
acadoWorkspace.state[0] = acadoVariables.x[lRun1 * 4];
acadoWorkspace.state[1] = acadoVariables.x[lRun1 * 4 + 1];
acadoWorkspace.state[2] = acadoVariables.x[lRun1 * 4 + 2];
acadoWorkspace.state[3] = acadoVariables.x[lRun1 * 4 + 3];

acadoWorkspace.state[24] = acadoVariables.u[lRun1];

ret = acado_integrate(acadoWorkspace.state, 1);

acadoWorkspace.d[lRun1 * 4] = acadoWorkspace.state[0] - acadoVariables.x[lRun1 * 4 + 4];
acadoWorkspace.d[lRun1 * 4 + 1] = acadoWorkspace.state[1] - acadoVariables.x[lRun1 * 4 + 5];
acadoWorkspace.d[lRun1 * 4 + 2] = acadoWorkspace.state[2] - acadoVariables.x[lRun1 * 4 + 6];
acadoWorkspace.d[lRun1 * 4 + 3] = acadoWorkspace.state[3] - acadoVariables.x[lRun1 * 4 + 7];

acadoWorkspace.evGx[lRun1 * 16] = acadoWorkspace.state[4];
acadoWorkspace.evGx[lRun1 * 16 + 1] = acadoWorkspace.state[5];
acadoWorkspace.evGx[lRun1 * 16 + 2] = acadoWorkspace.state[6];
acadoWorkspace.evGx[lRun1 * 16 + 3] = acadoWorkspace.state[7];
acadoWorkspace.evGx[lRun1 * 16 + 4] = acadoWorkspace.state[8];
acadoWorkspace.evGx[lRun1 * 16 + 5] = acadoWorkspace.state[9];
acadoWorkspace.evGx[lRun1 * 16 + 6] = acadoWorkspace.state[10];
acadoWorkspace.evGx[lRun1 * 16 + 7] = acadoWorkspace.state[11];
acadoWorkspace.evGx[lRun1 * 16 + 8] = acadoWorkspace.state[12];
acadoWorkspace.evGx[lRun1 * 16 + 9] = acadoWorkspace.state[13];
acadoWorkspace.evGx[lRun1 * 16 + 10] = acadoWorkspace.state[14];
acadoWorkspace.evGx[lRun1 * 16 + 11] = acadoWorkspace.state[15];
acadoWorkspace.evGx[lRun1 * 16 + 12] = acadoWorkspace.state[16];
acadoWorkspace.evGx[lRun1 * 16 + 13] = acadoWorkspace.state[17];
acadoWorkspace.evGx[lRun1 * 16 + 14] = acadoWorkspace.state[18];
acadoWorkspace.evGx[lRun1 * 16 + 15] = acadoWorkspace.state[19];

acadoWorkspace.evGu[lRun1 * 4] = acadoWorkspace.state[20];
acadoWorkspace.evGu[lRun1 * 4 + 1] = acadoWorkspace.state[21];
acadoWorkspace.evGu[lRun1 * 4 + 2] = acadoWorkspace.state[22];
acadoWorkspace.evGu[lRun1 * 4 + 3] = acadoWorkspace.state[23];
}
return ret;
}

void acado_evaluateLSQ(const real_t* in, real_t* out)
{
const real_t* xd = in;
const real_t* u = in + 4;
/* Vector of auxiliary variables; number of elements: 8. */
real_t* a = acadoWorkspace.objAuxVar;

/* Compute intermediate quantities: */
a[0] = (sqrt(((xd[0]*xd[0])+(xd[1]*xd[1]))));
a[1] = (1.0/sqrt(((xd[0]*xd[0])+(xd[1]*xd[1]))));
a[2] = (a[1]*(real_t)(5.0000000000000000e-01));
a[3] = ((xd[0]+xd[0])*a[2]);
a[4] = ((xd[1]+xd[1])*a[2]);
a[5] = (real_t)(0.0000000000000000e+00);
a[6] = (real_t)(0.0000000000000000e+00);
a[7] = (real_t)(0.0000000000000000e+00);

/* Compute outputs: */
out[0] = a[0];
out[1] = u[0];
out[2] = a[3];
out[3] = a[4];
out[4] = a[5];
out[5] = a[6];
out[6] = (real_t)(0.0000000000000000e+00);
out[7] = (real_t)(0.0000000000000000e+00);
out[8] = (real_t)(0.0000000000000000e+00);
out[9] = (real_t)(0.0000000000000000e+00);
out[10] = a[7];
out[11] = (real_t)(1.0000000000000000e+00);
}

void acado_evaluateLSQEndTerm(const real_t* in, real_t* out)
{
const real_t* xd = in;
/* Vector of auxiliary variables; number of elements: 7. */
real_t* a = acadoWorkspace.objAuxVar;

/* Compute intermediate quantities: */
a[0] = (sqrt(((xd[0]*xd[0])+(xd[1]*xd[1]))));
a[1] = (1.0/sqrt(((xd[0]*xd[0])+(xd[1]*xd[1]))));
a[2] = (a[1]*(real_t)(5.0000000000000000e-01));
a[3] = ((xd[0]+xd[0])*a[2]);
a[4] = ((xd[1]+xd[1])*a[2]);
a[5] = (real_t)(0.0000000000000000e+00);
a[6] = (real_t)(0.0000000000000000e+00);

/* Compute outputs: */
out[0] = a[0];
out[1] = a[3];
out[2] = a[4];
out[3] = a[5];
out[4] = a[6];
}

void acado_setObjQ1Q2( real_t* const tmpFx, real_t* const tmpObjS, real_t* const tmpQ1, real_t* const tmpQ2 )
{
tmpQ2[0] = + tmpFx[0]*tmpObjS[0] + tmpFx[4]*tmpObjS[2];
tmpQ2[1] = + tmpFx[0]*tmpObjS[1] + tmpFx[4]*tmpObjS[3];
tmpQ2[2] = + tmpFx[1]*tmpObjS[0] + tmpFx[5]*tmpObjS[2];
tmpQ2[3] = + tmpFx[1]*tmpObjS[1] + tmpFx[5]*tmpObjS[3];
tmpQ2[4] = + tmpFx[2]*tmpObjS[0] + tmpFx[6]*tmpObjS[2];
tmpQ2[5] = + tmpFx[2]*tmpObjS[1] + tmpFx[6]*tmpObjS[3];
tmpQ2[6] = + tmpFx[3]*tmpObjS[0] + tmpFx[7]*tmpObjS[2];
tmpQ2[7] = + tmpFx[3]*tmpObjS[1] + tmpFx[7]*tmpObjS[3];
tmpQ1[0] = + tmpQ2[0]*tmpFx[0] + tmpQ2[1]*tmpFx[4];
tmpQ1[1] = + tmpQ2[0]*tmpFx[1] + tmpQ2[1]*tmpFx[5];
tmpQ1[2] = + tmpQ2[0]*tmpFx[2] + tmpQ2[1]*tmpFx[6];
tmpQ1[3] = + tmpQ2[0]*tmpFx[3] + tmpQ2[1]*tmpFx[7];
tmpQ1[4] = + tmpQ2[2]*tmpFx[0] + tmpQ2[3]*tmpFx[4];
tmpQ1[5] = + tmpQ2[2]*tmpFx[1] + tmpQ2[3]*tmpFx[5];
tmpQ1[6] = + tmpQ2[2]*tmpFx[2] + tmpQ2[3]*tmpFx[6];
tmpQ1[7] = + tmpQ2[2]*tmpFx[3] + tmpQ2[3]*tmpFx[7];
tmpQ1[8] = + tmpQ2[4]*tmpFx[0] + tmpQ2[5]*tmpFx[4];
tmpQ1[9] = + tmpQ2[4]*tmpFx[1] + tmpQ2[5]*tmpFx[5];
tmpQ1[10] = + tmpQ2[4]*tmpFx[2] + tmpQ2[5]*tmpFx[6];
tmpQ1[11] = + tmpQ2[4]*tmpFx[3] + tmpQ2[5]*tmpFx[7];
tmpQ1[12] = + tmpQ2[6]*tmpFx[0] + tmpQ2[7]*tmpFx[4];
tmpQ1[13] = + tmpQ2[6]*tmpFx[1] + tmpQ2[7]*tmpFx[5];
tmpQ1[14] = + tmpQ2[6]*tmpFx[2] + tmpQ2[7]*tmpFx[6];
tmpQ1[15] = + tmpQ2[6]*tmpFx[3] + tmpQ2[7]*tmpFx[7];
}

void acado_setObjR1R2( real_t* const tmpFu, real_t* const tmpObjS, real_t* const tmpR1, real_t* const tmpR2 )
{
tmpR2[0] = + tmpFu[0]*tmpObjS[0] + tmpFu[1]*tmpObjS[2];
tmpR2[1] = + tmpFu[0]*tmpObjS[1] + tmpFu[1]*tmpObjS[3];
tmpR1[0] = + tmpR2[0]*tmpFu[0] + tmpR2[1]*tmpFu[1];
}

void acado_setObjQN1QN2( real_t* const tmpFx, real_t* const tmpObjSEndTerm, real_t* const tmpQN1, real_t* const tmpQN2 )
{
tmpQN2[0] = + tmpFx[0]*tmpObjSEndTerm[0];
tmpQN2[1] = + tmpFx[1]*tmpObjSEndTerm[0];
tmpQN2[2] = + tmpFx[2]*tmpObjSEndTerm[0];
tmpQN2[3] = + tmpFx[3]*tmpObjSEndTerm[0];
tmpQN1[0] = + tmpQN2[0]*tmpFx[0];
tmpQN1[1] = + tmpQN2[0]*tmpFx[1];
tmpQN1[2] = + tmpQN2[0]*tmpFx[2];
tmpQN1[3] = + tmpQN2[0]*tmpFx[3];
tmpQN1[4] = + tmpQN2[1]*tmpFx[0];
tmpQN1[5] = + tmpQN2[1]*tmpFx[1];
tmpQN1[6] = + tmpQN2[1]*tmpFx[2];
tmpQN1[7] = + tmpQN2[1]*tmpFx[3];
tmpQN1[8] = + tmpQN2[2]*tmpFx[0];
tmpQN1[9] = + tmpQN2[2]*tmpFx[1];
tmpQN1[10] = + tmpQN2[2]*tmpFx[2];
tmpQN1[11] = + tmpQN2[2]*tmpFx[3];
tmpQN1[12] = + tmpQN2[3]*tmpFx[0];
tmpQN1[13] = + tmpQN2[3]*tmpFx[1];
tmpQN1[14] = + tmpQN2[3]*tmpFx[2];
tmpQN1[15] = + tmpQN2[3]*tmpFx[3];
}

void acado_evaluateObjective(  )
{
int runObj;
for (runObj = 0; runObj < 60; ++runObj)
{
acadoWorkspace.objValueIn[0] = acadoVariables.x[runObj * 4];
acadoWorkspace.objValueIn[1] = acadoVariables.x[runObj * 4 + 1];
acadoWorkspace.objValueIn[2] = acadoVariables.x[runObj * 4 + 2];
acadoWorkspace.objValueIn[3] = acadoVariables.x[runObj * 4 + 3];
acadoWorkspace.objValueIn[4] = acadoVariables.u[runObj];

acado_evaluateLSQ( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.Dy[runObj * 2] = acadoWorkspace.objValueOut[0];
acadoWorkspace.Dy[runObj * 2 + 1] = acadoWorkspace.objValueOut[1];

acado_setObjQ1Q2( &(acadoWorkspace.objValueOut[ 2 ]), acadoVariables.W, &(acadoWorkspace.Q1[ runObj * 16 ]), &(acadoWorkspace.Q2[ runObj * 8 ]) );

acado_setObjR1R2( &(acadoWorkspace.objValueOut[ 10 ]), acadoVariables.W, &(acadoWorkspace.R1[ runObj ]), &(acadoWorkspace.R2[ runObj * 2 ]) );

}
acadoWorkspace.objValueIn[0] = acadoVariables.x[240];
acadoWorkspace.objValueIn[1] = acadoVariables.x[241];
acadoWorkspace.objValueIn[2] = acadoVariables.x[242];
acadoWorkspace.objValueIn[3] = acadoVariables.x[243];
acado_evaluateLSQEndTerm( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );

acadoWorkspace.DyN[0] = acadoWorkspace.objValueOut[0];

acado_setObjQN1QN2( &(acadoWorkspace.objValueOut[ 1 ]), acadoVariables.WN, acadoWorkspace.QN1, acadoWorkspace.QN2 );

}

void acado_multGxd( real_t* const dOld, real_t* const Gx1, real_t* const dNew )
{
dNew[0] += + Gx1[0]*dOld[0] + Gx1[1]*dOld[1] + Gx1[2]*dOld[2] + Gx1[3]*dOld[3];
dNew[1] += + Gx1[4]*dOld[0] + Gx1[5]*dOld[1] + Gx1[6]*dOld[2] + Gx1[7]*dOld[3];
dNew[2] += + Gx1[8]*dOld[0] + Gx1[9]*dOld[1] + Gx1[10]*dOld[2] + Gx1[11]*dOld[3];
dNew[3] += + Gx1[12]*dOld[0] + Gx1[13]*dOld[1] + Gx1[14]*dOld[2] + Gx1[15]*dOld[3];
}

void acado_moveGxT( real_t* const Gx1, real_t* const Gx2 )
{
Gx2[0] = Gx1[0];
Gx2[1] = Gx1[1];
Gx2[2] = Gx1[2];
Gx2[3] = Gx1[3];
Gx2[4] = Gx1[4];
Gx2[5] = Gx1[5];
Gx2[6] = Gx1[6];
Gx2[7] = Gx1[7];
Gx2[8] = Gx1[8];
Gx2[9] = Gx1[9];
Gx2[10] = Gx1[10];
Gx2[11] = Gx1[11];
Gx2[12] = Gx1[12];
Gx2[13] = Gx1[13];
Gx2[14] = Gx1[14];
Gx2[15] = Gx1[15];
}

void acado_multGxGx( real_t* const Gx1, real_t* const Gx2, real_t* const Gx3 )
{
Gx3[0] = + Gx1[0]*Gx2[0] + Gx1[1]*Gx2[4] + Gx1[2]*Gx2[8] + Gx1[3]*Gx2[12];
Gx3[1] = + Gx1[0]*Gx2[1] + Gx1[1]*Gx2[5] + Gx1[2]*Gx2[9] + Gx1[3]*Gx2[13];
Gx3[2] = + Gx1[0]*Gx2[2] + Gx1[1]*Gx2[6] + Gx1[2]*Gx2[10] + Gx1[3]*Gx2[14];
Gx3[3] = + Gx1[0]*Gx2[3] + Gx1[1]*Gx2[7] + Gx1[2]*Gx2[11] + Gx1[3]*Gx2[15];
Gx3[4] = + Gx1[4]*Gx2[0] + Gx1[5]*Gx2[4] + Gx1[6]*Gx2[8] + Gx1[7]*Gx2[12];
Gx3[5] = + Gx1[4]*Gx2[1] + Gx1[5]*Gx2[5] + Gx1[6]*Gx2[9] + Gx1[7]*Gx2[13];
Gx3[6] = + Gx1[4]*Gx2[2] + Gx1[5]*Gx2[6] + Gx1[6]*Gx2[10] + Gx1[7]*Gx2[14];
Gx3[7] = + Gx1[4]*Gx2[3] + Gx1[5]*Gx2[7] + Gx1[6]*Gx2[11] + Gx1[7]*Gx2[15];
Gx3[8] = + Gx1[8]*Gx2[0] + Gx1[9]*Gx2[4] + Gx1[10]*Gx2[8] + Gx1[11]*Gx2[12];
Gx3[9] = + Gx1[8]*Gx2[1] + Gx1[9]*Gx2[5] + Gx1[10]*Gx2[9] + Gx1[11]*Gx2[13];
Gx3[10] = + Gx1[8]*Gx2[2] + Gx1[9]*Gx2[6] + Gx1[10]*Gx2[10] + Gx1[11]*Gx2[14];
Gx3[11] = + Gx1[8]*Gx2[3] + Gx1[9]*Gx2[7] + Gx1[10]*Gx2[11] + Gx1[11]*Gx2[15];
Gx3[12] = + Gx1[12]*Gx2[0] + Gx1[13]*Gx2[4] + Gx1[14]*Gx2[8] + Gx1[15]*Gx2[12];
Gx3[13] = + Gx1[12]*Gx2[1] + Gx1[13]*Gx2[5] + Gx1[14]*Gx2[9] + Gx1[15]*Gx2[13];
Gx3[14] = + Gx1[12]*Gx2[2] + Gx1[13]*Gx2[6] + Gx1[14]*Gx2[10] + Gx1[15]*Gx2[14];
Gx3[15] = + Gx1[12]*Gx2[3] + Gx1[13]*Gx2[7] + Gx1[14]*Gx2[11] + Gx1[15]*Gx2[15];
}

void acado_multGxGu( real_t* const Gx1, real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = + Gx1[0]*Gu1[0] + Gx1[1]*Gu1[1] + Gx1[2]*Gu1[2] + Gx1[3]*Gu1[3];
Gu2[1] = + Gx1[4]*Gu1[0] + Gx1[5]*Gu1[1] + Gx1[6]*Gu1[2] + Gx1[7]*Gu1[3];
Gu2[2] = + Gx1[8]*Gu1[0] + Gx1[9]*Gu1[1] + Gx1[10]*Gu1[2] + Gx1[11]*Gu1[3];
Gu2[3] = + Gx1[12]*Gu1[0] + Gx1[13]*Gu1[1] + Gx1[14]*Gu1[2] + Gx1[15]*Gu1[3];
}

void acado_moveGuE( real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = Gu1[0];
Gu2[1] = Gu1[1];
Gu2[2] = Gu1[2];
Gu2[3] = Gu1[3];
}

void acado_setBlockH11( int iRow, int iCol, real_t* const Gu1, real_t* const Gu2 )
{
acadoWorkspace.H[(iRow * 60) + (iCol)] += + Gu1[0]*Gu2[0] + Gu1[1]*Gu2[1] + Gu1[2]*Gu2[2] + Gu1[3]*Gu2[3];
}

void acado_setBlockH11_R1( int iRow, int iCol, real_t* const R11 )
{
acadoWorkspace.H[(iRow * 60) + (iCol)] = R11[0];
}

void acado_zeroBlockH11( int iRow, int iCol )
{
acadoWorkspace.H[(iRow * 60) + (iCol)] = 0.0000000000000000e+00;
}

void acado_copyHTH( int iRow, int iCol )
{
acadoWorkspace.H[(iRow * 60) + (iCol)] = acadoWorkspace.H[(iCol * 60) + (iRow)];
}

void acado_multQ1d( real_t* const Gx1, real_t* const dOld, real_t* const dNew )
{
dNew[0] = + Gx1[0]*dOld[0] + Gx1[1]*dOld[1] + Gx1[2]*dOld[2] + Gx1[3]*dOld[3];
dNew[1] = + Gx1[4]*dOld[0] + Gx1[5]*dOld[1] + Gx1[6]*dOld[2] + Gx1[7]*dOld[3];
dNew[2] = + Gx1[8]*dOld[0] + Gx1[9]*dOld[1] + Gx1[10]*dOld[2] + Gx1[11]*dOld[3];
dNew[3] = + Gx1[12]*dOld[0] + Gx1[13]*dOld[1] + Gx1[14]*dOld[2] + Gx1[15]*dOld[3];
}

void acado_multQN1d( real_t* const QN1, real_t* const dOld, real_t* const dNew )
{
dNew[0] = + acadoWorkspace.QN1[0]*dOld[0] + acadoWorkspace.QN1[1]*dOld[1] + acadoWorkspace.QN1[2]*dOld[2] + acadoWorkspace.QN1[3]*dOld[3];
dNew[1] = + acadoWorkspace.QN1[4]*dOld[0] + acadoWorkspace.QN1[5]*dOld[1] + acadoWorkspace.QN1[6]*dOld[2] + acadoWorkspace.QN1[7]*dOld[3];
dNew[2] = + acadoWorkspace.QN1[8]*dOld[0] + acadoWorkspace.QN1[9]*dOld[1] + acadoWorkspace.QN1[10]*dOld[2] + acadoWorkspace.QN1[11]*dOld[3];
dNew[3] = + acadoWorkspace.QN1[12]*dOld[0] + acadoWorkspace.QN1[13]*dOld[1] + acadoWorkspace.QN1[14]*dOld[2] + acadoWorkspace.QN1[15]*dOld[3];
}

void acado_multRDy( real_t* const R2, real_t* const Dy1, real_t* const RDy1 )
{
RDy1[0] = + R2[0]*Dy1[0] + R2[1]*Dy1[1];
}

void acado_multQDy( real_t* const Q2, real_t* const Dy1, real_t* const QDy1 )
{
QDy1[0] = + Q2[0]*Dy1[0] + Q2[1]*Dy1[1];
QDy1[1] = + Q2[2]*Dy1[0] + Q2[3]*Dy1[1];
QDy1[2] = + Q2[4]*Dy1[0] + Q2[5]*Dy1[1];
QDy1[3] = + Q2[6]*Dy1[0] + Q2[7]*Dy1[1];
}

void acado_multEQDy( real_t* const E1, real_t* const QDy1, real_t* const U1 )
{
U1[0] += + E1[0]*QDy1[0] + E1[1]*QDy1[1] + E1[2]*QDy1[2] + E1[3]*QDy1[3];
}

void acado_multQETGx( real_t* const E1, real_t* const Gx1, real_t* const H101 )
{
H101[0] += + E1[0]*Gx1[0] + E1[1]*Gx1[4] + E1[2]*Gx1[8] + E1[3]*Gx1[12];
H101[1] += + E1[0]*Gx1[1] + E1[1]*Gx1[5] + E1[2]*Gx1[9] + E1[3]*Gx1[13];
H101[2] += + E1[0]*Gx1[2] + E1[1]*Gx1[6] + E1[2]*Gx1[10] + E1[3]*Gx1[14];
H101[3] += + E1[0]*Gx1[3] + E1[1]*Gx1[7] + E1[2]*Gx1[11] + E1[3]*Gx1[15];
}

void acado_zeroBlockH10( real_t* const H101 )
{
{ int lCopy; for (lCopy = 0; lCopy < 4; lCopy++) H101[ lCopy ] = 0; }
}

void acado_multEDu( real_t* const E1, real_t* const U1, real_t* const dNew )
{
dNew[0] += + E1[0]*U1[0];
dNew[1] += + E1[1]*U1[0];
dNew[2] += + E1[2]*U1[0];
dNew[3] += + E1[3]*U1[0];
}

void acado_macETSlu( real_t* const E0, real_t* const g1 )
{
g1[0] += 0.0;
;
}

void acado_condensePrep(  )
{
int lRun1;
int lRun2;
int lRun3;
int lRun4;
int lRun5;
/** Row vector of size: 60 */
static const int xBoundIndices[ 60 ] = 
{ 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99, 103, 107, 111, 115, 119, 123, 127, 131, 135, 139, 143, 147, 151, 155, 159, 163, 167, 171, 175, 179, 183, 187, 191, 195, 199, 203, 207, 211, 215, 219, 223, 227, 231, 235, 239, 243 };
acado_moveGuE( acadoWorkspace.evGu, acadoWorkspace.E );
for (lRun1 = 1; lRun1 < 60; ++lRun1)
{
acado_moveGxT( &(acadoWorkspace.evGx[ lRun1 * 16 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ lRun1 * 4-4 ]), &(acadoWorkspace.evGx[ lRun1 * 16 ]), &(acadoWorkspace.d[ lRun1 * 4 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ lRun1 * 16-16 ]), &(acadoWorkspace.evGx[ lRun1 * 16 ]) );
for (lRun2 = 0; lRun2 < lRun1; ++lRun2)
{
lRun4 = (((lRun1) * (lRun1-1)) / (2)) + (lRun2);
lRun3 = (((lRun1 + 1) * (lRun1)) / (2)) + (lRun2);
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ lRun4 * 4 ]), &(acadoWorkspace.E[ lRun3 * 4 ]) );
}
lRun3 = (((lRun1 + 1) * (lRun1)) / (2)) + (lRun2);
acado_moveGuE( &(acadoWorkspace.evGu[ lRun1 * 4 ]), &(acadoWorkspace.E[ lRun3 * 4 ]) );
}

for (lRun1 = 0; lRun1 < 59; ++lRun1)
{
for (lRun2 = 0; lRun2 < lRun1 + 1; ++lRun2)
{
lRun3 = (((lRun1 + 1) * (lRun1)) / (2)) + (lRun2);
acado_multGxGu( &(acadoWorkspace.Q1[ lRun1 * 16 + 16 ]), &(acadoWorkspace.E[ lRun3 * 4 ]), &(acadoWorkspace.QE[ lRun3 * 4 ]) );
}
}

for (lRun2 = 0; lRun2 < lRun1 + 1; ++lRun2)
{
lRun3 = (((lRun1 + 1) * (lRun1)) / (2)) + (lRun2);
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ lRun3 * 4 ]), &(acadoWorkspace.QE[ lRun3 * 4 ]) );
}

for (lRun1 = 0; lRun1 < 60; ++lRun1)
{
acado_zeroBlockH10( &(acadoWorkspace.H10[ lRun1 * 4 ]) );
for (lRun2 = lRun1; lRun2 < 60; ++lRun2)
{
lRun3 = (((lRun2 + 1) * (lRun2)) / (2)) + (lRun1);
acado_multQETGx( &(acadoWorkspace.QE[ lRun3 * 4 ]), &(acadoWorkspace.evGx[ lRun2 * 16 ]), &(acadoWorkspace.H10[ lRun1 * 4 ]) );
}
}

for (lRun1 = 0; lRun1 < 60; ++lRun1)
{
acado_setBlockH11_R1( lRun1, lRun1, &(acadoWorkspace.R1[ lRun1 ]) );
lRun2 = lRun1;
for (lRun3 = lRun1; lRun3 < 60; ++lRun3)
{
lRun4 = (((lRun3 + 1) * (lRun3)) / (2)) + (lRun1);
lRun5 = (((lRun3 + 1) * (lRun3)) / (2)) + (lRun2);
acado_setBlockH11( lRun1, lRun2, &(acadoWorkspace.E[ lRun4 * 4 ]), &(acadoWorkspace.QE[ lRun5 * 4 ]) );
}
for (lRun2 = lRun1 + 1; lRun2 < 60; ++lRun2)
{
acado_zeroBlockH11( lRun1, lRun2 );
for (lRun3 = lRun2; lRun3 < 60; ++lRun3)
{
lRun4 = (((lRun3 + 1) * (lRun3)) / (2)) + (lRun1);
lRun5 = (((lRun3 + 1) * (lRun3)) / (2)) + (lRun2);
acado_setBlockH11( lRun1, lRun2, &(acadoWorkspace.E[ lRun4 * 4 ]), &(acadoWorkspace.QE[ lRun5 * 4 ]) );
}
}
}

for (lRun1 = 0; lRun1 < 60; ++lRun1)
{
for (lRun2 = 0; lRun2 < lRun1; ++lRun2)
{
acado_copyHTH( lRun1, lRun2 );
}
}

acado_multQ1d( &(acadoWorkspace.Q1[ 16 ]), acadoWorkspace.d, acadoWorkspace.Qd );
acado_multQ1d( &(acadoWorkspace.Q1[ 32 ]), &(acadoWorkspace.d[ 4 ]), &(acadoWorkspace.Qd[ 4 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 48 ]), &(acadoWorkspace.d[ 8 ]), &(acadoWorkspace.Qd[ 8 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 64 ]), &(acadoWorkspace.d[ 12 ]), &(acadoWorkspace.Qd[ 12 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 80 ]), &(acadoWorkspace.d[ 16 ]), &(acadoWorkspace.Qd[ 16 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 96 ]), &(acadoWorkspace.d[ 20 ]), &(acadoWorkspace.Qd[ 20 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 112 ]), &(acadoWorkspace.d[ 24 ]), &(acadoWorkspace.Qd[ 24 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 128 ]), &(acadoWorkspace.d[ 28 ]), &(acadoWorkspace.Qd[ 28 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.d[ 32 ]), &(acadoWorkspace.Qd[ 32 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 160 ]), &(acadoWorkspace.d[ 36 ]), &(acadoWorkspace.Qd[ 36 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 176 ]), &(acadoWorkspace.d[ 40 ]), &(acadoWorkspace.Qd[ 40 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 192 ]), &(acadoWorkspace.d[ 44 ]), &(acadoWorkspace.Qd[ 44 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 208 ]), &(acadoWorkspace.d[ 48 ]), &(acadoWorkspace.Qd[ 48 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 224 ]), &(acadoWorkspace.d[ 52 ]), &(acadoWorkspace.Qd[ 52 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 240 ]), &(acadoWorkspace.d[ 56 ]), &(acadoWorkspace.Qd[ 56 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 256 ]), &(acadoWorkspace.d[ 60 ]), &(acadoWorkspace.Qd[ 60 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 272 ]), &(acadoWorkspace.d[ 64 ]), &(acadoWorkspace.Qd[ 64 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.d[ 68 ]), &(acadoWorkspace.Qd[ 68 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 304 ]), &(acadoWorkspace.d[ 72 ]), &(acadoWorkspace.Qd[ 72 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 320 ]), &(acadoWorkspace.d[ 76 ]), &(acadoWorkspace.Qd[ 76 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 336 ]), &(acadoWorkspace.d[ 80 ]), &(acadoWorkspace.Qd[ 80 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 352 ]), &(acadoWorkspace.d[ 84 ]), &(acadoWorkspace.Qd[ 84 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 368 ]), &(acadoWorkspace.d[ 88 ]), &(acadoWorkspace.Qd[ 88 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 384 ]), &(acadoWorkspace.d[ 92 ]), &(acadoWorkspace.Qd[ 92 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 400 ]), &(acadoWorkspace.d[ 96 ]), &(acadoWorkspace.Qd[ 96 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 416 ]), &(acadoWorkspace.d[ 100 ]), &(acadoWorkspace.Qd[ 100 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.d[ 104 ]), &(acadoWorkspace.Qd[ 104 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 448 ]), &(acadoWorkspace.d[ 108 ]), &(acadoWorkspace.Qd[ 108 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 464 ]), &(acadoWorkspace.d[ 112 ]), &(acadoWorkspace.Qd[ 112 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 480 ]), &(acadoWorkspace.d[ 116 ]), &(acadoWorkspace.Qd[ 116 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 496 ]), &(acadoWorkspace.d[ 120 ]), &(acadoWorkspace.Qd[ 120 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 512 ]), &(acadoWorkspace.d[ 124 ]), &(acadoWorkspace.Qd[ 124 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 528 ]), &(acadoWorkspace.d[ 128 ]), &(acadoWorkspace.Qd[ 128 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 544 ]), &(acadoWorkspace.d[ 132 ]), &(acadoWorkspace.Qd[ 132 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 560 ]), &(acadoWorkspace.d[ 136 ]), &(acadoWorkspace.Qd[ 136 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.d[ 140 ]), &(acadoWorkspace.Qd[ 140 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 592 ]), &(acadoWorkspace.d[ 144 ]), &(acadoWorkspace.Qd[ 144 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 608 ]), &(acadoWorkspace.d[ 148 ]), &(acadoWorkspace.Qd[ 148 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 624 ]), &(acadoWorkspace.d[ 152 ]), &(acadoWorkspace.Qd[ 152 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 640 ]), &(acadoWorkspace.d[ 156 ]), &(acadoWorkspace.Qd[ 156 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 656 ]), &(acadoWorkspace.d[ 160 ]), &(acadoWorkspace.Qd[ 160 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 672 ]), &(acadoWorkspace.d[ 164 ]), &(acadoWorkspace.Qd[ 164 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 688 ]), &(acadoWorkspace.d[ 168 ]), &(acadoWorkspace.Qd[ 168 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 704 ]), &(acadoWorkspace.d[ 172 ]), &(acadoWorkspace.Qd[ 172 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 720 ]), &(acadoWorkspace.d[ 176 ]), &(acadoWorkspace.Qd[ 176 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 736 ]), &(acadoWorkspace.d[ 180 ]), &(acadoWorkspace.Qd[ 180 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 752 ]), &(acadoWorkspace.d[ 184 ]), &(acadoWorkspace.Qd[ 184 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 768 ]), &(acadoWorkspace.d[ 188 ]), &(acadoWorkspace.Qd[ 188 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 784 ]), &(acadoWorkspace.d[ 192 ]), &(acadoWorkspace.Qd[ 192 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 800 ]), &(acadoWorkspace.d[ 196 ]), &(acadoWorkspace.Qd[ 196 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 816 ]), &(acadoWorkspace.d[ 200 ]), &(acadoWorkspace.Qd[ 200 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 832 ]), &(acadoWorkspace.d[ 204 ]), &(acadoWorkspace.Qd[ 204 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 848 ]), &(acadoWorkspace.d[ 208 ]), &(acadoWorkspace.Qd[ 208 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 864 ]), &(acadoWorkspace.d[ 212 ]), &(acadoWorkspace.Qd[ 212 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 880 ]), &(acadoWorkspace.d[ 216 ]), &(acadoWorkspace.Qd[ 216 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 896 ]), &(acadoWorkspace.d[ 220 ]), &(acadoWorkspace.Qd[ 220 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 912 ]), &(acadoWorkspace.d[ 224 ]), &(acadoWorkspace.Qd[ 224 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 928 ]), &(acadoWorkspace.d[ 228 ]), &(acadoWorkspace.Qd[ 228 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 944 ]), &(acadoWorkspace.d[ 232 ]), &(acadoWorkspace.Qd[ 232 ]) );
acado_multQN1d( acadoWorkspace.QN1, &(acadoWorkspace.d[ 236 ]), &(acadoWorkspace.Qd[ 236 ]) );

for (lRun1 = 0; lRun1 < 60; ++lRun1)
{
for (lRun2 = lRun1; lRun2 < 60; ++lRun2)
{
lRun3 = (((lRun2 + 1) * (lRun2)) / (2)) + (lRun1);
acado_macETSlu( &(acadoWorkspace.QE[ lRun3 * 4 ]), &(acadoWorkspace.g[ lRun1 ]) );
}
}
acadoWorkspace.lb[0] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[0];
acadoWorkspace.lb[1] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[1];
acadoWorkspace.lb[2] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[2];
acadoWorkspace.lb[3] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[3];
acadoWorkspace.lb[4] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[4];
acadoWorkspace.lb[5] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[5];
acadoWorkspace.lb[6] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[6];
acadoWorkspace.lb[7] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[7];
acadoWorkspace.lb[8] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[8];
acadoWorkspace.lb[9] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[9];
acadoWorkspace.lb[10] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[10];
acadoWorkspace.lb[11] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[11];
acadoWorkspace.lb[12] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[12];
acadoWorkspace.lb[13] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[13];
acadoWorkspace.lb[14] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[14];
acadoWorkspace.lb[15] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[15];
acadoWorkspace.lb[16] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[16];
acadoWorkspace.lb[17] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[17];
acadoWorkspace.lb[18] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[18];
acadoWorkspace.lb[19] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[19];
acadoWorkspace.lb[20] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[20];
acadoWorkspace.lb[21] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[21];
acadoWorkspace.lb[22] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[22];
acadoWorkspace.lb[23] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[23];
acadoWorkspace.lb[24] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[24];
acadoWorkspace.lb[25] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[25];
acadoWorkspace.lb[26] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[26];
acadoWorkspace.lb[27] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[27];
acadoWorkspace.lb[28] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[28];
acadoWorkspace.lb[29] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[29];
acadoWorkspace.lb[30] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[30];
acadoWorkspace.lb[31] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[31];
acadoWorkspace.lb[32] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[32];
acadoWorkspace.lb[33] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[33];
acadoWorkspace.lb[34] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[34];
acadoWorkspace.lb[35] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[35];
acadoWorkspace.lb[36] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[36];
acadoWorkspace.lb[37] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[37];
acadoWorkspace.lb[38] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[38];
acadoWorkspace.lb[39] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[39];
acadoWorkspace.lb[40] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[40];
acadoWorkspace.lb[41] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[41];
acadoWorkspace.lb[42] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[42];
acadoWorkspace.lb[43] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[43];
acadoWorkspace.lb[44] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[44];
acadoWorkspace.lb[45] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[45];
acadoWorkspace.lb[46] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[46];
acadoWorkspace.lb[47] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[47];
acadoWorkspace.lb[48] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[48];
acadoWorkspace.lb[49] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[49];
acadoWorkspace.lb[50] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[50];
acadoWorkspace.lb[51] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[51];
acadoWorkspace.lb[52] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[52];
acadoWorkspace.lb[53] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[53];
acadoWorkspace.lb[54] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[54];
acadoWorkspace.lb[55] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[55];
acadoWorkspace.lb[56] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[56];
acadoWorkspace.lb[57] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[57];
acadoWorkspace.lb[58] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[58];
acadoWorkspace.lb[59] = (real_t)-1.0471975511965976e+00 - acadoVariables.u[59];
acadoWorkspace.ub[0] = (real_t)1.0471975511965976e+00 - acadoVariables.u[0];
acadoWorkspace.ub[1] = (real_t)1.0471975511965976e+00 - acadoVariables.u[1];
acadoWorkspace.ub[2] = (real_t)1.0471975511965976e+00 - acadoVariables.u[2];
acadoWorkspace.ub[3] = (real_t)1.0471975511965976e+00 - acadoVariables.u[3];
acadoWorkspace.ub[4] = (real_t)1.0471975511965976e+00 - acadoVariables.u[4];
acadoWorkspace.ub[5] = (real_t)1.0471975511965976e+00 - acadoVariables.u[5];
acadoWorkspace.ub[6] = (real_t)1.0471975511965976e+00 - acadoVariables.u[6];
acadoWorkspace.ub[7] = (real_t)1.0471975511965976e+00 - acadoVariables.u[7];
acadoWorkspace.ub[8] = (real_t)1.0471975511965976e+00 - acadoVariables.u[8];
acadoWorkspace.ub[9] = (real_t)1.0471975511965976e+00 - acadoVariables.u[9];
acadoWorkspace.ub[10] = (real_t)1.0471975511965976e+00 - acadoVariables.u[10];
acadoWorkspace.ub[11] = (real_t)1.0471975511965976e+00 - acadoVariables.u[11];
acadoWorkspace.ub[12] = (real_t)1.0471975511965976e+00 - acadoVariables.u[12];
acadoWorkspace.ub[13] = (real_t)1.0471975511965976e+00 - acadoVariables.u[13];
acadoWorkspace.ub[14] = (real_t)1.0471975511965976e+00 - acadoVariables.u[14];
acadoWorkspace.ub[15] = (real_t)1.0471975511965976e+00 - acadoVariables.u[15];
acadoWorkspace.ub[16] = (real_t)1.0471975511965976e+00 - acadoVariables.u[16];
acadoWorkspace.ub[17] = (real_t)1.0471975511965976e+00 - acadoVariables.u[17];
acadoWorkspace.ub[18] = (real_t)1.0471975511965976e+00 - acadoVariables.u[18];
acadoWorkspace.ub[19] = (real_t)1.0471975511965976e+00 - acadoVariables.u[19];
acadoWorkspace.ub[20] = (real_t)1.0471975511965976e+00 - acadoVariables.u[20];
acadoWorkspace.ub[21] = (real_t)1.0471975511965976e+00 - acadoVariables.u[21];
acadoWorkspace.ub[22] = (real_t)1.0471975511965976e+00 - acadoVariables.u[22];
acadoWorkspace.ub[23] = (real_t)1.0471975511965976e+00 - acadoVariables.u[23];
acadoWorkspace.ub[24] = (real_t)1.0471975511965976e+00 - acadoVariables.u[24];
acadoWorkspace.ub[25] = (real_t)1.0471975511965976e+00 - acadoVariables.u[25];
acadoWorkspace.ub[26] = (real_t)1.0471975511965976e+00 - acadoVariables.u[26];
acadoWorkspace.ub[27] = (real_t)1.0471975511965976e+00 - acadoVariables.u[27];
acadoWorkspace.ub[28] = (real_t)1.0471975511965976e+00 - acadoVariables.u[28];
acadoWorkspace.ub[29] = (real_t)1.0471975511965976e+00 - acadoVariables.u[29];
acadoWorkspace.ub[30] = (real_t)1.0471975511965976e+00 - acadoVariables.u[30];
acadoWorkspace.ub[31] = (real_t)1.0471975511965976e+00 - acadoVariables.u[31];
acadoWorkspace.ub[32] = (real_t)1.0471975511965976e+00 - acadoVariables.u[32];
acadoWorkspace.ub[33] = (real_t)1.0471975511965976e+00 - acadoVariables.u[33];
acadoWorkspace.ub[34] = (real_t)1.0471975511965976e+00 - acadoVariables.u[34];
acadoWorkspace.ub[35] = (real_t)1.0471975511965976e+00 - acadoVariables.u[35];
acadoWorkspace.ub[36] = (real_t)1.0471975511965976e+00 - acadoVariables.u[36];
acadoWorkspace.ub[37] = (real_t)1.0471975511965976e+00 - acadoVariables.u[37];
acadoWorkspace.ub[38] = (real_t)1.0471975511965976e+00 - acadoVariables.u[38];
acadoWorkspace.ub[39] = (real_t)1.0471975511965976e+00 - acadoVariables.u[39];
acadoWorkspace.ub[40] = (real_t)1.0471975511965976e+00 - acadoVariables.u[40];
acadoWorkspace.ub[41] = (real_t)1.0471975511965976e+00 - acadoVariables.u[41];
acadoWorkspace.ub[42] = (real_t)1.0471975511965976e+00 - acadoVariables.u[42];
acadoWorkspace.ub[43] = (real_t)1.0471975511965976e+00 - acadoVariables.u[43];
acadoWorkspace.ub[44] = (real_t)1.0471975511965976e+00 - acadoVariables.u[44];
acadoWorkspace.ub[45] = (real_t)1.0471975511965976e+00 - acadoVariables.u[45];
acadoWorkspace.ub[46] = (real_t)1.0471975511965976e+00 - acadoVariables.u[46];
acadoWorkspace.ub[47] = (real_t)1.0471975511965976e+00 - acadoVariables.u[47];
acadoWorkspace.ub[48] = (real_t)1.0471975511965976e+00 - acadoVariables.u[48];
acadoWorkspace.ub[49] = (real_t)1.0471975511965976e+00 - acadoVariables.u[49];
acadoWorkspace.ub[50] = (real_t)1.0471975511965976e+00 - acadoVariables.u[50];
acadoWorkspace.ub[51] = (real_t)1.0471975511965976e+00 - acadoVariables.u[51];
acadoWorkspace.ub[52] = (real_t)1.0471975511965976e+00 - acadoVariables.u[52];
acadoWorkspace.ub[53] = (real_t)1.0471975511965976e+00 - acadoVariables.u[53];
acadoWorkspace.ub[54] = (real_t)1.0471975511965976e+00 - acadoVariables.u[54];
acadoWorkspace.ub[55] = (real_t)1.0471975511965976e+00 - acadoVariables.u[55];
acadoWorkspace.ub[56] = (real_t)1.0471975511965976e+00 - acadoVariables.u[56];
acadoWorkspace.ub[57] = (real_t)1.0471975511965976e+00 - acadoVariables.u[57];
acadoWorkspace.ub[58] = (real_t)1.0471975511965976e+00 - acadoVariables.u[58];
acadoWorkspace.ub[59] = (real_t)1.0471975511965976e+00 - acadoVariables.u[59];

for (lRun1 = 0; lRun1 < 60; ++lRun1)
{
lRun3 = xBoundIndices[ lRun1 ] - 4;
lRun4 = ((lRun3) / (4)) + (1);
for (lRun2 = 0; lRun2 < lRun4; ++lRun2)
{
lRun5 = (((((lRun4) * (lRun4-1)) / (2)) + (lRun2)) * (4)) + ((lRun3) % (4));
acadoWorkspace.A[(lRun1 * 60) + (lRun2)] = acadoWorkspace.E[lRun5];
}
}

}

void acado_condenseFdb(  )
{
int lRun1;
int lRun2;
int lRun3;
real_t tmp;

acadoWorkspace.Dx0[0] = acadoVariables.x0[0] - acadoVariables.x[0];
acadoWorkspace.Dx0[1] = acadoVariables.x0[1] - acadoVariables.x[1];
acadoWorkspace.Dx0[2] = acadoVariables.x0[2] - acadoVariables.x[2];
acadoWorkspace.Dx0[3] = acadoVariables.x0[3] - acadoVariables.x[3];

acadoWorkspace.Dy[0] -= acadoVariables.y[0];
acadoWorkspace.Dy[1] -= acadoVariables.y[1];
acadoWorkspace.Dy[2] -= acadoVariables.y[2];
acadoWorkspace.Dy[3] -= acadoVariables.y[3];
acadoWorkspace.Dy[4] -= acadoVariables.y[4];
acadoWorkspace.Dy[5] -= acadoVariables.y[5];
acadoWorkspace.Dy[6] -= acadoVariables.y[6];
acadoWorkspace.Dy[7] -= acadoVariables.y[7];
acadoWorkspace.Dy[8] -= acadoVariables.y[8];
acadoWorkspace.Dy[9] -= acadoVariables.y[9];
acadoWorkspace.Dy[10] -= acadoVariables.y[10];
acadoWorkspace.Dy[11] -= acadoVariables.y[11];
acadoWorkspace.Dy[12] -= acadoVariables.y[12];
acadoWorkspace.Dy[13] -= acadoVariables.y[13];
acadoWorkspace.Dy[14] -= acadoVariables.y[14];
acadoWorkspace.Dy[15] -= acadoVariables.y[15];
acadoWorkspace.Dy[16] -= acadoVariables.y[16];
acadoWorkspace.Dy[17] -= acadoVariables.y[17];
acadoWorkspace.Dy[18] -= acadoVariables.y[18];
acadoWorkspace.Dy[19] -= acadoVariables.y[19];
acadoWorkspace.Dy[20] -= acadoVariables.y[20];
acadoWorkspace.Dy[21] -= acadoVariables.y[21];
acadoWorkspace.Dy[22] -= acadoVariables.y[22];
acadoWorkspace.Dy[23] -= acadoVariables.y[23];
acadoWorkspace.Dy[24] -= acadoVariables.y[24];
acadoWorkspace.Dy[25] -= acadoVariables.y[25];
acadoWorkspace.Dy[26] -= acadoVariables.y[26];
acadoWorkspace.Dy[27] -= acadoVariables.y[27];
acadoWorkspace.Dy[28] -= acadoVariables.y[28];
acadoWorkspace.Dy[29] -= acadoVariables.y[29];
acadoWorkspace.Dy[30] -= acadoVariables.y[30];
acadoWorkspace.Dy[31] -= acadoVariables.y[31];
acadoWorkspace.Dy[32] -= acadoVariables.y[32];
acadoWorkspace.Dy[33] -= acadoVariables.y[33];
acadoWorkspace.Dy[34] -= acadoVariables.y[34];
acadoWorkspace.Dy[35] -= acadoVariables.y[35];
acadoWorkspace.Dy[36] -= acadoVariables.y[36];
acadoWorkspace.Dy[37] -= acadoVariables.y[37];
acadoWorkspace.Dy[38] -= acadoVariables.y[38];
acadoWorkspace.Dy[39] -= acadoVariables.y[39];
acadoWorkspace.Dy[40] -= acadoVariables.y[40];
acadoWorkspace.Dy[41] -= acadoVariables.y[41];
acadoWorkspace.Dy[42] -= acadoVariables.y[42];
acadoWorkspace.Dy[43] -= acadoVariables.y[43];
acadoWorkspace.Dy[44] -= acadoVariables.y[44];
acadoWorkspace.Dy[45] -= acadoVariables.y[45];
acadoWorkspace.Dy[46] -= acadoVariables.y[46];
acadoWorkspace.Dy[47] -= acadoVariables.y[47];
acadoWorkspace.Dy[48] -= acadoVariables.y[48];
acadoWorkspace.Dy[49] -= acadoVariables.y[49];
acadoWorkspace.Dy[50] -= acadoVariables.y[50];
acadoWorkspace.Dy[51] -= acadoVariables.y[51];
acadoWorkspace.Dy[52] -= acadoVariables.y[52];
acadoWorkspace.Dy[53] -= acadoVariables.y[53];
acadoWorkspace.Dy[54] -= acadoVariables.y[54];
acadoWorkspace.Dy[55] -= acadoVariables.y[55];
acadoWorkspace.Dy[56] -= acadoVariables.y[56];
acadoWorkspace.Dy[57] -= acadoVariables.y[57];
acadoWorkspace.Dy[58] -= acadoVariables.y[58];
acadoWorkspace.Dy[59] -= acadoVariables.y[59];
acadoWorkspace.Dy[60] -= acadoVariables.y[60];
acadoWorkspace.Dy[61] -= acadoVariables.y[61];
acadoWorkspace.Dy[62] -= acadoVariables.y[62];
acadoWorkspace.Dy[63] -= acadoVariables.y[63];
acadoWorkspace.Dy[64] -= acadoVariables.y[64];
acadoWorkspace.Dy[65] -= acadoVariables.y[65];
acadoWorkspace.Dy[66] -= acadoVariables.y[66];
acadoWorkspace.Dy[67] -= acadoVariables.y[67];
acadoWorkspace.Dy[68] -= acadoVariables.y[68];
acadoWorkspace.Dy[69] -= acadoVariables.y[69];
acadoWorkspace.Dy[70] -= acadoVariables.y[70];
acadoWorkspace.Dy[71] -= acadoVariables.y[71];
acadoWorkspace.Dy[72] -= acadoVariables.y[72];
acadoWorkspace.Dy[73] -= acadoVariables.y[73];
acadoWorkspace.Dy[74] -= acadoVariables.y[74];
acadoWorkspace.Dy[75] -= acadoVariables.y[75];
acadoWorkspace.Dy[76] -= acadoVariables.y[76];
acadoWorkspace.Dy[77] -= acadoVariables.y[77];
acadoWorkspace.Dy[78] -= acadoVariables.y[78];
acadoWorkspace.Dy[79] -= acadoVariables.y[79];
acadoWorkspace.Dy[80] -= acadoVariables.y[80];
acadoWorkspace.Dy[81] -= acadoVariables.y[81];
acadoWorkspace.Dy[82] -= acadoVariables.y[82];
acadoWorkspace.Dy[83] -= acadoVariables.y[83];
acadoWorkspace.Dy[84] -= acadoVariables.y[84];
acadoWorkspace.Dy[85] -= acadoVariables.y[85];
acadoWorkspace.Dy[86] -= acadoVariables.y[86];
acadoWorkspace.Dy[87] -= acadoVariables.y[87];
acadoWorkspace.Dy[88] -= acadoVariables.y[88];
acadoWorkspace.Dy[89] -= acadoVariables.y[89];
acadoWorkspace.Dy[90] -= acadoVariables.y[90];
acadoWorkspace.Dy[91] -= acadoVariables.y[91];
acadoWorkspace.Dy[92] -= acadoVariables.y[92];
acadoWorkspace.Dy[93] -= acadoVariables.y[93];
acadoWorkspace.Dy[94] -= acadoVariables.y[94];
acadoWorkspace.Dy[95] -= acadoVariables.y[95];
acadoWorkspace.Dy[96] -= acadoVariables.y[96];
acadoWorkspace.Dy[97] -= acadoVariables.y[97];
acadoWorkspace.Dy[98] -= acadoVariables.y[98];
acadoWorkspace.Dy[99] -= acadoVariables.y[99];
acadoWorkspace.Dy[100] -= acadoVariables.y[100];
acadoWorkspace.Dy[101] -= acadoVariables.y[101];
acadoWorkspace.Dy[102] -= acadoVariables.y[102];
acadoWorkspace.Dy[103] -= acadoVariables.y[103];
acadoWorkspace.Dy[104] -= acadoVariables.y[104];
acadoWorkspace.Dy[105] -= acadoVariables.y[105];
acadoWorkspace.Dy[106] -= acadoVariables.y[106];
acadoWorkspace.Dy[107] -= acadoVariables.y[107];
acadoWorkspace.Dy[108] -= acadoVariables.y[108];
acadoWorkspace.Dy[109] -= acadoVariables.y[109];
acadoWorkspace.Dy[110] -= acadoVariables.y[110];
acadoWorkspace.Dy[111] -= acadoVariables.y[111];
acadoWorkspace.Dy[112] -= acadoVariables.y[112];
acadoWorkspace.Dy[113] -= acadoVariables.y[113];
acadoWorkspace.Dy[114] -= acadoVariables.y[114];
acadoWorkspace.Dy[115] -= acadoVariables.y[115];
acadoWorkspace.Dy[116] -= acadoVariables.y[116];
acadoWorkspace.Dy[117] -= acadoVariables.y[117];
acadoWorkspace.Dy[118] -= acadoVariables.y[118];
acadoWorkspace.Dy[119] -= acadoVariables.y[119];
acadoWorkspace.DyN[0] -= acadoVariables.yN[0];

acado_multRDy( acadoWorkspace.R2, acadoWorkspace.Dy, acadoWorkspace.g );
acado_multRDy( &(acadoWorkspace.R2[ 2 ]), &(acadoWorkspace.Dy[ 2 ]), &(acadoWorkspace.g[ 1 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 4 ]), &(acadoWorkspace.Dy[ 4 ]), &(acadoWorkspace.g[ 2 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 6 ]), &(acadoWorkspace.Dy[ 6 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 8 ]), &(acadoWorkspace.Dy[ 8 ]), &(acadoWorkspace.g[ 4 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 10 ]), &(acadoWorkspace.Dy[ 10 ]), &(acadoWorkspace.g[ 5 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 12 ]), &(acadoWorkspace.Dy[ 12 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 14 ]), &(acadoWorkspace.Dy[ 14 ]), &(acadoWorkspace.g[ 7 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 16 ]), &(acadoWorkspace.Dy[ 16 ]), &(acadoWorkspace.g[ 8 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 18 ]), &(acadoWorkspace.Dy[ 18 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 20 ]), &(acadoWorkspace.Dy[ 20 ]), &(acadoWorkspace.g[ 10 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 22 ]), &(acadoWorkspace.Dy[ 22 ]), &(acadoWorkspace.g[ 11 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 24 ]), &(acadoWorkspace.Dy[ 24 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 26 ]), &(acadoWorkspace.Dy[ 26 ]), &(acadoWorkspace.g[ 13 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 28 ]), &(acadoWorkspace.Dy[ 28 ]), &(acadoWorkspace.g[ 14 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 30 ]), &(acadoWorkspace.Dy[ 30 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 32 ]), &(acadoWorkspace.Dy[ 32 ]), &(acadoWorkspace.g[ 16 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 34 ]), &(acadoWorkspace.Dy[ 34 ]), &(acadoWorkspace.g[ 17 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 36 ]), &(acadoWorkspace.Dy[ 36 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 38 ]), &(acadoWorkspace.Dy[ 38 ]), &(acadoWorkspace.g[ 19 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 40 ]), &(acadoWorkspace.Dy[ 40 ]), &(acadoWorkspace.g[ 20 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 42 ]), &(acadoWorkspace.Dy[ 42 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 44 ]), &(acadoWorkspace.Dy[ 44 ]), &(acadoWorkspace.g[ 22 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 46 ]), &(acadoWorkspace.Dy[ 46 ]), &(acadoWorkspace.g[ 23 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 48 ]), &(acadoWorkspace.Dy[ 48 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 50 ]), &(acadoWorkspace.Dy[ 50 ]), &(acadoWorkspace.g[ 25 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 52 ]), &(acadoWorkspace.Dy[ 52 ]), &(acadoWorkspace.g[ 26 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 54 ]), &(acadoWorkspace.Dy[ 54 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 56 ]), &(acadoWorkspace.Dy[ 56 ]), &(acadoWorkspace.g[ 28 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 58 ]), &(acadoWorkspace.Dy[ 58 ]), &(acadoWorkspace.g[ 29 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 60 ]), &(acadoWorkspace.Dy[ 60 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 62 ]), &(acadoWorkspace.Dy[ 62 ]), &(acadoWorkspace.g[ 31 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 64 ]), &(acadoWorkspace.Dy[ 64 ]), &(acadoWorkspace.g[ 32 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 66 ]), &(acadoWorkspace.Dy[ 66 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 68 ]), &(acadoWorkspace.Dy[ 68 ]), &(acadoWorkspace.g[ 34 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 70 ]), &(acadoWorkspace.Dy[ 70 ]), &(acadoWorkspace.g[ 35 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 72 ]), &(acadoWorkspace.Dy[ 72 ]), &(acadoWorkspace.g[ 36 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 74 ]), &(acadoWorkspace.Dy[ 74 ]), &(acadoWorkspace.g[ 37 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 76 ]), &(acadoWorkspace.Dy[ 76 ]), &(acadoWorkspace.g[ 38 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 78 ]), &(acadoWorkspace.Dy[ 78 ]), &(acadoWorkspace.g[ 39 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 80 ]), &(acadoWorkspace.Dy[ 80 ]), &(acadoWorkspace.g[ 40 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 82 ]), &(acadoWorkspace.Dy[ 82 ]), &(acadoWorkspace.g[ 41 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 84 ]), &(acadoWorkspace.Dy[ 84 ]), &(acadoWorkspace.g[ 42 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 86 ]), &(acadoWorkspace.Dy[ 86 ]), &(acadoWorkspace.g[ 43 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 88 ]), &(acadoWorkspace.Dy[ 88 ]), &(acadoWorkspace.g[ 44 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 90 ]), &(acadoWorkspace.Dy[ 90 ]), &(acadoWorkspace.g[ 45 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 92 ]), &(acadoWorkspace.Dy[ 92 ]), &(acadoWorkspace.g[ 46 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 94 ]), &(acadoWorkspace.Dy[ 94 ]), &(acadoWorkspace.g[ 47 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 96 ]), &(acadoWorkspace.Dy[ 96 ]), &(acadoWorkspace.g[ 48 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 98 ]), &(acadoWorkspace.Dy[ 98 ]), &(acadoWorkspace.g[ 49 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 100 ]), &(acadoWorkspace.Dy[ 100 ]), &(acadoWorkspace.g[ 50 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 102 ]), &(acadoWorkspace.Dy[ 102 ]), &(acadoWorkspace.g[ 51 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 104 ]), &(acadoWorkspace.Dy[ 104 ]), &(acadoWorkspace.g[ 52 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 106 ]), &(acadoWorkspace.Dy[ 106 ]), &(acadoWorkspace.g[ 53 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 108 ]), &(acadoWorkspace.Dy[ 108 ]), &(acadoWorkspace.g[ 54 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 110 ]), &(acadoWorkspace.Dy[ 110 ]), &(acadoWorkspace.g[ 55 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 112 ]), &(acadoWorkspace.Dy[ 112 ]), &(acadoWorkspace.g[ 56 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 114 ]), &(acadoWorkspace.Dy[ 114 ]), &(acadoWorkspace.g[ 57 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 116 ]), &(acadoWorkspace.Dy[ 116 ]), &(acadoWorkspace.g[ 58 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 118 ]), &(acadoWorkspace.Dy[ 118 ]), &(acadoWorkspace.g[ 59 ]) );

acado_multQDy( acadoWorkspace.Q2, acadoWorkspace.Dy, acadoWorkspace.QDy );
acado_multQDy( &(acadoWorkspace.Q2[ 8 ]), &(acadoWorkspace.Dy[ 2 ]), &(acadoWorkspace.QDy[ 4 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 16 ]), &(acadoWorkspace.Dy[ 4 ]), &(acadoWorkspace.QDy[ 8 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 24 ]), &(acadoWorkspace.Dy[ 6 ]), &(acadoWorkspace.QDy[ 12 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 32 ]), &(acadoWorkspace.Dy[ 8 ]), &(acadoWorkspace.QDy[ 16 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 40 ]), &(acadoWorkspace.Dy[ 10 ]), &(acadoWorkspace.QDy[ 20 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 48 ]), &(acadoWorkspace.Dy[ 12 ]), &(acadoWorkspace.QDy[ 24 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 56 ]), &(acadoWorkspace.Dy[ 14 ]), &(acadoWorkspace.QDy[ 28 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 64 ]), &(acadoWorkspace.Dy[ 16 ]), &(acadoWorkspace.QDy[ 32 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 72 ]), &(acadoWorkspace.Dy[ 18 ]), &(acadoWorkspace.QDy[ 36 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 80 ]), &(acadoWorkspace.Dy[ 20 ]), &(acadoWorkspace.QDy[ 40 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 88 ]), &(acadoWorkspace.Dy[ 22 ]), &(acadoWorkspace.QDy[ 44 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 96 ]), &(acadoWorkspace.Dy[ 24 ]), &(acadoWorkspace.QDy[ 48 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 104 ]), &(acadoWorkspace.Dy[ 26 ]), &(acadoWorkspace.QDy[ 52 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 112 ]), &(acadoWorkspace.Dy[ 28 ]), &(acadoWorkspace.QDy[ 56 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 120 ]), &(acadoWorkspace.Dy[ 30 ]), &(acadoWorkspace.QDy[ 60 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 128 ]), &(acadoWorkspace.Dy[ 32 ]), &(acadoWorkspace.QDy[ 64 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 136 ]), &(acadoWorkspace.Dy[ 34 ]), &(acadoWorkspace.QDy[ 68 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 144 ]), &(acadoWorkspace.Dy[ 36 ]), &(acadoWorkspace.QDy[ 72 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 152 ]), &(acadoWorkspace.Dy[ 38 ]), &(acadoWorkspace.QDy[ 76 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 160 ]), &(acadoWorkspace.Dy[ 40 ]), &(acadoWorkspace.QDy[ 80 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 168 ]), &(acadoWorkspace.Dy[ 42 ]), &(acadoWorkspace.QDy[ 84 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 176 ]), &(acadoWorkspace.Dy[ 44 ]), &(acadoWorkspace.QDy[ 88 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 184 ]), &(acadoWorkspace.Dy[ 46 ]), &(acadoWorkspace.QDy[ 92 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 192 ]), &(acadoWorkspace.Dy[ 48 ]), &(acadoWorkspace.QDy[ 96 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 200 ]), &(acadoWorkspace.Dy[ 50 ]), &(acadoWorkspace.QDy[ 100 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 208 ]), &(acadoWorkspace.Dy[ 52 ]), &(acadoWorkspace.QDy[ 104 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 216 ]), &(acadoWorkspace.Dy[ 54 ]), &(acadoWorkspace.QDy[ 108 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 224 ]), &(acadoWorkspace.Dy[ 56 ]), &(acadoWorkspace.QDy[ 112 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 232 ]), &(acadoWorkspace.Dy[ 58 ]), &(acadoWorkspace.QDy[ 116 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 240 ]), &(acadoWorkspace.Dy[ 60 ]), &(acadoWorkspace.QDy[ 120 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 248 ]), &(acadoWorkspace.Dy[ 62 ]), &(acadoWorkspace.QDy[ 124 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 256 ]), &(acadoWorkspace.Dy[ 64 ]), &(acadoWorkspace.QDy[ 128 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 264 ]), &(acadoWorkspace.Dy[ 66 ]), &(acadoWorkspace.QDy[ 132 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 272 ]), &(acadoWorkspace.Dy[ 68 ]), &(acadoWorkspace.QDy[ 136 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 280 ]), &(acadoWorkspace.Dy[ 70 ]), &(acadoWorkspace.QDy[ 140 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 288 ]), &(acadoWorkspace.Dy[ 72 ]), &(acadoWorkspace.QDy[ 144 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 296 ]), &(acadoWorkspace.Dy[ 74 ]), &(acadoWorkspace.QDy[ 148 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 304 ]), &(acadoWorkspace.Dy[ 76 ]), &(acadoWorkspace.QDy[ 152 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 312 ]), &(acadoWorkspace.Dy[ 78 ]), &(acadoWorkspace.QDy[ 156 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 320 ]), &(acadoWorkspace.Dy[ 80 ]), &(acadoWorkspace.QDy[ 160 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 328 ]), &(acadoWorkspace.Dy[ 82 ]), &(acadoWorkspace.QDy[ 164 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 336 ]), &(acadoWorkspace.Dy[ 84 ]), &(acadoWorkspace.QDy[ 168 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 344 ]), &(acadoWorkspace.Dy[ 86 ]), &(acadoWorkspace.QDy[ 172 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 352 ]), &(acadoWorkspace.Dy[ 88 ]), &(acadoWorkspace.QDy[ 176 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 360 ]), &(acadoWorkspace.Dy[ 90 ]), &(acadoWorkspace.QDy[ 180 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 368 ]), &(acadoWorkspace.Dy[ 92 ]), &(acadoWorkspace.QDy[ 184 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 376 ]), &(acadoWorkspace.Dy[ 94 ]), &(acadoWorkspace.QDy[ 188 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 384 ]), &(acadoWorkspace.Dy[ 96 ]), &(acadoWorkspace.QDy[ 192 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 392 ]), &(acadoWorkspace.Dy[ 98 ]), &(acadoWorkspace.QDy[ 196 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 400 ]), &(acadoWorkspace.Dy[ 100 ]), &(acadoWorkspace.QDy[ 200 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 408 ]), &(acadoWorkspace.Dy[ 102 ]), &(acadoWorkspace.QDy[ 204 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 416 ]), &(acadoWorkspace.Dy[ 104 ]), &(acadoWorkspace.QDy[ 208 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 424 ]), &(acadoWorkspace.Dy[ 106 ]), &(acadoWorkspace.QDy[ 212 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 432 ]), &(acadoWorkspace.Dy[ 108 ]), &(acadoWorkspace.QDy[ 216 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 440 ]), &(acadoWorkspace.Dy[ 110 ]), &(acadoWorkspace.QDy[ 220 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 448 ]), &(acadoWorkspace.Dy[ 112 ]), &(acadoWorkspace.QDy[ 224 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 456 ]), &(acadoWorkspace.Dy[ 114 ]), &(acadoWorkspace.QDy[ 228 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 464 ]), &(acadoWorkspace.Dy[ 116 ]), &(acadoWorkspace.QDy[ 232 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 472 ]), &(acadoWorkspace.Dy[ 118 ]), &(acadoWorkspace.QDy[ 236 ]) );

acadoWorkspace.QDy[240] = + acadoWorkspace.QN2[0]*acadoWorkspace.DyN[0];
acadoWorkspace.QDy[241] = + acadoWorkspace.QN2[1]*acadoWorkspace.DyN[0];
acadoWorkspace.QDy[242] = + acadoWorkspace.QN2[2]*acadoWorkspace.DyN[0];
acadoWorkspace.QDy[243] = + acadoWorkspace.QN2[3]*acadoWorkspace.DyN[0];

for (lRun2 = 0; lRun2 < 240; ++lRun2)
acadoWorkspace.QDy[lRun2 + 4] += acadoWorkspace.Qd[lRun2];


for (lRun1 = 0; lRun1 < 60; ++lRun1)
{
for (lRun2 = lRun1; lRun2 < 60; ++lRun2)
{
lRun3 = (((lRun2 + 1) * (lRun2)) / (2)) + (lRun1);
acado_multEQDy( &(acadoWorkspace.E[ lRun3 * 4 ]), &(acadoWorkspace.QDy[ lRun2 * 4 + 4 ]), &(acadoWorkspace.g[ lRun1 ]) );
}
}

acadoWorkspace.g[0] += + acadoWorkspace.H10[0]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[1]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[2]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[3]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[1] += + acadoWorkspace.H10[4]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[5]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[6]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[7]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[2] += + acadoWorkspace.H10[8]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[9]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[10]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[11]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[3] += + acadoWorkspace.H10[12]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[13]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[14]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[15]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[4] += + acadoWorkspace.H10[16]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[17]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[18]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[19]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[5] += + acadoWorkspace.H10[20]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[21]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[22]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[23]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[6] += + acadoWorkspace.H10[24]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[25]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[26]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[27]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[7] += + acadoWorkspace.H10[28]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[29]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[30]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[31]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[8] += + acadoWorkspace.H10[32]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[33]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[34]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[35]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[9] += + acadoWorkspace.H10[36]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[37]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[38]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[39]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[10] += + acadoWorkspace.H10[40]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[41]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[42]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[43]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[11] += + acadoWorkspace.H10[44]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[45]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[46]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[47]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[12] += + acadoWorkspace.H10[48]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[49]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[50]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[51]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[13] += + acadoWorkspace.H10[52]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[53]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[54]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[55]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[14] += + acadoWorkspace.H10[56]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[57]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[58]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[59]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[15] += + acadoWorkspace.H10[60]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[61]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[62]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[63]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[16] += + acadoWorkspace.H10[64]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[65]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[66]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[67]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[17] += + acadoWorkspace.H10[68]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[69]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[70]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[71]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[18] += + acadoWorkspace.H10[72]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[73]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[74]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[75]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[19] += + acadoWorkspace.H10[76]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[77]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[78]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[79]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[20] += + acadoWorkspace.H10[80]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[81]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[82]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[83]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[21] += + acadoWorkspace.H10[84]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[85]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[86]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[87]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[22] += + acadoWorkspace.H10[88]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[89]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[90]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[91]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[23] += + acadoWorkspace.H10[92]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[93]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[94]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[95]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[24] += + acadoWorkspace.H10[96]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[97]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[98]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[99]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[25] += + acadoWorkspace.H10[100]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[101]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[102]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[103]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[26] += + acadoWorkspace.H10[104]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[105]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[106]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[107]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[27] += + acadoWorkspace.H10[108]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[109]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[110]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[111]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[28] += + acadoWorkspace.H10[112]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[113]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[114]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[115]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[29] += + acadoWorkspace.H10[116]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[117]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[118]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[119]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[30] += + acadoWorkspace.H10[120]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[121]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[122]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[123]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[31] += + acadoWorkspace.H10[124]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[125]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[126]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[127]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[32] += + acadoWorkspace.H10[128]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[129]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[130]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[131]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[33] += + acadoWorkspace.H10[132]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[133]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[134]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[135]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[34] += + acadoWorkspace.H10[136]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[137]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[138]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[139]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[35] += + acadoWorkspace.H10[140]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[141]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[142]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[143]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[36] += + acadoWorkspace.H10[144]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[145]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[146]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[147]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[37] += + acadoWorkspace.H10[148]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[149]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[150]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[151]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[38] += + acadoWorkspace.H10[152]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[153]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[154]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[155]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[39] += + acadoWorkspace.H10[156]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[157]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[158]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[159]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[40] += + acadoWorkspace.H10[160]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[161]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[162]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[163]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[41] += + acadoWorkspace.H10[164]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[165]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[166]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[167]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[42] += + acadoWorkspace.H10[168]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[169]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[170]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[171]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[43] += + acadoWorkspace.H10[172]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[173]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[174]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[175]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[44] += + acadoWorkspace.H10[176]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[177]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[178]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[179]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[45] += + acadoWorkspace.H10[180]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[181]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[182]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[183]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[46] += + acadoWorkspace.H10[184]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[185]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[186]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[187]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[47] += + acadoWorkspace.H10[188]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[189]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[190]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[191]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[48] += + acadoWorkspace.H10[192]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[193]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[194]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[195]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[49] += + acadoWorkspace.H10[196]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[197]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[198]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[199]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[50] += + acadoWorkspace.H10[200]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[201]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[202]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[203]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[51] += + acadoWorkspace.H10[204]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[205]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[206]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[207]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[52] += + acadoWorkspace.H10[208]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[209]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[210]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[211]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[53] += + acadoWorkspace.H10[212]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[213]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[214]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[215]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[54] += + acadoWorkspace.H10[216]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[217]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[218]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[219]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[55] += + acadoWorkspace.H10[220]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[221]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[222]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[223]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[56] += + acadoWorkspace.H10[224]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[225]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[226]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[227]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[57] += + acadoWorkspace.H10[228]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[229]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[230]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[231]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[58] += + acadoWorkspace.H10[232]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[233]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[234]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[235]*acadoWorkspace.Dx0[3];
acadoWorkspace.g[59] += + acadoWorkspace.H10[236]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[237]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[238]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[239]*acadoWorkspace.Dx0[3];

tmp = + acadoWorkspace.evGx[12]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[13]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[14]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[15]*acadoWorkspace.Dx0[3] + acadoVariables.x[7];
tmp += acadoWorkspace.d[3];
acadoWorkspace.lbA[0] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[0] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[28]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[29]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[30]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[31]*acadoWorkspace.Dx0[3] + acadoVariables.x[11];
tmp += acadoWorkspace.d[7];
acadoWorkspace.lbA[1] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[1] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[44]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[45]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[46]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[47]*acadoWorkspace.Dx0[3] + acadoVariables.x[15];
tmp += acadoWorkspace.d[11];
acadoWorkspace.lbA[2] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[2] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[60]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[61]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[62]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[63]*acadoWorkspace.Dx0[3] + acadoVariables.x[19];
tmp += acadoWorkspace.d[15];
acadoWorkspace.lbA[3] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[3] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[76]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[77]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[78]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[79]*acadoWorkspace.Dx0[3] + acadoVariables.x[23];
tmp += acadoWorkspace.d[19];
acadoWorkspace.lbA[4] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[4] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[92]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[93]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[94]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[95]*acadoWorkspace.Dx0[3] + acadoVariables.x[27];
tmp += acadoWorkspace.d[23];
acadoWorkspace.lbA[5] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[5] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[108]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[109]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[110]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[111]*acadoWorkspace.Dx0[3] + acadoVariables.x[31];
tmp += acadoWorkspace.d[27];
acadoWorkspace.lbA[6] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[6] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[124]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[125]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[126]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[127]*acadoWorkspace.Dx0[3] + acadoVariables.x[35];
tmp += acadoWorkspace.d[31];
acadoWorkspace.lbA[7] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[7] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[140]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[141]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[142]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[143]*acadoWorkspace.Dx0[3] + acadoVariables.x[39];
tmp += acadoWorkspace.d[35];
acadoWorkspace.lbA[8] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[8] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[156]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[157]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[158]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[159]*acadoWorkspace.Dx0[3] + acadoVariables.x[43];
tmp += acadoWorkspace.d[39];
acadoWorkspace.lbA[9] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[9] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[172]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[173]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[174]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[175]*acadoWorkspace.Dx0[3] + acadoVariables.x[47];
tmp += acadoWorkspace.d[43];
acadoWorkspace.lbA[10] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[10] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[188]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[189]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[190]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[191]*acadoWorkspace.Dx0[3] + acadoVariables.x[51];
tmp += acadoWorkspace.d[47];
acadoWorkspace.lbA[11] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[11] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[204]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[205]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[206]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[207]*acadoWorkspace.Dx0[3] + acadoVariables.x[55];
tmp += acadoWorkspace.d[51];
acadoWorkspace.lbA[12] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[12] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[220]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[221]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[222]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[223]*acadoWorkspace.Dx0[3] + acadoVariables.x[59];
tmp += acadoWorkspace.d[55];
acadoWorkspace.lbA[13] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[13] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[236]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[237]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[238]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[239]*acadoWorkspace.Dx0[3] + acadoVariables.x[63];
tmp += acadoWorkspace.d[59];
acadoWorkspace.lbA[14] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[14] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[252]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[253]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[254]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[255]*acadoWorkspace.Dx0[3] + acadoVariables.x[67];
tmp += acadoWorkspace.d[63];
acadoWorkspace.lbA[15] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[15] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[268]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[269]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[270]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[271]*acadoWorkspace.Dx0[3] + acadoVariables.x[71];
tmp += acadoWorkspace.d[67];
acadoWorkspace.lbA[16] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[16] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[284]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[285]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[286]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[287]*acadoWorkspace.Dx0[3] + acadoVariables.x[75];
tmp += acadoWorkspace.d[71];
acadoWorkspace.lbA[17] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[17] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[300]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[301]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[302]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[303]*acadoWorkspace.Dx0[3] + acadoVariables.x[79];
tmp += acadoWorkspace.d[75];
acadoWorkspace.lbA[18] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[18] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[316]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[317]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[318]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[319]*acadoWorkspace.Dx0[3] + acadoVariables.x[83];
tmp += acadoWorkspace.d[79];
acadoWorkspace.lbA[19] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[19] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[332]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[333]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[334]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[335]*acadoWorkspace.Dx0[3] + acadoVariables.x[87];
tmp += acadoWorkspace.d[83];
acadoWorkspace.lbA[20] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[20] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[348]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[349]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[350]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[351]*acadoWorkspace.Dx0[3] + acadoVariables.x[91];
tmp += acadoWorkspace.d[87];
acadoWorkspace.lbA[21] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[21] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[364]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[365]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[366]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[367]*acadoWorkspace.Dx0[3] + acadoVariables.x[95];
tmp += acadoWorkspace.d[91];
acadoWorkspace.lbA[22] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[22] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[380]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[381]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[382]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[383]*acadoWorkspace.Dx0[3] + acadoVariables.x[99];
tmp += acadoWorkspace.d[95];
acadoWorkspace.lbA[23] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[23] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[396]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[397]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[398]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[399]*acadoWorkspace.Dx0[3] + acadoVariables.x[103];
tmp += acadoWorkspace.d[99];
acadoWorkspace.lbA[24] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[24] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[412]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[413]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[414]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[415]*acadoWorkspace.Dx0[3] + acadoVariables.x[107];
tmp += acadoWorkspace.d[103];
acadoWorkspace.lbA[25] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[25] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[428]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[429]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[430]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[431]*acadoWorkspace.Dx0[3] + acadoVariables.x[111];
tmp += acadoWorkspace.d[107];
acadoWorkspace.lbA[26] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[26] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[444]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[445]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[446]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[447]*acadoWorkspace.Dx0[3] + acadoVariables.x[115];
tmp += acadoWorkspace.d[111];
acadoWorkspace.lbA[27] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[27] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[460]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[461]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[462]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[463]*acadoWorkspace.Dx0[3] + acadoVariables.x[119];
tmp += acadoWorkspace.d[115];
acadoWorkspace.lbA[28] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[28] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[476]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[477]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[478]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[479]*acadoWorkspace.Dx0[3] + acadoVariables.x[123];
tmp += acadoWorkspace.d[119];
acadoWorkspace.lbA[29] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[29] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[492]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[493]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[494]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[495]*acadoWorkspace.Dx0[3] + acadoVariables.x[127];
tmp += acadoWorkspace.d[123];
acadoWorkspace.lbA[30] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[30] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[508]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[509]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[510]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[511]*acadoWorkspace.Dx0[3] + acadoVariables.x[131];
tmp += acadoWorkspace.d[127];
acadoWorkspace.lbA[31] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[31] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[524]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[525]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[526]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[527]*acadoWorkspace.Dx0[3] + acadoVariables.x[135];
tmp += acadoWorkspace.d[131];
acadoWorkspace.lbA[32] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[32] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[540]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[541]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[542]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[543]*acadoWorkspace.Dx0[3] + acadoVariables.x[139];
tmp += acadoWorkspace.d[135];
acadoWorkspace.lbA[33] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[33] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[556]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[557]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[558]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[559]*acadoWorkspace.Dx0[3] + acadoVariables.x[143];
tmp += acadoWorkspace.d[139];
acadoWorkspace.lbA[34] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[34] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[572]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[573]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[574]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[575]*acadoWorkspace.Dx0[3] + acadoVariables.x[147];
tmp += acadoWorkspace.d[143];
acadoWorkspace.lbA[35] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[35] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[588]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[589]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[590]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[591]*acadoWorkspace.Dx0[3] + acadoVariables.x[151];
tmp += acadoWorkspace.d[147];
acadoWorkspace.lbA[36] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[36] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[604]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[605]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[606]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[607]*acadoWorkspace.Dx0[3] + acadoVariables.x[155];
tmp += acadoWorkspace.d[151];
acadoWorkspace.lbA[37] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[37] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[620]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[621]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[622]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[623]*acadoWorkspace.Dx0[3] + acadoVariables.x[159];
tmp += acadoWorkspace.d[155];
acadoWorkspace.lbA[38] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[38] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[636]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[637]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[638]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[639]*acadoWorkspace.Dx0[3] + acadoVariables.x[163];
tmp += acadoWorkspace.d[159];
acadoWorkspace.lbA[39] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[39] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[652]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[653]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[654]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[655]*acadoWorkspace.Dx0[3] + acadoVariables.x[167];
tmp += acadoWorkspace.d[163];
acadoWorkspace.lbA[40] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[40] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[668]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[669]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[670]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[671]*acadoWorkspace.Dx0[3] + acadoVariables.x[171];
tmp += acadoWorkspace.d[167];
acadoWorkspace.lbA[41] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[41] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[684]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[685]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[686]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[687]*acadoWorkspace.Dx0[3] + acadoVariables.x[175];
tmp += acadoWorkspace.d[171];
acadoWorkspace.lbA[42] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[42] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[700]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[701]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[702]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[703]*acadoWorkspace.Dx0[3] + acadoVariables.x[179];
tmp += acadoWorkspace.d[175];
acadoWorkspace.lbA[43] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[43] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[716]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[717]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[718]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[719]*acadoWorkspace.Dx0[3] + acadoVariables.x[183];
tmp += acadoWorkspace.d[179];
acadoWorkspace.lbA[44] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[44] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[732]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[733]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[734]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[735]*acadoWorkspace.Dx0[3] + acadoVariables.x[187];
tmp += acadoWorkspace.d[183];
acadoWorkspace.lbA[45] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[45] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[748]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[749]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[750]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[751]*acadoWorkspace.Dx0[3] + acadoVariables.x[191];
tmp += acadoWorkspace.d[187];
acadoWorkspace.lbA[46] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[46] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[764]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[765]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[766]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[767]*acadoWorkspace.Dx0[3] + acadoVariables.x[195];
tmp += acadoWorkspace.d[191];
acadoWorkspace.lbA[47] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[47] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[780]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[781]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[782]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[783]*acadoWorkspace.Dx0[3] + acadoVariables.x[199];
tmp += acadoWorkspace.d[195];
acadoWorkspace.lbA[48] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[48] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[796]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[797]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[798]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[799]*acadoWorkspace.Dx0[3] + acadoVariables.x[203];
tmp += acadoWorkspace.d[199];
acadoWorkspace.lbA[49] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[49] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[812]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[813]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[814]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[815]*acadoWorkspace.Dx0[3] + acadoVariables.x[207];
tmp += acadoWorkspace.d[203];
acadoWorkspace.lbA[50] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[50] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[828]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[829]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[830]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[831]*acadoWorkspace.Dx0[3] + acadoVariables.x[211];
tmp += acadoWorkspace.d[207];
acadoWorkspace.lbA[51] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[51] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[844]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[845]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[846]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[847]*acadoWorkspace.Dx0[3] + acadoVariables.x[215];
tmp += acadoWorkspace.d[211];
acadoWorkspace.lbA[52] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[52] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[860]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[861]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[862]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[863]*acadoWorkspace.Dx0[3] + acadoVariables.x[219];
tmp += acadoWorkspace.d[215];
acadoWorkspace.lbA[53] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[53] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[876]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[877]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[878]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[879]*acadoWorkspace.Dx0[3] + acadoVariables.x[223];
tmp += acadoWorkspace.d[219];
acadoWorkspace.lbA[54] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[54] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[892]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[893]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[894]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[895]*acadoWorkspace.Dx0[3] + acadoVariables.x[227];
tmp += acadoWorkspace.d[223];
acadoWorkspace.lbA[55] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[55] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[908]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[909]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[910]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[911]*acadoWorkspace.Dx0[3] + acadoVariables.x[231];
tmp += acadoWorkspace.d[227];
acadoWorkspace.lbA[56] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[56] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[924]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[925]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[926]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[927]*acadoWorkspace.Dx0[3] + acadoVariables.x[235];
tmp += acadoWorkspace.d[231];
acadoWorkspace.lbA[57] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[57] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[940]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[941]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[942]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[943]*acadoWorkspace.Dx0[3] + acadoVariables.x[239];
tmp += acadoWorkspace.d[235];
acadoWorkspace.lbA[58] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[58] = (real_t)5.2359877559829882e-01 - tmp;
tmp = + acadoWorkspace.evGx[956]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[957]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[958]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[959]*acadoWorkspace.Dx0[3] + acadoVariables.x[243];
tmp += acadoWorkspace.d[239];
acadoWorkspace.lbA[59] = (real_t)-5.2359877559829882e-01 - tmp;
acadoWorkspace.ubA[59] = (real_t)5.2359877559829882e-01 - tmp;

}

void acado_expand(  )
{
int lRun1;
int lRun2;
int lRun3;
acadoVariables.u[0] += acadoWorkspace.x[0];
acadoVariables.u[1] += acadoWorkspace.x[1];
acadoVariables.u[2] += acadoWorkspace.x[2];
acadoVariables.u[3] += acadoWorkspace.x[3];
acadoVariables.u[4] += acadoWorkspace.x[4];
acadoVariables.u[5] += acadoWorkspace.x[5];
acadoVariables.u[6] += acadoWorkspace.x[6];
acadoVariables.u[7] += acadoWorkspace.x[7];
acadoVariables.u[8] += acadoWorkspace.x[8];
acadoVariables.u[9] += acadoWorkspace.x[9];
acadoVariables.u[10] += acadoWorkspace.x[10];
acadoVariables.u[11] += acadoWorkspace.x[11];
acadoVariables.u[12] += acadoWorkspace.x[12];
acadoVariables.u[13] += acadoWorkspace.x[13];
acadoVariables.u[14] += acadoWorkspace.x[14];
acadoVariables.u[15] += acadoWorkspace.x[15];
acadoVariables.u[16] += acadoWorkspace.x[16];
acadoVariables.u[17] += acadoWorkspace.x[17];
acadoVariables.u[18] += acadoWorkspace.x[18];
acadoVariables.u[19] += acadoWorkspace.x[19];
acadoVariables.u[20] += acadoWorkspace.x[20];
acadoVariables.u[21] += acadoWorkspace.x[21];
acadoVariables.u[22] += acadoWorkspace.x[22];
acadoVariables.u[23] += acadoWorkspace.x[23];
acadoVariables.u[24] += acadoWorkspace.x[24];
acadoVariables.u[25] += acadoWorkspace.x[25];
acadoVariables.u[26] += acadoWorkspace.x[26];
acadoVariables.u[27] += acadoWorkspace.x[27];
acadoVariables.u[28] += acadoWorkspace.x[28];
acadoVariables.u[29] += acadoWorkspace.x[29];
acadoVariables.u[30] += acadoWorkspace.x[30];
acadoVariables.u[31] += acadoWorkspace.x[31];
acadoVariables.u[32] += acadoWorkspace.x[32];
acadoVariables.u[33] += acadoWorkspace.x[33];
acadoVariables.u[34] += acadoWorkspace.x[34];
acadoVariables.u[35] += acadoWorkspace.x[35];
acadoVariables.u[36] += acadoWorkspace.x[36];
acadoVariables.u[37] += acadoWorkspace.x[37];
acadoVariables.u[38] += acadoWorkspace.x[38];
acadoVariables.u[39] += acadoWorkspace.x[39];
acadoVariables.u[40] += acadoWorkspace.x[40];
acadoVariables.u[41] += acadoWorkspace.x[41];
acadoVariables.u[42] += acadoWorkspace.x[42];
acadoVariables.u[43] += acadoWorkspace.x[43];
acadoVariables.u[44] += acadoWorkspace.x[44];
acadoVariables.u[45] += acadoWorkspace.x[45];
acadoVariables.u[46] += acadoWorkspace.x[46];
acadoVariables.u[47] += acadoWorkspace.x[47];
acadoVariables.u[48] += acadoWorkspace.x[48];
acadoVariables.u[49] += acadoWorkspace.x[49];
acadoVariables.u[50] += acadoWorkspace.x[50];
acadoVariables.u[51] += acadoWorkspace.x[51];
acadoVariables.u[52] += acadoWorkspace.x[52];
acadoVariables.u[53] += acadoWorkspace.x[53];
acadoVariables.u[54] += acadoWorkspace.x[54];
acadoVariables.u[55] += acadoWorkspace.x[55];
acadoVariables.u[56] += acadoWorkspace.x[56];
acadoVariables.u[57] += acadoWorkspace.x[57];
acadoVariables.u[58] += acadoWorkspace.x[58];
acadoVariables.u[59] += acadoWorkspace.x[59];

acadoVariables.x[0] += acadoWorkspace.Dx0[0];
acadoVariables.x[1] += acadoWorkspace.Dx0[1];
acadoVariables.x[2] += acadoWorkspace.Dx0[2];
acadoVariables.x[3] += acadoWorkspace.Dx0[3];

acadoVariables.x[4] += + acadoWorkspace.evGx[0]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[1]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[2]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[3]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[0];
acadoVariables.x[5] += + acadoWorkspace.evGx[4]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[5]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[6]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[7]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[1];
acadoVariables.x[6] += + acadoWorkspace.evGx[8]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[9]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[10]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[11]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[2];
acadoVariables.x[7] += + acadoWorkspace.evGx[12]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[13]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[14]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[15]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[3];
acadoVariables.x[8] += + acadoWorkspace.evGx[16]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[17]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[18]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[19]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[4];
acadoVariables.x[9] += + acadoWorkspace.evGx[20]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[21]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[22]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[23]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[5];
acadoVariables.x[10] += + acadoWorkspace.evGx[24]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[25]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[26]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[27]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[6];
acadoVariables.x[11] += + acadoWorkspace.evGx[28]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[29]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[30]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[31]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[7];
acadoVariables.x[12] += + acadoWorkspace.evGx[32]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[33]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[34]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[35]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[8];
acadoVariables.x[13] += + acadoWorkspace.evGx[36]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[37]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[38]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[39]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[9];
acadoVariables.x[14] += + acadoWorkspace.evGx[40]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[41]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[42]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[43]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[10];
acadoVariables.x[15] += + acadoWorkspace.evGx[44]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[45]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[46]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[47]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[11];
acadoVariables.x[16] += + acadoWorkspace.evGx[48]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[49]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[50]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[51]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[12];
acadoVariables.x[17] += + acadoWorkspace.evGx[52]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[53]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[54]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[55]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[13];
acadoVariables.x[18] += + acadoWorkspace.evGx[56]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[57]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[58]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[59]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[14];
acadoVariables.x[19] += + acadoWorkspace.evGx[60]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[61]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[62]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[63]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[15];
acadoVariables.x[20] += + acadoWorkspace.evGx[64]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[65]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[66]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[67]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[16];
acadoVariables.x[21] += + acadoWorkspace.evGx[68]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[69]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[70]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[71]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[17];
acadoVariables.x[22] += + acadoWorkspace.evGx[72]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[73]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[74]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[75]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[18];
acadoVariables.x[23] += + acadoWorkspace.evGx[76]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[77]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[78]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[79]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[19];
acadoVariables.x[24] += + acadoWorkspace.evGx[80]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[81]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[82]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[83]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[20];
acadoVariables.x[25] += + acadoWorkspace.evGx[84]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[85]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[86]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[87]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[21];
acadoVariables.x[26] += + acadoWorkspace.evGx[88]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[89]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[90]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[91]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[22];
acadoVariables.x[27] += + acadoWorkspace.evGx[92]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[93]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[94]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[95]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[23];
acadoVariables.x[28] += + acadoWorkspace.evGx[96]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[97]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[98]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[99]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[24];
acadoVariables.x[29] += + acadoWorkspace.evGx[100]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[101]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[102]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[103]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[25];
acadoVariables.x[30] += + acadoWorkspace.evGx[104]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[105]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[106]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[107]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[26];
acadoVariables.x[31] += + acadoWorkspace.evGx[108]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[109]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[110]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[111]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[27];
acadoVariables.x[32] += + acadoWorkspace.evGx[112]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[113]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[114]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[115]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[28];
acadoVariables.x[33] += + acadoWorkspace.evGx[116]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[117]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[118]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[119]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[29];
acadoVariables.x[34] += + acadoWorkspace.evGx[120]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[121]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[122]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[123]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[30];
acadoVariables.x[35] += + acadoWorkspace.evGx[124]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[125]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[126]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[127]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[31];
acadoVariables.x[36] += + acadoWorkspace.evGx[128]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[129]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[130]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[131]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[32];
acadoVariables.x[37] += + acadoWorkspace.evGx[132]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[133]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[134]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[135]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[33];
acadoVariables.x[38] += + acadoWorkspace.evGx[136]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[137]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[138]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[139]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[34];
acadoVariables.x[39] += + acadoWorkspace.evGx[140]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[141]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[142]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[143]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[35];
acadoVariables.x[40] += + acadoWorkspace.evGx[144]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[145]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[146]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[147]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[36];
acadoVariables.x[41] += + acadoWorkspace.evGx[148]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[149]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[150]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[151]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[37];
acadoVariables.x[42] += + acadoWorkspace.evGx[152]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[153]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[154]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[155]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[38];
acadoVariables.x[43] += + acadoWorkspace.evGx[156]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[157]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[158]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[159]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[39];
acadoVariables.x[44] += + acadoWorkspace.evGx[160]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[161]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[162]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[163]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[40];
acadoVariables.x[45] += + acadoWorkspace.evGx[164]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[165]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[166]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[167]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[41];
acadoVariables.x[46] += + acadoWorkspace.evGx[168]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[169]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[170]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[171]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[42];
acadoVariables.x[47] += + acadoWorkspace.evGx[172]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[173]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[174]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[175]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[43];
acadoVariables.x[48] += + acadoWorkspace.evGx[176]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[177]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[178]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[179]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[44];
acadoVariables.x[49] += + acadoWorkspace.evGx[180]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[181]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[182]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[183]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[45];
acadoVariables.x[50] += + acadoWorkspace.evGx[184]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[185]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[186]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[187]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[46];
acadoVariables.x[51] += + acadoWorkspace.evGx[188]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[189]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[190]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[191]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[47];
acadoVariables.x[52] += + acadoWorkspace.evGx[192]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[193]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[194]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[195]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[48];
acadoVariables.x[53] += + acadoWorkspace.evGx[196]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[197]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[198]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[199]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[49];
acadoVariables.x[54] += + acadoWorkspace.evGx[200]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[201]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[202]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[203]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[50];
acadoVariables.x[55] += + acadoWorkspace.evGx[204]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[205]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[206]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[207]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[51];
acadoVariables.x[56] += + acadoWorkspace.evGx[208]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[209]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[210]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[211]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[52];
acadoVariables.x[57] += + acadoWorkspace.evGx[212]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[213]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[214]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[215]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[53];
acadoVariables.x[58] += + acadoWorkspace.evGx[216]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[217]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[218]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[219]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[54];
acadoVariables.x[59] += + acadoWorkspace.evGx[220]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[221]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[222]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[223]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[55];
acadoVariables.x[60] += + acadoWorkspace.evGx[224]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[225]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[226]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[227]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[56];
acadoVariables.x[61] += + acadoWorkspace.evGx[228]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[229]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[230]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[231]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[57];
acadoVariables.x[62] += + acadoWorkspace.evGx[232]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[233]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[234]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[235]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[58];
acadoVariables.x[63] += + acadoWorkspace.evGx[236]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[237]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[238]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[239]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[59];
acadoVariables.x[64] += + acadoWorkspace.evGx[240]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[241]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[242]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[243]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[60];
acadoVariables.x[65] += + acadoWorkspace.evGx[244]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[245]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[246]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[247]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[61];
acadoVariables.x[66] += + acadoWorkspace.evGx[248]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[249]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[250]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[251]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[62];
acadoVariables.x[67] += + acadoWorkspace.evGx[252]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[253]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[254]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[255]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[63];
acadoVariables.x[68] += + acadoWorkspace.evGx[256]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[257]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[258]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[259]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[64];
acadoVariables.x[69] += + acadoWorkspace.evGx[260]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[261]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[262]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[263]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[65];
acadoVariables.x[70] += + acadoWorkspace.evGx[264]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[265]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[266]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[267]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[66];
acadoVariables.x[71] += + acadoWorkspace.evGx[268]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[269]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[270]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[271]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[67];
acadoVariables.x[72] += + acadoWorkspace.evGx[272]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[273]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[274]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[275]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[68];
acadoVariables.x[73] += + acadoWorkspace.evGx[276]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[277]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[278]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[279]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[69];
acadoVariables.x[74] += + acadoWorkspace.evGx[280]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[281]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[282]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[283]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[70];
acadoVariables.x[75] += + acadoWorkspace.evGx[284]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[285]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[286]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[287]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[71];
acadoVariables.x[76] += + acadoWorkspace.evGx[288]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[289]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[290]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[291]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[72];
acadoVariables.x[77] += + acadoWorkspace.evGx[292]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[293]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[294]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[295]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[73];
acadoVariables.x[78] += + acadoWorkspace.evGx[296]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[297]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[298]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[299]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[74];
acadoVariables.x[79] += + acadoWorkspace.evGx[300]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[301]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[302]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[303]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[75];
acadoVariables.x[80] += + acadoWorkspace.evGx[304]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[305]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[306]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[307]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[76];
acadoVariables.x[81] += + acadoWorkspace.evGx[308]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[309]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[310]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[311]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[77];
acadoVariables.x[82] += + acadoWorkspace.evGx[312]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[313]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[314]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[315]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[78];
acadoVariables.x[83] += + acadoWorkspace.evGx[316]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[317]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[318]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[319]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[79];
acadoVariables.x[84] += + acadoWorkspace.evGx[320]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[321]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[322]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[323]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[80];
acadoVariables.x[85] += + acadoWorkspace.evGx[324]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[325]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[326]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[327]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[81];
acadoVariables.x[86] += + acadoWorkspace.evGx[328]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[329]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[330]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[331]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[82];
acadoVariables.x[87] += + acadoWorkspace.evGx[332]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[333]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[334]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[335]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[83];
acadoVariables.x[88] += + acadoWorkspace.evGx[336]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[337]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[338]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[339]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[84];
acadoVariables.x[89] += + acadoWorkspace.evGx[340]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[341]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[342]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[343]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[85];
acadoVariables.x[90] += + acadoWorkspace.evGx[344]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[345]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[346]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[347]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[86];
acadoVariables.x[91] += + acadoWorkspace.evGx[348]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[349]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[350]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[351]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[87];
acadoVariables.x[92] += + acadoWorkspace.evGx[352]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[353]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[354]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[355]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[88];
acadoVariables.x[93] += + acadoWorkspace.evGx[356]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[357]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[358]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[359]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[89];
acadoVariables.x[94] += + acadoWorkspace.evGx[360]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[361]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[362]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[363]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[90];
acadoVariables.x[95] += + acadoWorkspace.evGx[364]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[365]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[366]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[367]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[91];
acadoVariables.x[96] += + acadoWorkspace.evGx[368]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[369]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[370]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[371]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[92];
acadoVariables.x[97] += + acadoWorkspace.evGx[372]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[373]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[374]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[375]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[93];
acadoVariables.x[98] += + acadoWorkspace.evGx[376]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[377]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[378]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[379]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[94];
acadoVariables.x[99] += + acadoWorkspace.evGx[380]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[381]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[382]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[383]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[95];
acadoVariables.x[100] += + acadoWorkspace.evGx[384]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[385]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[386]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[387]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[96];
acadoVariables.x[101] += + acadoWorkspace.evGx[388]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[389]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[390]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[391]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[97];
acadoVariables.x[102] += + acadoWorkspace.evGx[392]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[393]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[394]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[395]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[98];
acadoVariables.x[103] += + acadoWorkspace.evGx[396]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[397]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[398]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[399]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[99];
acadoVariables.x[104] += + acadoWorkspace.evGx[400]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[401]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[402]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[403]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[100];
acadoVariables.x[105] += + acadoWorkspace.evGx[404]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[405]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[406]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[407]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[101];
acadoVariables.x[106] += + acadoWorkspace.evGx[408]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[409]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[410]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[411]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[102];
acadoVariables.x[107] += + acadoWorkspace.evGx[412]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[413]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[414]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[415]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[103];
acadoVariables.x[108] += + acadoWorkspace.evGx[416]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[417]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[418]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[419]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[104];
acadoVariables.x[109] += + acadoWorkspace.evGx[420]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[421]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[422]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[423]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[105];
acadoVariables.x[110] += + acadoWorkspace.evGx[424]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[425]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[426]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[427]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[106];
acadoVariables.x[111] += + acadoWorkspace.evGx[428]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[429]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[430]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[431]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[107];
acadoVariables.x[112] += + acadoWorkspace.evGx[432]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[433]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[434]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[435]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[108];
acadoVariables.x[113] += + acadoWorkspace.evGx[436]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[437]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[438]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[439]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[109];
acadoVariables.x[114] += + acadoWorkspace.evGx[440]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[441]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[442]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[443]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[110];
acadoVariables.x[115] += + acadoWorkspace.evGx[444]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[445]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[446]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[447]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[111];
acadoVariables.x[116] += + acadoWorkspace.evGx[448]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[449]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[450]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[451]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[112];
acadoVariables.x[117] += + acadoWorkspace.evGx[452]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[453]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[454]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[455]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[113];
acadoVariables.x[118] += + acadoWorkspace.evGx[456]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[457]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[458]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[459]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[114];
acadoVariables.x[119] += + acadoWorkspace.evGx[460]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[461]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[462]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[463]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[115];
acadoVariables.x[120] += + acadoWorkspace.evGx[464]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[465]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[466]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[467]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[116];
acadoVariables.x[121] += + acadoWorkspace.evGx[468]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[469]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[470]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[471]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[117];
acadoVariables.x[122] += + acadoWorkspace.evGx[472]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[473]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[474]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[475]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[118];
acadoVariables.x[123] += + acadoWorkspace.evGx[476]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[477]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[478]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[479]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[119];
acadoVariables.x[124] += + acadoWorkspace.evGx[480]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[481]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[482]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[483]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[120];
acadoVariables.x[125] += + acadoWorkspace.evGx[484]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[485]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[486]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[487]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[121];
acadoVariables.x[126] += + acadoWorkspace.evGx[488]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[489]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[490]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[491]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[122];
acadoVariables.x[127] += + acadoWorkspace.evGx[492]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[493]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[494]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[495]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[123];
acadoVariables.x[128] += + acadoWorkspace.evGx[496]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[497]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[498]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[499]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[124];
acadoVariables.x[129] += + acadoWorkspace.evGx[500]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[501]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[502]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[503]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[125];
acadoVariables.x[130] += + acadoWorkspace.evGx[504]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[505]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[506]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[507]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[126];
acadoVariables.x[131] += + acadoWorkspace.evGx[508]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[509]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[510]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[511]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[127];
acadoVariables.x[132] += + acadoWorkspace.evGx[512]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[513]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[514]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[515]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[128];
acadoVariables.x[133] += + acadoWorkspace.evGx[516]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[517]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[518]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[519]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[129];
acadoVariables.x[134] += + acadoWorkspace.evGx[520]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[521]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[522]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[523]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[130];
acadoVariables.x[135] += + acadoWorkspace.evGx[524]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[525]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[526]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[527]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[131];
acadoVariables.x[136] += + acadoWorkspace.evGx[528]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[529]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[530]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[531]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[132];
acadoVariables.x[137] += + acadoWorkspace.evGx[532]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[533]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[534]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[535]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[133];
acadoVariables.x[138] += + acadoWorkspace.evGx[536]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[537]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[538]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[539]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[134];
acadoVariables.x[139] += + acadoWorkspace.evGx[540]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[541]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[542]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[543]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[135];
acadoVariables.x[140] += + acadoWorkspace.evGx[544]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[545]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[546]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[547]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[136];
acadoVariables.x[141] += + acadoWorkspace.evGx[548]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[549]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[550]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[551]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[137];
acadoVariables.x[142] += + acadoWorkspace.evGx[552]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[553]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[554]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[555]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[138];
acadoVariables.x[143] += + acadoWorkspace.evGx[556]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[557]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[558]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[559]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[139];
acadoVariables.x[144] += + acadoWorkspace.evGx[560]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[561]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[562]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[563]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[140];
acadoVariables.x[145] += + acadoWorkspace.evGx[564]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[565]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[566]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[567]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[141];
acadoVariables.x[146] += + acadoWorkspace.evGx[568]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[569]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[570]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[571]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[142];
acadoVariables.x[147] += + acadoWorkspace.evGx[572]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[573]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[574]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[575]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[143];
acadoVariables.x[148] += + acadoWorkspace.evGx[576]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[577]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[578]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[579]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[144];
acadoVariables.x[149] += + acadoWorkspace.evGx[580]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[581]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[582]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[583]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[145];
acadoVariables.x[150] += + acadoWorkspace.evGx[584]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[585]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[586]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[587]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[146];
acadoVariables.x[151] += + acadoWorkspace.evGx[588]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[589]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[590]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[591]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[147];
acadoVariables.x[152] += + acadoWorkspace.evGx[592]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[593]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[594]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[595]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[148];
acadoVariables.x[153] += + acadoWorkspace.evGx[596]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[597]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[598]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[599]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[149];
acadoVariables.x[154] += + acadoWorkspace.evGx[600]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[601]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[602]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[603]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[150];
acadoVariables.x[155] += + acadoWorkspace.evGx[604]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[605]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[606]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[607]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[151];
acadoVariables.x[156] += + acadoWorkspace.evGx[608]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[609]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[610]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[611]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[152];
acadoVariables.x[157] += + acadoWorkspace.evGx[612]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[613]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[614]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[615]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[153];
acadoVariables.x[158] += + acadoWorkspace.evGx[616]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[617]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[618]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[619]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[154];
acadoVariables.x[159] += + acadoWorkspace.evGx[620]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[621]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[622]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[623]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[155];
acadoVariables.x[160] += + acadoWorkspace.evGx[624]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[625]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[626]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[627]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[156];
acadoVariables.x[161] += + acadoWorkspace.evGx[628]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[629]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[630]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[631]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[157];
acadoVariables.x[162] += + acadoWorkspace.evGx[632]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[633]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[634]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[635]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[158];
acadoVariables.x[163] += + acadoWorkspace.evGx[636]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[637]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[638]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[639]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[159];
acadoVariables.x[164] += + acadoWorkspace.evGx[640]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[641]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[642]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[643]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[160];
acadoVariables.x[165] += + acadoWorkspace.evGx[644]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[645]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[646]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[647]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[161];
acadoVariables.x[166] += + acadoWorkspace.evGx[648]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[649]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[650]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[651]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[162];
acadoVariables.x[167] += + acadoWorkspace.evGx[652]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[653]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[654]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[655]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[163];
acadoVariables.x[168] += + acadoWorkspace.evGx[656]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[657]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[658]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[659]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[164];
acadoVariables.x[169] += + acadoWorkspace.evGx[660]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[661]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[662]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[663]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[165];
acadoVariables.x[170] += + acadoWorkspace.evGx[664]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[665]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[666]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[667]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[166];
acadoVariables.x[171] += + acadoWorkspace.evGx[668]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[669]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[670]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[671]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[167];
acadoVariables.x[172] += + acadoWorkspace.evGx[672]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[673]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[674]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[675]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[168];
acadoVariables.x[173] += + acadoWorkspace.evGx[676]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[677]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[678]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[679]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[169];
acadoVariables.x[174] += + acadoWorkspace.evGx[680]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[681]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[682]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[683]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[170];
acadoVariables.x[175] += + acadoWorkspace.evGx[684]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[685]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[686]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[687]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[171];
acadoVariables.x[176] += + acadoWorkspace.evGx[688]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[689]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[690]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[691]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[172];
acadoVariables.x[177] += + acadoWorkspace.evGx[692]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[693]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[694]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[695]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[173];
acadoVariables.x[178] += + acadoWorkspace.evGx[696]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[697]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[698]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[699]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[174];
acadoVariables.x[179] += + acadoWorkspace.evGx[700]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[701]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[702]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[703]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[175];
acadoVariables.x[180] += + acadoWorkspace.evGx[704]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[705]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[706]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[707]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[176];
acadoVariables.x[181] += + acadoWorkspace.evGx[708]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[709]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[710]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[711]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[177];
acadoVariables.x[182] += + acadoWorkspace.evGx[712]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[713]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[714]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[715]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[178];
acadoVariables.x[183] += + acadoWorkspace.evGx[716]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[717]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[718]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[719]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[179];
acadoVariables.x[184] += + acadoWorkspace.evGx[720]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[721]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[722]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[723]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[180];
acadoVariables.x[185] += + acadoWorkspace.evGx[724]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[725]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[726]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[727]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[181];
acadoVariables.x[186] += + acadoWorkspace.evGx[728]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[729]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[730]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[731]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[182];
acadoVariables.x[187] += + acadoWorkspace.evGx[732]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[733]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[734]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[735]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[183];
acadoVariables.x[188] += + acadoWorkspace.evGx[736]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[737]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[738]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[739]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[184];
acadoVariables.x[189] += + acadoWorkspace.evGx[740]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[741]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[742]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[743]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[185];
acadoVariables.x[190] += + acadoWorkspace.evGx[744]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[745]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[746]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[747]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[186];
acadoVariables.x[191] += + acadoWorkspace.evGx[748]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[749]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[750]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[751]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[187];
acadoVariables.x[192] += + acadoWorkspace.evGx[752]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[753]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[754]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[755]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[188];
acadoVariables.x[193] += + acadoWorkspace.evGx[756]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[757]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[758]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[759]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[189];
acadoVariables.x[194] += + acadoWorkspace.evGx[760]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[761]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[762]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[763]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[190];
acadoVariables.x[195] += + acadoWorkspace.evGx[764]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[765]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[766]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[767]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[191];
acadoVariables.x[196] += + acadoWorkspace.evGx[768]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[769]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[770]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[771]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[192];
acadoVariables.x[197] += + acadoWorkspace.evGx[772]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[773]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[774]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[775]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[193];
acadoVariables.x[198] += + acadoWorkspace.evGx[776]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[777]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[778]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[779]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[194];
acadoVariables.x[199] += + acadoWorkspace.evGx[780]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[781]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[782]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[783]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[195];
acadoVariables.x[200] += + acadoWorkspace.evGx[784]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[785]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[786]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[787]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[196];
acadoVariables.x[201] += + acadoWorkspace.evGx[788]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[789]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[790]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[791]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[197];
acadoVariables.x[202] += + acadoWorkspace.evGx[792]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[793]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[794]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[795]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[198];
acadoVariables.x[203] += + acadoWorkspace.evGx[796]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[797]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[798]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[799]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[199];
acadoVariables.x[204] += + acadoWorkspace.evGx[800]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[801]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[802]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[803]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[200];
acadoVariables.x[205] += + acadoWorkspace.evGx[804]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[805]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[806]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[807]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[201];
acadoVariables.x[206] += + acadoWorkspace.evGx[808]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[809]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[810]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[811]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[202];
acadoVariables.x[207] += + acadoWorkspace.evGx[812]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[813]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[814]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[815]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[203];
acadoVariables.x[208] += + acadoWorkspace.evGx[816]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[817]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[818]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[819]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[204];
acadoVariables.x[209] += + acadoWorkspace.evGx[820]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[821]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[822]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[823]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[205];
acadoVariables.x[210] += + acadoWorkspace.evGx[824]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[825]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[826]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[827]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[206];
acadoVariables.x[211] += + acadoWorkspace.evGx[828]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[829]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[830]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[831]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[207];
acadoVariables.x[212] += + acadoWorkspace.evGx[832]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[833]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[834]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[835]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[208];
acadoVariables.x[213] += + acadoWorkspace.evGx[836]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[837]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[838]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[839]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[209];
acadoVariables.x[214] += + acadoWorkspace.evGx[840]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[841]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[842]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[843]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[210];
acadoVariables.x[215] += + acadoWorkspace.evGx[844]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[845]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[846]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[847]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[211];
acadoVariables.x[216] += + acadoWorkspace.evGx[848]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[849]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[850]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[851]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[212];
acadoVariables.x[217] += + acadoWorkspace.evGx[852]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[853]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[854]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[855]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[213];
acadoVariables.x[218] += + acadoWorkspace.evGx[856]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[857]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[858]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[859]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[214];
acadoVariables.x[219] += + acadoWorkspace.evGx[860]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[861]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[862]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[863]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[215];
acadoVariables.x[220] += + acadoWorkspace.evGx[864]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[865]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[866]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[867]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[216];
acadoVariables.x[221] += + acadoWorkspace.evGx[868]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[869]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[870]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[871]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[217];
acadoVariables.x[222] += + acadoWorkspace.evGx[872]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[873]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[874]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[875]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[218];
acadoVariables.x[223] += + acadoWorkspace.evGx[876]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[877]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[878]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[879]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[219];
acadoVariables.x[224] += + acadoWorkspace.evGx[880]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[881]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[882]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[883]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[220];
acadoVariables.x[225] += + acadoWorkspace.evGx[884]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[885]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[886]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[887]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[221];
acadoVariables.x[226] += + acadoWorkspace.evGx[888]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[889]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[890]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[891]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[222];
acadoVariables.x[227] += + acadoWorkspace.evGx[892]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[893]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[894]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[895]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[223];
acadoVariables.x[228] += + acadoWorkspace.evGx[896]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[897]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[898]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[899]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[224];
acadoVariables.x[229] += + acadoWorkspace.evGx[900]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[901]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[902]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[903]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[225];
acadoVariables.x[230] += + acadoWorkspace.evGx[904]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[905]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[906]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[907]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[226];
acadoVariables.x[231] += + acadoWorkspace.evGx[908]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[909]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[910]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[911]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[227];
acadoVariables.x[232] += + acadoWorkspace.evGx[912]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[913]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[914]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[915]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[228];
acadoVariables.x[233] += + acadoWorkspace.evGx[916]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[917]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[918]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[919]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[229];
acadoVariables.x[234] += + acadoWorkspace.evGx[920]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[921]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[922]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[923]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[230];
acadoVariables.x[235] += + acadoWorkspace.evGx[924]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[925]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[926]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[927]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[231];
acadoVariables.x[236] += + acadoWorkspace.evGx[928]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[929]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[930]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[931]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[232];
acadoVariables.x[237] += + acadoWorkspace.evGx[932]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[933]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[934]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[935]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[233];
acadoVariables.x[238] += + acadoWorkspace.evGx[936]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[937]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[938]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[939]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[234];
acadoVariables.x[239] += + acadoWorkspace.evGx[940]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[941]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[942]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[943]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[235];
acadoVariables.x[240] += + acadoWorkspace.evGx[944]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[945]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[946]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[947]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[236];
acadoVariables.x[241] += + acadoWorkspace.evGx[948]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[949]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[950]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[951]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[237];
acadoVariables.x[242] += + acadoWorkspace.evGx[952]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[953]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[954]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[955]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[238];
acadoVariables.x[243] += + acadoWorkspace.evGx[956]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[957]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[958]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[959]*acadoWorkspace.Dx0[3] + acadoWorkspace.d[239];

for (lRun1 = 0; lRun1 < 60; ++lRun1)
{
for (lRun2 = 0; lRun2 < lRun1 + 1; ++lRun2)
{
lRun3 = (((lRun1 + 1) * (lRun1)) / (2)) + (lRun2);
acado_multEDu( &(acadoWorkspace.E[ lRun3 * 4 ]), &(acadoWorkspace.x[ lRun2 ]), &(acadoVariables.x[ lRun1 * 4 + 4 ]) );
}
}
}

int acado_preparationStep(  )
{
int ret;

ret = acado_modelSimulation();
acado_evaluateObjective(  );
acado_condensePrep(  );
return ret;
}

int acado_feedbackStep(  )
{
int tmp;

acado_condenseFdb(  );

tmp = acado_solve( );

acado_expand(  );
return tmp;
}

int acado_initializeSolver(  )
{
int ret;

/* This is a function which must be called once before any other function call! */


ret = 0;

memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
return ret;
}

void acado_initializeNodesByForwardSimulation(  )
{
int index;
for (index = 0; index < 60; ++index)
{
acadoWorkspace.state[0] = acadoVariables.x[index * 4];
acadoWorkspace.state[1] = acadoVariables.x[index * 4 + 1];
acadoWorkspace.state[2] = acadoVariables.x[index * 4 + 2];
acadoWorkspace.state[3] = acadoVariables.x[index * 4 + 3];
acadoWorkspace.state[24] = acadoVariables.u[index];

acado_integrate(acadoWorkspace.state, index == 0);

acadoVariables.x[index * 4 + 4] = acadoWorkspace.state[0];
acadoVariables.x[index * 4 + 5] = acadoWorkspace.state[1];
acadoVariables.x[index * 4 + 6] = acadoWorkspace.state[2];
acadoVariables.x[index * 4 + 7] = acadoWorkspace.state[3];
}
}

void acado_shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd )
{
int index;
for (index = 0; index < 60; ++index)
{
acadoVariables.x[index * 4] = acadoVariables.x[index * 4 + 4];
acadoVariables.x[index * 4 + 1] = acadoVariables.x[index * 4 + 5];
acadoVariables.x[index * 4 + 2] = acadoVariables.x[index * 4 + 6];
acadoVariables.x[index * 4 + 3] = acadoVariables.x[index * 4 + 7];
}

if (strategy == 1 && xEnd != 0)
{
acadoVariables.x[240] = xEnd[0];
acadoVariables.x[241] = xEnd[1];
acadoVariables.x[242] = xEnd[2];
acadoVariables.x[243] = xEnd[3];
}
else if (strategy == 2) 
{
acadoWorkspace.state[0] = acadoVariables.x[240];
acadoWorkspace.state[1] = acadoVariables.x[241];
acadoWorkspace.state[2] = acadoVariables.x[242];
acadoWorkspace.state[3] = acadoVariables.x[243];
if (uEnd != 0)
{
acadoWorkspace.state[24] = uEnd[0];
}
else
{
acadoWorkspace.state[24] = acadoVariables.u[59];
}

acado_integrate(acadoWorkspace.state, 1);

acadoVariables.x[240] = acadoWorkspace.state[0];
acadoVariables.x[241] = acadoWorkspace.state[1];
acadoVariables.x[242] = acadoWorkspace.state[2];
acadoVariables.x[243] = acadoWorkspace.state[3];
}
}

void acado_shiftControls( real_t* const uEnd )
{
int index;
for (index = 0; index < 59; ++index)
{
acadoVariables.u[index] = acadoVariables.u[index + 1];
}

if (uEnd != 0)
{
acadoVariables.u[59] = uEnd[0];
}
}

real_t acado_getKKT(  )
{
real_t kkt;

int index;
real_t prd;

kkt = + acadoWorkspace.g[0]*acadoWorkspace.x[0] + acadoWorkspace.g[1]*acadoWorkspace.x[1] + acadoWorkspace.g[2]*acadoWorkspace.x[2] + acadoWorkspace.g[3]*acadoWorkspace.x[3] + acadoWorkspace.g[4]*acadoWorkspace.x[4] + acadoWorkspace.g[5]*acadoWorkspace.x[5] + acadoWorkspace.g[6]*acadoWorkspace.x[6] + acadoWorkspace.g[7]*acadoWorkspace.x[7] + acadoWorkspace.g[8]*acadoWorkspace.x[8] + acadoWorkspace.g[9]*acadoWorkspace.x[9] + acadoWorkspace.g[10]*acadoWorkspace.x[10] + acadoWorkspace.g[11]*acadoWorkspace.x[11] + acadoWorkspace.g[12]*acadoWorkspace.x[12] + acadoWorkspace.g[13]*acadoWorkspace.x[13] + acadoWorkspace.g[14]*acadoWorkspace.x[14] + acadoWorkspace.g[15]*acadoWorkspace.x[15] + acadoWorkspace.g[16]*acadoWorkspace.x[16] + acadoWorkspace.g[17]*acadoWorkspace.x[17] + acadoWorkspace.g[18]*acadoWorkspace.x[18] + acadoWorkspace.g[19]*acadoWorkspace.x[19] + acadoWorkspace.g[20]*acadoWorkspace.x[20] + acadoWorkspace.g[21]*acadoWorkspace.x[21] + acadoWorkspace.g[22]*acadoWorkspace.x[22] + acadoWorkspace.g[23]*acadoWorkspace.x[23] + acadoWorkspace.g[24]*acadoWorkspace.x[24] + acadoWorkspace.g[25]*acadoWorkspace.x[25] + acadoWorkspace.g[26]*acadoWorkspace.x[26] + acadoWorkspace.g[27]*acadoWorkspace.x[27] + acadoWorkspace.g[28]*acadoWorkspace.x[28] + acadoWorkspace.g[29]*acadoWorkspace.x[29] + acadoWorkspace.g[30]*acadoWorkspace.x[30] + acadoWorkspace.g[31]*acadoWorkspace.x[31] + acadoWorkspace.g[32]*acadoWorkspace.x[32] + acadoWorkspace.g[33]*acadoWorkspace.x[33] + acadoWorkspace.g[34]*acadoWorkspace.x[34] + acadoWorkspace.g[35]*acadoWorkspace.x[35] + acadoWorkspace.g[36]*acadoWorkspace.x[36] + acadoWorkspace.g[37]*acadoWorkspace.x[37] + acadoWorkspace.g[38]*acadoWorkspace.x[38] + acadoWorkspace.g[39]*acadoWorkspace.x[39] + acadoWorkspace.g[40]*acadoWorkspace.x[40] + acadoWorkspace.g[41]*acadoWorkspace.x[41] + acadoWorkspace.g[42]*acadoWorkspace.x[42] + acadoWorkspace.g[43]*acadoWorkspace.x[43] + acadoWorkspace.g[44]*acadoWorkspace.x[44] + acadoWorkspace.g[45]*acadoWorkspace.x[45] + acadoWorkspace.g[46]*acadoWorkspace.x[46] + acadoWorkspace.g[47]*acadoWorkspace.x[47] + acadoWorkspace.g[48]*acadoWorkspace.x[48] + acadoWorkspace.g[49]*acadoWorkspace.x[49] + acadoWorkspace.g[50]*acadoWorkspace.x[50] + acadoWorkspace.g[51]*acadoWorkspace.x[51] + acadoWorkspace.g[52]*acadoWorkspace.x[52] + acadoWorkspace.g[53]*acadoWorkspace.x[53] + acadoWorkspace.g[54]*acadoWorkspace.x[54] + acadoWorkspace.g[55]*acadoWorkspace.x[55] + acadoWorkspace.g[56]*acadoWorkspace.x[56] + acadoWorkspace.g[57]*acadoWorkspace.x[57] + acadoWorkspace.g[58]*acadoWorkspace.x[58] + acadoWorkspace.g[59]*acadoWorkspace.x[59];
kkt = fabs( kkt );
for (index = 0; index < 60; ++index)
{
prd = acadoWorkspace.y[index];
if (prd > 1e-12)
kkt += fabs(acadoWorkspace.lb[index] * prd);
else if (prd < -1e-12)
kkt += fabs(acadoWorkspace.ub[index] * prd);
}
for (index = 0; index < 60; ++index)
{
prd = acadoWorkspace.y[index + 60];
if (prd > 1e-12)
kkt += fabs(acadoWorkspace.lbA[index] * prd);
else if (prd < -1e-12)
kkt += fabs(acadoWorkspace.ubA[index] * prd);
}
return kkt;
}

real_t acado_getObjective(  )
{
real_t objVal;

int lRun1;
/** Row vector of size: 2 */
real_t tmpDy[ 2 ];

/** Column vector of size: 1 */
real_t tmpDyN[ 1 ];

for (lRun1 = 0; lRun1 < 60; ++lRun1)
{
acadoWorkspace.objValueIn[0] = acadoVariables.x[lRun1 * 4];
acadoWorkspace.objValueIn[1] = acadoVariables.x[lRun1 * 4 + 1];
acadoWorkspace.objValueIn[2] = acadoVariables.x[lRun1 * 4 + 2];
acadoWorkspace.objValueIn[3] = acadoVariables.x[lRun1 * 4 + 3];
acadoWorkspace.objValueIn[4] = acadoVariables.u[lRun1];

acado_evaluateLSQ( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.Dy[lRun1 * 2] = acadoWorkspace.objValueOut[0] - acadoVariables.y[lRun1 * 2];
acadoWorkspace.Dy[lRun1 * 2 + 1] = acadoWorkspace.objValueOut[1] - acadoVariables.y[lRun1 * 2 + 1];
}
acadoWorkspace.objValueIn[0] = acadoVariables.x[240];
acadoWorkspace.objValueIn[1] = acadoVariables.x[241];
acadoWorkspace.objValueIn[2] = acadoVariables.x[242];
acadoWorkspace.objValueIn[3] = acadoVariables.x[243];
acado_evaluateLSQEndTerm( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.DyN[0] = acadoWorkspace.objValueOut[0] - acadoVariables.yN[0];
objVal = 0.0000000000000000e+00;
for (lRun1 = 0; lRun1 < 60; ++lRun1)
{
tmpDy[0] = + acadoWorkspace.Dy[lRun1 * 2]*acadoVariables.W[0];
tmpDy[1] = + acadoWorkspace.Dy[lRun1 * 2 + 1]*acadoVariables.W[3];
objVal += + acadoWorkspace.Dy[lRun1 * 2]*tmpDy[0] + acadoWorkspace.Dy[lRun1 * 2 + 1]*tmpDy[1];
}

tmpDyN[0] = + acadoWorkspace.DyN[0]*acadoVariables.WN[0];
objVal += + acadoWorkspace.DyN[0]*tmpDyN[0];

objVal *= 0.5;
return objVal;
}

