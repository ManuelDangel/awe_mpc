/*
*    This file is part of ACADO Toolkit.
*
*    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
*    Copyright (C) 2008-2009 by Boris Houska and Hans Joachim Ferreau, K.U.Leuven.
*    Developed within the Optimization in Engineering Center (OPTEC) under
*    supervision of Moritz Diehl. All rights reserved.
*
*    ACADO Toolkit is free software; you can redistribute it and/or
*    modify it under the terms of the GNU Lesser General Public
*    License as published by the Free Software Foundation; either
*    version 3 of the License, or (at your option) any later version.
*
*    ACADO Toolkit is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public
*    License along with ACADO Toolkit; if not, write to the Free Software
*    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*
*/


/**
*    Author David Ariens, Rien Quirynen
*    Date 2009-2013
*    http://www.acadotoolkit.org/matlab 
*/

#include <acado_optimal_control.hpp>
#include <acado_toolkit.hpp>
#include <acado/utils/matlab_acado_utils.hpp>

USING_NAMESPACE_ACADO

#include <mex.h>


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
 { 
 
    MatlabConsoleStreamBuf mybuf;
    RedirectStream redirect(std::cout, mybuf);
    clearAllStaticCounters( ); 
 
    mexPrintf("\nACADO Toolkit for Matlab - Developed by David Ariens and Rien Quirynen, 2009-2013 \n"); 
    mexPrintf("Support available at http://www.acadotoolkit.org/matlab \n \n"); 

    if (nrhs != 0){ 
      mexErrMsgTxt("This problem expects 0 right hand side argument(s) since you have defined 0 MexInput(s)");
    } 
 
    TIME autotime;
    DifferentialState psi;
    DifferentialState theta;
    DifferentialState gamma;
    DifferentialState phi;
    DifferentialState vt;
    Control dphi;
    OnlineData vw; 
    OnlineData r; 
    OnlineData circle_azimut; 
    OnlineData circle_elevation; 
    OnlineData init; 
    BMatrix acadodata_M1;
    acadodata_M1.read( "awe_mpc_data_acadodata_M1.txt" );
    Function acadodata_f2;
    acadodata_f2 << init*sqrt(((-circle_azimut+psi)*(-circle_azimut+psi)+(-circle_elevation+theta)*(-circle_elevation+theta)));
    acadodata_f2 << 1.97310749999999979210e-01*cos((-circle_azimut+psi))*cos(phi)*cos(theta)*sqrt(((-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)*(-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)+(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)*(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)+(-(-cos(theta))*sin(psi)*vw)*(-(-cos(theta))*sin(psi)*vw)))*sqrt(((-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)*(-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)+(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)*(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)+(-(-cos(theta))*sin(psi)*vw)*(-(-cos(theta))*sin(psi)*vw)));
    acadodata_f2 << dphi;
    BMatrix acadodata_M2;
    acadodata_M2.read( "awe_mpc_data_acadodata_M2.txt" );
    Function acadodata_f3;
    acadodata_f3 << init*sqrt(((-circle_azimut+psi)*(-circle_azimut+psi)+(-circle_elevation+theta)*(-circle_elevation+theta)));
    acadodata_f3 << 1.97310749999999979210e-01*cos((-circle_azimut+psi))*cos(phi)*cos(theta)*sqrt(((-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)*(-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)+(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)*(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)+(-(-cos(theta))*sin(psi)*vw)*(-(-cos(theta))*sin(psi)*vw)))*sqrt(((-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)*(-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)+(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)*(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)+(-(-cos(theta))*sin(psi)*vw)*(-(-cos(theta))*sin(psi)*vw)));
    OCP ocp1(0, 3, 30);
    ocp1.minimizeLSQ(acadodata_M1, acadodata_f2);
    ocp1.minimizeLSQEndTerm(acadodata_M2, acadodata_f3);
    ocp1.subjectTo((-1.22173047639603060688e+00) <= phi <= 1.22173047639603060688e+00);
    ocp1.subjectTo((-2.44346095279206121376e+00) <= dphi <= 2.44346095279206121376e+00);
    ocp1.subjectTo(1.74532925199432947805e-01 <= theta <= 1.22173047639603060688e+00);
    ocp1.subjectTo(1.00000000000000000000e+01 <= vt <= 2.00000000000000000000e+02);
    DifferentialEquation acadodata_f1;
    acadodata_f1 << dot(psi) == 1/cos(theta)/r*sin(gamma)*vt;
    acadodata_f1 << dot(theta) == cos(gamma)/r*vt;
    acadodata_f1 << dot(gamma) == (1.97310749999999979210e-01/2.64999999999999991118e+00*sin(phi)*sqrt(((-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)*(-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)+(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)*(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)+(-(-cos(theta))*sin(psi)*vw)*(-(-cos(theta))*sin(psi)*vw)))*sqrt(((-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)*(-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)+(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)*(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)+(-(-cos(theta))*sin(psi)*vw)*(-(-cos(theta))*sin(psi)*vw)))+9.81000000000000049738e+00*cos(theta)*sin(gamma))/vt;
    acadodata_f1 << dot(phi) == dphi;
    acadodata_f1 << dot(vt) == ((1.97310749999999979210e-01*cos(phi)*sin(asin((-(-cos(theta))*sin(psi)*vw)/sqrt(((-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)*(-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)+(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)*(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)+(-(-cos(theta))*sin(psi)*vw)*(-(-cos(theta))*sin(psi)*vw)))))-4.44785250000000048631e-02*cos(asin((-(-cos(theta))*sin(psi)*vw)/sqrt(((-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)*(-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)+(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)*(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)+(-(-cos(theta))*sin(psi)*vw)*(-(-cos(theta))*sin(psi)*vw))))))/2.64999999999999991118e+00*sqrt(((-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)*(-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)+(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)*(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)+(-(-cos(theta))*sin(psi)*vw)*(-(-cos(theta))*sin(psi)*vw)))*sqrt(((-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)*(-((-sin(gamma))*(-sin(theta))*sin(psi)+cos(gamma)*cos(psi))*vw)+(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)*(-((-sin(theta))*sin(psi)*cos(gamma)+cos(psi)*sin(gamma))*vw+vt)+(-(-cos(theta))*sin(psi)*vw)*(-(-cos(theta))*sin(psi)*vw)))-9.81000000000000049738e+00*cos(gamma)*cos(theta));

    ocp1.setModel( acadodata_f1 );


    ocp1.setNU( 1 );
    ocp1.setNP( 0 );
    ocp1.setNOD( 5 );
    OCPexport ExportModule1( ocp1 );
    ExportModule1.set( GENERATE_MATLAB_INTERFACE, 1 );
    uint options_flag;
    options_flag = ExportModule1.set( HESSIAN_APPROXIMATION, GAUSS_NEWTON );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: HESSIAN_APPROXIMATION");
    options_flag = ExportModule1.set( DISCRETIZATION_TYPE, MULTIPLE_SHOOTING );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: DISCRETIZATION_TYPE");
    options_flag = ExportModule1.set( INTEGRATOR_TYPE, INT_IRK_GL4 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: INTEGRATOR_TYPE");
    options_flag = ExportModule1.set( DYNAMIC_SENSITIVITY, FORWARD );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: DYNAMIC_SENSITIVITY");
    options_flag = ExportModule1.set( LINEAR_ALGEBRA_SOLVER, GAUSS_LU );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: LINEAR_ALGEBRA_SOLVER");
    options_flag = ExportModule1.set( NUM_INTEGRATOR_STEPS, 30 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: NUM_INTEGRATOR_STEPS");
    options_flag = ExportModule1.set( QP_SOLVER, QP_QPOASES );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: QP_SOLVER");
    options_flag = ExportModule1.set( HOTSTART_QP, NO );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: HOTSTART_QP");
    options_flag = ExportModule1.set( SPARSE_QP_SOLUTION, FULL_CONDENSING );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: SPARSE_QP_SOLUTION");
    uint export_flag;
    export_flag = ExportModule1.exportCode( "export_awe_mpc" );
    if(export_flag != 0) mexErrMsgTxt("ACADO export failed because of the above error(s)!");


    clearAllStaticCounters( ); 
 
} 

