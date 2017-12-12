#ifndef __c2_awe_mpc_sim_h__
#define __c2_awe_mpc_sim_h__

/* Type Definitions */
#ifndef struct_tag_suhmJzAUluz5Cp5Ic1cYFe
#define struct_tag_suhmJzAUluz5Cp5Ic1cYFe

struct tag_suhmJzAUluz5Cp5Ic1cYFe
{
  real_T status;
  real_T cpuTime;
  real_T kktValue;
  real_T objValue;
  real_T nIterations;
};

#endif                                 /*struct_tag_suhmJzAUluz5Cp5Ic1cYFe*/

#ifndef typedef_c2_suhmJzAUluz5Cp5Ic1cYFe
#define typedef_c2_suhmJzAUluz5Cp5Ic1cYFe

typedef struct tag_suhmJzAUluz5Cp5Ic1cYFe c2_suhmJzAUluz5Cp5Ic1cYFe;

#endif                                 /*typedef_c2_suhmJzAUluz5Cp5Ic1cYFe*/

#ifndef struct_tag_sXq4UNdTyVJzP7dB2SfY0TH
#define struct_tag_sXq4UNdTyVJzP7dB2SfY0TH

struct tag_sXq4UNdTyVJzP7dB2SfY0TH
{
  real_T x0[5];
  real_T x[155];
  real_T od[155];
  real_T u[30];
  real_T y[90];
  real_T yN[2];
  real_T W[9];
  real_T WN[4];
};

#endif                                 /*struct_tag_sXq4UNdTyVJzP7dB2SfY0TH*/

#ifndef typedef_c2_sXq4UNdTyVJzP7dB2SfY0TH
#define typedef_c2_sXq4UNdTyVJzP7dB2SfY0TH

typedef struct tag_sXq4UNdTyVJzP7dB2SfY0TH c2_sXq4UNdTyVJzP7dB2SfY0TH;

#endif                                 /*typedef_c2_sXq4UNdTyVJzP7dB2SfY0TH*/

#ifndef struct_tag_sYQ3DxPN71J245DEMZxJdLC
#define struct_tag_sYQ3DxPN71J245DEMZxJdLC

struct tag_sYQ3DxPN71J245DEMZxJdLC
{
  real_T x[155];
  real_T u[30];
  c2_suhmJzAUluz5Cp5Ic1cYFe info;
};

#endif                                 /*struct_tag_sYQ3DxPN71J245DEMZxJdLC*/

#ifndef typedef_c2_sYQ3DxPN71J245DEMZxJdLC
#define typedef_c2_sYQ3DxPN71J245DEMZxJdLC

typedef struct tag_sYQ3DxPN71J245DEMZxJdLC c2_sYQ3DxPN71J245DEMZxJdLC;

#endif                                 /*typedef_c2_sYQ3DxPN71J245DEMZxJdLC*/

#ifndef typedef_SFc2_awe_mpc_simInstanceStruct
#define typedef_SFc2_awe_mpc_simInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c2_sfEvent;
  boolean_T c2_doneDoubleBufferReInit;
  uint8_T c2_is_active_c2_awe_mpc_sim;
  void *c2_fEmlrtCtx;
  real_T (*c2_X0)[5];
  real_T (*c2_x_out)[155];
  real_T (*c2_x_in)[155];
  real_T (*c2_u_in)[30];
  real_T (*c2_u_out)[30];
  real_T *c2_roll;
  real_T (*c2_x_next)[5];
} SFc2_awe_mpc_simInstanceStruct;

#endif                                 /*typedef_SFc2_awe_mpc_simInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c2_awe_mpc_sim_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c2_awe_mpc_sim_get_check_sum(mxArray *plhs[]);
extern void c2_awe_mpc_sim_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
