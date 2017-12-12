/* Include files */

#include "awe_mpc_sim_sfun.h"
#include "c2_awe_mpc_sim.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "awe_mpc_sim_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c_with_debugger(S, sfGlobalDebugInstanceStruct);

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization);
static void chart_debug_initialize_data_addresses(SimStruct *S);
static const mxArray* sf_opaque_get_hover_data_for_msg(void *chartInstance,
  int32_T msgSSID);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c2_debug_family_names[20] = { "vw", "r", "circle_azimut",
  "circle_elevation", "R", "init", "N", "input", "Uref", "Yref", "output",
  "nargin", "nargout", "X0", "x_in", "u_in", "x_out", "u_out", "roll", "x_next"
};

/* Function Declarations */
static void initialize_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance);
static void initialize_params_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance);
static void enable_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct *chartInstance);
static void disable_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct *chartInstance);
static void c2_update_debugger_state_c2_awe_mpc_sim
  (SFc2_awe_mpc_simInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c2_awe_mpc_sim
  (SFc2_awe_mpc_simInstanceStruct *chartInstance);
static void set_sim_state_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_st);
static void finalize_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance);
static void sf_gateway_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance);
static void mdl_start_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance);
static void c2_chartstep_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance);
static void initSimStructsc2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber, uint32_T c2_instanceNumber);
static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData);
static void c2_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_b_x_next, const char_T *c2_identifier, real_T c2_y[5]);
static void c2_b_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[5]);
static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static real_T c2_c_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_b_roll, const char_T *c2_identifier);
static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_d_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_b_u_out, const char_T *c2_identifier, real_T c2_y[30]);
static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_e_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_b_x_out, const char_T *c2_identifier, real_T c2_y[155]);
static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_e_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_f_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_f_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[60]);
static void c2_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_g_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_g_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  c2_sXq4UNdTyVJzP7dB2SfY0TH *c2_y);
static void c2_h_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[90]);
static void c2_i_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[2]);
static void c2_j_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[9]);
static void c2_k_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[4]);
static void c2_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static void c2_repmat(SFc2_awe_mpc_simInstanceStruct *chartInstance, real_T
                      c2_b[155]);
static void c2_diag(SFc2_awe_mpc_simInstanceStruct *chartInstance, real_T c2_v[3],
                    real_T c2_d[9]);
static void c2_b_diag(SFc2_awe_mpc_simInstanceStruct *chartInstance, real_T
                      c2_v[2], real_T c2_d[4]);
static const mxArray *c2_emlrt_marshallOut(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const c2_sXq4UNdTyVJzP7dB2SfY0TH *c2_u);
static void c2_l_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_awe_MPCstep, const char_T *c2_identifier,
  c2_sYQ3DxPN71J245DEMZxJdLC *c2_y);
static void c2_m_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  c2_sYQ3DxPN71J245DEMZxJdLC *c2_y);
static void c2_n_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[155]);
static void c2_o_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[30]);
static c2_suhmJzAUluz5Cp5Ic1cYFe c2_p_emlrt_marshallIn
  (SFc2_awe_mpc_simInstanceStruct *chartInstance, const mxArray *c2_u, const
   emlrtMsgIdentifier *c2_parentId);
static real_T c2_q_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static const mxArray *c2_h_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static int32_T c2_r_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static uint8_T c2_s_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_b_is_active_c2_awe_mpc_sim, const char_T
  *c2_identifier);
static uint8_T c2_t_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void init_dsm_address_info(SFc2_awe_mpc_simInstanceStruct *chartInstance);
static void init_simulink_io_address(SFc2_awe_mpc_simInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance)
{
  if (sf_is_first_init_cond(chartInstance->S)) {
    initSimStructsc2_awe_mpc_sim(chartInstance);
    chart_debug_initialize_data_addresses(chartInstance->S);
  }

  chartInstance->c2_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c2_is_active_c2_awe_mpc_sim = 0U;
}

static void initialize_params_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c2_update_debugger_state_c2_awe_mpc_sim
  (SFc2_awe_mpc_simInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c2_awe_mpc_sim
  (SFc2_awe_mpc_simInstanceStruct *chartInstance)
{
  const mxArray *c2_st;
  const mxArray *c2_y = NULL;
  real_T c2_hoistedGlobal;
  const mxArray *c2_b_y = NULL;
  const mxArray *c2_c_y = NULL;
  const mxArray *c2_d_y = NULL;
  const mxArray *c2_e_y = NULL;
  uint8_T c2_b_hoistedGlobal;
  const mxArray *c2_f_y = NULL;
  c2_st = NULL;
  c2_st = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createcellmatrix(5, 1), false);
  c2_hoistedGlobal = *chartInstance->c2_roll;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_hoistedGlobal, 0, 0U, 0U, 0U, 0),
                false);
  sf_mex_setcell(c2_y, 0, c2_b_y);
  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", *chartInstance->c2_u_out, 0, 0U, 1U,
    0U, 1, 30), false);
  sf_mex_setcell(c2_y, 1, c2_c_y);
  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_create("y", *chartInstance->c2_x_next, 0, 0U, 1U,
    0U, 2, 1, 5), false);
  sf_mex_setcell(c2_y, 2, c2_d_y);
  c2_e_y = NULL;
  sf_mex_assign(&c2_e_y, sf_mex_create("y", *chartInstance->c2_x_out, 0, 0U, 1U,
    0U, 2, 31, 5), false);
  sf_mex_setcell(c2_y, 3, c2_e_y);
  c2_b_hoistedGlobal = chartInstance->c2_is_active_c2_awe_mpc_sim;
  c2_f_y = NULL;
  sf_mex_assign(&c2_f_y, sf_mex_create("y", &c2_b_hoistedGlobal, 3, 0U, 0U, 0U,
    0), false);
  sf_mex_setcell(c2_y, 4, c2_f_y);
  sf_mex_assign(&c2_st, c2_y, false);
  return c2_st;
}

static void set_sim_state_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_st)
{
  const mxArray *c2_u;
  real_T c2_dv0[30];
  int32_T c2_i0;
  real_T c2_dv1[5];
  int32_T c2_i1;
  real_T c2_dv2[155];
  int32_T c2_i2;
  chartInstance->c2_doneDoubleBufferReInit = true;
  c2_u = sf_mex_dup(c2_st);
  *chartInstance->c2_roll = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c2_u, 0)), "roll");
  c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 1)),
                        "u_out", c2_dv0);
  for (c2_i0 = 0; c2_i0 < 30; c2_i0++) {
    (*chartInstance->c2_u_out)[c2_i0] = c2_dv0[c2_i0];
  }

  c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 2)),
                      "x_next", c2_dv1);
  for (c2_i1 = 0; c2_i1 < 5; c2_i1++) {
    (*chartInstance->c2_x_next)[c2_i1] = c2_dv1[c2_i1];
  }

  c2_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 3)),
                        "x_out", c2_dv2);
  for (c2_i2 = 0; c2_i2 < 155; c2_i2++) {
    (*chartInstance->c2_x_out)[c2_i2] = c2_dv2[c2_i2];
  }

  chartInstance->c2_is_active_c2_awe_mpc_sim = c2_s_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 4)),
     "is_active_c2_awe_mpc_sim");
  sf_mex_destroy(&c2_u);
  c2_update_debugger_state_c2_awe_mpc_sim(chartInstance);
  sf_mex_destroy(&c2_st);
}

static void finalize_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance)
{
  int32_T c2_i3;
  int32_T c2_i4;
  int32_T c2_i5;
  int32_T c2_i6;
  int32_T c2_i7;
  int32_T c2_i8;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  for (c2_i3 = 0; c2_i3 < 30; c2_i3++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_u_in)[c2_i3], 2U);
  }

  for (c2_i4 = 0; c2_i4 < 155; c2_i4++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_x_in)[c2_i4], 1U);
  }

  for (c2_i5 = 0; c2_i5 < 5; c2_i5++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_X0)[c2_i5], 0U);
  }

  chartInstance->c2_sfEvent = CALL_EVENT;
  c2_chartstep_c2_awe_mpc_sim(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_awe_mpc_simMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c2_i6 = 0; c2_i6 < 155; c2_i6++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_x_out)[c2_i6], 3U);
  }

  for (c2_i7 = 0; c2_i7 < 30; c2_i7++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_u_out)[c2_i7], 4U);
  }

  _SFD_DATA_RANGE_CHECK(*chartInstance->c2_roll, 5U);
  for (c2_i8 = 0; c2_i8 < 5; c2_i8++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_x_next)[c2_i8], 6U);
  }
}

static void mdl_start_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance)
{
  sim_mode_is_external(chartInstance->S);
}

static void c2_chartstep_c2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance)
{
  int32_T c2_i9;
  int32_T c2_i10;
  real_T c2_b_X0[5];
  int32_T c2_i11;
  real_T c2_b_x_in[155];
  uint32_T c2_debug_family_var_map[20];
  real_T c2_b_u_in[30];
  real_T c2_vw;
  real_T c2_r;
  real_T c2_circle_azimut;
  real_T c2_circle_elevation;
  real_T c2_R;
  real_T c2_init;
  real_T c2_N;
  c2_sXq4UNdTyVJzP7dB2SfY0TH c2_input;
  real_T c2_Uref[30];
  real_T c2_Yref[60];
  c2_sYQ3DxPN71J245DEMZxJdLC c2_output;
  real_T c2_nargin = 3.0;
  real_T c2_nargout = 4.0;
  real_T c2_b_x_out[155];
  real_T c2_b_u_out[30];
  real_T c2_b_roll;
  real_T c2_b_x_next[5];
  int32_T c2_i12;
  int32_T c2_i13;
  int32_T c2_i14;
  int32_T c2_i15;
  int32_T c2_i16;
  int32_T c2_i17;
  int32_T c2_i18;
  real_T c2_c_x_in[150];
  int32_T c2_i19;
  int32_T c2_i20;
  int32_T c2_i21;
  real_T c2_c_u_in[29];
  int32_T c2_i22;
  int32_T c2_i23;
  real_T c2_dv3[155];
  int32_T c2_i24;
  int32_T c2_i25;
  int32_T c2_i26;
  int32_T c2_i27;
  static real_T c2_dv4[60] = { 0.388318718172466, 0.388318718172466,
    0.388318718172466, 0.388318718172466, 0.388318718172466, 0.388318718172466,
    0.388318718172466, 0.388318718172466, 0.388318718172466, 0.388318718172466,
    0.388318718172466, 0.388318718172466, 0.388318718172466, 0.388318718172466,
    0.388318718172466, 0.388318718172466, 0.388318718172466, 0.388318718172466,
    0.388318718172466, 0.388318718172466, 0.388318718172466, 0.388318718172466,
    0.388318718172466, 0.388318718172466, 0.388318718172466, 0.388318718172466,
    0.388318718172466, 0.388318718172466, 0.388318718172466, 0.388318718172466,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

  int32_T c2_i28;
  int32_T c2_i29;
  int32_T c2_i30;
  int32_T c2_i31;
  int32_T c2_i32;
  real_T c2_dv5[3];
  real_T c2_dv6[9];
  static real_T c2_dv7[3] = { 1.0, 0.0, 0.01 };

  int32_T c2_i33;
  int32_T c2_i34;
  real_T c2_dv8[2];
  real_T c2_dv9[4];
  int32_T c2_i35;
  int32_T c2_i36;
  int32_T c2_i37;
  c2_sYQ3DxPN71J245DEMZxJdLC c2_r0;
  int32_T c2_i38;
  int32_T c2_i39;
  int32_T c2_i40;
  int32_T c2_i41;
  int32_T c2_i42;
  int32_T c2_i43;
  int32_T c2_i44;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  for (c2_i9 = 0; c2_i9 < 5; c2_i9++) {
    c2_b_X0[c2_i9] = (*chartInstance->c2_X0)[c2_i9];
  }

  for (c2_i10 = 0; c2_i10 < 155; c2_i10++) {
    c2_b_x_in[c2_i10] = (*chartInstance->c2_x_in)[c2_i10];
  }

  for (c2_i11 = 0; c2_i11 < 30; c2_i11++) {
    c2_b_u_in[c2_i11] = (*chartInstance->c2_u_in)[c2_i11];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 20U, 20U, c2_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_vw, 0U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_r, 1U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_circle_azimut, 2U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_circle_elevation, 3U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_R, 4U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_init, 5U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_N, 6U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_input, 7U, c2_g_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_Uref, 8U, c2_c_sf_marshallOut,
    c2_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_Yref, 9U, c2_f_sf_marshallOut,
    c2_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_output, 10U, c2_e_sf_marshallOut,
    c2_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 11U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 12U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_b_X0, 13U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_b_x_in, 14U, c2_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_b_u_in, 15U, c2_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_b_x_out, 16U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_b_u_out, 17U, c2_c_sf_marshallOut,
    c2_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b_roll, 18U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_b_x_next, 19U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 5);
  c2_vw = -5.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 6);
  c2_r = 220.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 7);
  c2_circle_azimut = -1.5707963267948966;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 8);
  c2_circle_elevation = 0.3490658503988659;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 9);
  c2_R = 0.388318718172466;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 10);
  c2_init = 1.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 12);
  c2_N = 30.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 15);
  CV_EML_IF(0, 1, 0, false);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 18);
  CV_EML_IF(0, 1, 1, false);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 21);
  CV_EML_IF(0, 1, 2, false);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 29);
  c2_i12 = 0;
  c2_i13 = 0;
  for (c2_i14 = 0; c2_i14 < 5; c2_i14++) {
    for (c2_i16 = 0; c2_i16 < 30; c2_i16++) {
      c2_c_x_in[c2_i16 + c2_i12] = c2_b_x_in[(c2_i16 + c2_i13) + 1];
    }

    c2_i12 += 30;
    c2_i13 += 31;
  }

  c2_i15 = 0;
  c2_i17 = 0;
  for (c2_i18 = 0; c2_i18 < 5; c2_i18++) {
    for (c2_i19 = 0; c2_i19 < 30; c2_i19++) {
      c2_b_x_in[c2_i19 + c2_i15] = c2_c_x_in[c2_i19 + c2_i17];
    }

    c2_i15 += 31;
    c2_i17 += 30;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 30);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 31);
  for (c2_i20 = 0; c2_i20 < 29; c2_i20++) {
    c2_c_u_in[c2_i20] = c2_b_u_in[c2_i20 + 1];
  }

  for (c2_i21 = 0; c2_i21 < 29; c2_i21++) {
    c2_b_u_in[c2_i21] = c2_c_u_in[c2_i21];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 32);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 39);
  for (c2_i22 = 0; c2_i22 < 5; c2_i22++) {
    c2_input.x0[c2_i22] = c2_b_X0[c2_i22];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 40);
  for (c2_i23 = 0; c2_i23 < 155; c2_i23++) {
    c2_input.x[c2_i23] = c2_b_x_in[c2_i23];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 41);
  c2_repmat(chartInstance, c2_dv3);
  for (c2_i24 = 0; c2_i24 < 155; c2_i24++) {
    c2_input.od[c2_i24] = c2_dv3[c2_i24];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 43);
  for (c2_i25 = 0; c2_i25 < 30; c2_i25++) {
    c2_Uref[c2_i25] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 44);
  for (c2_i26 = 0; c2_i26 < 30; c2_i26++) {
    c2_input.u[c2_i26] = c2_b_u_in[c2_i26];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 46);
  for (c2_i27 = 0; c2_i27 < 60; c2_i27++) {
    c2_Yref[c2_i27] = c2_dv4[c2_i27];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 48);
  c2_i28 = 0;
  for (c2_i29 = 0; c2_i29 < 2; c2_i29++) {
    for (c2_i31 = 0; c2_i31 < 30; c2_i31++) {
      c2_input.y[c2_i31 + c2_i28] = c2_Yref[c2_i31 + c2_i28];
    }

    c2_i28 += 30;
  }

  for (c2_i30 = 0; c2_i30 < 30; c2_i30++) {
    c2_input.y[c2_i30 + 60] = c2_Uref[c2_i30];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 49);
  c2_input.yN[0] = c2_R;
  c2_input.yN[1] = 1.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 51);
  for (c2_i32 = 0; c2_i32 < 3; c2_i32++) {
    c2_dv5[c2_i32] = c2_dv7[c2_i32];
  }

  c2_diag(chartInstance, c2_dv5, c2_dv6);
  for (c2_i33 = 0; c2_i33 < 9; c2_i33++) {
    c2_input.W[c2_i33] = c2_dv6[c2_i33];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 52);
  for (c2_i34 = 0; c2_i34 < 2; c2_i34++) {
    c2_dv8[c2_i34] = 1.0 - (real_T)c2_i34;
  }

  c2_b_diag(chartInstance, c2_dv8, c2_dv9);
  for (c2_i35 = 0; c2_i35 < 4; c2_i35++) {
    c2_input.WN[c2_i35] = c2_dv9[c2_i35];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 55);
  for (c2_i36 = 0; c2_i36 < 155; c2_i36++) {
    c2_output.x[c2_i36] = c2_input.x[c2_i36];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 56);
  for (c2_i37 = 0; c2_i37 < 30; c2_i37++) {
    c2_output.u[c2_i37] = c2_input.u[c2_i37];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 57);
  c2_output.info.status = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 58);
  c2_output.info.cpuTime = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 59);
  c2_output.info.kktValue = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 60);
  c2_output.info.objValue = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 61);
  c2_output.info.nIterations = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 63);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 64);
  c2_l_emlrt_marshallIn(chartInstance, sf_mex_call_debug
                        (sfGlobalDebugInstanceStruct, "awe_MPCstep", 1U, 1U, 14,
    c2_emlrt_marshallOut(chartInstance, &c2_input)), "awe_MPCstep", &c2_r0);
  c2_output = c2_r0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 70);
  for (c2_i38 = 0; c2_i38 < 155; c2_i38++) {
    c2_b_x_out[c2_i38] = c2_output.x[c2_i38];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 71);
  for (c2_i39 = 0; c2_i39 < 30; c2_i39++) {
    c2_b_u_out[c2_i39] = c2_output.u[c2_i39];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 72);
  c2_b_roll = c2_output.x[93];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 73);
  c2_i40 = 0;
  for (c2_i41 = 0; c2_i41 < 5; c2_i41++) {
    c2_b_x_next[c2_i41] = c2_output.x[c2_i40 + 1];
    c2_i40 += 31;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, -73);
  _SFD_SYMBOL_SCOPE_POP();
  for (c2_i42 = 0; c2_i42 < 155; c2_i42++) {
    (*chartInstance->c2_x_out)[c2_i42] = c2_b_x_out[c2_i42];
  }

  for (c2_i43 = 0; c2_i43 < 30; c2_i43++) {
    (*chartInstance->c2_u_out)[c2_i43] = c2_b_u_out[c2_i43];
  }

  *chartInstance->c2_roll = c2_b_roll;
  for (c2_i44 = 0; c2_i44 < 5; c2_i44++) {
    (*chartInstance->c2_x_next)[c2_i44] = c2_b_x_next[c2_i44];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
}

static void initSimStructsc2_awe_mpc_sim(SFc2_awe_mpc_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber, uint32_T c2_instanceNumber)
{
  (void)(c2_machineNumber);
  (void)(c2_chartNumber);
  (void)(c2_instanceNumber);
}

static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData)
{
  const mxArray *c2_mxArrayOutData;
  int32_T c2_i45;
  const mxArray *c2_y = NULL;
  real_T c2_u[5];
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_mxArrayOutData = NULL;
  for (c2_i45 = 0; c2_i45 < 5; c2_i45++) {
    c2_u[c2_i45] = (*(real_T (*)[5])c2_inData)[c2_i45];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 1, 5), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_b_x_next, const char_T *c2_identifier, real_T c2_y[5])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = (const char *)c2_identifier;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_x_next), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_x_next);
}

static void c2_b_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[5])
{
  real_T c2_dv10[5];
  int32_T c2_i46;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv10, 1, 0, 0U, 1, 0U, 2, 1, 5);
  for (c2_i46 = 0; c2_i46 < 5; c2_i46++) {
    c2_y[c2_i46] = c2_dv10[c2_i46];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_x_next;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[5];
  int32_T c2_i47;
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_b_x_next = sf_mex_dup(c2_mxArrayInData);
  c2_thisId.fIdentifier = (const char *)c2_varName;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_x_next), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_x_next);
  for (c2_i47 = 0; c2_i47 < 5; c2_i47++) {
    (*(real_T (*)[5])c2_outData)[c2_i47] = c2_y[c2_i47];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData;
  real_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_mxArrayOutData = NULL;
  c2_u = *(real_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static real_T c2_c_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_b_roll, const char_T *c2_identifier)
{
  real_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = (const char *)c2_identifier;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_y = c2_q_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_roll), &c2_thisId);
  sf_mex_destroy(&c2_b_roll);
  return c2_y;
}

static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_roll;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y;
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_b_roll = sf_mex_dup(c2_mxArrayInData);
  c2_thisId.fIdentifier = (const char *)c2_varName;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_y = c2_q_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_roll), &c2_thisId);
  sf_mex_destroy(&c2_b_roll);
  *(real_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData;
  int32_T c2_i48;
  const mxArray *c2_y = NULL;
  real_T c2_u[30];
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_mxArrayOutData = NULL;
  for (c2_i48 = 0; c2_i48 < 30; c2_i48++) {
    c2_u[c2_i48] = (*(real_T (*)[30])c2_inData)[c2_i48];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 30), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_d_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_b_u_out, const char_T *c2_identifier, real_T c2_y[30])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = (const char *)c2_identifier;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_o_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_u_out), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_u_out);
}

static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_u_out;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[30];
  int32_T c2_i49;
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_b_u_out = sf_mex_dup(c2_mxArrayInData);
  c2_thisId.fIdentifier = (const char *)c2_varName;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_o_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_u_out), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_u_out);
  for (c2_i49 = 0; c2_i49 < 30; c2_i49++) {
    (*(real_T (*)[30])c2_outData)[c2_i49] = c2_y[c2_i49];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData;
  int32_T c2_i50;
  int32_T c2_i51;
  const mxArray *c2_y = NULL;
  int32_T c2_i52;
  real_T c2_u[155];
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_mxArrayOutData = NULL;
  c2_i50 = 0;
  for (c2_i51 = 0; c2_i51 < 5; c2_i51++) {
    for (c2_i52 = 0; c2_i52 < 31; c2_i52++) {
      c2_u[c2_i52 + c2_i50] = (*(real_T (*)[155])c2_inData)[c2_i52 + c2_i50];
    }

    c2_i50 += 31;
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 31, 5), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_e_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_b_x_out, const char_T *c2_identifier, real_T c2_y[155])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = (const char *)c2_identifier;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_n_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_x_out), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_x_out);
}

static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_x_out;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[155];
  int32_T c2_i53;
  int32_T c2_i54;
  int32_T c2_i55;
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_b_x_out = sf_mex_dup(c2_mxArrayInData);
  c2_thisId.fIdentifier = (const char *)c2_varName;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_n_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_x_out), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_b_x_out);
  c2_i53 = 0;
  for (c2_i54 = 0; c2_i54 < 5; c2_i54++) {
    for (c2_i55 = 0; c2_i55 < 31; c2_i55++) {
      (*(real_T (*)[155])c2_outData)[c2_i55 + c2_i53] = c2_y[c2_i55 + c2_i53];
    }

    c2_i53 += 31;
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_e_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData;
  c2_sYQ3DxPN71J245DEMZxJdLC c2_u;
  const mxArray *c2_y = NULL;
  int32_T c2_i56;
  const mxArray *c2_b_y = NULL;
  real_T c2_b_u[155];
  int32_T c2_i57;
  const mxArray *c2_c_y = NULL;
  real_T c2_c_u[30];
  const mxArray *c2_d_y = NULL;
  real_T c2_d_u;
  const mxArray *c2_e_y = NULL;
  real_T c2_e_u;
  const mxArray *c2_f_y = NULL;
  real_T c2_f_u;
  const mxArray *c2_g_y = NULL;
  real_T c2_g_u;
  const mxArray *c2_h_y = NULL;
  real_T c2_h_u;
  const mxArray *c2_i_y = NULL;
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_mxArrayOutData = NULL;
  c2_u = *(c2_sYQ3DxPN71J245DEMZxJdLC *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  for (c2_i56 = 0; c2_i56 < 155; c2_i56++) {
    c2_b_u[c2_i56] = c2_u.x[c2_i56];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_b_u, 0, 0U, 1U, 0U, 2, 31, 5),
                false);
  sf_mex_addfield(c2_y, c2_b_y, "x", "x", 0);
  for (c2_i57 = 0; c2_i57 < 30; c2_i57++) {
    c2_c_u[c2_i57] = c2_u.u[c2_i57];
  }

  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", c2_c_u, 0, 0U, 1U, 0U, 1, 30), false);
  sf_mex_addfield(c2_y, c2_c_y, "u", "u", 0);
  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  c2_d_u = c2_u.info.status;
  c2_e_y = NULL;
  sf_mex_assign(&c2_e_y, sf_mex_create("y", &c2_d_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_d_y, c2_e_y, "status", "status", 0);
  c2_e_u = c2_u.info.cpuTime;
  c2_f_y = NULL;
  sf_mex_assign(&c2_f_y, sf_mex_create("y", &c2_e_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_d_y, c2_f_y, "cpuTime", "cpuTime", 0);
  c2_f_u = c2_u.info.kktValue;
  c2_g_y = NULL;
  sf_mex_assign(&c2_g_y, sf_mex_create("y", &c2_f_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_d_y, c2_g_y, "kktValue", "kktValue", 0);
  c2_g_u = c2_u.info.objValue;
  c2_h_y = NULL;
  sf_mex_assign(&c2_h_y, sf_mex_create("y", &c2_g_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_d_y, c2_h_y, "objValue", "objValue", 0);
  c2_h_u = c2_u.info.nIterations;
  c2_i_y = NULL;
  sf_mex_assign(&c2_i_y, sf_mex_create("y", &c2_h_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_d_y, c2_i_y, "nIterations", "nIterations", 0);
  sf_mex_addfield(c2_y, c2_d_y, "info", "info", 0);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_awe_MPCstep;
  emlrtMsgIdentifier c2_thisId;
  c2_sYQ3DxPN71J245DEMZxJdLC c2_y;
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_awe_MPCstep = sf_mex_dup(c2_mxArrayInData);
  c2_thisId.fIdentifier = (const char *)c2_varName;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_m_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_awe_MPCstep), &c2_thisId,
                        &c2_y);
  sf_mex_destroy(&c2_awe_MPCstep);
  *(c2_sYQ3DxPN71J245DEMZxJdLC *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_f_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData;
  int32_T c2_i58;
  int32_T c2_i59;
  const mxArray *c2_y = NULL;
  int32_T c2_i60;
  real_T c2_u[60];
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_mxArrayOutData = NULL;
  c2_i58 = 0;
  for (c2_i59 = 0; c2_i59 < 2; c2_i59++) {
    for (c2_i60 = 0; c2_i60 < 30; c2_i60++) {
      c2_u[c2_i60 + c2_i58] = (*(real_T (*)[60])c2_inData)[c2_i60 + c2_i58];
    }

    c2_i58 += 30;
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 30, 2), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_f_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[60])
{
  real_T c2_dv11[60];
  int32_T c2_i61;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv11, 1, 0, 0U, 1, 0U, 2, 30,
                2);
  for (c2_i61 = 0; c2_i61 < 60; c2_i61++) {
    c2_y[c2_i61] = c2_dv11[c2_i61];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_Yref;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[60];
  int32_T c2_i62;
  int32_T c2_i63;
  int32_T c2_i64;
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_Yref = sf_mex_dup(c2_mxArrayInData);
  c2_thisId.fIdentifier = (const char *)c2_varName;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_Yref), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_Yref);
  c2_i62 = 0;
  for (c2_i63 = 0; c2_i63 < 2; c2_i63++) {
    for (c2_i64 = 0; c2_i64 < 30; c2_i64++) {
      (*(real_T (*)[60])c2_outData)[c2_i64 + c2_i62] = c2_y[c2_i64 + c2_i62];
    }

    c2_i62 += 30;
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_g_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData;
  c2_sXq4UNdTyVJzP7dB2SfY0TH c2_u;
  const mxArray *c2_y = NULL;
  int32_T c2_i65;
  const mxArray *c2_b_y = NULL;
  real_T c2_b_u[5];
  int32_T c2_i66;
  const mxArray *c2_c_y = NULL;
  real_T c2_c_u[155];
  int32_T c2_i67;
  const mxArray *c2_d_y = NULL;
  int32_T c2_i68;
  const mxArray *c2_e_y = NULL;
  real_T c2_d_u[30];
  int32_T c2_i69;
  const mxArray *c2_f_y = NULL;
  real_T c2_e_u[90];
  int32_T c2_i70;
  const mxArray *c2_g_y = NULL;
  real_T c2_f_u[2];
  int32_T c2_i71;
  const mxArray *c2_h_y = NULL;
  real_T c2_g_u[9];
  int32_T c2_i72;
  const mxArray *c2_i_y = NULL;
  real_T c2_h_u[4];
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_mxArrayOutData = NULL;
  c2_u = *(c2_sXq4UNdTyVJzP7dB2SfY0TH *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  for (c2_i65 = 0; c2_i65 < 5; c2_i65++) {
    c2_b_u[c2_i65] = c2_u.x0[c2_i65];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_b_u, 0, 0U, 1U, 0U, 2, 1, 5),
                false);
  sf_mex_addfield(c2_y, c2_b_y, "x0", "x0", 0);
  for (c2_i66 = 0; c2_i66 < 155; c2_i66++) {
    c2_c_u[c2_i66] = c2_u.x[c2_i66];
  }

  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", c2_c_u, 0, 0U, 1U, 0U, 2, 31, 5),
                false);
  sf_mex_addfield(c2_y, c2_c_y, "x", "x", 0);
  for (c2_i67 = 0; c2_i67 < 155; c2_i67++) {
    c2_c_u[c2_i67] = c2_u.od[c2_i67];
  }

  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_create("y", c2_c_u, 0, 0U, 1U, 0U, 2, 31, 5),
                false);
  sf_mex_addfield(c2_y, c2_d_y, "od", "od", 0);
  for (c2_i68 = 0; c2_i68 < 30; c2_i68++) {
    c2_d_u[c2_i68] = c2_u.u[c2_i68];
  }

  c2_e_y = NULL;
  sf_mex_assign(&c2_e_y, sf_mex_create("y", c2_d_u, 0, 0U, 1U, 0U, 1, 30), false);
  sf_mex_addfield(c2_y, c2_e_y, "u", "u", 0);
  for (c2_i69 = 0; c2_i69 < 90; c2_i69++) {
    c2_e_u[c2_i69] = c2_u.y[c2_i69];
  }

  c2_f_y = NULL;
  sf_mex_assign(&c2_f_y, sf_mex_create("y", c2_e_u, 0, 0U, 1U, 0U, 2, 30, 3),
                false);
  sf_mex_addfield(c2_y, c2_f_y, "y", "y", 0);
  for (c2_i70 = 0; c2_i70 < 2; c2_i70++) {
    c2_f_u[c2_i70] = c2_u.yN[c2_i70];
  }

  c2_g_y = NULL;
  sf_mex_assign(&c2_g_y, sf_mex_create("y", c2_f_u, 0, 0U, 1U, 0U, 2, 1, 2),
                false);
  sf_mex_addfield(c2_y, c2_g_y, "yN", "yN", 0);
  for (c2_i71 = 0; c2_i71 < 9; c2_i71++) {
    c2_g_u[c2_i71] = c2_u.W[c2_i71];
  }

  c2_h_y = NULL;
  sf_mex_assign(&c2_h_y, sf_mex_create("y", c2_g_u, 0, 0U, 1U, 0U, 2, 3, 3),
                false);
  sf_mex_addfield(c2_y, c2_h_y, "W", "W", 0);
  for (c2_i72 = 0; c2_i72 < 4; c2_i72++) {
    c2_h_u[c2_i72] = c2_u.WN[c2_i72];
  }

  c2_i_y = NULL;
  sf_mex_assign(&c2_i_y, sf_mex_create("y", c2_h_u, 0, 0U, 1U, 0U, 2, 2, 2),
                false);
  sf_mex_addfield(c2_y, c2_i_y, "WN", "WN", 0);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_g_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  c2_sXq4UNdTyVJzP7dB2SfY0TH *c2_y)
{
  emlrtMsgIdentifier c2_thisId;
  static const char * c2_fieldNames[8] = { "x0", "x", "od", "u", "y", "yN", "W",
    "WN" };

  c2_thisId.fParent = c2_parentId;
  c2_thisId.bParentIsCell = false;
  sf_mex_check_struct(c2_parentId, c2_u, 8, c2_fieldNames, 0U, NULL);
  c2_thisId.fIdentifier = "x0";
  c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_u, "x0",
    "x0", 0)), &c2_thisId, c2_y->x0);
  c2_thisId.fIdentifier = "x";
  c2_n_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_u, "x", "x",
    0)), &c2_thisId, c2_y->x);
  c2_thisId.fIdentifier = "od";
  c2_n_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_u, "od",
    "od", 0)), &c2_thisId, c2_y->od);
  c2_thisId.fIdentifier = "u";
  c2_o_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_u, "u", "u",
    0)), &c2_thisId, c2_y->u);
  c2_thisId.fIdentifier = "y";
  c2_h_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_u, "y", "y",
    0)), &c2_thisId, c2_y->y);
  c2_thisId.fIdentifier = "yN";
  c2_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_u, "yN",
    "yN", 0)), &c2_thisId, c2_y->yN);
  c2_thisId.fIdentifier = "W";
  c2_j_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_u, "W", "W",
    0)), &c2_thisId, c2_y->W);
  c2_thisId.fIdentifier = "WN";
  c2_k_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_u, "WN",
    "WN", 0)), &c2_thisId, c2_y->WN);
  sf_mex_destroy(&c2_u);
}

static void c2_h_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[90])
{
  real_T c2_dv12[90];
  int32_T c2_i73;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv12, 1, 0, 0U, 1, 0U, 2, 30,
                3);
  for (c2_i73 = 0; c2_i73 < 90; c2_i73++) {
    c2_y[c2_i73] = c2_dv12[c2_i73];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_i_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[2])
{
  real_T c2_dv13[2];
  int32_T c2_i74;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv13, 1, 0, 0U, 1, 0U, 2, 1, 2);
  for (c2_i74 = 0; c2_i74 < 2; c2_i74++) {
    c2_y[c2_i74] = c2_dv13[c2_i74];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_j_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[9])
{
  real_T c2_dv14[9];
  int32_T c2_i75;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv14, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c2_i75 = 0; c2_i75 < 9; c2_i75++) {
    c2_y[c2_i75] = c2_dv14[c2_i75];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_k_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[4])
{
  real_T c2_dv15[4];
  int32_T c2_i76;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv15, 1, 0, 0U, 1, 0U, 2, 2, 2);
  for (c2_i76 = 0; c2_i76 < 4; c2_i76++) {
    c2_y[c2_i76] = c2_dv15[c2_i76];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_input;
  emlrtMsgIdentifier c2_thisId;
  c2_sXq4UNdTyVJzP7dB2SfY0TH c2_y;
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_input = sf_mex_dup(c2_mxArrayInData);
  c2_thisId.fIdentifier = (const char *)c2_varName;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_input), &c2_thisId, &c2_y);
  sf_mex_destroy(&c2_input);
  *(c2_sXq4UNdTyVJzP7dB2SfY0TH *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

const mxArray *sf_c2_awe_mpc_sim_get_eml_resolved_functions_info(void)
{
  const mxArray *c2_nameCaptureInfo = NULL;
  c2_nameCaptureInfo = NULL;
  sf_mex_assign(&c2_nameCaptureInfo, sf_mex_create("nameCaptureInfo", NULL, 0,
    0U, 1U, 0U, 2, 0, 1), false);
  return c2_nameCaptureInfo;
}

static void c2_repmat(SFc2_awe_mpc_simInstanceStruct *chartInstance, real_T
                      c2_b[155])
{
  int32_T c2_jcol;
  int32_T c2_ibmat;
  int32_T c2_itilerow;
  static real_T c2_dv16[5] = { -5.0, 220.0, -1.5707963267948966,
    0.3490658503988659, 1.0 };

  (void)chartInstance;
  for (c2_jcol = 0; c2_jcol < 5; c2_jcol++) {
    c2_ibmat = c2_jcol * 31;
    for (c2_itilerow = 0; c2_itilerow < 31; c2_itilerow++) {
      c2_b[c2_ibmat + c2_itilerow] = c2_dv16[c2_jcol];
    }
  }
}

static void c2_diag(SFc2_awe_mpc_simInstanceStruct *chartInstance, real_T c2_v[3],
                    real_T c2_d[9])
{
  int32_T c2_i77;
  int32_T c2_j;
  int32_T c2_b_j;
  (void)chartInstance;
  for (c2_i77 = 0; c2_i77 < 9; c2_i77++) {
    c2_d[c2_i77] = 0.0;
  }

  c2_j = 0;
  for (c2_b_j = 0; c2_b_j < 3; c2_b_j++) {
    c2_d[c2_j] = c2_v[c2_b_j];
    c2_j += 4;
  }
}

static void c2_b_diag(SFc2_awe_mpc_simInstanceStruct *chartInstance, real_T
                      c2_v[2], real_T c2_d[4])
{
  int32_T c2_i78;
  int32_T c2_j;
  int32_T c2_b_j;
  (void)chartInstance;
  for (c2_i78 = 0; c2_i78 < 4; c2_i78++) {
    c2_d[c2_i78] = 0.0;
  }

  c2_j = 0;
  for (c2_b_j = 0; c2_b_j < 2; c2_b_j++) {
    c2_d[c2_j] = c2_v[c2_b_j];
    c2_j += 3;
  }
}

static const mxArray *c2_emlrt_marshallOut(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const c2_sXq4UNdTyVJzP7dB2SfY0TH *c2_u)
{
  const mxArray *c2_y;
  int32_T c2_i79;
  const mxArray *c2_b_y = NULL;
  real_T c2_b_u[5];
  int32_T c2_i80;
  const mxArray *c2_c_y = NULL;
  real_T c2_c_u[155];
  int32_T c2_i81;
  const mxArray *c2_d_y = NULL;
  int32_T c2_i82;
  const mxArray *c2_e_y = NULL;
  real_T c2_d_u[30];
  int32_T c2_i83;
  const mxArray *c2_f_y = NULL;
  real_T c2_e_u[90];
  int32_T c2_i84;
  const mxArray *c2_g_y = NULL;
  real_T c2_f_u[2];
  int32_T c2_i85;
  const mxArray *c2_h_y = NULL;
  real_T c2_g_u[9];
  int32_T c2_i86;
  const mxArray *c2_i_y = NULL;
  real_T c2_h_u[4];
  (void)chartInstance;
  c2_y = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  for (c2_i79 = 0; c2_i79 < 5; c2_i79++) {
    c2_b_u[c2_i79] = c2_u->x0[c2_i79];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_b_u, 0, 0U, 1U, 0U, 2, 1, 5),
                false);
  sf_mex_addfield(c2_y, c2_b_y, "x0", "x0", 0);
  for (c2_i80 = 0; c2_i80 < 155; c2_i80++) {
    c2_c_u[c2_i80] = c2_u->x[c2_i80];
  }

  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", c2_c_u, 0, 0U, 1U, 0U, 2, 31, 5),
                false);
  sf_mex_addfield(c2_y, c2_c_y, "x", "x", 0);
  for (c2_i81 = 0; c2_i81 < 155; c2_i81++) {
    c2_c_u[c2_i81] = c2_u->od[c2_i81];
  }

  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_create("y", c2_c_u, 0, 0U, 1U, 0U, 2, 31, 5),
                false);
  sf_mex_addfield(c2_y, c2_d_y, "od", "od", 0);
  for (c2_i82 = 0; c2_i82 < 30; c2_i82++) {
    c2_d_u[c2_i82] = c2_u->u[c2_i82];
  }

  c2_e_y = NULL;
  sf_mex_assign(&c2_e_y, sf_mex_create("y", c2_d_u, 0, 0U, 1U, 0U, 1, 30), false);
  sf_mex_addfield(c2_y, c2_e_y, "u", "u", 0);
  for (c2_i83 = 0; c2_i83 < 90; c2_i83++) {
    c2_e_u[c2_i83] = c2_u->y[c2_i83];
  }

  c2_f_y = NULL;
  sf_mex_assign(&c2_f_y, sf_mex_create("y", c2_e_u, 0, 0U, 1U, 0U, 2, 30, 3),
                false);
  sf_mex_addfield(c2_y, c2_f_y, "y", "y", 0);
  for (c2_i84 = 0; c2_i84 < 2; c2_i84++) {
    c2_f_u[c2_i84] = c2_u->yN[c2_i84];
  }

  c2_g_y = NULL;
  sf_mex_assign(&c2_g_y, sf_mex_create("y", c2_f_u, 0, 0U, 1U, 0U, 2, 1, 2),
                false);
  sf_mex_addfield(c2_y, c2_g_y, "yN", "yN", 0);
  for (c2_i85 = 0; c2_i85 < 9; c2_i85++) {
    c2_g_u[c2_i85] = c2_u->W[c2_i85];
  }

  c2_h_y = NULL;
  sf_mex_assign(&c2_h_y, sf_mex_create("y", c2_g_u, 0, 0U, 1U, 0U, 2, 3, 3),
                false);
  sf_mex_addfield(c2_y, c2_h_y, "W", "W", 0);
  for (c2_i86 = 0; c2_i86 < 4; c2_i86++) {
    c2_h_u[c2_i86] = c2_u->WN[c2_i86];
  }

  c2_i_y = NULL;
  sf_mex_assign(&c2_i_y, sf_mex_create("y", c2_h_u, 0, 0U, 1U, 0U, 2, 2, 2),
                false);
  sf_mex_addfield(c2_y, c2_i_y, "WN", "WN", 0);
  return c2_y;
}

static void c2_l_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_awe_MPCstep, const char_T *c2_identifier,
  c2_sYQ3DxPN71J245DEMZxJdLC *c2_y)
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = (const char *)c2_identifier;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_m_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_awe_MPCstep), &c2_thisId,
                        c2_y);
  sf_mex_destroy(&c2_awe_MPCstep);
}

static void c2_m_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  c2_sYQ3DxPN71J245DEMZxJdLC *c2_y)
{
  emlrtMsgIdentifier c2_thisId;
  static const char * c2_fieldNames[3] = { "x", "u", "info" };

  c2_thisId.fParent = c2_parentId;
  c2_thisId.bParentIsCell = false;
  sf_mex_check_struct(c2_parentId, c2_u, 3, c2_fieldNames, 0U, NULL);
  c2_thisId.fIdentifier = "x";
  c2_n_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_u, "x", "x",
    0)), &c2_thisId, c2_y->x);
  c2_thisId.fIdentifier = "u";
  c2_o_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield(c2_u, "u", "u",
    0)), &c2_thisId, c2_y->u);
  c2_thisId.fIdentifier = "info";
  c2_y->info = c2_p_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "info", "info", 0)), &c2_thisId);
  sf_mex_destroy(&c2_u);
}

static void c2_n_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[155])
{
  real_T c2_dv17[155];
  int32_T c2_i87;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv17, 1, 0, 0U, 1, 0U, 2, 31,
                5);
  for (c2_i87 = 0; c2_i87 < 155; c2_i87++) {
    c2_y[c2_i87] = c2_dv17[c2_i87];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_o_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[30])
{
  real_T c2_dv18[30];
  int32_T c2_i88;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv18, 1, 0, 0U, 1, 0U, 1, 30);
  for (c2_i88 = 0; c2_i88 < 30; c2_i88++) {
    c2_y[c2_i88] = c2_dv18[c2_i88];
  }

  sf_mex_destroy(&c2_u);
}

static c2_suhmJzAUluz5Cp5Ic1cYFe c2_p_emlrt_marshallIn
  (SFc2_awe_mpc_simInstanceStruct *chartInstance, const mxArray *c2_u, const
   emlrtMsgIdentifier *c2_parentId)
{
  c2_suhmJzAUluz5Cp5Ic1cYFe c2_y;
  emlrtMsgIdentifier c2_thisId;
  static const char * c2_fieldNames[5] = { "status", "cpuTime", "kktValue",
    "objValue", "nIterations" };

  c2_thisId.fParent = c2_parentId;
  c2_thisId.bParentIsCell = false;
  sf_mex_check_struct(c2_parentId, c2_u, 5, c2_fieldNames, 0U, NULL);
  c2_thisId.fIdentifier = "status";
  c2_y.status = c2_q_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "status", "status", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "cpuTime";
  c2_y.cpuTime = c2_q_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "cpuTime", "cpuTime", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "kktValue";
  c2_y.kktValue = c2_q_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c2_u, "kktValue", "kktValue", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "objValue";
  c2_y.objValue = c2_q_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c2_u, "objValue", "objValue", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "nIterations";
  c2_y.nIterations = c2_q_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getfield(c2_u, "nIterations", "nIterations", 0)), &c2_thisId);
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static real_T c2_q_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d0;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d0, 1, 0, 0U, 0, 0U, 0);
  c2_y = c2_d0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static const mxArray *c2_h_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData;
  int32_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_mxArrayOutData = NULL;
  c2_u = *(int32_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static int32_T c2_r_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  int32_T c2_y;
  int32_T c2_i89;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_i89, 1, 6, 0U, 0, 0U, 0);
  c2_y = c2_i89;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_sfEvent;
  emlrtMsgIdentifier c2_thisId;
  int32_T c2_y;
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)chartInstanceVoid;
  c2_b_sfEvent = sf_mex_dup(c2_mxArrayInData);
  c2_thisId.fIdentifier = (const char *)c2_varName;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_y = c2_r_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_sfEvent),
    &c2_thisId);
  sf_mex_destroy(&c2_b_sfEvent);
  *(int32_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static uint8_T c2_s_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_b_is_active_c2_awe_mpc_sim, const char_T
  *c2_identifier)
{
  uint8_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = (const char *)c2_identifier;
  c2_thisId.fParent = NULL;
  c2_thisId.bParentIsCell = false;
  c2_y = c2_t_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c2_b_is_active_c2_awe_mpc_sim), &c2_thisId);
  sf_mex_destroy(&c2_b_is_active_c2_awe_mpc_sim);
  return c2_y;
}

static uint8_T c2_t_emlrt_marshallIn(SFc2_awe_mpc_simInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  uint8_T c2_y;
  uint8_T c2_u0;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_u0, 1, 3, 0U, 0, 0U, 0);
  c2_y = c2_u0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void init_dsm_address_info(SFc2_awe_mpc_simInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc2_awe_mpc_simInstanceStruct
  *chartInstance)
{
  chartInstance->c2_fEmlrtCtx = (void *)sfrtGetEmlrtCtx(chartInstance->S);
  chartInstance->c2_X0 = (real_T (*)[5])ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c2_x_out = (real_T (*)[155])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c2_x_in = (real_T (*)[155])ssGetInputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c2_u_in = (real_T (*)[30])ssGetInputPortSignal_wrapper
    (chartInstance->S, 2);
  chartInstance->c2_u_out = (real_T (*)[30])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 2);
  chartInstance->c2_roll = (real_T *)ssGetOutputPortSignal_wrapper
    (chartInstance->S, 3);
  chartInstance->c2_x_next = (real_T (*)[5])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 4);
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c2_awe_mpc_sim_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2000624560U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1378670350U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1760992601U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1870544348U);
}

mxArray* sf_c2_awe_mpc_sim_get_post_codegen_info(void);
mxArray *sf_c2_awe_mpc_sim_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("90hjXyzypor5Q1oDxQnhNE");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(5);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt", "isFixedPointType" };

      mxArray *mxType = mxCreateStructMatrix(1,1,sizeof(typeFields)/sizeof
        (typeFields[0]),typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxType,0,"isFixedPointType",mxCreateDoubleScalar(0));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(31);
      pr[1] = (double)(5);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt", "isFixedPointType" };

      mxArray *mxType = mxCreateStructMatrix(1,1,sizeof(typeFields)/sizeof
        (typeFields[0]),typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxType,0,"isFixedPointType",mxCreateDoubleScalar(0));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,1,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(30);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt", "isFixedPointType" };

      mxArray *mxType = mxCreateStructMatrix(1,1,sizeof(typeFields)/sizeof
        (typeFields[0]),typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxType,0,"isFixedPointType",mxCreateDoubleScalar(0));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(31);
      pr[1] = (double)(5);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt", "isFixedPointType" };

      mxArray *mxType = mxCreateStructMatrix(1,1,sizeof(typeFields)/sizeof
        (typeFields[0]),typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxType,0,"isFixedPointType",mxCreateDoubleScalar(0));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,1,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(30);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt", "isFixedPointType" };

      mxArray *mxType = mxCreateStructMatrix(1,1,sizeof(typeFields)/sizeof
        (typeFields[0]),typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxType,0,"isFixedPointType",mxCreateDoubleScalar(0));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,0,mxREAL);
      double *pr = mxGetPr(mxSize);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt", "isFixedPointType" };

      mxArray *mxType = mxCreateStructMatrix(1,1,sizeof(typeFields)/sizeof
        (typeFields[0]),typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxType,0,"isFixedPointType",mxCreateDoubleScalar(0));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(5);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt", "isFixedPointType" };

      mxArray *mxType = mxCreateStructMatrix(1,1,sizeof(typeFields)/sizeof
        (typeFields[0]),typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxType,0,"isFixedPointType",mxCreateDoubleScalar(0));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    mxArray* mxPostCodegenInfo = sf_c2_awe_mpc_sim_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c2_awe_mpc_sim_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c2_awe_mpc_sim_jit_fallback_info(void)
{
  const char *infoFields[] = { "fallbackType", "fallbackReason",
    "hiddenFallbackType", "hiddenFallbackReason", "incompatibleSymbol" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 5, infoFields);
  mxArray *fallbackType = mxCreateString("pre");
  mxArray *fallbackReason = mxCreateString("hasBreakpoints");
  mxArray *hiddenFallbackType = mxCreateString("none");
  mxArray *hiddenFallbackReason = mxCreateString("");
  mxArray *incompatibleSymbol = mxCreateString("");
  mxSetField(mxInfo, 0, infoFields[0], fallbackType);
  mxSetField(mxInfo, 0, infoFields[1], fallbackReason);
  mxSetField(mxInfo, 0, infoFields[2], hiddenFallbackType);
  mxSetField(mxInfo, 0, infoFields[3], hiddenFallbackReason);
  mxSetField(mxInfo, 0, infoFields[4], incompatibleSymbol);
  return mxInfo;
}

mxArray *sf_c2_awe_mpc_sim_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c2_awe_mpc_sim_get_post_codegen_info(void)
{
  const char* fieldNames[] = { "exportedFunctionsUsedByThisChart",
    "exportedFunctionsChecksum" };

  mwSize dims[2] = { 1, 1 };

  mxArray* mxPostCodegenInfo = mxCreateStructArray(2, dims, sizeof(fieldNames)/
    sizeof(fieldNames[0]), fieldNames);

  {
    mxArray* mxExportedFunctionsChecksum = mxCreateString("");
    mwSize exp_dims[2] = { 0, 1 };

    mxArray* mxExportedFunctionsUsedByThisChart = mxCreateCellArray(2, exp_dims);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsUsedByThisChart",
               mxExportedFunctionsUsedByThisChart);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsChecksum",
               mxExportedFunctionsChecksum);
  }

  return mxPostCodegenInfo;
}

static const mxArray *sf_get_sim_state_info_c2_awe_mpc_sim(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x5'type','srcId','name','auxInfo'{{M[1],M[9],T\"roll\",},{M[1],M[8],T\"u_out\",},{M[1],M[10],T\"x_next\",},{M[1],M[5],T\"x_out\",},{M[8],M[0],T\"is_active_c2_awe_mpc_sim\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 5, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c2_awe_mpc_sim_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc2_awe_mpc_simInstanceStruct *chartInstance =
      (SFc2_awe_mpc_simInstanceStruct *)sf_get_chart_instance_ptr(S);
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _awe_mpc_simMachineNumber_,
           2,
           1,
           1,
           0,
           7,
           0,
           0,
           0,
           0,
           0,
           &chartInstance->chartNumber,
           &chartInstance->instanceNumber,
           (void *)S);

        /* Each instance must initialize its own list of scripts */
        init_script_number_translation(_awe_mpc_simMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_awe_mpc_simMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _awe_mpc_simMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"X0");
          _SFD_SET_DATA_PROPS(1,1,1,0,"x_in");
          _SFD_SET_DATA_PROPS(2,1,1,0,"u_in");
          _SFD_SET_DATA_PROPS(3,2,0,1,"x_out");
          _SFD_SET_DATA_PROPS(4,2,0,1,"u_out");
          _SFD_SET_DATA_PROPS(5,2,0,1,"roll");
          _SFD_SET_DATA_PROPS(6,2,0,1,"x_next");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,3,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,1469);
        _SFD_CV_INIT_EML_IF(0,1,0,314,328,-1,387);
        _SFD_CV_INIT_EML_IF(0,1,1,388,404,-1,437);
        _SFD_CV_INIT_EML_IF(0,1,2,438,454,-1,481);

        {
          unsigned int dimVector[2];
          dimVector[0]= 1U;
          dimVector[1]= 5U;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 31U;
          dimVector[1]= 5U;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 30U;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 31U;
          dimVector[1]= 5U;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_d_sf_marshallOut,(MexInFcnForType)
            c2_d_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 30U;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_c_sf_marshallOut,(MexInFcnForType)
            c2_c_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)c2_b_sf_marshallIn);

        {
          unsigned int dimVector[2];
          dimVector[0]= 1U;
          dimVector[1]= 5U;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)
            c2_sf_marshallIn);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _awe_mpc_simMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static void chart_debug_initialize_data_addresses(SimStruct *S)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc2_awe_mpc_simInstanceStruct *chartInstance =
      (SFc2_awe_mpc_simInstanceStruct *)sf_get_chart_instance_ptr(S);
    if (ssIsFirstInitCond(S)) {
      /* do this only if simulation is starting and after we know the addresses of all data */
      {
        _SFD_SET_DATA_VALUE_PTR(0U, (void *)chartInstance->c2_X0);
        _SFD_SET_DATA_VALUE_PTR(3U, (void *)chartInstance->c2_x_out);
        _SFD_SET_DATA_VALUE_PTR(1U, (void *)chartInstance->c2_x_in);
        _SFD_SET_DATA_VALUE_PTR(2U, (void *)chartInstance->c2_u_in);
        _SFD_SET_DATA_VALUE_PTR(4U, (void *)chartInstance->c2_u_out);
        _SFD_SET_DATA_VALUE_PTR(5U, (void *)chartInstance->c2_roll);
        _SFD_SET_DATA_VALUE_PTR(6U, (void *)chartInstance->c2_x_next);
      }
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "sscGCk0gflXHLAayenlGPeE";
}

static void sf_opaque_initialize_c2_awe_mpc_sim(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc2_awe_mpc_simInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c2_awe_mpc_sim((SFc2_awe_mpc_simInstanceStruct*)
    chartInstanceVar);
  initialize_c2_awe_mpc_sim((SFc2_awe_mpc_simInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c2_awe_mpc_sim(void *chartInstanceVar)
{
  enable_c2_awe_mpc_sim((SFc2_awe_mpc_simInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c2_awe_mpc_sim(void *chartInstanceVar)
{
  disable_c2_awe_mpc_sim((SFc2_awe_mpc_simInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c2_awe_mpc_sim(void *chartInstanceVar)
{
  sf_gateway_c2_awe_mpc_sim((SFc2_awe_mpc_simInstanceStruct*) chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c2_awe_mpc_sim(SimStruct* S)
{
  return get_sim_state_c2_awe_mpc_sim((SFc2_awe_mpc_simInstanceStruct *)
    sf_get_chart_instance_ptr(S));     /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c2_awe_mpc_sim(SimStruct* S, const mxArray
  *st)
{
  set_sim_state_c2_awe_mpc_sim((SFc2_awe_mpc_simInstanceStruct*)
    sf_get_chart_instance_ptr(S), st);
}

static void sf_opaque_terminate_c2_awe_mpc_sim(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc2_awe_mpc_simInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_awe_mpc_sim_optimization_info();
    }

    finalize_c2_awe_mpc_sim((SFc2_awe_mpc_simInstanceStruct*) chartInstanceVar);
    utFree(chartInstanceVar);
    if (ssGetUserData(S)!= NULL) {
      sf_free_ChartRunTimeInfo(S);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc2_awe_mpc_sim((SFc2_awe_mpc_simInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c2_awe_mpc_sim(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c2_awe_mpc_sim((SFc2_awe_mpc_simInstanceStruct*)
      sf_get_chart_instance_ptr(S));
  }
}

static void mdlSetWorkWidths_c2_awe_mpc_sim(SimStruct *S)
{
  /* Set overwritable ports for inplace optimization */
  ssSetInputPortDirectFeedThrough(S, 0, 1);
  ssSetInputPortDirectFeedThrough(S, 1, 1);
  ssSetInputPortDirectFeedThrough(S, 2, 1);
  ssSetStatesModifiedOnlyInUpdate(S, 1);
  ssSetBlockIsPurelyCombinatorial_wrapper(S, 1);
  ssMdlUpdateIsEmpty(S, 1);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_awe_mpc_sim_optimization_info(sim_mode_is_rtw_gen
      (S), sim_mode_is_modelref_sim(S), sim_mode_is_external(S));
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,2);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,1);
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,2,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_set_chart_accesses_machine_info(S, sf_get_instance_specialization(),
      infoStruct, 2);
    sf_update_buildInfo(S, sf_get_instance_specialization(),infoStruct,2);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,2,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,2,4);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=4; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 3; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,2);
    sf_register_codegen_names_for_scoped_functions_defined_by_chart(S);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(780678167U));
  ssSetChecksum1(S,(3528562378U));
  ssSetChecksum2(S,(2451483719U));
  ssSetChecksum3(S,(1987758974U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSetStateSemanticsClassicAndSynchronous(S, true);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c2_awe_mpc_sim(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c2_awe_mpc_sim(SimStruct *S)
{
  SFc2_awe_mpc_simInstanceStruct *chartInstance;
  chartInstance = (SFc2_awe_mpc_simInstanceStruct *)utMalloc(sizeof
    (SFc2_awe_mpc_simInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  memset(chartInstance, 0, sizeof(SFc2_awe_mpc_simInstanceStruct));
  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c2_awe_mpc_sim;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c2_awe_mpc_sim;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c2_awe_mpc_sim;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c2_awe_mpc_sim;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c2_awe_mpc_sim;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c2_awe_mpc_sim;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c2_awe_mpc_sim;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c2_awe_mpc_sim;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c2_awe_mpc_sim;
  chartInstance->chartInfo.mdlStart = mdlStart_c2_awe_mpc_sim;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c2_awe_mpc_sim;
  chartInstance->chartInfo.callGetHoverDataForMsg = NULL;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.callAtomicSubchartUserFcn = NULL;
  chartInstance->chartInfo.callAtomicSubchartAutoFcn = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  sf_init_ChartRunTimeInfo(S, &(chartInstance->chartInfo), false, 0);
  init_dsm_address_info(chartInstance);
  init_simulink_io_address(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  chart_debug_initialization(S,1);
  mdl_start_c2_awe_mpc_sim(chartInstance);
}

void c2_awe_mpc_sim_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c2_awe_mpc_sim(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c2_awe_mpc_sim(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c2_awe_mpc_sim(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c2_awe_mpc_sim_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
