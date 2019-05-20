#define  PHYSICS                 MHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              VECTOR
#define  COOLING                 TABULATED
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           CHARACTERISTIC_TRACING
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 2
#define  USER_DEF_PARAMETERS     23
#define  USER_DEF_CONSTANTS      9

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          YES
#define  MHD_FORMULATION         DIV_CLEANING
#define  BACKGROUND_FIELD        NO
#define  RESISTIVE_MHD           NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  GAMMA_GAS               0
#define  CLOUD_RADIUS            1
#define  CLOUD_CENTREX1          2
#define  CLOUD_CENTREX2          3
#define  CLOUD_CENTREX3          4
#define  RHO_W                   5
#define  RHO_C                   6
#define  PRS_C                   7
#define  VX1_W                   8
#define  VX1_C                   9
#define  BETA_W                  10
#define  BETA_C                  11
#define  DENSITY_STEEPNESS       12
#define  CORE_RADIUS             13
#define  CUT_RADIUS              14
#define  Z_W                     15
#define  Z_C                     16
#define  TEMP_C                  17
#define  TEMP_W                  18
#define  MU_C                    19
#define  GRAV_ACC                20
#define  VEL_THRESH              21
#define  N_0                     22

/* -- user-defined symbolic constants -- */

#define  UNIT_LENGTH             1000*CONST_pc
#define  FIELD_NONE              0
#define  FIELD_ORDERED           1
#define  FIELD_TANGLED           2
#define  FIELD_TRANSVERSE        1
#define  FIELD_PARALLEL          2
#define  FIELD_OBLIQUE           3
#define  CLOUD_FIELD             FIELD_TRANSVERSE
#define  BG_FIELD                FIELD_TRANSVERSE

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          NO
#define  PRINT_TO_FILE             YES
#define  INTERNAL_BOUNDARY         YES
#define  SHOCK_FLATTENING          MULTID
#define  CHOMBO_EN_SWITCH          YES
#define  CHOMBO_REF_VAR            RHO
#define  CHOMBO_LOGR               NO
#define  ARTIFICIAL_VISCOSITY      NO
#define  CHAR_LIMITING             YES
#define  LIMITER                   MINMOD_LIM
#define  ASSIGN_VECTOR_POTENTIAL   NO
