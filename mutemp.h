#ifndef HEADER_FT
#define HEADER_FT

//double g_CenterOfMassZ;

double *g_Tmu_tab, *g_mu_tab;
int g_ntabmu;
int g_oof;

typedef struct temperature_params temperature_params;

struct temperature_params
{
  double rho;
  double prs;
  double mu;
};


void LoadMuTable();
double GetMuFromTable(double);
double fT(double, void*);

#endif //HEADER_FT
