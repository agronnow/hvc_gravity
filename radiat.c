#include "pluto.h"
#include "mutemp.h"

//#define WOLFIRE95_CIICOOL
//#define WOLFIRE95_UVHEATING
//#define ANALYTIC_HEATING

#define frac_Z   1.e-3   /* = N(Z) / N(H), fractional number density of 
                              metals (Z) with respect to hydrogen (H) */ 
#define frac_He  0.082   /* = N(Z) / N(H), fractional number density of 
                              helium (He) with respect to hydrogen (H) */ 

#define TRC_Z TRC+1

//double *g_Tmu_tab, *g_mu_tab;
//int g_ntabmu = 0;
//int g_oof = 0;

void LoadMuTable()
{
  FILE *fmu;
  double dummy1, dummy2;

  if (g_mu_tab == NULL) //Read mu table
  {
    print1 (" > Reading mu table from disk...\n");
    fmu = fopen("mutable.dat","r");
    if (fmu == NULL){
      print1 ("! mutemp: mutable.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    g_mu_tab = ARRAY_1D(20000, double);
    g_Tmu_tab = ARRAY_1D(20000, double);

    g_ntabmu = 0;
    while (fscanf(fmu, "%lf %lf %lf %lf\n", g_Tmu_tab + g_ntabmu, &dummy1, &dummy2,
                  g_mu_tab + g_ntabmu)!=EOF) {
      g_ntabmu++;
    }
  }
}

double GetMuFromTable(double T)
{
  int    klo, khi, kmid;
  double  Tmid, dT;

  if (T > 2.e5) {return (double)g_mu_tab[g_ntabmu-1];} //Above this temperature mu is approximately constant
  if (T < 1.e4) return (double)g_mu_tab[0]; //Below this temperature mu is approximately constant

  klo = 0;
  khi = g_ntabmu - 1;

  if (T > g_Tmu_tab[khi])
  {
    g_oof++;//print (" ! T out of range   %12.6e\n",T);
    //    QUIT_PLUTO(1);
    return g_mu_tab[khi];
  }
  else if (T < g_Tmu_tab[klo])
  {
    g_oof++;//print (" ! T out of range   %12.6e\n",T);
    //    QUIT_PLUTO(1);
    return g_mu_tab[klo];
  }
  else
  {
    /* ----------------------------------------------
              Table lookup by binary search
       ---------------------------------------------- */
    while (klo != (khi - 1))
    {
      kmid = (klo + khi)/2;
      Tmid = g_Tmu_tab[kmid];
      if (T <= Tmid){
        khi = kmid;
      }else if (T > Tmid){
        klo = kmid;
      }
    }
    dT = g_Tmu_tab[khi] - g_Tmu_tab[klo];
    return g_mu_tab[klo]*(g_Tmu_tab[khi] - T)/dT + g_mu_tab[khi]*(T - g_Tmu_tab[klo])/dT;
  }
}

double fT(double T, void *par)
{
  struct temperature_params *p = (struct temperature_params *) par;
  p->mu = GetMuFromTable(T);
  //if (T < 2e5) {print("T: %f mu: %f\n",T,*mu);}
  //  print("T, mu, T_prs, T_mu: %f %f %f %f\n",T,*mu, prs*KELVIN/rho, T/(*mu));
  //ncall++;
  //if (ncall % 10000 == 0) {print("Process %d: ncall: %d\n", prank, ncall);}
  return p->prs*KELVIN/p->rho*p->mu - T;
}


/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 * 
 ******************************************************************* */
{
#ifdef WOLFIRE95_CIICOOL
  g_minCoolingTemp = 1.0;
#else
  g_minCoolingTemp = 1.e4;
#endif
  int    klo, khi, kmid;
  static int ntab;
  double  mu, T, Tmid, scrh, dT, prs;
  static double *L_tab_z0, *L_tab_z05, *L_tab_z1, *T_tab, E_cost;
#ifdef WOLFIRE95_CIICOOL
  static int ntab_CII;
  static double *L_tab_CII, *T_tab_CII;
#endif
#ifdef WOLFIRE95_UVHEATING
  static int ntab_heat;
  static double *T_tab_heat, *Gamma_tab;
#endif
  //static int niter=0;  
  FILE *fcool;
  int OutOfBounds_low = 0;
  int OutOfBounds_hi = 0;

  //Temperatures below/above which mu=const. Should match the table used for mu(T).
  const double T_neutral = 10000;
  const double T_ionized = 15000;

/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

//Modified to read multiple tables with different metallcities
//IMPORTANT: All tables must have the exact same number of lines
//and the same temperatures

  if (T_tab == NULL){
    double dummy = 0.0;
    print1 (" > Reading tables from disk...\n");
    fcool = fopen("cooltable-z0.dat","r");
    if (fcool == NULL){
      print1 ("! Radiat: cooltable-z0.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab_z0 = ARRAY_1D(7000, double);
    T_tab = ARRAY_1D(7000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab, 
                                       L_tab_z0 + ntab)!=EOF) {
      ntab++;
    }
    fclose(fcool);

    fcool = fopen("cooltable-z05.dat","r");
    if (fcool == NULL){
      print1 ("! Radiat: cooltable-z05.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab_z05 = ARRAY_1D(7000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", &dummy, 
                                       L_tab_z05 + ntab)!=EOF) {
      ntab++;
    }
    fclose(fcool);

    fcool = fopen("cooltable-z1.dat","r");
    if (fcool == NULL){
      print1 ("! Radiat: cooltable-z1.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab_z1 = ARRAY_1D(7000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", &dummy, 
                                       L_tab_z1 + ntab)!=EOF) {
      ntab++;
    }
    fclose(fcool);
    E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);

#ifdef WOLFIRE95_CIICOOL
//For T < 10^4 K cooling, generally dominated by CII. Load tabulated cooling curve based on Wolfire et al. (1995).
//This does not take metallicity into account 
    fcool = fopen("cooltable_wolfire95.dat","r");
    if (fcool == NULL){
      print1 ("! Radiat: cooltable-wolfire95.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    T_tab_CII = ARRAY_1D(1000, double);
    L_tab_CII = ARRAY_1D(1000, double);

    ntab_CII = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab_CII + ntab_CII, 
                                       L_tab_CII + ntab_CII)!=EOF) {
      ntab_CII++;
    }
    fclose(fcool);
    //Join boundaries of tables together
    L_tab_z0[0] = L_tab_CII[ntab_CII-1];
    L_tab_z05[0] = L_tab_CII[ntab_CII-1];
    L_tab_z1[0] = L_tab_CII[ntab_CII-1];
#endif
#if defined(WOLFIRE95_UVHEATING) && !defined(ANALYTIC_HEATING)
    fcool = fopen("heattable_w95.dat","r");
    if (fcool == NULL){
      print1 ("! Radiat: heattable-w95.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    T_tab_heat = ARRAY_1D(1000, double);
    Gamma_tab = ARRAY_1D(1000, double);

    ntab_heat = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab_heat + ntab_heat,
                                       Gamma_tab + ntab_heat)!=EOF) {
      ntab_heat++;
    }
    fclose(fcool);
#endif
  }


/* ---------------------------------------------
            Get pressure and temperature 
   --------------------------------------------- */

  prs = v[RHOE]*(g_gamma-1.0);
  if (prs < 0.0) {
    prs     = g_smallPressure;
    v[RHOE] = prs/(g_gamma - 1.0);
  }



  double a, b, c;
  int count, maxiter;
  double tolerance;
  count = 0;
  tolerance = 0.05;

  double rho = v[RHO];

  a = prs*((double)g_mu_tab[0]+tolerance)*KELVIN/rho;
  b = prs*((double)g_mu_tab[g_ntabmu-1]-tolerance)*KELVIN/rho;

  if ((a > T_ionized) && (b > T_ionized)) //Temperature is guaranteed to be sufficient for the medium to be fully ionized
  {
    mu = g_mu_tab[g_ntabmu-1];
    T = prs*mu*KELVIN/rho;
  }
  else if ((a < T_neutral) && (b < T_neutral)) //Temperature is guaranteed to be low enough for the medium to be fully neutral
  {
    mu = g_mu_tab[0];
    T = prs*mu*KELVIN/rho;
  }
  else
  {
    //Temperature might be in between fully neutral and fully ionized
    //Find temperature and mean molecuar weight as the root of mu(T)*P/(k_B*rho) - T=0 using Brent's method
    struct temperature_params par;
    par.mu = -1;
    par.rho = rho;
    par.prs = prs;
    int status = Brent(fT, &par, a, b, -1, 1.e-12, &T);
    if (status != 0)
    {
      print("! radiat: Failed to find root of temperature function on proc %d\nStatus: %d\n", prank, status);
      QUIT_PLUTO(1);
    }
    mu = par.mu;
  }

  if (T != T){
    printf (" ! Nan found in radiat \n");
    printf (" ! rho = %12.6e, prs = %12.6e\n",v[RHO], prs);
    QUIT_PLUTO(1);
  }

  rhs[RHOE] = 0.0;

  double ndens = v[RHO]*UNIT_DENSITY/(mu*CONST_mp);

#ifdef WOLFIRE95_UVHEATING
#ifdef ANALYTIC_HEATING //Use analytical heating rate calculation based on Wolfire et al. (1995)
  const double G0 = 0.463282; //G0=1.1/(1+(z/8.53 kpc)^2) with z = 10 kpc following Heitsch & Putman (2009)
  double eff = 0.0;
  //if (T < 1.1e4) eff = 0.05
  //else
  //{
    double electron_frac = 0.519049;
    double Teff = T;
    if (T < 1.e4) Teff = 1.e4;
    if (T < 7.e4) electron_frac = 0.494062/(1.0+exp(-0.0008*(Teff-1.5e4))); //Sigmoid fit to n_e/n from the model of Sutherland & Dopita (1993) when not fully ionized
    double n_e = ndens*electron_frac;
    double phi = G0*sqrt(Teff)/n_e; //Grain charging parameter
    eff = 4.9e-2/(1.0+pow(phi/1925.0,0.73))+3.7e-2*pow(Teff/1.e4,0.7)/(1.0+(phi/5000.0)); //Heating efficiency for dust grains from Bakes & Tielens (1994)
  //}
  scrh = 1.e-24 * eff * G0;
#else //Use tabulated heating rates
  if (T < T_tab_heat[0]) scrh = Gamma_tab[0];
  else
  {
    klo = 0;
    khi = ntab_heat - 1;
    while (klo != (khi - 1))
    {
      kmid = (klo + khi)/2;
      Tmid = T_tab_heat[kmid];
      if (T <= Tmid)
      {
       	khi = kmid;
      }
      else if (T > Tmid)
      {
       	klo = kmid;
      }
    }
    dT       = T_tab_heat[khi] - T_tab_heat[klo];
    scrh = Gamma_tab[klo]*(T_tab_heat[khi] - T)/dT + Gamma_tab[khi]*(T - T_tab_heat[klo])/dT;
  }
#endif //ANALYTIC_HEATING
  rhs[RHOE] = scrh*ndens*E_cost;
//  print1("T: %e, rho: %e, Gamma_cgs: %e, Gamma_low: %e, Gamma_hi: %e, Gamma: %e, rhs_cgs: %e, rhs: %e\n", T, v[RHO], scrh, Gamma_tab[klo], Gamma_tab[khi], scrh*E_cost, rhs[RHOE]/E_cost, rhs[RHOE]);
#else
  if (T < g_minCoolingTemp)
  { 
    rhs[RHOE] = 0.0;
    return;
  }
#endif

#ifdef WOLFIRE95_CIICOOL
  if (T < T_tab[0])
  {
    klo = 0;
    khi = ntab_CII - 1;
    while (klo != (khi - 1))
    {
      kmid = (klo + khi)/2;
      Tmid = T_tab_CII[kmid];
      if (T <= Tmid)
      {
        khi = kmid;
      }
      else if (T > Tmid)
      {
        klo = kmid;
      }
    }
    dT       = T_tab_CII[khi] - T_tab_CII[klo];
    scrh = L_tab_CII[klo]*(T_tab_CII[khi] - T)/dT + L_tab_CII[khi]*(T - T_tab_CII[klo])/dT;
    rhs[RHOE] -= scrh*ndens*ndens*E_cost;
    return;
  }
#endif

  int Zstatus = 0;
  double *L_tab_loZ = NULL;
  double *L_tab_hiZ = NULL;
  double Zhi = 0.0;
  double Zlo = 0.0;
  if (v[TRC_Z] > 1.0) Zstatus = 2; //Metallicity greater than max table
  else if (v[TRC_Z] > 0.3162)
  {
    Zhi = 1.0;
    Zlo = 0.3162;
    L_tab_loZ = L_tab_z05;
    L_tab_hiZ = L_tab_z0;
  }
  else if (v[TRC_Z] > 0.1)
  {
    Zhi = 0.3162;
    Zlo = 0.1;
    L_tab_loZ = L_tab_z1;
    L_tab_hiZ = L_tab_z05;
  }
  else Zstatus = 1; //Metallicity lower than min table
  double dZ = Zhi - Zlo;

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

  klo = 0;
  khi = ntab - 1;

  if (T > T_tab[khi]){
    OutOfBounds_hi++;
    //    print (" ! T out of range   %12.6e\n",T);
//    QUIT_PLUTO(1);
    if (Zstatus == 1) scrh = L_tab_z1[khi];
    else if (Zstatus == 2) scrh = L_tab_z0[khi];
    else scrh = L_tab_loZ[khi]*(Zhi - v[TRC_Z])/dZ + L_tab_hiZ[khi]*(v[TRC_Z] - Zlo)/dZ;

  }
  else if (T < T_tab[klo]){
    OutOfBounds_low++;
    //    print (" ! T out of range   %12.6e\n",T);
//    QUIT_PLUTO(1);
    if (Zstatus == 1) scrh = L_tab_z1[klo];
    else if (Zstatus == 2) scrh = L_tab_z0[klo];
    else scrh = L_tab_loZ[klo]*(Zhi - v[TRC_Z])/dZ + L_tab_hiZ[klo]*(v[TRC_Z] - Zlo)/dZ;
  }
  else
  {
    while (klo != (khi - 1))
    {
      kmid = (klo + khi)/2;
      Tmid = T_tab[kmid];
      if (T <= Tmid)
      {
        khi = kmid;
      }
      else if (T > Tmid)
      {
        klo = kmid;
      }
    }
    dT       = T_tab[khi] - T_tab[klo];
    if (Zstatus == 1) scrh = L_tab_z1[klo]*(T_tab[khi] - T)/dT + L_tab_z1[khi]*(T - T_tab[klo])/dT;
    else if (Zstatus == 2) scrh = L_tab_z0[klo]*(T_tab[khi] - T)/dT + L_tab_z0[khi]*(T - T_tab[klo])/dT;
    else
    {
      double L_loZ = L_tab_loZ[klo]*(T_tab[khi] - T)/dT + L_tab_loZ[khi]*(T - T_tab[klo])/dT;
      double L_hiZ = L_tab_hiZ[klo]*(T_tab[khi] - T)/dT + L_tab_hiZ[khi]*(T - T_tab[klo])/dT;
      scrh     = L_loZ*(Zhi - v[TRC_Z])/dZ + L_hiZ*(v[TRC_Z] - Zlo)/dZ;
    }
  }

  rhs[RHOE] -= scrh*ndens*ndens*E_cost;

//  print("T, Z, L, rhs: %6.2e %f %6.2e %6.2e\n", T, v[TRC_Z], scrh, rhs[RHOE]);


#ifdef CELLINFO_VERBOSE
  if (OutOfBounds_low > 0)  print("! radiat: %d cells on proc %d had temperatures below the minimum in table\n", OutOfBounds_low, prank);
  if (OutOfBounds_hi > 0)  print("! radiat: %d cells on proc %d had temperatures above the maximum in table\n", OutOfBounds_hi, prank);
#endif
}
#undef T_MIN
/* ******************************************************************* */
double MeanMolecularWeight (double *V)
/*
 *
 *
 *
 ********************************************************************* */
{
  return (0.59);//v[TRC_MU]);
/*
  return  ( (A_H + frac_He*A_He + frac_Z*A_Z) /
            (2.0 + frac_He + 2.0*frac_Z - 0.0));
*/
}



