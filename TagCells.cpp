#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"
extern "C"
{
#include "mutemp.h"
}

static void computeRefVar(double ***UU[], double ***q, double, RBox *Ubox);

#if (EOS != ISOTHERMAL) && (CHOMBO_EN_SWITCH == NO)
 #ifndef CHOMBO_REF_VAR  
  #define CHOMBO_REF_VAR ENG 
 #endif
#else
 #ifndef CHOMBO_REF_VAR  
  #define CHOMBO_REF_VAR RHO
 #endif
#endif

#define TRC_Z TRC+1
#define REF_CRIT 2   /* 1 == first derivative, 2 == second derivative */

/* ************************************************************************* */
void PatchPluto::computeRefGradient(FArrayBox& gFab, FArrayBox& UFab, 
                                    const FArrayBox& a_dV, const Box& b)
/*!
 * Tag zones for refinement using gradient of the conservative 
 * variables.
 * The gradient is computed by standard finite differences using
 *
 * - REF_CRIT equal to 1 --> compute (normalized) gradient using 1st 
 *                           derivative of the solution;
 * - REF_CRIT equal to 2 --> compute (normalized) gradient using 2nd 
 *                           derivative of the solution (default);
 *                           This approach is based on Lohner (1987).
 *
 * Zones will be flagged for refinement whenever grad[k][j][i] exceeds 
 * the threshold value specified by the 'Refine_thresh' parameter read in
 * pluto.ini.
 *
 * Derivatives are computed using the conserved variable
 * U[CHOMBO_REF_VAR] 
 * where CHOMBO_REF_VAR is taken to be energy density (default).
 * However, by setting CHOMBO_REF_VAR = -1, you can provide your own 
 * physical variable through the function computeRefVar().
 * 
 * \authors C. Zanni   (zanni@oato.inaf.it)\n
 *          A. Mignone (mignone@ph.unito.it)
 * \date    Oct 11, 2012
 *************************************************************************** */
{
  CH_assert(m_isDefined);

  int nv, i, j, k;
  double x1, dqx_p, dqx_m, dqx, d2qx, den_x;
  double x2, dqy_p, dqy_m, dqy, d2qy, den_y;
  double x3, dqz_p, dqz_m, dqz, d2qz, den_z;
  double gr1, gr2, eps = 0.01;
 #if CHOMBO_REF_VAR == -1 || REF_CRIT == 3
  double ***UU[NVAR];
 #endif
  double ***q, ***grad, ***bmag;
  RBox  Ubox, Gbox;
  static double RefThresh = atof(ParamFileGet("Refine_thresh",1));
  int nrefgrad = 0;

#if REF_CRIT == 3
  int    klo, khi, kmid;
  static int ntab;
  double  mu, T, Tmid, scrh, dT, prs, Lambda;
  static double *L_tab, *T_tab, E_cost;

  FILE *fcool;
  static int niter = 0;
  int nrefcool = 0;

  /* -------------------------------------------                                                                                                                                           
        Read tabulated cooling function                                                                                                                                                  
  ------------------------------------------- */

  if (T_tab == NULL){
    print1 (" > Reading table from disk...\n");
    fcool = fopen("cooltable-z05.dat","r");
    if (fcool == NULL){
      print1 ("! Radiat: cooltable-z05.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab,
		  L_tab + ntab)!=EOF) {
      ntab++;
    }
    E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
  }
  int OutOfBounds_low = 0;
  int OutOfBounds_hi = 0;
#endif //REF_CRIT == 3

/* -- check ref criterion -- */

  #if REF_CRIT != 1 && REF_CRIT != 2 && REF_CRIT != 3
   print ("! TagCells.cpp: Refinement criterion not valid\n");
   QUIT_PLUTO(1);
  #endif

/* -----------------------------------------------------
   1. The solution array U is defined on the box 
      [Uib, Uie] x [Ujb, Uje] x [Ukb, Uke], which 
      differs from that of gFab ([Gib,...Gke]), 
      typically one point larger in each direction. 
   ----------------------------------------------------- */
    
  Ubox.jb = Ubox.je = Ubox.kb = Ubox.ke = 0;
  Gbox.jb = Gbox.je = Gbox.kb = Gbox.ke = 0;

  D_EXPAND(Ubox.ib = UFab.loVect()[IDIR]; Ubox.ie = UFab.hiVect()[IDIR]; ,
           Ubox.jb = UFab.loVect()[JDIR]; Ubox.je = UFab.hiVect()[JDIR]; ,
           Ubox.kb = UFab.loVect()[KDIR]; Ubox.ke = UFab.hiVect()[KDIR]; );

  D_EXPAND(Gbox.ib = gFab.loVect()[IDIR]; Gbox.ie = gFab.hiVect()[IDIR]; ,
           Gbox.jb = gFab.loVect()[JDIR]; Gbox.je = gFab.hiVect()[JDIR]; ,
           Gbox.kb = gFab.loVect()[KDIR]; Gbox.ke = gFab.hiVect()[KDIR]; );

/* --------------------------------------------------------
   2. Input solution array (UFab.dataPtr(nv)) is defined 
      as dV*U/dx^3, where U is an array of conservative 
      variables and dV is the zone volume. 
      To obtain U we must divide by volume.
   -------------------------------------------------------- */

  #if CHOMBO_REF_VAR == -1 || REF_CRIT == 3
   FArrayBox tmpU(UFab.box(),NVAR);
   tmpU.copy(UFab);
  #else
   FArrayBox tmpU(UFab.box(),1);
   tmpU.copy(UFab,CHOMBO_REF_VAR,0);
  #endif 

  #if GEOMETRY != CARTESIAN

   #if CHOMBO_REF_VAR == -1

    for (nv = 0; nv < NVAR; nv++) tmpU.divide(a_dV,0,nv);

    #if CHOMBO_CONS_AM == YES
     #if ROTATING_FRAME == YES
      Box curBox = UFab.box();
      for(BoxIterator bit(curBox); bit.ok(); ++bit) {
        const IntVect& iv = bit();
        tmpU(iv,iMPHI) /= a_dV(iv,1);
        tmpU(iv,iMPHI) -= tmpU(iv,RHO)*a_dV(iv,1)*g_OmegaZ;
      }
     #else
      tmpU.divide(a_dV,1,iMPHI);
     #endif
    #endif

   #else

    tmpU.divide(a_dV,0,0);

   #endif

  #else

   if (g_stretch_fact != 1.) tmpU /= g_stretch_fact;

  #endif

/* ---------------------------------------------
   3. Set refinement variable
   --------------------------------------------- */

  #if CHOMBO_REF_VAR == -1
   for (nv = 0; nv < NVAR; nv++)
     UU[nv] = ArrayBoxMap(Ubox.kb, Ubox.ke,
                          Ubox.jb, Ubox.je,
                          Ubox.ib, Ubox.ie, tmpU.dataPtr(nv));
   q = ArrayBox(Ubox.kb, Ubox.ke, Ubox.jb, Ubox.je, Ubox.ib, Ubox.ie);
   computeRefVar(UU, q, m_dx, &Ubox);
  #elif REF_CRIT == 3
   for (nv = 0; nv < NVAR; nv++)
     UU[nv] = ArrayBoxMap(Ubox.kb, Ubox.ke,
                          Ubox.jb, Ubox.je,
                          Ubox.ib, Ubox.ie, tmpU.dataPtr(nv));


   /* ----------------------------------------------
     convert from conservative to primitive 
     ---------------------------------------------- */

   static unsigned char *flag;
   static double **u, **v;

   if (u == NULL){
     u    = ARRAY_2D(NMAX_POINT, NVAR, double);
     v    = ARRAY_2D(NMAX_POINT, NVAR, double);
     flag = ARRAY_1D(NMAX_POINT, unsigned char);
   }

   /* --------------------------------------------------------
        Recover total energy if necessary.
	-------------------------------------------------------- */

  #if CHOMBO_EN_SWITCH == YES
   totEnergySwitch (UU, Ubox.ib, Ubox.ie, 
		    Ubox.jb, Ubox.je, 
		    Ubox.kb, Ubox.ke, +1);
  #endif

   for (k = Ubox.kb; k <= Ubox.ke; k++){
     for (j = Ubox.jb; j <= Ubox.je; j++){
       for (i = Ubox.ib; i <= Ubox.ie; i++){
	 for (nv=0; nv < NVAR; nv++){
	   u[i-Ubox.ib][nv] = UU[nv][k][j][i];
	 }}
       ConsToPrim (u, v, 0, Ubox.ie-Ubox.ib, flag);
       for (i = Ubox.ib; i <= Ubox.ie; i++){
	 for (nv=0 ; nv < NVAR; nv++){
	   UU[nv][k][j][i] = v[i-Ubox.ib][nv];
	 }}
     }}

  q = ArrayBox(Ubox.kb, Ubox.ke, Ubox.jb, Ubox.je, Ubox.ib, Ubox.ie);
  BOX_LOOP(&Ubox, k, j, i) {
    q[k][j][i] = UU[RHO][k][j][i];
  }
  #else
   q = ArrayBoxMap(Ubox.kb, Ubox.ke,
                   Ubox.jb, Ubox.je,
                   Ubox.ib, Ubox.ie, tmpU.dataPtr(0));
  #endif
  
  grad = ArrayBoxMap(Gbox.kb, Gbox.ke, 
                     Gbox.jb, Gbox.je, 
                     Gbox.ib, Gbox.ie, gFab.dataPtr(0));

/* ----------------------------------------------------------------
   4. Main spatial loop for zone tagging based on 1st 
     (REF_CRIT = 1) or 2nd (REF_CRIT = 2) derivative error norm. 
   ---------------------------------------------------------------- */

  BOX_LOOP(&Gbox, k, j, i){
    x3 = (k + 0.5)*m_dx*g_x3stretch + g_domBeg[KDIR];
    x2 = (j + 0.5)*m_dx*g_x2stretch + g_domBeg[JDIR];
   #if CHOMBO_LOGR == NO
    x1 = (i + 0.5)*m_dx          + g_domBeg[IDIR];
   #else
    double xl = g_domBeg[IDIR] + i*m_dx;
    double xr = xl + m_dx;
    x1 = g_domBeg[IDIR]*0.5*(exp(xr)+exp(xl));
   #endif 

    D_EXPAND(dqx_p =    q[k][j][i+1] - q[k][j][i];
             dqx_m = - (q[k][j][i-1] - q[k][j][i]);  ,
             dqy_p =    q[k][j+1][i] - q[k][j][i];
             dqy_m = - (q[k][j-1][i] - q[k][j][i]);  ,
             dqz_p =    q[k+1][j][i] - q[k][j][i];
             dqz_m = - (q[k-1][j][i] - q[k][j][i]);)

  /* --------------------------------------------------------------
      Physical boundary values are not up to date and should be 
      excluded from gradient computation. 
      In this case, left and right derivatives are set equal to 
      each other. This will not trigger refinement in the leftmost 
      and rightmost internal zones (using 2nd derivative) but we 
      really don't care since buffer size will do the job.
     -------------------------------------------------------------- */
      
    D_EXPAND(if (i == 0) dqx_m = dqx_p;  ,
             if (j == 0) dqy_m = dqy_p;  ,
             if (k == 0) dqz_m = dqz_p;)

    D_EXPAND(if (i == m_domain.size(IDIR)-1) dqx_p = dqx_m;  ,
             if (j == m_domain.size(JDIR)-1) dqy_p = dqy_m;  ,
             if (k == m_domain.size(KDIR)-1) dqz_p = dqz_m;)

  /* -----------------------------------------------
         Compute gradient using 1st derivative 
      ---------------------------------------------- */

    #if REF_CRIT == 1
     D_EXPAND(dqx = dqx_p + dqx_m;  ,
              dqy = dqy_p + dqy_m;  ,
              dqz = dqz_p + dqz_m;)

     D_EXPAND(den_x = fabs(q[k][j][i+1]) + fabs(q[k][j][i-1]);  ,
              den_y = fabs(q[k][j+1][i]) + fabs(q[k][j-1][i]);  ,
              den_z = fabs(q[k+1][j][i]) + fabs(q[k-1][j][i]);)

     gr1  = D_EXPAND(dqx*dqx, + dqy*dqy, + dqz*dqz);
     gr1 /= D_EXPAND(den_x*den_x, + den_y*den_y, + den_z*den_z);

     grad[k][j][i] = sqrt(gr1);
    #endif

  /* -----------------------------------------------
         Compute gradient using 2nd derivative 
      ---------------------------------------------- */

    #if REF_CRIT == 2 || REF_CRIT == 3
     D_EXPAND(d2qx = dqx_p - dqx_m;  ,
              d2qy = dqy_p - dqy_m;  ,
              d2qz = dqz_p - dqz_m;)

     D_EXPAND(
       den_x = 2.0*fabs(q[k][j][i]) + fabs(q[k][j][i+1]) + fabs(q[k][j][i-1]);
       den_x = fabs(dqx_p) + fabs(dqx_m) + eps*den_x;    ,

       den_y = 2.0*fabs(q[k][j][i]) + fabs(q[k][j+1][i]) + fabs(q[k][j-1][i]);
       den_y = fabs(dqy_p) + fabs(dqy_m) + eps*den_y;    ,

       den_z = 2.0*fabs(q[k][j][i]) + fabs(q[k+1][j][i]) + fabs(q[k-1][j][i]);
       den_z = fabs(dqz_p) + fabs(dqz_m) + eps*den_z;
     )

     gr2  = D_EXPAND(d2qx*d2qx,   + d2qy*d2qy,   + d2qz*d2qz);
     gr2 /= D_EXPAND(den_x*den_x, + den_y*den_y, + den_z*den_z);

     grad[k][j][i] = sqrt(gr2);
     
    if (grad[k][j][i] > RefThresh) {nrefgrad++;}
    #endif

    #if REF_CRIT == 3
    if (grad[k][j][i] < RefThresh) //Don't bother if the cell is going to be flagged for refinement anyway
    {
      /* ---------------------------------------------                                                                                                                                         
                Get pressure and temperature                                                                                                                                                 
      --------------------------------------------- */

      /*      static int ntabmu;
      static double *mu_tab, *Tmu_tab;
      FILE *fmu;
      double cmu;*/
      double a, b, c;
      int count, maxiter;
      double tolerance = 0.05;
/*      
      double dummy1, dummy2;

      if (Tmu_tab == NULL) //Read mu table
      {
	  print1 (" > Reading mu table from disk...\n");
	  fmu = fopen("mutable.dat","r");
	  if (fmu == NULL){
	    print1 ("! Radiat: mutable.dat could not be found.\n");
	    QUIT_PLUTO(1);
	  }
	  mu_tab = ARRAY_1D(20000, double);
	  Tmu_tab = ARRAY_1D(20000, double);

	  ntabmu = 0;
	  while (fscanf(fmu, "%lf %lf %lf %lf\n", Tmu_tab + ntabmu, &dummy1, &dummy2,
			mu_tab + ntabmu)!=EOF) {
	    ntabmu++;
	  }
	  }*/

      //  print("Ntab: %d\n",ntab);
        count = 0;
      tolerance = 0.05;
      maxiter = 100;
/*      cmu = 0.0;
      int niter;
      niter = 0;
*/

      double rho = UU[RHO][k][j][i];
      double prs = UU[PRS][k][j][i];
      //Find temperature and mu as root of mu(T)*P/(k_B*rho) - T=0 using secant method
      a = prs*((double)g_mu_tab[0]+tolerance)*KELVIN/rho;
      b = prs*((double)g_mu_tab[g_ntabmu-1]-tolerance)*KELVIN/rho;
      //      do
      //{
	  /*
	    if(f(a)==f(b))
	    {
	    printf("\nSolution cannot be found as the values of a and b are same.\n");
	    return;
	    }
	  */
      /*	  c=(a*fT(b, &cmu, prs, rho, T_tab, mu_tab, ntab)-b*fT(a, &cmu, prs, rho, T_tab, mu_tab, ntab))/(fT(b, &cmu, prs, rho, T_tab, mu_tab, ntab)-fT(a, &cmu, prs, rho, T_tab, mu_tab, ntab));

	  a=b;
	  b=c;
	  //  print("Iteration No-%d    a,b,c=%f %f %f\n",count,a,b,c);
	  count++;
	  if(count == maxiter)
	    {
	      printf("T did not converge within %d steps!\n", maxiter);
	      break;
	    }
      } while(fabs(c-b)/b > tolerance);  //fabs(fT(c, &cmu, d, T_tab, mu_tab, ntab, i, j, k)) > tolerance);

      fT(c, &cmu, prs, rho, T_tab, mu_tab, ntab); //Update mu
      T = c;
      mu = cmu;*/
      temperature_params par;
      par.mu = -1;
      par.rho = rho;
      par.prs = prs;
      int status = Brent(fT, &par, a, b, -1, 1.e-12, &T);
      if (status != 0)
      {
	print("! TagCells: Failed to find root of temperature function on proc %d at %f, %f, %f (%d, %d, %d)\nStatus: %d\n", prank, x1,x2,x3,i,j,k,status);
	QUIT_PLUTO(1);
      }
      mu = par.mu;
      //    if (niter==1000) {print("i,j,k: %d %d %d T: %f mu: %f rho: %f\n",i,j,k,T[k][j][i],mu[k][j][i],d->Vc[RHO][k][j][i]);niter=0;}
      niter++;


      //  mu  = g_inputParam[MU_C];
      //  T   = UU[PRS][k][j][i]/UU[RHO][k][j][i]*KELVIN*mu;
      //  print1("T, prs, mu, gamma: %12.6e %6.4e %12.6e\n", T, UU[PRS][k][j][i], mu, g_gamma);

      if (T != T){
	printf (" ! Nan found in TagCells on proc %d at %f, %f, %f (%d, %d, %d)\n", prank, x1, x2, x3, i, j, k);
	printf (" ! rho = %12.6e, prs = %12.6e\n",UU[RHO][k][j][i], UU[PRS][k][j][i]);
	QUIT_PLUTO(1);
      }
      if (T < g_minCoolingTemp) {
	Lambda = 0.0;
      }
      else
      {
	/* ----------------------------------------------                                                                                                                                        
	     Table lookup by binary search                                                                                                                                                    
         ---------------------------------------------- */
	klo = 0;
	khi = ntab - 1;
	if (T > T_tab[khi]){
          OutOfBounds_low++;
	  //print (" ! T out of range   %12.6e\n",T);
	  //    QUIT_PLUTO(1);
	  scrh     = UU[TRC_Z][k][j][i]*L_tab[khi]/0.3162;
	}
	else if (T < T_tab[klo]){
          OutOfBounds_hi++;
	  //print (" ! T out of range   %12.6e\n",T);
	  //    QUIT_PLUTO(1);
	  scrh     = UU[TRC_Z][k][j][i]*L_tab[0]/0.3162;
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
	  scrh     = (UU[TRC_Z][k][j][i]/0.3162)*(L_tab[klo]*(T_tab[khi] - T)/dT + L_tab[khi]*(T - T_tab[klo])/dT);
	}

	Lambda = scrh*UU[RHO][k][j][i]*UU[RHO][k][j][i];
	Lambda *= UNIT_DENSITY*UNIT_DENSITY/(mu*mu*CONST_mp*CONST_mp);
	double t_cool = 3*CONST_kB*T*mu/(2*UU[RHO][k][j][i]*Lambda);

	double c_s = sqrt(g_gamma*UU[PRS][k][j][i]/UU[RHO][k][j][i])*UNIT_VELOCITY;
	double t_sc = m_dx*UNIT_LENGTH/c_s;
	//      if (niter == 100000) {print("T, n, c_s, delta_x, L, t_cool, t_sc: %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n", T, UU[RHO][k][j][i]/mu, c_s, m_dx*UNIT_LENGTH, Lambda, t_cool, t_sc);niter=0;}
	//    niter++;
	if (t_cool < t_sc)
	{
          //Tag cell for refinement
	  //print("Cooling refinement cell at: %d %d %d\n", i, j , k);
	  grad[k][j][i] = 100.0;
          nrefcool++;
	}
      }
    }
    #endif //REF_CRIT == 3
  }
  #ifdef CELLINFO_VERBOSE
  if (nrefgrad > 0) print("%d cells refined based on gradient on proc %d\n", nrefgrad, prank);
  #if REF_CRIT == 3
  if (nrefcool > 0) print("%d cells refined based on cooling time on proc %d\n", nrefcool, prank);
  if (OutOfBounds_low > 0)  print("! TagCells: %d cells on proc %d had temperatures below the minimum in table\n", OutOfBounds_low, prank);
  if (OutOfBounds_hi > 0)  print("! TagCells: %d cells on proc %d had temperatures above the maximum in table\n", OutOfBounds_hi, prank);
  #endif
  #endif

/* --------------------------------------------------------------
   6. Free array
   -------------------------------------------------------------- */
   
  FreeArrayBoxMap(grad, Gbox.kb, Gbox.ke, Gbox.jb, Gbox.je, Gbox.ib, Gbox.ie);

  #if CHOMBO_REF_VAR == -1
   for (nv = 0; nv < NVAR; nv++){
     FreeArrayBoxMap(UU[nv], Ubox.kb, Ubox.ke,
                             Ubox.jb, Ubox.je, Ubox.ib, Ubox.ie);
   }
   FreeArrayBox(q, Ubox.kb, Ubox.jb, Ubox.ib);
  #elif REF_CRIT == 3
   for (nv = 0; nv < NVAR; nv++){
     FreeArrayBoxMap(UU[nv], Ubox.kb, Ubox.ke,
                             Ubox.jb, Ubox.je, Ubox.ib, Ubox.ie);
   }
   FreeArrayBox(q, Ubox.kb, Ubox.jb, Ubox.ib);
  #else
   FreeArrayBoxMap(q, Ubox.kb, Ubox.ke, Ubox.jb, Ubox.je, Ubox.ib, Ubox.ie);
  #endif

}

/* ********************************************************************* */
void computeRefVar(double ***UU[], double ***q, double dx, RBox *Ubox)
/*!
 * Compute a user-defined array q(U) function of the conserved
 * variables.
 *
 *
 *********************************************************************** */
{
  int nv, i, j, k;
  double us[NVAR], vs[NVAR];

  BOX_LOOP(Ubox, k, j, i) {
    VAR_LOOP(nv) us[nv] = UU[nv][k][j][i]; 
    q[k][j][i] = us[RHO];
  }

}
#undef REF_CRIT
