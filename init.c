/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author Asger Gronnow (asger.gronnow@usyd.edu.au)
  \date   Nov 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "mutemp.h"

/*struct temperature_params
{
  double rho;
  double prs;
  double mu;
  };*/


#define PRS_FROM_TEMP
//#define VELOCITY_MACH
#define TRC_CLOUD TRC
#define TRC_Z TRC+1
//#define TRC_MU TRC+4

#define UNIT_MICROGAUSS 2.181
#define BINIT_MICROGAUSS
#define USE_TANH_DENSPROFILE
#define VAR_MU

 //#define WOLFIRE95_CIICOOL

 //#define USE_HP09_PROFILE

#define HYDROSTATIC_HALO
#define UNIT_TIME (UNIT_LENGTH/UNIT_VELOCITY)
#define UNIT_CMPERS2 (UNIT_LENGTH/(UNIT_TIME*UNIT_TIME))
#define UNIT_KPCPERMYR2 ((UNIT_TIME/3.154e13)*(UNIT_TIME/3.154e13))
#define UNIT_KPCPERMYR_TO_CMPERS2 (CONST_pc*1000.0/(3.154e13*3.154e13))

#if defined(VARIABLE_INFLOW) || defined(HYDROSTATIC_HALO)
double g_tabulardistance[1024];
double g_tabulardensity[1024];
double g_tabularbfield[1024];
int g_curidx = 0;

/* ********************************************************************* */
int ReadTable(char TableFile[256])
/*!
 * Read table of z-distance, halo density and magnetic field strength
 * used when calculating variable inflow.
 *
 * \param [in] TableFile     Path to file containing table with distance
 *                           in kpc in the first column, halo number
 *                           density per cubic cm in the second column
 *                           and halo magnetic field magnitude in
 *                           microGauss in the third column.
 *
 * \return 0 on success or one of the following negative integers in
 *         case of an error:
 *         -1 Table file could not be opened
 *         -2 Table has too many data points (> 1024)
 *********************************************************************** */
{
  int i = 0;
  FILE* table = fopen(TableFile, "r");
  if (table == NULL) {return -1;}
  else
    {
      char line[256];
      double dist;
      double dens;
      double bfield;
      while (fgets(line, sizeof(line), table) != NULL) 
	{
	  if (line[0] != '#')
	    {
	      sscanf(line, "%lf %lf %lf", &dist, &dens, &bfield);
	      g_tabulardistance[i] = dist;
	      g_tabulardensity[i] = dens;
	      g_tabularbfield[i] = bfield;
	      i++;
	    }
	  if (i > 1023) {return -2;}
	}
    }
  fclose(table);
  g_curidx = i - 1;
  return 0;
}

#endif

double CenterOfMassX(const Data *d, Grid *grid, int dir, int trc, int DoReduce, double* numerator, double* denom)
{
  int i, j, k;
  int var = NVAR;

  double LocalCenterOfMassX = 0.0;
  double LocalMass = 0.0;
  double GlobalCenterOfMassX = 0.0;
  double GlobalMass = 0.0;

//  print("x,y,z(0): %f,%f,%f\n",grid[IDIR].x[0],grid[JDIR].x[0],grid[KDIR].x[0]);

//  print("prank %d(%d): 1\n",prank,grid[KDIR].level);

  DOM_LOOP(k,j,i)
  {
    //print("k,j,i: %d,%d,%d\n",k,j,i);
    double dens = d->Vc[RHO][k][j][i];
    //print("dens: %f\n",dens);
    if (trc) dens *= d->Vc[trc][k][j][i];
    //print("trcdens: %f\n",dens);
    LocalMass += dens;// * dV; //For a uniform grid the volume term is constant cancels out when computing the global center of mass
    int idx = -10000000;
    switch (dir)
      {
      case IDIR:
	{
	  idx = i;
	  break;
	}
      case JDIR:
	{
	  idx = j;
	  break;
	}
      case KDIR:
	{
	  idx = k;
	  break;
	}
      }
    //int GlobalIndex = idx + (grid[dir].beg - grid[dir].nghost);
    //print("idx: %d\n",GlobalIndex);
    LocalCenterOfMassX += grid[dir].x[idx] * dens;//(grid[dir].x_glob[GlobalIndex] + grid[dir].dx[idx]/2.0) * dens;// * dV;
  }
//  print("prank %d(%d): 2\n",prank,grid[KDIR].level);

  if (DoReduce)
  {
    MPI_Allreduce(&LocalCenterOfMassX, &GlobalCenterOfMassX, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&LocalMass, &GlobalMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  print("prank %d(%d): 3\n",prank,grid[KDIR].level);
    GlobalCenterOfMassX /= GlobalMass;
    MPI_Barrier(MPI_COMM_WORLD);
    print("XCM(%d): %d %12.8e, GM %12.8e, LCM %12.8e, LM %12.8e\n",grid[KDIR].level,prank, GlobalCenterOfMassX,GlobalMass,LocalCenterOfMassX,LocalMass);
    return GlobalCenterOfMassX;
  }
  else return LocalCenterOfMassX/LocalMass;
}

#ifdef USE_TANH_DENSPROFILE
double TanhProfile(double CurRad, double CloudRad, double Width)
{
    return 1-tanh(1+(CurRad-(CloudRad+Width))/Width);
}

//Same profile as above but in terms of steepness rather than transition width
double TanhProfileSteepness(double CurRad, double CloudRad, double Steepness)
{
    return 1-tanh(Steepness*(CurRad/CloudRad-1));
}
#endif



/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rd dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  static int firstcall = 1;
  double CloudRad = g_inputParam[CLOUD_RADIUS];
  double CurRad;
  double TCloud = g_inputParam[TEMP_C];
  double nWind = g_inputParam[RHO_W];
  double nCloud = g_inputParam[RHO_C];
  #ifdef PRS_FROM_TEMP
  double PrsCloud = nCloud*TCloud/KELVIN;
  #else
  double PrsCloud = g_inputParam[PRS_C];
  #endif
  double ZCloud = g_inputParam[Z_C];
  double ZWind = g_inputParam[Z_W];
  double TWind = g_inputParam[TEMP_W];
  double PrsWind = PrsCloud;

  const double n0 = g_inputParam[N_0];

//  print1("1\n");

  g_gamma = g_inputParam[GAMMA_GAS];

  #if PHYSICS == MHD
  #ifdef BINIT_MICROGAUSS
  double B0 = g_inputParam[BETA_W]*UNIT_MICROGAUSS;
  #else
  double B0 = sqrt(2*PrsCloud/g_inputParam[BETA_W]);
  #endif
  #endif

  g_smallDensity   = 1.0000e-10;  // included to avoid negative densities at cloud's front
  g_smallPressure  = 1.0000e-07;  // included to avoid negative pressures at cloud's front


  #ifdef VAR_MU
  if (firstcall)
  {
    LoadMuTable();
  }
  #endif

  #ifdef HYDROSTATIC_HALO
//  if (firstcall)
//  {
//    int tableret = ReadTable("distdensbfield-sun10.txt");
//    if (tableret != 0) {print("Error reading table on proc %d: %d\n",prank,tableret);}
//  }
//  print1("2\n");

  #if DIMENSIONS == 2
  double dist = x2;
  #else
  double dist = x3;
  #endif
  //dist += g_inputParam[INIT_ZDIST];
  #if PHYSICS == MHD
  //  while (g_tabulardistance[g_curidx] < g_inputParam[INIT_ZDIST]) {g_curidx++;}
  //if (fabs(g_tabulardistance[g_curidx]-g_inputParam[INIT_ZDIST]) > fabs(g_tabulardistance[g_curidx-1]-g_inputParam[INIT_ZDIST])) g_curidx--;
  //B0 = g_tabularbfield[g_curidx] * UNIT_MICROGAUSS;
  double zprime = (dist - 1.5)/4.0;
  //Azimuthal field from Sun & Reich (2010) in cylindrical phi coordinate
  double Bphi = 2.0/(1+zprime*zprime);
  B0 = Bphi * UNIT_MICROGAUSS;
  #endif
  #ifdef VAR_MU
  double T = TWind;//PrsWind*KELVIN/Numdens;
//  print1("P, N, T: %f, %f, %f, %f, %f\n", PrsWind, Numdens, T, nWind, nCloud);
  double mu = GetMuFromTable(T);
  #else
  double mu = g_inputParam[MU_C];
  #endif

  nWind = n0*exp(-g_inputParam[GRAV_ACC]*UNIT_KPCPERMYR_TO_CMPERS2*mu*CONST_mp*dist*1000.0*CONST_pc/(CONST_kB*TWind));
//  print1("g, dist, TWind: %f, %f, %f, %f\n", g_inputParam[GRAV_ACC], dist, TWind, UNIT_KPCPERMYR_TO_CMPERS2);
  PrsWind = nWind*TWind/KELVIN;
  #endif

//  print1("3\n");


  double x = x1 - g_inputParam[CLOUD_CENTREX1];
  double y = x2 - g_inputParam[CLOUD_CENTREX2];
  double z = x3 - g_inputParam[CLOUD_CENTREX3];

  CurRad = EXPAND(x*x , + y*y , + z*z);
  CurRad = sqrt(CurRad);

  v[TRC_CLOUD] = 0;

  v[TRC_Z] = ZWind;

  #ifdef USE_TANH_DENSPROFILE
  #ifdef USE_HP09_PROFILE
  if ((CloudRad < 0.085) || (CloudRad > 0.095)) {print1("!HP09 profile requires CloudRad to be approx. 0.09 kpc. Quitting.\n");QUIT_PLUTO(1);}
  double Numdens = 0.0;
  if (CurRad < 0.1) Numdens = exp(-CurRad/0.05)*(nWind + (nCloud - nWind)*0.5*TanhProfileSteepness(CurRad,0.085,7.0));
  else Numdens = nWind + (nCloud - nWind)*0.5*TanhProfileSteepness(CurRad,0.0752,7.0);
  #else
  double Numdens = nWind + (nCloud - nWind)*0.5*TanhProfile(CurRad,CloudRad,CloudRad/g_inputParam[DENSITY_STEEPNESS]);
  #endif
  double ZRad = (CloudRad/g_inputParam[DENSITY_STEEPNESS])*atanh((nCloud-3*nWind)/(nCloud-nWind)) + CloudRad;  
  #else
  double Numdens = nWind + (nCloud - nWind)/(1+pow(CurRad/CloudRad,g_inputParam[DENSITY_STEEPNESS]));
  double ZRad = CloudRad*pow(nCloud/nWind-2,1.0/g_inputParam[DENSITY_STEEPNESS]);
  #endif

  //Use radius where the density is twice the halo density as cloud-wind metallicity boundary
  if (CurRad < ZRad) {v[TRC_Z] = ZCloud;} else {v[TRC_Z] = ZWind;}

//  print1("4\n");

//  print1("4b\n");

  v[PRS] = PrsWind;

//  printf("n, prs, mu: %f, %f, %f\n",Numdens, PrsCloud, mu);
  T = KELVIN*PrsWind/Numdens;
  mu = GetMuFromTable(T);
  v[RHO] = Numdens*mu;

  if (CurRad < CloudRad*(1.0 + 3.0/g_inputParam[DENSITY_STEEPNESS]))
  {
//    v[RHO] = RhoCloud;//RhoWind + (RhoCloud - RhoWind)/(1+pow(CurRad/CoreRad,g_inputParam[DENSITY_STEEPNESS]));

    v[VX1] = 0.0;
    #if DIMENSIONS == 2
    v[VX2] = g_inputParam[VX1_C];
    v[VX3] = 0.0;
    #else
    v[VX2] = 0.0;
    v[VX3] = g_inputParam[VX1_C];
    #endif

    if (CurRad < CloudRad) v[TRC_CLOUD] = 1.0;
  }
  else
  {
    #ifdef VELOCITY_MACH
    double VelWind = g_inputParam[VX1_W] * sqrt(g_gamma*PrsCloud/v[RHO]);
    #else
    double VelWind = g_inputParam[VX1_W];
    #endif
  //  v[RHO] = RhoWind;
    v[VX1] = 0.0;
    #if	DIMENSIONS == 2
    v[VX2] = VelWind;
    v[VX3] = 0.0; 
    #else
    v[VX2] = 0.0;
    v[VX3] = VelWind;
    #endif
//    v[TRC_Z] = ZWind;
  }

//  print1("5\n");

  #if PHYSICS == MHD
  double Bx = 0;
  double By = 0;
  double Bz = 0;
  double Ax = 0;
  double Ay = 0;
  double Az = 0;

#if BG_FIELD == FIELD_TRANSVERSE /* Background field is B=(0, B0, 0) */
  By = B0;
  Az = -x*B0;
#elif BG_FIELD == FIELD_PARALLEL /* Background field is B=(0, 0, B0) */
  Bz = B0;
//  Az = y*B0;
#elif BG_FIELD == FIELD_OBLIQUE  /* Background field is B=(B0/sqrt(3), B0/sqrt(3), B0/sqrt(3)) */
  Bx = B0/sqrt(3);
  By = B0/sqrt(3);
  Bz = B0/sqrt(3);
  Ay = B0*x/sqrt(3);
  Az = B0*(y-x)/sqrt(3);
#endif

  v[BX1] = Bx;
  v[BX2] = By;
  v[BX3] = Bz;

  v[AX1] = Ax;
  v[AX2] = Ay;
  v[AX3] = Az;
#endif // PHYSICS == MHD


  if (firstcall) firstcall = 0;
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
#ifdef OUTP_BOV
  double dt = 0.1;
  int filenumber = g_time/dt;
  if (prank == 0)
  {
    char* dir;
    char filename[512];
    FILE* BOVout;
    dir = GetOutputDir();
    Output vnames;
    Output* varnames = &vnames;
    varnames->var_name = ARRAY_2D(64,128,char);
    SetDefaultVarNames(varnames);

    int nv = NVAR;
    #ifdef STAGGERED_MHD
     D_EXPAND(
       varnames->var_name[nv]   = "bx1s";
       nv++; ,

       varnames->var_name[nv]   = "bx2s"; 
       nv++; ,

       varnames->var_name[nv]   = "bx3s"; 
       nv++;
     )
    #endif
    #if UPDATE_VECTOR_POTENTIAL == YES
     #if DIMENSIONS == 3   
      varnames->var_name[nv]   = "Ax1";
      nv++;
      varnames->var_name[nv]   = "Ax2";
      nv++;
     #endif
     varnames->var_name[nv]   = "Ax3";
     nv++;
    #endif

    int i;
    for (i = 0; i < nv; i++)
    {
      sprintf(filename, "%s/%s.%04d%s", dir, varnames->var_name[i], filenumber, ".bov");
      BOVout = fopen(filename, "w");
      fprintf(BOVout, "TIME: %f\nDATA_FILE: %s.%04d%s\nDATA_SIZE: %d %d %d\nDATA_FORMAT: DOUBLE\n"
	      "VARIABLE: %s\nDATA_ENDIAN: LITTLE\nCENTERING: zonal\n"
	      "BRICK_ORIGIN: %f %f %f\nBRICK_SIZE: %f %f %f",
	      g_time, varnames->var_name[i], filenumber, ".dbl", grid[IDIR].np_int_glob, grid[JDIR].np_int_glob, grid[KDIR].np_int_glob,
	      varnames->var_name[i], g_inputParam[CLOUD_CENTREX1], g_inputParam[CLOUD_CENTREX2], g_inputParam[CLOUD_CENTREX3],
	      g_domEnd[IDIR]-g_domBeg[IDIR], g_domEnd[JDIR]-g_domBeg[JDIR], g_domEnd[KDIR]-g_domBeg[KDIR]);
      fclose(BOVout);
    }

//    filenumber++;
  }
#endif //OUTP_BOV


#ifdef OUTP_CENTER_OF_MASS
  int i, j, k;
    int var = NVAR;
    //  print1("var %f", d->Vc[TRC][middlez][middley][middlex+1]);

    double LocalCenterOfMassX = 0.0;
    double LocalMass = 0.0;
    double GlobalCenterOfMassX = 0.0;
    double GlobalMass = 0.0;
    //double dV = grid[IDIR].dx[0] * grid[JDIR].dx[0] * grid[KDIR].dx[0];

    //    print1("dV %f", grid[IDIR].dx[0] * UNIT_LENGTH);
    //print("rank %i, gbeg %i, gend %i ,ghost %i, globbeg %f globend %f ",prank, grid[IDIR].beg,grid[IDIR].end,grid[IDIR].nghost,grid[IDIR].x_glob[grid[IDIR].beg],grid[IDIR].x_glob[grid[IDIR].end]);
    //    print("rank %i, x1 %f, x2 %f ",prank,(grid[IDIR].x_glob[(grid[IDIR].beg - grid[IDIR].nghost)] + grid[IDIR].dx[0]/2.0), (grid[IDIR].x_glob[IEND + (grid[IDIR].beg - grid[IDIR].nghost)] + grid[IDIR].dx[IEND]/2.0));

    TOT_LOOP(k,j,i)
    {
      LocalMass += d->Vc[RHO][k][j][i];// * dV; //For a uniform grid the volume term is constant cancels out when computing the global center of mass
      int GlobalIndex = i + (grid[IDIR].beg - grid[IDIR].nghost);
      LocalCenterOfMassX += (grid[IDIR].x_glob[GlobalIndex] + grid[IDIR].dx[i]/2.0) * d->Vc[RHO][k][j][i];// * dV;
    }

    /* TODO: MPI STUFF SO THAT THIS WORKS IN PARALLEL
             Hopefully this won't slow the processes down too much. If it does, maybe find a way to only do this ocassionaly rather than at every timestep*/\

    MPI_Allreduce(&LocalCenterOfMassX, &GlobalCenterOfMassX, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&LocalMass, &GlobalMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    GlobalCenterOfMassX /= GlobalMass;

    double maxrho = -1.0;
    double maxrhoX = 0.0;
    TOT_LOOP(k, j, i)
    {
      if (d->Vc[RHO][k][j][i] > maxrho)
      {
	int GlobalIndex = i + (grid[IDIR].beg - grid[IDIR].nghost);
	maxrho = d->Vc[RHO][k][j][i];
	maxrhoX = grid[IDIR].x_glob[GlobalIndex];
      }
    }

    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    double *GlobalmaxrhoX = (double *)malloc(sizeof(double) * numprocs);
    double *Globalmaxrho = (double *)malloc(sizeof(double) * numprocs);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(&maxrhoX, 1, MPI_DOUBLE, GlobalmaxrhoX, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&maxrho, 1, MPI_DOUBLE, Globalmaxrho, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    int maxelem = -1.0;
    int elem;
    for (elem=0; elem < numprocs; elem++)
    {
      print("Glob %f ", Globalmaxrho[elem]);
	 if (Globalmaxrho[elem] > Globalmaxrho[maxelem])
	 {
	    maxelem = elem;
	 }
    }

    maxrhoX = GlobalmaxrhoX[maxelem];
    
    free(GlobalmaxrhoX);
    free(Globalmaxrho);

    if (prank == 0)
      {
	char *dir, fname[512];
	FILE *fp;
	dir = GetOutputDir();
	sprintf(fname, "%s/centerofmassx.txt", dir);
	fp = fopen(fname, "a");
	fprintf(fp, "%12.6e  %12.6e %12.6e\n", g_time, GlobalCenterOfMassX, maxrhoX);
	fclose(fp);
      }
#endif //OUTP_CENTER_OF_MASS

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int i, j, k;
  int nv;
  static int firstcall = 1;
  static int GridCanMove = 1;
  static double TotalNum0 = 0.0;
  static double TotalDenom0 = 0.0;
  static double TotalNum1 = 0.0;
  static double TotalDenom1 = 0.0;
  static double TotalNum2 = 0.0;
  static double TotalDenom2 = 0.0;
  static double TotalNum3 = 0.0;
  static double TotalDenom3 = 0.0;
  static double TotalNum4 = 0.0;
  static double TotalDenom4 = 0.0;
  static double CoM_prev = -1.0;
  static double vel_prev = 0.0;
  static double s_RefFrameVelocity = 1.e10;
  static double t_prev = 1.e30;
  static double s_dist = 0.0;
  static int maxlevel = 0;
  static int level_prev = 0;
  double VelThresh = g_inputParam[VEL_THRESH];

//  print("rank: %d, level: %d\n",prank,grid[KDIR].level);

//  print("rho: %12.8e, prs: %12.8e, vx1: %12.8e, vx2 %12.8e, vx3: %12.8e\n",d->Vc[RHO][10][10][10], d->Vc[PRS][10][10][10], d->Vc[VX1][10][10][10], d->Vc[VX2][10][10][10], d->Vc[VX3][10][10][10]);
//  print("B1\n");
  if ((firstcall) && (g_time > 0.0) && (0==1))
  {
//    print("B1b\n");
    FILE *file;
    if (file = fopen("framevel.dat", "r"))
    {
      //Read reference frame velocity from file
      double time = 0.0;
      double ptime = 0.0;
      double vel = 0.0;
      double pvel = 0.0;
      double dist = 0.0;
      while (fscanf(file, "%lf %lf %lf", &time, &vel, &dist) != EOF)
      {
	if (time > g_time)
	{
	  s_dist += pvel*(g_time - ptime);
	  break;
	}
	s_RefFrameVelocity = vel;
	s_dist += pvel*(time - ptime);
	ptime = time;
	pvel = vel;
      }
      if (s_RefFrameVelocity == 1.e10) s_RefFrameVelocity = 0.0;
      fclose(file);
      print1("time: %f, ftime: %f, vel: %f, dist: %f\n",g_time, time, s_RefFrameVelocity, s_dist);
    }
    else s_RefFrameVelocity = 0.0;
    t_prev = g_time;
    firstcall = 0;
  }
  static int cstep = 0;


  const double n0 = g_inputParam[N_0];
  double TWind = g_inputParam[TEMP_W];

  if ((side == 0) && (!firstcall))    /* -- check solution inside domain -- */
  {
    TOT_LOOP(k,j,i){
      if (d->Vc[RHO][k][j][i] < g_smallDensity) {d->Vc[RHO][k][j][i] = g_smallDensity;}
      if (d->Vc[PRS][k][j][i] < g_smallPressure) {d->Vc[PRS][k][j][i] = g_smallPressure;}
      if (d->Vc[TRC][k][j][i] < 0.0) {d->Vc[TRC][k][j][i] = 0.0;}     
      if (d->Vc[TRC_Z][k][j][i] < 0.0) {d->Vc[TRC_Z][k][j][i] = 0.0;}
      /*      if (grid[IDIR].x[i] < -1.8)
      {
	d->Vc[BX1][k][j][i] = sqrt(2*PrsWind/(3*g_inputParam[BETA_W]));
	d->Vc[BX2][k][j][i] = sqrt(2*PrsWind/(3*g_inputParam[BETA_W]));
	d->Vc[BX3][k][j][i] = sqrt(2*PrsWind/(3*g_inputParam[BETA_W]));
	d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
      }*/
    }
  }

  #if DIMENSIONS == 2
  if (((side == X2_BEG) || (side == X2_END)))// && (!firstcall))
  #else
    if (((side == X3_BEG) || (side == X3_END)))// && (!firstcall))
  #endif
  {
    if (box->vpos == CENTER)
    {
//      double mu = GetMuFromTable(TWind);

      double a, b, c;
      double rho = 0.0;
      double prs = 0.0;
      double vx1 = 0.0;
      double vx2 = 0.0;
      double vx3 = 0.0;
      double T = 0.0;
      double mu = 0.0;
      double zbound = 0.0;
      int count, maxiter;
      double tolerance;
      count = 0;
      tolerance = 0.05;

      int nz = grid[KDIR].np_tot;
      int nghost = grid[KDIR].nghost;

      if (side == X3_BEG)
      {
	rho = d->Vc[RHO][nghost][0][0];
	prs = d->Vc[PRS][nghost][0][0];
	zbound = grid[KDIR].x[nghost];
	vx1 = d->Vc[VX1][nghost][0][0];
	vx2 = d->Vc[VX2][nghost][0][0];
	vx3 = d->Vc[VX3][nghost][0][0];
      }
      else
      {
	rho = d->Vc[RHO][nz-nghost-1][0][0];
	prs = d->Vc[PRS][nz-nghost-1][0][0];
	zbound = grid[KDIR].x[nz-nghost-1];
	vx1 = d->Vc[VX1][nz-nghost-1][0][0];
	vx2 = d->Vc[VX2][nz-nghost-1][0][0];
	vx3 = d->Vc[VX3][nz-nghost-1][0][0];
      }
      //Find temperature and mu as root of mu(T)*P/(k_B*rho) - T=0
      a = prs*((double)g_mu_tab[0]+tolerance)*KELVIN/rho;
      b = prs*((double)g_mu_tab[g_ntabmu-1]-tolerance)*KELVIN/rho;

      struct temperature_params par;
      par.mu = -1;
      par.rho = rho;
      par.prs = prs;
      int status = Brent(fT, &par, a, b, -1, 1.e-12, &T);
      if (status != 0)
      {
	print("! z boundary: Failed to find root of temperature function on proc %d\nStatus: %d\n", prank, status);
	QUIT_PLUTO(1);
      }
      mu = par.mu;

      if (T != T){
	printf (" ! Nan found in z boundary \n");
	printf (" ! rho = %12.6e, prs = %12.6e\n",rho, prs);
	QUIT_PLUTO(1);
      }


      BOX_LOOP(box, k, j, i)
      {
	double c_dist = grid[KDIR].CoM;
#if DIMENSIONS == 2
	double dist = grid[JDIR].x[j] + c_dist;//s_dist;
#else
	double dist = grid[KDIR].x[k] + c_dist;//s_dist;
	double delta_dist = zbound - grid[KDIR].x[k];
#endif
//	double nWind = n0*exp(-g_inputParam[GRAV_ACC]*UNIT_KPCPERMYR_TO_CMPERS2*mu*CONST_mp*dist*1000.0*CONST_pc/(CONST_kB*TWind));

	double nWind = (rho/mu)*exp(g_inputParam[GRAV_ACC]*UNIT_KPCPERMYR_TO_CMPERS2*mu*CONST_mp*delta_dist*1000.0*CONST_pc/(CONST_kB*T)); //Extrapolate hydrostatic density profile
	//print1("i, j, k, nwind, rhowind, rho, mu, T, dist, delta_dist, s_dist, zbound: %d, %d, %d, %lf, %lf, %lf, %f, %f, %f, %f, %f, %f\n",i,j,k,nWind,nWind*mu,rho,mu,T,dist,delta_dist,s_dist,zbound);


	d->Vc[RHO][k][j][i] = mu*nWind;
	d->Vc[PRS][k][j][i] = nWind*T/KELVIN; //nWind*Twind/KELVIN
	d->Vc[VX1][k][j][i] = vx1;
	d->Vc[VX2][k][j][i] = vx2;
#if DIMENSIONS == 3
	d->Vc[VX3][k][j][i] = vx3;
#endif
	d->Vc[TRC][k][j][i] = 0.0;
	d->Vc[TRC_Z][k][j][i] = g_inputParam[Z_W];
#if PHYSICS == MHD
	double zprime = (dist - 1.5)/4.0;
	//Azimuthal field from Sun & Reich (2010) in cylindrical phi coordinate
	double Bphi = 2.0/(1+zprime*zprime);
	double B0 = Bphi * UNIT_MICROGAUSS;
	#if BG_FIELD == FIELD_TRANSVERSE
	d->Vc[BX1][k][j][i] = 0.0;
	d->Vc[BX2][k][j][i] = B0;
        #if DIMENSIONS == 3
	d->Vc[BX3][k][j][i] = 0.0;
        #endif
        #endif //FIELD_TRANSVERSE
	#if BG_FIELD == FIELD_PARALLEL
	d->Vc[BX1][k][j][i] = 0.0;
	d->Vc[BX2][k][j][i] = 0.0;
        #if DIMENSIONS == 2
	d->Vc[BX2][k][j][i] = B0;
        #else
	d->Vc[BX3][k][j][i] = B0;
        #endif
        #endif //FIELD_PARALLEL
#endif //PHYSICS == MHD
      }
    }
    else if (box->vpos == X1FACE)
    {
//  print("B6\n");

      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i)
      {
        #if (DIMENSIONS == 2) && (BG_FIELD == FIELD_TRANSVERSE)
	double dist = grid[KDIR].x[k] + grid[KDIR].CoM;//s_dist;
	double zprime = (dist - 1.5)/4.0;
	//Azimuthal field from Sun & Reich (2010) in cylindrical phi coordinate
	double Bphi = 2.0/(1+zprime*zprime);
	double B0 = Bphi * UNIT_MICROGAUSS;
	d->Vs[BX1s][k][j][i] = B0;
	#else
	d->Vs[BX1s][k][j][i] = 0.0;
	#endif

      }
      #endif
    }
    else if (box->vpos == X2FACE)
    {

#if defined(STAGGERED_MHD) && (DIMENSIONS == 3) && (BG_FIELD == FIELD_TRANSVERSE)
      BOX_LOOP(box,k,j,i)
      {
	double dist = grid[KDIR].x[k] + grid[KDIR].CoM;//s_dist;
	double zprime = (dist - 1.5)/4.0;
	//Azimuthal field from Sun & Reich (2010) in cylindrical phi coordinate
	double Bphi = 2.0/(1+zprime*zprime);
	double B0 = Bphi * UNIT_MICROGAUSS;
	d->Vs[BX2s][k][j][i] = B0;
      }
#endif //STAGGERED_MHD
    }
    else if (box->vpos == X3FACE) {}
  }
//print("B8\n");
//  cstep++;
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  #if DIMENSIONS == 2
  g[IDIR] = 0.0;
  g[JDIR] = -g_inputParam[GRAV_ACC]*UNIT_KPCPERMYR2;
  g[KDIR] = 0.0;
  #else
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = -g_inputParam[GRAV_ACC]*UNIT_KPCPERMYR2;
  #endif
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif

