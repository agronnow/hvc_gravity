#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>

#include "LevelPluto.H"
#include "BoxIterator.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "LoHiSide.H"

#include "CH_Timer.H"

#include "NamespaceHeader.H"

#define RHO_CLOUD TRC

// Constructor - set up some defaults
LevelPluto::LevelPluto()
{
  m_RefFrameVelocity = 0.0;
  m_dist = 0.0;
  m_vel = 0.0;
  m_VelUpdated = false;
  m_dx           = 0.0;
  m_dl_min       = 1.e30;
  m_refineCoarse = 0;
  m_patchPluto   = NULL;
  m_isDefined    = false;
}

// Destructor - free up storage
LevelPluto::~LevelPluto()
{
  if (m_patchPluto != NULL)
    {
      delete m_patchPluto;
    }
}

// Define the object so that time stepping can begin
void LevelPluto::define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                        const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                        const ProblemDomain&      a_domain,
                        const int&                a_refineCoarse,
                        const int&                a_level,
                        const Real&               a_dx,
                        const PatchPluto*         a_patchPlutoFactory,
                        const bool&               a_hasCoarser,
                        const bool&               a_hasFiner)
{
  CH_TIME("LevelPluto::define");

  // Sanity checks
  CH_assert(a_refineCoarse > 0);
  CH_assert(a_dx > 0.0);

  // Make a copy of the current grids
  m_grids  = a_thisDisjointBoxLayout;

  // Cache data
  m_dx = a_dx;
  m_level = a_level;
  m_domain = a_domain;
  m_refineCoarse = a_refineCoarse;
  m_hasCoarser = a_hasCoarser;
  m_hasFiner = a_hasFiner;

  // Remove old patch integrator (if any), create a new one, and initialize
  if (m_patchPluto != NULL)
    {
      delete m_patchPluto;
    }

 // Determing the number of ghost cells necessary here
  m_numGhost = GetNghost(NULL);

  m_patchPluto = a_patchPlutoFactory->new_patchPluto();
  m_patchPluto->define(m_domain,m_dx,m_level,m_numGhost);
 
  // Set the grid for the entire level
  setGridLevel();

  // Get the number of conserved variable and face centered fluxes
  m_numCons   = m_patchPluto->numConserved();
  m_numFluxes = m_patchPluto->numFluxes();

  m_exchangeCopier.exchangeDefine(a_thisDisjointBoxLayout,
                                  m_numGhost*IntVect::Unit);

  // Setup an interval corresponding to the conserved variables
  Interval UInterval(0,m_numCons-1);

 #if (TIME_STEPPING == RK2)
  // Create temporary storage with a layer of "m_numGhost" ghost cells
  // for the flags passing from predictor to corrector (RK2 only)
  m_Flags.define(m_grids,1,m_numGhost*IntVect::Unit);
  m_Utmp.define(m_grids,m_numCons,m_numGhost*IntVect::Unit);
 #endif

  // Create temporary storage with a layer of "m_numGhost" ghost cells
  {
    CH_TIME("setup::Udefine");
    m_U.define(m_grids,m_numCons,m_numGhost*IntVect::Unit);
  }

  // Set up the interpolator if there is a coarser level
  if (m_hasCoarser)
    {
      m_patcher.define(a_thisDisjointBoxLayout,
                       a_coarserDisjointBoxLayout,
                       m_numCons,
                       coarsen(a_domain,a_refineCoarse),
                       a_refineCoarse,
                       m_dx,
                       m_numGhost);
    }

  // Everything is defined
  m_isDefined = true;
}

// Advance the solution by "a_dt" by using an unsplit method.
// "a_finerFluxRegister" is the flux register with the next finer level.
// "a_coarseFluxRegister" is flux register with the next coarser level.
// If source terms do not exist, "a_S" should be null constructed and not
// defined (i.e. its define() should not be called).
Real LevelPluto::step(LevelData<FArrayBox>&       a_U,
                      LevelData<FArrayBox>        a_flux[CH_SPACEDIM],
                      LevelFluxRegister&          a_finerFluxRegister,
                      LevelFluxRegister&          a_coarserFluxRegister,
                      LevelData<FArrayBox>&       a_split_tags,
                      const LevelData<FArrayBox>& a_UCoarseOld,
                      const Real&                 a_TCoarseOld,
                      const LevelData<FArrayBox>& a_UCoarseNew,
                      const Real&                 a_TCoarseNew,
                      const Real&                 a_time,
                      const Real&                 a_dt,
                      const Real&                 a_cfl)
{
  CH_TIMERS("LevelPluto::step");

  CH_TIMER("LevelPluto::step::setup"   ,timeSetup);
  CH_TIMER("LevelPluto::step::update"  ,timeUpdate);
  CH_TIMER("LevelPluto::step::reflux"  ,timeReflux);
  CH_TIMER("LevelPluto::step::conclude",timeConclude);

  // Make sure everything is defined
  CH_assert(m_isDefined);

  CH_START(timeSetup);

  // Clear flux registers with next finer level
  if (m_hasFiner && (g_intStage == 1))
    {
      a_finerFluxRegister.setToZero();
    }

  // Setup an interval corresponding to the conserved variables
  Interval UInterval(0,m_numCons-1);

  {
    CH_TIME("setup::localU");
    for (DataIterator dit = m_U.dataIterator(); dit.ok(); ++dit)
    {
      m_U[dit].setVal(0.0); // Gets rid of denormalized crap.
      m_U[dit].copy(a_U[dit]);
    }

    m_U.exchange(m_exchangeCopier);
  }

  // Fill m_U's ghost cells using fillInterp
  if (m_hasCoarser)
    {
      // Fraction "a_time" falls between the old and the new coarse times
      Real alpha = (a_time - a_TCoarseOld) / (a_TCoarseNew - a_TCoarseOld);

    // Truncate the fraction to the range [0,1] to remove floating-point
    // subtraction roundoff effects
      Real eps = 0.04 * a_dt / m_refineCoarse;

      if (Abs(alpha) < eps)     alpha = 0.0;
      if (Abs(1.0-alpha) < eps) alpha = 1.0;

      // Current time before old coarse time
      if (alpha < 0.0)
        {
          MayDay::Error( "LevelPluto::step: alpha < 0.0");
        }

      // Current time after new coarse time
      if (alpha > 1.0)
        {
          MayDay::Error( "LevelPluto::step: alpha > 1.0");
        }

      // Interpolate ghost cells from next coarser level using both space
      // and time interpolation
      m_patcher.fillInterp(m_U,
                           a_UCoarseOld,
                           a_UCoarseNew,
                           alpha,
                           0,0,m_numCons);
    }

  // Potentially used in boundary conditions
  m_patchPluto->setCurrentTime(a_time);

  // Use to restrict maximum wave speed away from zero
  Real maxWaveSpeed = 1.e-12;
  Real minDtCool    = 1.e38;

  int lev = -1;

  // The grid structure
  Grid *grid;
  static Time_Step Dts;
  Real inv_dt;
  
  #ifdef GLM_MHD
   glm_ch = g_coeff_dl_min*m_dx/(a_dt + 1.e-16)*a_cfl;
//   glm_ch = g_coeff_dl_min/(a_dt + 1.e-16)*a_cfl; /* If subcycling is turned off */
   glm_ch = MIN(glm_ch,glm_ch_max*g_coeff_dl_min);
  #endif

  CH_STOP(timeSetup);
  g_level_dx = m_dx;


  // Beginning of loop through patches/grids.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    CH_START(timeUpdate);

    // The current box
    Box curBox = m_grids.get(dit());

    // The current grid of conserved variables
    FArrayBox& curU = m_U[dit];

    // The current grid of volumes
    #if GEOMETRY != CARTESIAN
     const FArrayBox& curdV = m_dV[dit()];
    #else
     const FArrayBox  curdV;
    #endif
 
   #ifdef SKIP_SPLIT_CELLS
    // The current grid of split/unsplit tags
    FArrayBox& split_tags = a_split_tags[dit];
   #else
    FArrayBox split_tags;
   #endif

   #if (TIME_STEPPING == RK2)
    // The current storage for flags (RK2 only)
    BaseFab<unsigned char>& flags = m_Flags[dit];
    // Local temporary storage for conserved variables
    FArrayBox& curUtmp = m_Utmp[dit];
   #else
    BaseFab<unsigned char> flags;
    FArrayBox curUtmp;
   #endif

    // The fluxes computed for this grid - used for refluxing and returning
    // other face centered quantities
    FluxBox flux;

    // Set the current box for the patch integrator
    m_patchPluto->setCurrentBox(curBox);

    Real minDtCoolGrid;
    
    grid = m_structs_grid[dit].getGrid();
 
    IBEG = grid[IDIR].lbeg; IEND = grid[IDIR].lend;
    JBEG = grid[JDIR].lbeg; JEND = grid[JDIR].lend;
    KBEG = grid[KDIR].lbeg; KEND = grid[KDIR].lend;

    NX1 = grid[IDIR].np_int;
    NX2 = grid[JDIR].np_int;
    NX3 = grid[KDIR].np_int;

    NX1_TOT = grid[IDIR].np_tot;
    NX2_TOT = grid[JDIR].np_tot;
    NX3_TOT = grid[KDIR].np_tot;
 
    g_dt   = a_dt;
    g_time = a_time;
    g_maxRiemannIter = 0;
    PLM_CoefficientsSet (grid);  /* -- these may be needed by
                                       shock flattening algorithms */
    #if INTERPOLATION == PARABOLIC
     PPM_CoefficientsSet (grid);  
    #endif
    
    // reset time step coefficients 
    if (Dts.cmax == NULL) Dts.cmax = ARRAY_1D(NMAX_POINT, double);
    int id;
    Dts.inv_dta = 1.e-18;
    Dts.inv_dtp = 1.e-18;
    Dts.dt_cool = 1.e18;
    Dts.cfl     = a_cfl;
    Where(-1, grid); /* -- store grid for subsequent calls -- */

    grid[KDIR].CoM = m_dist;//GetDist(g_time);

    if (!m_VelUpdated)
    {
    double*** UU[NVAR];
    int k,j,i,nv;
    for (nv = 0; nv < NVAR; nv++)
    {
      UU[nv] = ArrayMap(NX3_TOT, NX2_TOT, NX1_TOT, curU.dataPtr(nv));
    }
    TOT_LOOP(k,j,i)
    {
      UU[MX3][k][j][i] -= UU[RHO][k][j][i]*m_vel;
    }
    for (nv = 0; nv < NVAR; nv++) FreeArrayMap(UU[nv]);
    }

    // updateSolution!!!!! Finally...
    m_patchPluto->updateSolution(curU, curUtmp, curdV, split_tags, flags, flux,
                                 &Dts, curBox, grid);

    inv_dt = Dts.inv_dta + 2.0*Dts.inv_dtp;
    maxWaveSpeed = Max(maxWaveSpeed, inv_dt); // Now the inverse of the timestep

    minDtCool = Min(minDtCool, Dts.dt_cool/a_cfl);

    CH_STOP(timeUpdate);

    CH_START(timeReflux);

    // Do flux register updates
    for (int idir = 0; idir < SpaceDim; idir++) {
    // Increment coarse flux register between this level and the next
    // finer level - this level is the next coarser level with respect
    // to the next finer level
      if (m_hasFiner) {
        a_finerFluxRegister.incrementCoarse(flux[idir],a_dt,dit(),
                                            UInterval, UInterval,idir);
/*
        const Vector<Box>& intersectlo =
	  a_finerFluxRegister.getCoarseLocations(idir,Side::Lo)[dit()];
	const Vector<Box>& intersecthi =
	  a_finerFluxRegister.getCoarseLocations(idir,Side::Hi)[dit()];
	intersects += intersectlo.size() + intersecthi.size();
*/
      }

      // Increment fine flux registers between this level and the next
      // coarser level - this level is the next finer level with respect
      // to the next coarser level
       if (m_hasCoarser) {
         a_coarserFluxRegister.incrementFine(flux[idir],a_dt,dit(),
                                             UInterval, UInterval,idir);
       }
    }

    CH_STOP(timeReflux);
  }

  m_VelUpdated = true;
/*  #ifdef COMOVING_FRAME
//  pout() << "t, Level, intersects: " << g_time << ", " << lev << ", " << intersects << endl;

  //Add the rho*z and cloud masses of all patches on each processor together to find the global center of mass for the cloud
  double GlobalRhoz = 0.0;
  double GlobalMass = 0.0;
  LevelRhoz[lev] += TotRhoz;
  LevelMass[lev] += TotMass;
  if (intersects == 0) //Finest level with existing patches in the domain of this processor
  {
    pout() << "!!!" << endl;
    double MassAllLevels = 0.0;
    double RhozAllLevels = 0.0;
    for (int l = 0 ; l <= lev; l++)
    {
      RhozAllLevels += LevelRhoz[l];
      MassAllLevels += LevelMass[l];
      LevelRhoz[l] = 0.0;
      LevelMass[l] = 0.0;
    }
*/
/*    MPI_Barrier(Chombo_MPI::comm);
    int red = MPI_Allreduce(&RhozAllLevels, &GlobalRhoz, 1, MPI_DOUBLE, MPI_SUM, Chombo_MPI::comm);
    if (red != MPI_SUCCESS) MayDay::Error("sorry, but I had a communcation error on calculating integral rho*z");
    red = MPI_Allreduce(&MassAllLevels, &GlobalMass, 1, MPI_DOUBLE, MPI_SUM, Chombo_MPI::comm);
    if (red != MPI_SUCCESS) MayDay::Error("sorry, but I had a communcation error on calculating total cloud mass");*/
/*    GlobalMass = 1.0;
    pout () << "Global: Num: " << GlobalRhoz << " GlobalMass: " << GlobalMass << " CoM: " << GlobalRhoz/GlobalMass << endl;
    CoM_prev = CenterOfMassZ;
    CenterOfMassZ = GlobalRhoz/GlobalMass;
    UpdateVelocities = true;
  }
  if ((CoM_prev != -1.0) && (g_time > 0.0) && (lev == 0) && (UpdateVelocities))
  {
    s_dist += s_RefFrameVelocity*(g_time - t_prev);
    double CoM_vel = (CenterOfMassZ - CoM_prev)/(g_time - t_prev);
    pout().precision(12);
    pout() << std::scientific;
    pout() << "t, dt, CoM, CoM_prev, CoM vel: " << g_time << ", " << g_dt << ", " << CenterOfMassZ << ", " << CoM_prev << ", " << CoM_vel << endl;
    int sgn = DSIGN(CoM_vel);
    t_prev = g_time;
  }
//  }
  #endif //COMOVING_FRAME
*/
  CH_START(timeConclude);

  {
    CH_TIME("conclude::copyU");
    // Now that we have completed the updates of all the patches, we copy the
    // contents of temporary storage, U, into the permanent storage, a_U.
    for(DataIterator dit = m_U.dataIterator(); dit.ok(); ++dit){
      a_U[dit].copy(m_U[dit]);
    }
   }

  // Find the minimum of dt's over this level
  Real local_dtNew = 1. / maxWaveSpeed;
  local_dtNew = Min(local_dtNew,minDtCool);
  Real dtNew;

  {
    CH_TIME("conclude::getDt");
 #ifdef CH_MPI
  #if (TIME_STEPPING == RK2) && (COOLING == NO)
  if (g_intStage == 1) {
  #endif
   int result = MPI_Allreduce(&local_dtNew, &dtNew, 1, MPI_CH_REAL,
                                  MPI_MIN, Chombo_MPI::comm);
   if(result != MPI_SUCCESS){ //bark!!!
      MayDay::Error("sorry, but I had a communcation error on new dt");
   }
  #if (TIME_STEPPING == RK2) && (COOLING == NO)
  } else {
   dtNew = local_dtNew;
  }
  #endif
 #else
   dtNew = local_dtNew;
 #endif
  }

  CH_STOP(timeConclude);

  // Return the maximum stable time step
  return dtNew;
}


void LevelPluto::GetCenterOfMassVel(double& Rhov, double& Mass, LevelData<FArrayBox>& a_split_tags)
{
  double TotVel = 0.0;
  double TotMass = 0.0;

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    int k, j, i, nv;
    double*** UU[NVAR];
    FArrayBox& curU = m_U[dit];
    Grid* grid = m_structs_grid[dit].getGrid();
 
    IBEG = grid[IDIR].lbeg; IEND = grid[IDIR].lend;
    JBEG = grid[JDIR].lbeg; JEND = grid[JDIR].lend;
    KBEG = grid[KDIR].lbeg; KEND = grid[KDIR].lend;

    NX1 = grid[IDIR].np_int;
    NX2 = grid[JDIR].np_int;
    NX3 = grid[KDIR].np_int;

    NX1_TOT = grid[IDIR].np_tot;
    NX2_TOT = grid[JDIR].np_tot;
    NX3_TOT = grid[KDIR].np_tot;

    for (nv = 0; nv < NVAR; nv++)
    {
      UU[nv] = ArrayMap(NX3_TOT, NX2_TOT, NX1_TOT, curU.dataPtr(nv));
    }
/*    RBox Ubox;
    D_EXPAND(Ubox.ib = curU.loVect()[IDIR]; Ubox.ie = curU.hiVect()[IDIR]; ,
	     Ubox.jb = curU.loVect()[JDIR]; Ubox.je = curU.hiVect()[JDIR]; ,
	     Ubox.kb = curU.loVect()[KDIR]; Ubox.ke = curU.hiVect()[KDIR]; );*/
	  
    FArrayBox& split_tags = a_split_tags[dit];
    double ***splitcells = ArrayBoxMap(KBEG, KEND, JBEG, JEND, IBEG, IEND, split_tags.dataPtr(0));
    DOM_LOOP(k,j,i)
    {
      if (splitcells[k][j][i] >= 0.5) //Cell is not split
      {
	    //UU holds conserved variables and so tracer mass densities rather than tracers. Thus, we do not have to multiply the tracer with rho.
	    double curmass = UU[RHO_CLOUD][k][j][i]*m_dx*m_dx*m_dx;
        double vel = UU[MX3][k][j][i] / UU[RHO][k][j][i];
	    TotMass += curmass;
	    TotVel += vel*curmass;
      }
    }
    for (nv = 0; nv < NVAR; nv++) FreeArrayMap(UU[nv]);
    FreeArrayBoxMap(splitcells, KBEG, KEND, JBEG, JEND, IBEG, IEND);
  }
  Rhov = TotVel;
  Mass = TotMass;
}

void LevelPluto::SubtractVelocity(double vel, double CoM, double t_prev)
{
/*  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    int k, j, i, nv;
    double*** UU[NVAR];
    FArrayBox& curU = m_U[dit];
    Grid* grid = m_structs_grid[dit].getGrid();
 
    IBEG = grid[IDIR].lbeg; IEND = grid[IDIR].lend;
    JBEG = grid[JDIR].lbeg; JEND = grid[JDIR].lend;
    KBEG = grid[KDIR].lbeg; KEND = grid[KDIR].lend;

    NX1 = grid[IDIR].np_int;
    NX2 = grid[JDIR].np_int;
    NX3 = grid[KDIR].np_int;

    NX1_TOT = grid[IDIR].np_tot;
    NX2_TOT = grid[JDIR].np_tot;
    NX3_TOT = grid[KDIR].np_tot;

    for (nv = 0; nv < NVAR; nv++)
    {
      UU[nv] = ArrayMap(NX3_TOT, NX2_TOT, NX1_TOT, curU.dataPtr(nv));
    }
    TOT_LOOP(k,j,i)
    {
      UU[MX3][k][j][i] -= UU[RHO][k][j][i]*vel;
    }
    for (nv = 0; nv < NVAR; nv++) FreeArrayMap(UU[nv]);
    }*/
  m_VelUpdated = false;
  m_vel = vel;
  m_dist += m_RefFrameVelocity*(g_time - t_prev);
  //  pout() << "t, tp, dist " << g_time << " " << t_prev << " " << m_dist << endl;
  m_RefFrameVelocity += vel;
  if ((prank == 0) && (m_level == 0))
  {
    FILE* f = fopen("framevel.dat", "a");
    if (f == NULL)
    {
      pout() << "Could not open framevel.dat file for writing." << endl;
      QUIT_PLUTO(1);
    }
    fprintf(f, "%.17g %.17g %.17g %.17g\n", g_time, m_RefFrameVelocity, m_dist);
    fclose(f);
  }
}

void LevelPluto::SetComovingFrame(double curtime, double& t_prev)
{
//  pout() << "SetFrame" << endl;
  m_RefFrameVelocity = 0.0;
  FILE *file;
  bool DefaultInit = false;
  if (curtime > 0.0)
  {
    if (file = fopen("framevel.dat", "r"))
    {
      //Read reference frame velocity from file                                                                                                                                          
      double time = 0.0;
      double ptime = 0.0;
      double vel = 0.0;
      double pvel = 0.0;
      double dist = 0.0;
      double pdist = 0.0;
      double CoM = 0.0;
      while (fscanf(file, "%lf %lf %lf %lf", &time, &vel, &dist) != EOF)
	{
	  //      pout() << "time, vel " << time << " " << vel << endl;
	  if (time > curtime)
	    {
	      m_RefFrameVelocity = pvel;
	      m_dist = pdist;
	      t_prev = ptime;
	      break;
	    }
	  ptime = time;
	  pvel = vel;
	  pdist = dist;
	}
      fclose(file);
      //pout() << "fdist: " << m_dist << endl;
          if (prank == 0) pout() << "time, t_prev, framevel, framedist " << curtime << ", " << t_prev << ", " << m_RefFrameVelocity << ", " << m_dist << endl;
    }
    else 
    {
      pout() << "WARNING: File framevel.dat not found for restart." << endl;
      DefaultInit = true;
    }
  }
  else DefaultInit = true;
  if (DefaultInit)
  {
    m_RefFrameVelocity = 0.0;
    m_dist = 0.0;
    t_prev = 0.0;
  }
}

void LevelPluto::setGridLevel()
{

 CH_TIME("LevelPluto::setGrid");

 m_structs_grid.define(m_grids);

 #if GEOMETRY != CARTESIAN
  m_dV.define(m_grids,CHOMBO_NDV,m_numGhost*IntVect::Unit);
 #endif

 Real dlMinLoc = 1.e30;

 for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
   {
      // The current box
      Box curBox = m_grids.get(dit());
      struct GRID* grid = m_structs_grid[dit].getGrid();

      #if GEOMETRY != CARTESIAN 
       FArrayBox& curdV = m_dV[dit()];    
      #else
       FArrayBox  curdV; 
      #endif
       
      m_patchPluto->setGrid(curBox, grid, curdV);           
       
      for (int idir = 0; idir < SpaceDim; idir++) dlMinLoc = Min(dlMinLoc,grid[idir].dl_min);   
   }

#if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)

   D_EXPAND(m_dl_min = m_dx; ,
            m_dl_min = MIN(m_dl_min,m_dx*g_x2stretch); ,
            m_dl_min = MIN(m_dl_min,m_dx*g_x3stretch); )
#else

 #ifdef CH_MPI
  Real dlMin;
  int result = MPI_Allreduce(&dlMinLoc, &dlMin, 1, MPI_CH_REAL,
                             MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){ //bark!!!
   MayDay::Error("sorry, but I had a communcation error on dlMin");
  }
  m_dl_min = dlMin;
 #else
  m_dl_min = dlMinLoc;
 #endif

#endif

}

#if GEOMETRY != CARTESIAN
const LevelData<FArrayBox>& LevelPluto::getdV() const
{
  return m_dV;
}
#endif

Real LevelPluto::getDlMin()
{

 CH_TIME("LevelPluto::getDlMin");

 return m_dl_min / m_dx;
// return m_dl_min; /* If subcycling is turned off */


}

#include "NamespaceFooter.H"
