#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iomanip>

#include "parstream.H"
#include "ParmParse.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "computeNorm.H"
#include "CH_HDF5.H"
#include "AMRIO.H"
#include "AMRLevel.H"
#include "LEINDEX.H"

#include "AMRLevelLinElast.H"

#include "LinElastPhysicsF_F.H"
#include "GodunovUtilitiesF_F.H"

// Constructor
AMRLevelLinElast::AMRLevelLinElast()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast default constructor" << endl;
  }

  m_linElastPhysics = NULL;
  m_paramsDefined = false;
  m_bdryFaceBox.resize(2*CH_SPACEDIM);
  m_bdryGrids.resize(2*CH_SPACEDIM);
  m_BNew.resize(2*CH_SPACEDIM,NULL);
  m_BOld.resize(2*CH_SPACEDIM,NULL);
  m_relateUB.resize(2*CH_SPACEDIM,NULL);
  m_bdryCoarseAverage.resize(2*CH_SPACEDIM,NULL);
  m_bdryFineInterp.resize(2*CH_SPACEDIM,NULL);
}

// Destructor
AMRLevelLinElast::~AMRLevelLinElast()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast destructor" << endl;
  }

  if (m_linElastPhysics != NULL)
  {
    delete m_linElastPhysics;
    m_linElastPhysics = NULL;
  }

  m_paramsDefined = false;
}

void AMRLevelLinElast::defineParams(const Real&                 a_cfl,
    const Real&                 a_domainLength,
    const int&                  a_verbosity,
    const Real&                 a_refineThresh,
    const Real&                 a_plasticThresh,
    const int&                  a_tagBufferSize,
    const Real&                 a_initialDtMultiplier,
    const LinElastPhysics* const a_linElastPhysics,
    const int&                  a_normalPredOrder,
    const bool&                 a_useFourthOrderSlopes,
    const bool&                 a_usePrimLimiting,
    const bool&                 a_useCharLimiting,
    const bool&                 a_useFlattening,
    const bool&                 a_useArtificialViscosity,
    const Real&                 a_artificialViscosity,
    const bool&                 a_useSourceTerm,
    const Real&                 a_sourceTermScaling,
    const bool&                 a_highOrderLimiter,
    const Vector<Real>&         a_xCoarsen,
    const Vector<Real>&         a_yCoarsen,
    const Vector<Real>&         a_zCoarsen,
    const Real&                 a_slopeCoarsen,
    const Real&                 a_widthCoarsen,
    const Vector<Real>&         a_xFaultStations,
    const Vector<Real>&         a_zFaultStations,
    const Vector<Real>&         a_xBodyStations,
    const Vector<Real>&         a_yBodyStations,
    const Vector<Real>&         a_zBodyStations,
    const Vector<Real>&         a_domainCenter,
    const string&               a_dataPrefix,
    const int&                  a_plotInterval,
    const bool&                 a_usePlasticity)
{
  // Set the CFL number
  m_cfl = a_cfl;

  // Set the physical dimension of the longest side of the domain
  m_domainLength = a_domainLength;

  verbosity(a_verbosity);

  // Set the refinement threshold
  m_refineThresh = a_refineThresh;
  m_plasticThresh = a_plasticThresh;

  // Set the tag buffer size
  m_tagBufferSize = a_tagBufferSize;

  initialDtMultiplier(a_initialDtMultiplier);

  if (m_linElastPhysics != NULL)
  {
    delete m_linElastPhysics;
    m_linElastPhysics = NULL;
  }

  m_linElastPhysics = (LinElastPhysics*) a_linElastPhysics->new_godunovPhysics();

  //BD m_bdryUseData = ((LEPhysIBC*) m_linElastPhysics->getPhysIBC())->hasBdryData();

  m_normalPredOrder = a_normalPredOrder;

  // Store the slope computation parameters
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_usePrimLimiting      = a_usePrimLimiting;
  m_useCharLimiting      = a_useCharLimiting;
  m_useFlattening        = a_useFlattening;

  // Artificial viscosity coefficient must be greater than zero
  CH_assert(!a_useArtificialViscosity || (a_artificialViscosity >= 0.0));

  // Store the artificial viscosity flag and coefficient
  m_useArtificialViscosity = a_useArtificialViscosity;
  m_artificialViscosity    = a_artificialViscosity;

  // Supply a source term to the computation
  m_useSourceTerm     = a_useSourceTerm;
  m_sourceTermScaling = a_sourceTermScaling;

  // Use a high-order limiter?
  m_highOrderLimiter = a_highOrderLimiter;

  // define the fault stations
  m_xFaultStations   = a_xFaultStations;
  m_zFaultStations   = a_zFaultStations;


  // define the body stations
  m_xBodyStations   = a_xBodyStations;
  m_yBodyStations   = a_yBodyStations;
  m_zBodyStations   = a_zBodyStations;

  // Set the center of the domain
  m_domainCenter = a_domainCenter;

  m_dataPrefix = a_dataPrefix;

  m_plotInterval = a_plotInterval;

  m_usePlasticity = a_usePlasticity;

  m_paramsDefined = true;

  m_xCoarsen     = a_xCoarsen;
  m_yCoarsen     = a_yCoarsen;
  m_zCoarsen     = a_zCoarsen;
  m_slopeCoarsen = a_slopeCoarsen;
  m_widthCoarsen = a_widthCoarsen;
}

// This instance should never get called - historical
void AMRLevelLinElast::define(AMRLevel*  a_coarserLevelPtr,
    const Box& a_problemDomain,
    int        a_level,
    int        a_refRatio)
{
  ProblemDomain physdomain(a_problemDomain);

  MayDay::Error("AMRLevelLinElast::define:\n\tShould never be called with a Box for a problem domain");
}

// Define new AMR level
void AMRLevelLinElast::define(AMRLevel*            a_coarserLevelPtr,
    const ProblemDomain& a_problemDomain,
    int                  a_level,
    int                  a_refRatio)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::define " << a_level << endl;
  }

  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
      a_problemDomain,
      a_level,
      a_refRatio);

  // Get setup information from the next coarser level
  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelLinElast* amrGodPtr = dynamic_cast<AMRLevelLinElast*>(a_coarserLevelPtr);

    if (amrGodPtr != NULL)
    {
      m_cfl           = amrGodPtr->m_cfl;
      m_domainLength  = amrGodPtr->m_domainLength;
      m_refineThresh  = amrGodPtr->m_refineThresh;
      m_tagBufferSize = amrGodPtr->m_tagBufferSize;
    }
    else
    {
      MayDay::Error("AMRLevelLinElast::define: a_coarserLevelPtr is not castable to AMRLevelLinElast*");
    }
  }

  // Compute the grid spacing
  m_dx = m_domainLength / a_problemDomain.domainBox().size(0);

  // Nominally, one layer of ghost cells is maintained permanently and
  // individual computations may create local data with more
  m_numGhost = 1;

  CH_assert(m_linElastPhysics != NULL);
  CH_assert(isDefined());
  m_linElastPhysics->define(m_problem_domain,m_dx);

  // Number and names of conserved states
  m_numStates  = m_linElastPhysics->numConserved();
  m_stateNames = m_linElastPhysics->stateNames();
  m_numBdryVars = ((LEPhysIBC*) m_linElastPhysics->getPhysIBC())->numBdryVars();
  for(int ix = 0; ix < 2*CH_SPACEDIM; ix++)
  {
    m_bdryNames.push_back(((LEPhysIBC*) m_linElastPhysics->getPhysIBC())->bdryNames(ix));
  }

  // Setup the boundary Face box
  for(int idim = 0;idim < CH_SPACEDIM;idim++)
  {
    m_bdryFaceBox[2*idim  ] = bdryLo(a_problemDomain.domainBox(),idim,1);
    m_bdryFaceBox[2*idim+1] = bdryHi(a_problemDomain.domainBox(),idim,1);
  }

  // have to do this here becase m_dx isn't defined earlier

  // Define the fault stations

  // Determine the station locations and interpolation parameters
  for(int itor = 0; itor < m_xFaultStations.size(); itor++)
  {
    // find closest station below (will do intersection between this and
    // nearst neighbor below)
    IntVect tmpLoc = 
      BASISV(0)*floor((m_xFaultStations[itor]+m_domainCenter[0])/m_dx-0.5)
      +
      BASISV(2)*floor((m_zFaultStations[itor]+m_domainCenter[2])/m_dx-0.5);
    m_ivFaultStations.push_back(tmpLoc);
    // distance to the station (for intepolation purposes)
    m_dxFaultStation.push_back(-(tmpLoc[0]+0.5)*m_dx + (m_xFaultStations[itor]+m_domainCenter[0]));
    if(SpaceDim > 2)
    {
      m_dzFaultStation.push_back(-(tmpLoc[2]+0.5)*m_dx + (m_zFaultStations[itor]+m_domainCenter[2]));
    }
    else
    {
      m_dzFaultStation.push_back(0);
    }
  }


  // Define the station name for this level
  for(int itor = 0; itor < m_xFaultStations.size(); itor++)
  {
    char tmpName[128];
    sprintf(tmpName,"%sFx%5.4fz%5.4f.L%d.dat",m_dataPrefix.c_str(),m_xFaultStations[itor],m_zFaultStations[itor],m_level);
    // pout() << tmpName << endl;
    m_faultSationNames.push_back(tmpName);
    if(procID() == 0)
    {
      FILE * stationData;
      stationData = fopen(m_faultSationNames[itor].c_str(),"w");
      // fprintf(stationData,"%% dx :: %E   IV :: (%d %d)  delta :: (%E %E)\n",
      //     m_dx,
      //     m_ivFaultStations[itor][0],
      //     m_ivFaultStations[itor][2],
      //     m_dxFaultStation[itor],
      //     m_dzFaultStation[itor]);

      fprintf(stationData,"#number of values = %d\n",m_numBdryVars[2]);
      fprintf(stationData,"#size of real = %d\n",sizeof(Real));
      fprintf(stationData,"#time");
      for(int i_c = 0; i_c < m_numBdryVars[2]; i_c++)
      {
        fprintf(stationData,"\t%s",m_bdryNames[2][i_c].c_str());
      }
      fprintf(stationData,"\tlevel\tprocID\n");
      fclose(stationData);
    }
  }

  // Define the Body stations

  // Determine the station locations and interpolation parameters
  for(int itor = 0; itor < m_xBodyStations.size(); itor++)
  {
    // find closest station above (will do intersection between this and
    // nearst neighbor bellow)
    IntVect tmpLoc = 
      BASISV(0)*floor((m_xBodyStations[itor]+m_domainCenter[0])/m_dx-0.5)
      +
      BASISV(1)*floor((m_yBodyStations[itor]+m_domainCenter[1])/m_dx-0.5)
      +
      BASISV(2)*floor((m_zBodyStations[itor]+m_domainCenter[2])/m_dx-0.5);
    m_ivBodyStations.push_back(tmpLoc);

    // distance to the station (for intepolation purposes)
    m_dxBodyStation.push_back(-(tmpLoc[0]+0.5)*m_dx + (m_xBodyStations[itor]+m_domainCenter[0]));
    m_dyBodyStation.push_back(-(tmpLoc[1]+0.5)*m_dx + (m_yBodyStations[itor]+m_domainCenter[1]));
    if(SpaceDim > 2)
      m_dzBodyStation.push_back(-(tmpLoc[2]+0.5)*m_dx + (m_zBodyStations[itor]+m_domainCenter[2]));
    else
      m_dzBodyStation.push_back(0);
  }


  // Define the station name for this level
  for(int itor = 0; itor < m_xBodyStations.size(); itor++)
  {
    char tmpName[128];
    sprintf(tmpName,"%sBx%5.4fy%5.4fz%5.4f.L%d.dat",m_dataPrefix.c_str(),m_xBodyStations[itor],m_yBodyStations[itor],m_zBodyStations[itor],m_level);
    // pout() << tmpName << endl;
    m_bodySationNames.push_back(tmpName);
    if(procID() == 0)
    {
      FILE * stationData;
      stationData = fopen(m_bodySationNames[itor].c_str(),"w");
      // fprintf(stationData,"%% dx :: %E   IV :: (%d %d %d)   delta :: (%E %E %E)\n",
      //     m_dx,
      //     m_ivBodyStations[itor][0],
      //     m_ivBodyStations[itor][1],
      //     m_ivBodyStations[itor][2],
      //     m_dxBodyStation[itor],
      //     m_dyBodyStation[itor],
      //     m_dzBodyStation[itor]);

      fprintf(stationData,"#number of values = %d\n",m_numStates);
      fprintf(stationData,"#size of real = %d\n",sizeof(Real));
      fprintf(stationData,"#time");
      for(int i_c = 0; i_c < m_numStates; i_c++)
      {
        fprintf(stationData,"\t%s",m_stateNames[i_c].c_str());
      }
      fprintf(stationData,"\tlevel\tprocID\n");
      fclose(stationData);
    }
  }


}

// Advance by one timestep
Real AMRLevelLinElast::advance()
{
  CH_assert(allDefined());

  if (s_verbosity >= 1 && procID() == 0)
  {
    for(int i = 0; i < m_level;i++) pout() << "   ";
    pout() << "AMRLevelLinElast::advance level " << m_level << " to time " << scientific << m_time + m_dt << " by dt " << m_dt << endl;
  }

  if (s_verbosity >= 2)
  {
    for(int i = 0; i < m_level;i++) pout() << "   ";
    pout() << "AMRLevelLinElast::advance level " << m_level << " to time " << scientific << m_time + m_dt << " by dt " << m_dt << endl;
  }

  // Copy the new to the old
  for (DataIterator dit = m_UNew.dataIterator(); dit.ok(); ++dit)
  {
    m_UOld[dit()].copy(m_UNew[dit()]);
  }

  for(int ix = 0;ix < 2*CH_SPACEDIM; ix++)
    for (DataIterator dit = (*m_BNew[ix]).dataIterator(); dit.ok(); ++dit)
      (*m_BOld[ix])[dit()].copy((*m_BNew[ix])[dit()]);

  Real newDt = 0.0;

  // Set up arguments to LinElastLevelGodunov::step based on whether there are
  // coarser and finer levels

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;

  // Undefined leveldata in case we need it
  const LevelData<FArrayBox> dummyData;

  // Set arguments to dummy values and then fix if real values are available
  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  const LevelData<FArrayBox>* coarserDataOld = &dummyData;
  const LevelData<FArrayBox>* coarserDataNew = &dummyData;

  Vector<LevelData<FArrayBox>* > coarserBdryOld;
  coarserBdryOld.resize(2*CH_SPACEDIM,NULL);
  Vector<LevelData<FArrayBox>* > coarserBdryNew;
  coarserBdryNew.resize(2*CH_SPACEDIM,NULL);

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  // A coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelLinElast* coarserPtr = getCoarserLevel();

    // Recall that my flux register goes between my level and the next
    // finer level
    coarserFR = &coarserPtr->m_fluxRegister;

    coarserDataOld = &coarserPtr->m_UOld;
    coarserDataNew = &coarserPtr->m_UNew;

    coarserBdryOld = coarserPtr->m_BOld;
    coarserBdryNew = coarserPtr->m_BNew;

    tCoarserNew = coarserPtr->m_time;
    tCoarserOld = tCoarserNew - coarserPtr->m_dt;
  }

  // A finer level exists
  if (m_hasFiner)
  {
    // Recall that my flux register goes between my level and the next
    // finer level
    finerFR = &m_fluxRegister;
  }

  // Source term leveldata
  LevelData<FArrayBox> sourceData;

  // we don't need the flux in the simple hyperbolic case...
  LevelData<FArrayBox> flux[SpaceDim];

  // Advance the solve one timestep
  newDt = m_LElevelGodunov.step(m_UNew,
      m_BNew,
      m_relateUB,
      flux,
      *finerFR,
      *coarserFR,
      sourceData,
      *coarserDataOld,
      coarserBdryOld,
      tCoarserOld,
      *coarserDataNew,
      coarserBdryNew,
      tCoarserNew,
      m_time,
      m_dt);

  if(m_usePlasticity) m_LElevelGodunov.plasticUpdate(m_UNew,m_dt);

  // Update the time and store the new timestep
  m_time += m_dt;
  Real returnDt = m_cfl * newDt;

  m_dtNew = returnDt;

  return returnDt;
}

// Things to do after a timestep
void AMRLevelLinElast::postTimeStep()
{
  CH_assert(allDefined());

  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::postTimeStep " << m_level << endl;
  }

  if (m_hasFiner)
  {
    // Reflux
    Real scale = -1.0/m_dx;
    m_fluxRegister.reflux(m_UNew,scale);

    // Average from finer level data
    AMRLevelLinElast* amrGodFinerPtr = getFinerLevel();

    amrGodFinerPtr->m_coarseAverage.averageToCoarse(m_UNew,
        amrGodFinerPtr->m_UNew);

    for(int ix = 0;ix < 2*CH_SPACEDIM;ix++)
      (*amrGodFinerPtr->m_bdryCoarseAverage[ix]).averageToCoarse((*m_BNew[ix]),
          *(amrGodFinerPtr->m_BNew[ix]));
  }

  if (s_verbosity >= 2 && m_level == 0)
  {
    int nRefFine = 1;

    pout() << "AMRLevelLinElast::postTimeStep:" << endl;
    pout() << "  Sums:" << endl;
    for (int comp = 0; comp < m_numStates; comp++)
    {
      Interval curComp(comp,comp);
      // Real integral = computeSum(m_UNew,NULL,nRefFine,m_dx,curComp);
      Real integral = computeNorm(m_UNew,NULL,nRefFine,m_dx,curComp);

      pout() << "    " << setw(23)
        << setprecision(16)
        << setiosflags(ios::showpoint)
        << setiosflags(ios::scientific)
        << integral
        << " --- " << m_stateNames[comp];

      if (comp == 0 && !first)
      {
        pout() << " (" << setw(23)
          << setprecision(16)
          << setiosflags(ios::showpoint)
          << setiosflags(ios::scientific)
          << (integral-last_integral)
          << " " << setw(23)
          << setprecision(16)
          << setiosflags(ios::showpoint)
          << setiosflags(ios::scientific)
          << (integral-orig_integral)
          << ")";
      }

      pout() << endl;

      if (comp == 0)
      {
        if (first)
        {
          orig_integral = integral;
          first = false;
        }

        last_integral = integral;
      }
    }
  }

  writeStationLevel();

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::postTimeStep " << m_level << " finished" << endl;
  }
}

// Create tags for regridding
void AMRLevelLinElast::tagCells(IntVectSet& a_tags)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::tagCells " << m_level << endl;
  }

  // Create tags based on undivided gradient of density
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  IntVectSet localTags;

  // If there is a coarser level interpolate undefined ghost cells
  if (m_hasCoarser)
  {
    const AMRLevelLinElast* amrGodCoarserPtr = getCoarserLevel();

    PiecewiseLinearFillPatch pwl(levelDomain,
        amrGodCoarserPtr->m_UNew.disjointBoxLayout(),
        m_numStates,
        amrGodCoarserPtr->m_problem_domain,
        amrGodCoarserPtr->m_ref_ratio,
        1);

    pwl.fillInterp(m_UNew,
        amrGodCoarserPtr->m_UNew,
        amrGodCoarserPtr->m_UNew,
        1.0,
        0,
        0,
        m_numStates);
  }
  m_UNew.exchange(Interval(0,m_numStates-1));

  // Compute relative gradient
  DataIterator dit = levelDomain.dataIterator();
  LEPhysIBC* lephysIBCPtr = (LEPhysIBC*) m_linElastPhysics->getPhysIBC();

  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& b = levelDomain[dit()];
    FArrayBox gradFab(b,SpaceDim);
    const FArrayBox& UFab = m_UNew[dit()];

    for (int dir = 0; dir < SpaceDim; ++dir)
    {
      const Box bCenter = b & grow(m_problem_domain,-BASISV(dir));

      const Box bLo     = b & adjCellLo(bCenter,dir);
      const int hasLo = ! bLo.isEmpty();

      const Box bHi     = b & adjCellHi(bCenter,dir);
      const int hasHi = ! bHi.isEmpty();

      //JK Chombo sample code all uses Relative Gradients, but this doesn't
      //JK work since we have near zero values

      //JK FORT_GETGRADF(
      //JK     CHF_FRA1(gradFab,dir),
      //JK     CHF_CONST_FRA1(UFab,8),
      //JK     CHF_CONST_INT(dir),
      //JK     CHF_BOX(bLo),
      //JK     CHF_CONST_INT(hasLo),
      //JK     CHF_BOX(bHi),
      //JK     CHF_CONST_INT(hasHi),
      //JK     CHF_BOX(bCenter));

      FORT_GETFULLGRADF(
          CHF_FRA1(gradFab,dir),
          CHF_CONST_FRA(UFab),
          CHF_CONST_INT(dir),
          CHF_BOX(bLo),
          CHF_CONST_INT(hasLo),
          CHF_BOX(bHi),
          CHF_CONST_INT(hasHi),
          CHF_BOX(bCenter));
    }

    FArrayBox gradMagFab(b,1);
    FORT_MAXF(CHF_FRA1(gradMagFab,0),
        CHF_CONST_FRA(gradFab),
        CHF_BOX(b));

    lephysIBCPtr->tagCells(gradMagFab,m_time,m_refineThresh);

    // Tag where gradient exceeds threshold
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();


      Real r = 0;
#if(CH_SPACEDIM >= 1)
      Real x = iv[0]*m_dx;
      r = max(r,max(m_xCoarsen[0]-x,x-m_xCoarsen[1]));
#endif
#if(CH_SPACEDIM >= 2)
      Real y = iv[1]*m_dx;
      r = max(r,max(m_yCoarsen[0]-y,y-m_yCoarsen[1]));
#endif
#if(CH_SPACEDIM >= 3)
      Real z = iv[1]*m_dx;
      r = max(r,max(m_zCoarsen[0]-z,z-m_zCoarsen[1]));
#endif

      Real refineThresh = m_refineThresh*(1+r*m_slopeCoarsen/m_widthCoarsen);

      if (gradMagFab(iv) >= refineThresh || UFab.get(iv,IX1_GAM) >= m_plasticThresh || UFab.get(iv,IX2_GAM) >= m_plasticThresh)
      {
        localTags |= iv;
        //pout() << gradFab(iv,0) << "     " << gradFab(iv,1) << "     " << gradMagFab(iv) << "    " << m_refineThresh << endl;
      }
    }
  }

  localTags.grow(m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;

  a_tags = localTags;
}

// Create tags at initialization
void AMRLevelLinElast::tagCellsInit(IntVectSet& a_tags)
{
  CH_assert(allDefined());
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::tagCellsInit " << m_level << endl;
  }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.

  // tagCells(a_tags);


  //JK This is a silly way to do things, but it works...
  CH_assert(allDefined());


  // Create tags based on undivided gradient of density
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  IntVectSet localTags;

  // Compute relative gradient
  DataIterator dit = levelDomain.dataIterator();

  LEPhysIBC* lephysIBCPtr = (LEPhysIBC*) m_linElastPhysics->getPhysIBC();

  for (dit.begin(); dit.ok(); ++dit)
  {
    // Let the Physics IBC handle the initial griding
    const Box& b = levelDomain[dit()];
    FArrayBox markFAB(b,1);
    markFAB.setVal(0.0);
    if(lephysIBCPtr->tagCellsInit(markFAB,1))
    {
      // Tag where gradient exceeds threshold
      BoxIterator bit(b);
      for (bit.begin(); bit.ok(); ++bit)
      {
        const IntVect& iv = bit();

        if (markFAB(iv) > 0)
        {
          // pout() << iv << endl;
          localTags |= iv;
        }
      }
    }
  }

  localTags.grow(m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;

  a_tags = localTags;
}

// Set up data on this level after regridding
void AMRLevelLinElast::regrid(const Vector<Box>& a_newGrids)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::regrid " << m_level << endl;
  }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);

  /////////
  /////////
  /////////
  /////////
  for(int idim = 0;idim < CH_SPACEDIM;idim++)
  {
    // Indicate/guarantee that the indexing below is only for reading
    // otherwise an error/assertion failure occurs
    const DisjointBoxLayout& constGrids = m_grids;
    Vector<Box> vectBndBoxLo;
    Vector<int> vectBndPIDLo;

    Vector<Box> vectBndBoxHi;
    Vector<int> vectBndPIDHi;

    if (s_verbosity >= 4)
    {
      pout() << "new grids: " << endl;
    }
    for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
    {
      if (s_verbosity >= 4)
      {
        pout() << "grid:          " << constGrids[lit()] << endl;
      }

      // We set the boundary box for this level to the intersection of the
      // box edge with the boundary. I believe that this should let the
      // boundary grid have the same layout structure as the underlying
      // grid
      const Box tmpBndBoxLo = (m_bdryFaceBox[2*idim] & bdryLo(constGrids[lit()],idim,1));
      if(!tmpBndBoxLo.isEmpty())
      {
        vectBndBoxLo.push_back(tmpBndBoxLo);
        vectBndPIDLo.push_back(constGrids.procID(lit()));
        if (s_verbosity >= 4)
        {
          pout() << "boundary grid["<< 2*idim<<"]: " << tmpBndBoxLo << " for: " << constGrids[lit()] << endl;
        }
      }

      const Box tmpBndBoxHi = (m_bdryFaceBox[2*idim+1] & bdryHi(constGrids[lit()],idim,1));
      if(!tmpBndBoxHi.isEmpty())
      {
        vectBndBoxHi.push_back(tmpBndBoxHi);
        vectBndPIDHi.push_back(constGrids.procID(lit()));
        if (s_verbosity >= 4)
        {
          pout() << "boundary grid["<< 2*idim+1<<"]: " << tmpBndBoxHi << " for: " << constGrids[lit()] << endl;
        }
      }
    }

    DisjointBoxLayout tmpBndGridsLo(vectBndBoxLo,vectBndPIDLo);
    m_bdryGrids[2*idim] = tmpBndGridsLo;

    DisjointBoxLayout tmpBndGridsHi(vectBndBoxHi,vectBndPIDHi);
    m_bdryGrids[2*idim+1] = tmpBndGridsHi;
  }
  /////////
  /////////
  /////////
  /////////

  // Save data for later
  for(DataIterator dit = m_UNew.dataIterator(); dit.ok(); ++dit)
  {
    m_UOld[dit()].copy(m_UNew[dit()]);
  }

  for(int ix = 0;ix < 2*CH_SPACEDIM; ix++)
    for (DataIterator dit = (*m_BNew[ix]).dataIterator(); dit.ok(); ++dit)
      (*m_BOld[ix])[dit()].copy((*m_BNew[ix])[dit()]);

  // Reshape state with new grids
  IntVect ivGhost = m_numGhost * IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);
  m_UNew.define(m_grids,m_numStates,ivGhost);
  for(int idim = 0;idim < CH_SPACEDIM;idim++)
  {
    IntVect ivBGhost = m_numGhost * (IntVect::Unit - BASISV(idim));
    for(int ix = 2*idim;ix < 2*(idim+1); ix++)
    {
      (*m_BNew[ix]).define(m_bdryGrids[ix],m_numBdryVars[ix],ivBGhost);
      (*m_relateUB[ix]).define(m_grids);
    }
  }

  setupRelateUB();

  // Reshape the boundary grid
  LEPhysIBC* lephysIBCPtr = (LEPhysIBC*) m_linElastPhysics->getPhysIBC();
  lephysIBCPtr->initialize(m_UNew);
  //BD if(m_bdryUseData)
  //BD {
  for(int ix = 0; ix < 2*CH_SPACEDIM; ix++) lephysIBCPtr->initializeBdry((*m_BNew[ix]), ix);
  //BD }

  // Set up data structures
  levelSetup();

  // Interpolate from coarser level
  if (m_hasCoarser)
  {
    AMRLevelLinElast* amrGodCoarserPtr = getCoarserLevel();
    m_fineInterp.interpToFine(m_UNew,amrGodCoarserPtr->m_UNew);
    for(int ix = 0; ix < 2*CH_SPACEDIM; ix++)
      (*m_bdryFineInterp[ix]).interpToFine((*m_BNew[ix]),*(amrGodCoarserPtr->m_BNew[ix]));
  }

  // Copy from old state
  m_UOld.copyTo(m_UOld.interval(),
      m_UNew,
      m_UNew.interval());

  m_UOld.define(m_grids,m_numStates,ivGhost);

  for(int idim = 0;idim < CH_SPACEDIM;idim++)
  {
    IntVect ivBGhost = m_numGhost * (IntVect::Unit - BASISV(idim));
    for(int ix = 2*idim;ix < 2*(idim+1); ix++)
    {
      (*m_BOld[ix]).copyTo((*m_BOld[ix]).interval(),(*m_BNew[ix]),(*m_BNew[ix]).interval());
      (*m_BOld[ix]).define(m_bdryGrids[ix],m_numBdryVars[ix],ivBGhost);
    }
  }
}

// Initialize grids
void AMRLevelLinElast::initialGrid(const Vector<Box>& a_newGrids)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::initialGrid " << m_level << endl;
  }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);

  /////////
  /////////
  /////////
  /////////
  for(int idim = 0;idim < CH_SPACEDIM;idim++)
  {
    // Indicate/guarantee that the indexing below is only for reading
    // otherwise an error/assertion failure occurs
    const DisjointBoxLayout& constGrids = m_grids;
    Vector<Box> vectBndBoxLo;
    Vector<int> vectBndPIDLo;

    Vector<Box> vectBndBoxHi;
    Vector<int> vectBndPIDHi;

    if (s_verbosity >= 4)
    {
      pout() << "new grids: " << endl;
    }
    for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
    {
      if (s_verbosity >= 4)
      {
        pout() << "grid:          " << constGrids[lit()] << endl;
      }

      // We set the boundary box for this level to the intersection of the
      // box edge with the boundary. I believe that this should let the
      // boundary grid have the same layout structure as the underlying
      // grid
      const Box tmpBndBoxLo = (m_bdryFaceBox[2*idim] & bdryLo(constGrids[lit()],idim,1));
      if(!tmpBndBoxLo.isEmpty())
      {
        vectBndBoxLo.push_back(tmpBndBoxLo);
        vectBndPIDLo.push_back(constGrids.procID(lit()));
        if (s_verbosity >= 4)
        {
          pout() << "boundary grid["<<2*idim<<"]: " << tmpBndBoxLo << " for: " << constGrids[lit()] << endl;
        }
      }

      const Box tmpBndBoxHi = (m_bdryFaceBox[2*idim+1] & bdryHi(constGrids[lit()],idim,1));
      if(!tmpBndBoxHi.isEmpty())
      {
        vectBndBoxHi.push_back(tmpBndBoxHi);
        vectBndPIDHi.push_back(constGrids.procID(lit()));
        if (s_verbosity >= 4)
        {
          pout() << "boundary grid["<<2*idim+1<<"]: " << tmpBndBoxHi << " for: " << constGrids[lit()] << endl;
        }
      }
    }


    DisjointBoxLayout tmpBndGridsLo(vectBndBoxLo,vectBndPIDLo);
    m_bdryGrids[2*idim] = tmpBndGridsLo;

    DisjointBoxLayout tmpBndGridsHi(vectBndBoxHi,vectBndPIDHi);
    m_bdryGrids[2*idim+1] = tmpBndGridsHi;
  }
  /////////
  /////////
  /////////
  /////////

  // Define old and new state data structures
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);
  m_UOld.define(m_grids,m_numStates,ivGhost);


  for(int idim = 0;idim < CH_SPACEDIM;idim++)
  {
    // Define the old and new boundary data structures
    IntVect ivBGhost = m_numGhost*(IntVect::Unit-BASISV(idim));
    for(int ix = 2*idim;ix < 2*(idim+1);ix++)
    {
      m_BNew[ix] = new LevelData<FArrayBox>(m_bdryGrids[ix],m_numBdryVars[ix],ivBGhost);
      m_BOld[ix] = new LevelData<FArrayBox>(m_bdryGrids[ix],m_numBdryVars[ix],ivBGhost);

      // Define the data holder to relate the grid data to the boundary data
      m_relateUB[ix] = new LayoutData<GridData>();
      (*m_relateUB[ix]).define(m_grids);
    }
  }

  setupRelateUB();

  // Set up data structures
  levelSetup();

}

// Initialize data
void AMRLevelLinElast::initialData()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::initialData " << m_level << endl;
  }

  LEPhysIBC* lephysIBCPtr = (LEPhysIBC*) m_linElastPhysics->getPhysIBC();
  lephysIBCPtr->initialize(m_UNew);
  //BD if(m_bdryUseData)
  //BD {
  for(int ix = 0; ix < 2*CH_SPACEDIM; ix++) lephysIBCPtr->initializeBdry(*m_BNew[ix], ix);
  //BD }
}

// Things to do after initialization
void AMRLevelLinElast::postInitialize()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::postInitialize " << m_level << endl;
  }

  if (m_hasFiner)
  {
    // Volume weighted average from finer level data
    AMRLevelLinElast* amrGodFinerPtr = getFinerLevel();

    amrGodFinerPtr->m_coarseAverage.averageToCoarse(m_UNew,
        amrGodFinerPtr->m_UNew);

    for(int ix = 0;ix < 2*CH_SPACEDIM;ix++)
      (*amrGodFinerPtr->m_bdryCoarseAverage[ix]).averageToCoarse(*m_BNew[ix],
          *(amrGodFinerPtr->m_BNew[ix]));
  }
}

#ifdef CH_USE_HDF5

// Write checkpoint header
void AMRLevelLinElast::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::writeCheckpointHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = m_stateNames[comp];
  }

  header.m_int["num_boundary_data"] = 2*CH_SPACEDIM;
  
  for(int ix = 0; ix < 2*CH_SPACEDIM; ix++)
  {
    sprintf(compStr,"num_boundary_%d_variables",ix);
    header.m_int[compStr] = m_numBdryVars[ix];
    for (int comp = 0; comp < m_numBdryVars[ix]; ++comp)
    {
      sprintf(compStr,"boundary_%d_variable_%d",ix,comp);
      header.m_string[compStr] = m_bdryNames[ix][comp];
    }
  }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

// Write checkpoint data for this level
void AMRLevelLinElast::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::writeCheckpointLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_int ["tag_buffer_size"] = m_tagBufferSize;
  header.m_real["dx"]              = m_dx;
  header.m_real["dt"]              = m_dt;
  header.m_real["time"]            = m_time;
  header.m_box ["prob_domain"]     = m_problem_domain.domainBox();

  // Setup the periodicity info
  D_TERM(if (m_problem_domain.isPeriodic(0))
      {
      header.m_int ["is_periodic_0"] = 1;
      }
      else
      {
      header.m_int ["is_periodic_0"] = 0;
      } ,

      if (m_problem_domain.isPeriodic(1))
      {
      header.m_int ["is_periodic_1"] = 1;
      }
      else
      {
      header.m_int ["is_periodic_1"] = 0;
      } ,

      if (m_problem_domain.isPeriodic(2))
      {
      header.m_int ["is_periodic_2"] = 1;
      }
      else
      {
        header.m_int ["is_periodic_2"] = 0;
      } );

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  // Write the data for this level
  write(a_handle,m_UNew.boxLayout());
  write(a_handle,m_UNew,"data");

  // Write the data for this level
  for(int ix = 0; ix < 2*CH_SPACEDIM; ix++)
  {
    char tmpName[60];
    sprintf(tmpName,"boundary_%d_data",ix);
    pout() << "START: "<<  tmpName<<endl;
    if((*m_BNew[ix]).boxLayout().numCells() > 0)
      write(a_handle,*m_BNew[ix],tmpName);
    pout() << "DONE: "<<  tmpName<<endl;
  }
}

// Read checkpoint header
void AMRLevelLinElast::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::readCheckpointHeader" << endl;
  }

  // Reader the header
  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << "hdf5 header data:" << endl;
    pout() << header << endl;
  }

  // Get the number of components
  if (header.m_int.find("num_components") == header.m_int.end())
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointHeader: checkpoint file does not have num_components");
  }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates)
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointHeader: num_components in checkpoint file does not match solver");
  }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    if (header.m_string.find(compStr) == header.m_string.end())
    {
      MayDay::Error("AMRLevelLinElast::readCheckpointHeader: checkpoint file does not have enough component names");
    }

    stateName = header.m_string[compStr];
    if (stateName != m_stateNames[comp])
    {
      MayDay::Error("AMRLevelLinElast::readCheckpointHeader: state_name in checkpoint does not match solver");
    }
  }

  //
  // Get the number of components
  if (header.m_int.find("num_boundary_data") == header.m_int.end())
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointHeader: checkpoint file does not have num_boundary_data");
  }

  int numBdryData = header.m_int["num_boundary_data"];
  if (numBdryData != 2*CH_SPACEDIM)
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointHeader: num_boundary_data in checkpoint file does not match solver");
  }

  //
  // Get the number of components
  for(int ix = 0; ix < 2*CH_SPACEDIM; ix++)
  {
    sprintf(compStr,"num_boundary_%d_variables",ix);
    if (header.m_int.find(compStr) == header.m_int.end())
    {
      sprintf(compStr,"AMRLevelLinElast::readCheckpointHeader: checkpoint file does not have %s",compStr);
      MayDay::Error(compStr);
    }

    int numBdryVars = header.m_int[compStr];
    if (numBdryVars != m_numBdryVars[ix])
    {
      sprintf(compStr,"AMRLevelLinElast::readCheckpointHeader: %s in checkpoint file does not match solver",compStr);
      MayDay::Error(compStr);
    }

    // Get the component names
    std::string bdryStateName;
    // char compStr[60];
    for (int comp = 0; comp < m_numBdryVars[ix]; ++comp)
    {
      sprintf(compStr,"boundary_%d_variable_%d",ix,comp);
      if (header.m_string.find(compStr) == header.m_string.end())
      {
        MayDay::Error("AMRLevelLinElast::readCheckpointHeader: checkpoint file does not have enough component names");
      }

      bdryStateName = header.m_string[compStr];
      if (bdryStateName != m_bdryNames[ix][comp])
      {
        sprintf(compStr,"AMRLevelLinElast::readCheckpointHeader: %s in checkpoint does not match solver",compStr);
        MayDay::Error(compStr);
      }
    }
  }
}

// Read checkpoint data for this level
void AMRLevelLinElast::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::readCheckpointLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << "hdf5 header data:" << endl;
    pout() << header << endl;
  }

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointLevel: file does not contain ref_ratio");
  }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
  {
    pout() << "read ref_ratio = " << m_ref_ratio << endl;
  }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointLevel: file does not contain tag_buffer_size");
  }
  m_tagBufferSize = header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
  }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointLevel: file does not contain dx");
  }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
  {
    pout() << "read dx = " << m_dx << endl;
  }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointLevel: file does not contain dt");
  }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
  {
    pout() << "read dt = " << m_dt << endl;
  }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointLevel: file does not contain time");
  }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
  {
    pout() << "read time = " << m_time << endl;
  }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointLevel: file does not contain prob_domain");
  }

  Box domainBox = header.m_box["prob_domain"];

  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  bool isPeriodic[SpaceDim];
  D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
      {
      isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
      }
      else
      {
      isPeriodic[0] = false;
      } ,

      if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
      {
      isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
      }
      else
      {
      isPeriodic[1] = false;
      } ,

      if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
      {
      isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
      }
      else
      {
        isPeriodic[2] = false;
      } );

  m_problem_domain = ProblemDomain(domainBox,isPeriodic);

  // Get the grids
  Vector<Box> grids;
  const int gridStatus = read(a_handle,grids);

  if (gridStatus != 0)
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointLevel: file does not contain a Vector<Box>");
  }

  // Create level domain
  m_grids = loadBalance(grids);

  // Indicate/guarantee that the indexing below is only for reading
  // otherwise an error/assertion failure occurs
  const DisjointBoxLayout& constGrids = m_grids;

  LayoutIterator lit = constGrids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
  {
    const Box& b = constGrids[lit()];
    m_level_grids.push_back(b);
  }

  if (s_verbosity >= 4)
  {
    pout() << "read level domain: " << endl;
    LayoutIterator lit = m_grids.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = m_grids[lit()];
      pout() << lit().intCode() << ": " << b << endl;
    }
    pout() << endl;
  }

  // Reshape state with new grids
  m_UNew.define(m_grids,m_numStates);
  const int dataStatus = read<FArrayBox>(a_handle,
      m_UNew,
      "data",
      m_grids);

  if (dataStatus != 0)
  {
    MayDay::Error("AMRLevelLinElast::readCheckpointLevel: file does not contain state data");
  }
  m_UOld.define(m_grids,m_numStates);

  for(int idim = 0;idim < CH_SPACEDIM;idim++)
  {
    // Indicate/guarantee that the indexing below is only for reading
    // otherwise an error/assertion failure occurs
    const DisjointBoxLayout& constGrids = m_grids;
    Vector<Box> vectBndBoxLo;
    Vector<int> vectBndPIDLo;

    Vector<Box> vectBndBoxHi;
    Vector<int> vectBndPIDHi;

    if (s_verbosity >= 4)
    {
      pout() << "new grids: " << endl;
    }
    for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
    {
      if (s_verbosity >= 4)
      {
        pout() << "grid:          " << constGrids[lit()] << endl;
      }

      // We set the boundary box for this level to the intersection of the
      // box edge with the boundary. I believe that this should let the
      // boundary grid have the same layout structure as the underlying
      // grid
      const Box tmpBndBoxLo = (m_bdryFaceBox[2*idim] & bdryLo(constGrids[lit()],idim,1));
      if(!tmpBndBoxLo.isEmpty())
      {
        vectBndBoxLo.push_back(tmpBndBoxLo);
        vectBndPIDLo.push_back(constGrids.procID(lit()));
        if (s_verbosity >= 4)
        {
          pout() << "boundary grid["<< 2*idim<<"]: " << tmpBndBoxLo << " for: " << constGrids[lit()] << endl;
        }
      }

      const Box tmpBndBoxHi = (m_bdryFaceBox[2*idim+1] & bdryHi(constGrids[lit()],idim,1));
      if(!tmpBndBoxHi.isEmpty())
      {
        vectBndBoxHi.push_back(tmpBndBoxHi);
        vectBndPIDHi.push_back(constGrids.procID(lit()));
        if (s_verbosity >= 4)
        {
          pout() << "boundary grid["<< 2*idim+1<<"]: " << tmpBndBoxHi << " for: " << constGrids[lit()] << endl;
        }
      }
    }

    DisjointBoxLayout tmpBndGridsLo(vectBndBoxLo,vectBndPIDLo);
    m_bdryGrids[2*idim] = tmpBndGridsLo;

    DisjointBoxLayout tmpBndGridsHi(vectBndBoxHi,vectBndPIDHi);
    m_bdryGrids[2*idim+1] = tmpBndGridsHi;
  }

  for(int idim = 0;idim < CH_SPACEDIM;idim++)
  {
    for(int ix = 2*idim; ix < 2*(idim+1); ix++)
    {

      IntVect ivBGhost = m_numGhost*(IntVect::Unit-BASISV(idim));
      m_BNew[ix] = new LevelData<FArrayBox>(m_bdryGrids[ix],m_numBdryVars[ix],ivBGhost);
      m_BOld[ix] = new LevelData<FArrayBox>(m_bdryGrids[ix],m_numBdryVars[ix],ivBGhost);

      char tmpName[60];
      sprintf(tmpName,"boundary_%d_data",ix);

      // only query the checkpoint file if we have data on this level!
      if(m_bdryGrids[ix].size() != 0)
      {
        const int bdryDataStatus = read<FArrayBox>(a_handle,
            *m_BNew[ix], //TODO: BOUNDARY
            tmpName,
            m_bdryGrids[ix]);

        if (bdryDataStatus != 0)
        {
          sprintf(tmpName,"AMRLevelLinElast::readCheckpointLevel: file does not contain %s",tmpName);
          MayDay::Error(tmpName);
        }
      }

      // Define the data holder to relate the grid data to the boundary data
      m_relateUB[ix] = new LayoutData<GridData>();
      (*m_relateUB[ix]).define(m_grids);
    }
  }

  setupRelateUB();

  // Set up data structures
  levelSetup();
}

// Write plotfile header
void AMRLevelLinElast::writePlotHeader(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::writePlotHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = m_stateNames[comp];
  }

  // Write the header
  header.writeToFile(a_handle);
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  m_LElevelGodunov.getGodunovPhysicsPtrConst()->expressions(expressions);
  expressions.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

// Write plotfile data for this level
void AMRLevelLinElast::writePlotLevel(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::writePlotLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx;
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  // Write the data for this level
  write(a_handle,m_UNew.boxLayout());
  write(a_handle,m_UNew,"data", IntVect::Unit);
}


int AMRLevelLinElast::numBdryLevels(int ix)
{
  int numLevel = 0;
  if((*m_BNew[ix]).boxLayout().numCells()>0) numLevel = m_level;
  if(m_hasFiner) numLevel = max(numLevel, getFinerLevel()->numBdryLevels(ix));
  return numLevel;
}

// Write boundary plotfile data for this level
void AMRLevelLinElast::writeCustomPlotFile(const std::string& a_prefix,
    const int& a_max_level,
    const int& a_finest_level,
    const int& a_cur_step,
    const Real& a_cur_time,
    const bool& a_conclude)
{
  CH_TIME("AMR::writeCustomPlotFile");

  // writeStationLevel();

  if(((m_plotInterval > 0)
        && (a_cur_step % m_plotInterval == 0))
      ||
      ((m_plotInterval >= 0) && a_conclude))
  {

    CH_assert(m_isDefined);

    if (s_verbosity >= 3)
    {
      pout() << "AMRLevelLinElast::writeCustomPlotFile" << endl;
    }
    
    for(int ix = 0; ix < 2*CH_SPACEDIM; ix++)
    {

      char iter_str[80];

      sprintf(iter_str,
          "%sboundary%d.%06d.%dd.hdf5",
          a_prefix.c_str(), ix, a_cur_step, SpaceDim);

      if (s_verbosity >= 2)
      {
        pout() << "plot file name = " << iter_str << endl;
      }

      HDF5Handle handle(iter_str, HDF5Handle::CREATE);

      // write amr data
      HDF5HeaderData header;
      header.m_int ["max_level"]  = a_max_level;
      header.m_int ["num_levels"] = numBdryLevels(ix) + 1;
      header.m_int ["iteration"]  = a_cur_step;
      header.m_real["time"]       = a_cur_time;

      // write the boundary physics class header data

      // Setup the number of components
      header.m_int["num_components"] = m_numBdryVars[ix];

      // Setup the component names
      char compStr[30];
      for (int comp = 0; comp < m_numBdryVars[ix]; ++comp)
      {
        sprintf(compStr,"component_%d",comp);
        header.m_string[compStr] = m_bdryNames[ix][comp];
      }

      if (s_verbosity >= 3)
      {
        pout() << header << endl;
      }

      header.writeToFile(handle);

      // handle.setGroup("/Expressions");
      // HDF5HeaderData expressions;
      // m_LElevelGodunov.getGodunovPhysicsPtrConst()->expressions(expressions);
      // expressions.writeToFile(handle);

      if (s_verbosity >= 3)
      {
        pout() << header << endl;
      }

      // write physics class per-level data
      writeThisBdryLevel(handle,ix);

      handle.close();
    }
  }
}

void AMRLevelLinElast::writeThisBdryLevel(HDF5Handle& a_handle, int ix)
{
  //JK Need to fix this so that it only puts data on the boundary
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::writeThisBdryLevel" << endl;
  }
  if((*m_BNew[ix]).boxLayout().numCells()>0)
  {

    // Setup the level string
    char levelStr[20];
    sprintf(levelStr,"%d",m_level);
    const std::string label = std::string("level_") + levelStr;

    a_handle.setGroup(label);

    // Setup the level header information
    HDF5HeaderData header;

    header.m_int ["ref_ratio"]   = m_ref_ratio;
    header.m_real["dx"]          = m_dx;
    header.m_real["dt"]          = m_dt;
    header.m_real["time"]        = m_time;
    //header.m_box ["prob_domain"] = m_problem_domain.domainBox();
    header.m_box ["prob_domain"] = m_bdryFaceBox[ix];

    // Write the header for this level
    header.writeToFile(a_handle);

    if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }

    // Write the data for this level
    // write(a_handle,(*m_BNew[ix]).boxLayout());
    // write(a_handle,(*m_BNew[ix]),"data", IntVect::Unit);
    write(a_handle,(*m_BNew[ix]).boxLayout());
    write(a_handle,(*m_BNew[ix]),"data", IntVect::Unit);
  }

  if (m_hasFiner)
  {
    getFinerLevel()->writeThisBdryLevel(a_handle,ix);
  }
}

void AMRLevelLinElast::writeStationLevel()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::writeStationLevel" << endl;
  }

  // First we do that fault stations
  if((*m_BNew[2]).boxLayout().numCells()>0)
  {
    for (DataIterator dit = (*m_BNew[2]).dataIterator(); dit.ok(); ++dit)
    {
      for(int itor = 0; itor < m_ivFaultStations.size();itor++)
      {
        // We get the box from the disjoint box layout since the that
        // data has a layer of ghost cells. Since we have one layer of
        // ghost cells, we know that if we have the m_ivFaultStations
        // data in the disjoint box we also have this plus  one in each
        // direction (which might be needed for interpolation)
        const Box tmpBndBox = m_bdryGrids[2][dit()];
        IntVect tmpLoc = m_ivFaultStations[itor];
        if(tmpBndBox.contains(tmpLoc))
        {
          FILE * stationData;
          stationData = fopen(m_faultSationNames[itor].c_str(),"a");

          fwrite (&m_time, sizeof(Real), 1, stationData );
          for(int i_c = 0; i_c < m_numBdryVars[2];i_c++)
          {
            if(SpaceDim < 3)
            {
              Real value = (*m_BNew[2])[dit()].get(tmpLoc,i_c)*(1-m_dxFaultStation[itor])
                +(*m_BNew[2])[dit()].get(tmpLoc+BASISV(0),i_c)*m_dxFaultStation[itor];
              fwrite(&value,sizeof(Real), 1, stationData );
            }
            else
            {
              Real value = (*m_BNew[2])[dit()].get(tmpLoc,i_c)*(1-m_dxFaultStation[itor])*(1-m_dzFaultStation[itor])
                +(*m_BNew[2])[dit()].get(tmpLoc+BASISV(0),i_c)*m_dxFaultStation[itor]*(1-m_dzFaultStation[itor])
                +(*m_BNew[2])[dit()].get(tmpLoc+BASISV(2),i_c)*(1-m_dxFaultStation[itor])*m_dzFaultStation[itor]
                +(*m_BNew[2])[dit()].get(tmpLoc+BASISV(0)+BASISV(2),i_c)*m_dxFaultStation[itor]*m_dzFaultStation[itor];
              fwrite(&value,sizeof(Real), 1, stationData );
            }
          }
          Real t_level = m_level;
          fwrite (&t_level, sizeof(Real), 1, stationData );
          Real t_procID = procID();
          fwrite (&t_procID, sizeof(Real), 1, stationData );
          fclose(stationData);
        }
      }
    }
  }

  // Now the Body stations
  if(m_UNew.boxLayout().numCells()>0)
  {
    for (DataIterator dit = m_UNew.dataIterator(); dit.ok(); ++dit)
    {
      for(int itor = 0; itor < m_ivBodyStations.size();itor++)
      {
        // We get the box from the disjoint box layout since the that
        // data has a layer of ghost cells. Since we have one layer of
        // ghost cells, we know that if we have the m_ivBodyStations
        // data in the disjoint box we also have this plus  one in each
        // direction (which might be needed for interpolation)
        const Box tmpBndBox = m_grids[dit()];
        IntVect tmpLoc = m_ivBodyStations[itor];
        if(tmpBndBox.contains(tmpLoc))
        {
          FILE * stationData;
          stationData = fopen(m_bodySationNames[itor].c_str(),"a");

          fwrite (&m_time, sizeof(Real), 1, stationData );
          for(int i_c = 0; i_c < m_numStates;i_c++)
          {
            if(SpaceDim < 3)
            {
              Real value = m_UNew[dit()].get(tmpLoc,i_c)*(1-m_dxBodyStation[itor])*(1-m_dyBodyStation[itor])
                +m_UNew[dit()].get(tmpLoc+BASISV(0),i_c)*m_dxBodyStation[itor]*(1-m_dyBodyStation[itor])
                +m_UNew[dit()].get(tmpLoc+BASISV(1),i_c)*(1-m_dxBodyStation[itor])*m_dyBodyStation[itor]
                +m_UNew[dit()].get(tmpLoc+BASISV(0)+BASISV(1),i_c)*m_dxBodyStation[itor]*m_dyBodyStation[itor];
              fwrite(&value,sizeof(Real), 1, stationData );
            }
            else
            {
              Real value = m_UNew[dit()].get(tmpLoc,i_c)*(1-m_dxBodyStation[itor])*(1-m_dyBodyStation[itor])*(1-m_dzBodyStation[itor])
                +m_UNew[dit()].get(tmpLoc+BASISV(0),i_c)*m_dxBodyStation[itor]*(1-m_dyBodyStation[itor])*(1-m_dzBodyStation[itor])
                +m_UNew[dit()].get(tmpLoc+BASISV(1),i_c)*(1-m_dxBodyStation[itor])*m_dyBodyStation[itor]*(1-m_dzBodyStation[itor])
                +m_UNew[dit()].get(tmpLoc+BASISV(0)+BASISV(1),i_c)*m_dxBodyStation[itor]*m_dyBodyStation[itor]*(1-m_dzBodyStation[itor])
                +m_UNew[dit()].get(tmpLoc+BASISV(2),i_c)*(1-m_dxBodyStation[itor])*(1-m_dyBodyStation[itor])*m_dzBodyStation[itor]
                +m_UNew[dit()].get(tmpLoc+BASISV(2)+BASISV(0),i_c)*m_dxBodyStation[itor]*(1-m_dyBodyStation[itor])*m_dzBodyStation[itor]
                +m_UNew[dit()].get(tmpLoc+BASISV(2)+BASISV(1),i_c)*(1-m_dxBodyStation[itor])*m_dyBodyStation[itor]*m_dzBodyStation[itor]
                +m_UNew[dit()].get(tmpLoc+BASISV(2)+BASISV(0)+BASISV(1),i_c)*m_dxBodyStation[itor]*m_dyBodyStation[itor]*m_dzBodyStation[itor];
              fwrite(&value,sizeof(Real), 1, stationData );
            }
          }
          Real t_level = m_level;
          fwrite (&t_level, sizeof(Real), 1, stationData );
          Real t_procID = procID();
          fwrite (&t_procID, sizeof(Real), 1, stationData );
          fclose(stationData);
        }
      }
    }
  }

  // if (m_hasFiner)
  // {
  //     getFinerLevel()->writeStationLevel();
  // }
}

#endif

// Returns the dt computed earlier for this level
Real AMRLevelLinElast::computeDt()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::computeDt " << m_level << endl;
  }

  Real newDt;
  newDt = m_dtNew;

  return newDt;
}

// Compute dt using initial data
Real AMRLevelLinElast::computeInitialDt()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::computeInitialDt " << m_level << endl;
  }

  Real newDT = m_initial_dt_multiplier * m_dx / m_LElevelGodunov.getMaxWaveSpeed(m_UNew);

  return newDT;
}

const LevelData<FArrayBox>& AMRLevelLinElast::getStateNew() const
{
  CH_assert(allDefined());

  return m_UNew;
}

const LevelData<FArrayBox>& AMRLevelLinElast::getStateOld() const
{
  CH_assert(allDefined());

  return m_UOld;
}

bool AMRLevelLinElast::allDefined() const
{
  return isDefined()     &&
    m_paramsDefined ;
}

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelLinElast::loadBalance(const Vector<Box>& a_grids)
{
  CH_assert(allDefined());

  // Load balance and create boxlayout
  Vector<int> procMap;

  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4)
  {
    pout() << "AMRLevelLinElast::loadBalance: procesor map: " << endl;

    for (int igrid = 0; igrid < a_grids.size(); ++igrid)
    {
      pout() << igrid << ": " << procMap[igrid] << "  " << endl;
    }

    pout() << endl;
  }

  DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
  dbl.close();

  return dbl;
}

// Setup menagerie of data structures
void AMRLevelLinElast::levelSetup()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelLinElast::levelSetup " << m_level << endl;
  }

  AMRLevelLinElast* amrGodCoarserPtr = getCoarserLevel();
  AMRLevelLinElast* amrGodFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrGodCoarserPtr != NULL);
  m_hasFiner   = (amrGodFinerPtr   != NULL);

  if (m_hasCoarser)
  {
    int nRefCrse = m_coarser_level_ptr->refRatio();

    m_coarseAverage.define(m_grids,
        m_numStates,
        nRefCrse);

    m_fineInterp.define(m_grids,
        m_numStates,
        nRefCrse,
        m_problem_domain);

    for(int idim = 0;idim < CH_SPACEDIM;idim++)
    {
      for(int ix = 2*idim; ix < 2*(idim+1); ix++)
      {
        m_bdryCoarseAverage[ix] = new D1CoarseAverage();
        (*m_bdryCoarseAverage[ix]).define(m_bdryGrids[ix],
            m_numBdryVars[ix],
            nRefCrse,
            idim);

        m_bdryFineInterp[ix] = new D1FineInterp();
        (*m_bdryFineInterp[ix]).define(m_bdryGrids[ix],
            m_numBdryVars[ix],
            nRefCrse,
            idim,
            m_bdryFaceBox[ix]);
      }
    }

    const DisjointBoxLayout& coarserLevelDomain = amrGodCoarserPtr->m_grids;
    const Vector<DisjointBoxLayout>& bdryCoarserLevelDomain = amrGodCoarserPtr->m_bdryGrids;

    // Maintain levelGodunov
    m_LElevelGodunov.define(m_grids,
        m_bdryGrids,
        coarserLevelDomain,
        bdryCoarserLevelDomain,
        m_problem_domain,
        m_bdryFaceBox,
        nRefCrse,
        m_dx,
        m_linElastPhysics,
        m_normalPredOrder,
        m_useFourthOrderSlopes,
        m_usePrimLimiting,
        m_useCharLimiting,
        m_useFlattening,
        m_useArtificialViscosity,
        m_artificialViscosity,
        m_hasCoarser,
        m_hasFiner);
    m_LElevelGodunov.highOrderLimiter(m_highOrderLimiter);

    // This may look twisted but you have to do this this way because the
    // coarser levels get setup before the finer levels so, since a flux
    // register lives between this level and the next FINER level, the finer
    // level has to do the setup because it is the only one with the
    // information at the time of construction.

    // Maintain flux registers
    amrGodCoarserPtr->m_fluxRegister.define(m_grids,
        amrGodCoarserPtr->m_grids,
        m_problem_domain,
        amrGodCoarserPtr->m_ref_ratio,
        m_numStates);
    amrGodCoarserPtr->m_fluxRegister.setToZero();
  }
  else
  {
    m_LElevelGodunov.define(m_grids,
        m_bdryGrids,
        DisjointBoxLayout(),
        Vector<DisjointBoxLayout>(2*CH_SPACEDIM,DisjointBoxLayout()),
        m_problem_domain,
        m_bdryFaceBox,
        m_ref_ratio,
        m_dx,
        m_linElastPhysics,
        m_normalPredOrder,
        m_useFourthOrderSlopes,
        m_usePrimLimiting,
        m_useCharLimiting,
        m_useFlattening,
        m_useArtificialViscosity,
        m_artificialViscosity,
        m_hasCoarser,
        m_hasFiner);
    m_LElevelGodunov.highOrderLimiter(m_highOrderLimiter);
  }
}

// Get the next coarser level
AMRLevelLinElast* AMRLevelLinElast::getCoarserLevel() const
{
  CH_assert(allDefined());

  AMRLevelLinElast* amrGodCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
  {
    amrGodCoarserPtr = dynamic_cast<AMRLevelLinElast*>(m_coarser_level_ptr);

    if (amrGodCoarserPtr == NULL)
    {
      MayDay::Error("AMRLevelLinElast::getCoarserLevel: dynamic cast failed");
    }
  }

  return amrGodCoarserPtr;
}

// Get the next finer level
AMRLevelLinElast* AMRLevelLinElast::getFinerLevel() const
{
  CH_assert(allDefined());

  AMRLevelLinElast* amrGodFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
  {
    amrGodFinerPtr = dynamic_cast<AMRLevelLinElast*>(m_finer_level_ptr);

    if (amrGodFinerPtr == NULL)
    {
      MayDay::Error("AMRLevelLinElast::getFinerLevel: dynamic cast failed");
    }
  }

  return amrGodFinerPtr;
}

// Setup the relateUB object
void AMRLevelLinElast::setupRelateUB()
{
  for(int idim = 0; idim < CH_SPACEDIM; idim++)
  {
    for (DataIterator dit = m_UNew.dataIterator(); dit.ok(); ++dit)
    {
      const Box tmpBndBoxLo = (m_bdryFaceBox[2*idim] & bdryLo(m_UNew[dit()].box(),idim,1+m_numGhost));
      if(!tmpBndBoxLo.isEmpty())
      {
        for (DataIterator bit = (*m_BNew[2*idim]).dataIterator(); bit.ok(); ++bit)
        {
          if((*m_BNew[2*idim])[bit()].box().contains(tmpBndBoxLo))
          {
            (*m_relateUB[2*idim])[dit()].setBdryIndex(bit());
            break;
          }
        }
      }

      const Box tmpBndBoxHi = (m_bdryFaceBox[2*idim+1] & bdryHi(m_UNew[dit()].box(),idim,1+m_numGhost));
      if(!tmpBndBoxHi.isEmpty())
      {
        for (DataIterator bit = (*m_BNew[2*idim+1]).dataIterator(); bit.ok(); ++bit)
        {
          if((*m_BNew[2*idim+1])[bit()].box().contains(tmpBndBoxHi))
          {
            (*m_relateUB[2*idim+1])[dit()].setBdryIndex(bit());
            break;
          }
        }
      }
    }
  }
}

void AMRLevelLinElast::dumpBdryData()
{
  return;
  if(m_level == 0)
  {
    char prefix[30];
    sprintf(prefix,"boundary_data.%d",procID());
    FILE * boundaryData;
    boundaryData = fopen(prefix,"w");
    if(boundaryData != NULL)
    {
      LEPhysIBC* lephysIBCPtr = (LEPhysIBC*) m_linElastPhysics->getPhysIBC();
      for (DataIterator dit = (*m_BNew[2]).dataIterator(); dit.ok(); ++dit) //TODO: BOUNDARY
      {
        lephysIBCPtr->setBdryData(&(*m_BNew[2])[dit()],2); //TODO: BOUNDARY
        lephysIBCPtr->dumpBdryData(boundaryData);
      }
      fclose(boundaryData);
    }
  }
}
