#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMRLEVELLINELAST_H_
#define _AMRLEVELLINELAST_H_

#include "FArrayBox.H"
#include "GridData.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "D1CoarseAverage.H"
#include "FineInterp.H"
#include "D1FineInterp.H"
#include "LevelFluxRegister.H"
#include "Box.H"
#include "IntVectSet.H"
#include "IVSFAB.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"

#include "LinElastLevelGodunov.H"
#include "PhysIBC.H"
#include "LEPhysIBC.H"
#include "LinElastPhysics.H"

#include "UsingNamespace.H"

/// AMR Godunov
/**
 */
class AMRLevelLinElast: public AMRLevel
{
public:
  /// Constructor
  /**
   */
  AMRLevelLinElast();

  /// Destructor
  /**
   */
  virtual ~AMRLevelLinElast();

  /// Define the parameters the object needs
  /**
   */
  void defineParams(const Real&                 a_cfl,
                    const Real&                 a_domainLength,
                    const int&                  a_verbosity,
                    const Real&                 a_refineThresh,
                    const Real&                 a_plasticThresh,
                    const int&                  a_tagBufferSize,
                    const Real&                 a_initialDtMultiplier,
                    const LinElastPhysics* const a_godunovPhysics,
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
                    const Vector<Real>&         a_xCoarsen2,
                    const Vector<Real>&         a_yCoarsen2,
                    const Vector<Real>&         a_zCoarsen2,
                    const int&                  a_coarsen2MaxLevel,
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
                    const bool&                  a_usePlasticity);

  /// This instance should never get called - historical
  /**
   */
  virtual void define(AMRLevel*  a_coarserLevelPtr,
                      const Box& a_problemDomain,
                      int        a_level,
                      int        a_refRatio);

  /// Define new AMR level
  /**
   */
  virtual void define(AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);

  /// Advance by one timestep
  /**
   */
  virtual Real advance();

  /// Things to do after a timestep
  /**
   */
  virtual void postTimeStep();

  /// Create tags for regridding
  /**
   */
  virtual void tagCells(IntVectSet& a_tags);

  /// Create tags at initialization
  /**
   */
  virtual void tagCellsInit(IntVectSet& a_tags);

  /// Set up data on this level after regridding
  /**
   */
  virtual void regrid(const Vector<Box>& a_newGrids);

  /// Initialize grids
  /**
   */
  virtual void initialGrid(const Vector<Box>& a_newGrids);

  /// Initialize data
  /**
   */
  virtual void initialData();

  /// Things to do after initialization
  /**
   */
  virtual void postInitialize();

#ifdef CH_USE_HDF5
  /// Write checkpoint header
  /**
   */
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// Write checkpoint data for this level
  /**
   */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// Read checkpoint header
  /**
   */
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// Read checkpoint data for this level
  /**
   */
  virtual void readCheckpointLevel(HDF5Handle& a_handle);

  /// Write plotfile header
  /**
   */
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  /// Write plotfile data for this level
  /**
   */
  virtual void writePlotLevel(HDF5Handle& a_handle) const;

  //// pass the user stuff to do their own plot files
  /**
   * Allows doing your own plot files. Only called on the first level, user must
   * handle all other functionality.
   */
  virtual void writeCustomPlotFile(const std::string& a_prefix,
      const int& a_max_level,
      const int& a_finest_level,
      const int& a_step,
      const Real& a_time,
      const bool& a_conclude);

#endif

  /// Returns the dt computed earlier for this level
  /**
   */
  virtual Real computeDt();

  /// Compute dt using initial data
  /**
   */
  virtual Real computeInitialDt();

  ///
  const LevelData<FArrayBox>& getStateNew() const;

  ///
  const LevelData<FArrayBox>& getStateOld() const;

  ///
  bool allDefined() const;


  ///
  int level(){return m_level;}

  /// Calls the IBC routine on the boundary data do dump it to a file
  /**
   * Primarily used to save the rupture and seismic data
   */
  virtual void dumpBdryData();

protected:
  // Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  // setup the pointers in the relateUB layoutData object
  void setupRelateUB();

  // Setup menagerie of data structures
  void levelSetup();

  // Get the next coarser level
  AMRLevelLinElast* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelLinElast* getFinerLevel() const;

  // Conserved state, U, at old and new time
  LevelData<FArrayBox> m_UOld,m_UNew;
  LevelData<FArrayBox> m_BOld,m_BNew;
  LevelData<FArrayBox> m_BOld_old,m_BNew_old;
  LayoutData<GridData> m_relateUB;

  // CFL number
  Real m_cfl;

  // Grid spacing
  Real m_dx;

  // Interpolation from coarse to fine level
  FineInterp m_fineInterp;

  // Averaging from fine to coarse level
  CoarseAverage m_coarseAverage;

  // Interpolation from coarse to fine level
  D1FineInterp m_bdryFineInterp;

  // Averaging from fine to coarse level
  D1CoarseAverage m_bdryCoarseAverage;

  // New time step
  Real m_dtNew;

  // Number of conserved states
  int m_numStates;

  // Number of boundary varialbes
  int m_numBdryVars;

  // Names of conserved states
  Vector<string> m_stateNames;
  Vector<string> m_bdryNames;

  // Number of ghost cells (in each direction)
  int m_numGhost;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // Level integrator
  LinElastLevelGodunov m_LElevelGodunov;

  // Flux register
  LevelFluxRegister m_fluxRegister;

  // Pointer to the class defining the physics of the problem
  LinElastPhysics* m_linElastPhysics;

  // Order of the normal predictor (1 -> PLM, 2-> PPM)
  int m_normalPredOrder;

  // Use 4th order slope computations (otherwise, use 2nd order)
  bool m_useFourthOrderSlopes;

  // Do slope limiting in the primitive or characteristic variables, respect.
  bool m_usePrimLimiting;
  bool m_useCharLimiting;

  // Do slope flattening - MUST BE USING 4th order slopes
  bool m_useFlattening;

  // Apply artificial viscosity of a set value
  bool m_useArtificialViscosity;
  Real m_artificialViscosity;

  // Supply a source term to the computation
  bool m_useSourceTerm;
  Real m_sourceTermScaling;

  // Use a high-order limiter?
  bool m_highOrderLimiter;
  
  // Refinement threshold for gradient
  Real m_refineThresh;
  Real m_plasticThresh;

  // Tag buffer size
  int m_tagBufferSize;

  // Flag coarser and finer levels
  bool m_hasCoarser;
  bool m_hasFiner;

  // Grid layout for this level
  DisjointBoxLayout m_grids;

  // True if all the parameters for this object are defined
  bool m_paramsDefined;

  // boolean to keep track of boundary data usage
  bool m_bdryUseData;

  // Boundary Box, face centered
  Box m_bdryFaceBox;

  // Boundary Grid Layout for this level
  DisjointBoxLayout m_bdryGrids;

  // Save the fault stations we will need
  Vector<Real>     m_xFaultStations;
  Vector<Real>     m_zFaultStations;
  Vector<string>   m_faultSationNames;
  Vector<IntVect> m_ivFaultStations;
  Vector<Real>    m_dxFaultStation;
  Vector<Real>    m_dzFaultStation;

  // Save the fault stations we will need
  Vector<Real>    m_xBodyStations;
  Vector<Real>    m_yBodyStations;
  Vector<Real>    m_zBodyStations;
  Vector<string>  m_bodySationNames;
  Vector<IntVect> m_ivBodyStations;
  Vector<Real>    m_dxBodyStation;
  Vector<Real>    m_dyBodyStation;
  Vector<Real>    m_dzBodyStation;

  // Need to know where the center is for moving the stations
  Vector<Real> m_domainCenter;

  // Used for save the files
  string m_dataPrefix;

  int m_plotInterval;


private:
  Vector<Real> m_xCoarsen2;
  Vector<Real> m_yCoarsen2;
  Vector<Real> m_zCoarsen2;
  int m_coarsen2MaxLevel;
  Vector<Real> m_xCoarsen;
  Vector<Real> m_yCoarsen;
  Vector<Real> m_zCoarsen;
  Real m_slopeCoarsen;
  Real m_widthCoarsen;

  bool m_usePlasticity;

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelLinElast& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelLinElast(const AMRLevelLinElast& a_input)
  {
    MayDay::Error("invalid operator");
  }
#ifdef CH_USE_HDF5
  void writeThisBdryLevel(HDF5Handle& a_handle);
  void writeStationLevel();
#endif
};

#endif
