#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMRLEVELLINELASTFACTORY_H_
#define _AMRLEVELLINELASTFACTORY_H_

#include "AMRLevelFactory.H"

#include "AMRLevelLinElast.H"

/// AMR Godunov factory
/**
 */
class AMRLevelLinElastFactory: public AMRLevelFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelLinElastFactory();

  /// Destructor
  /**
   */
  virtual ~AMRLevelLinElastFactory();

  /// Define the parameters of the factory
  /**
   */
  virtual void define(const Real&                 a_cfl,
                      const Real&                 a_domainLength,
                      const int&                  a_verbosity,
                      const Real&                 a_refineThresh,
                      const int&                  a_tagBufferSize,
                      const Real&                 a_initialDtMultiplier,
                      const GodunovPhysics* const a_godunovPhysics,
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
                      const Vector<Real>&         a_xFaultStations,
                      const Vector<Real>&         a_zFaultStations,
                      const Vector<Real>&         a_xBodyStations,
                      const Vector<Real>&         a_yBodyStations,
                      const Vector<Real>&         a_zBodyStations,
                      const Vector<Real>&         a_domainCenter);

  /// Create a new AMRLevel (for polytropic gases)
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Check that this object has been defined
  /**
   */
  bool isDefined() const;

  /// Set the prefix for output files
  /**
   */
  void dataPrefix(const string& a_dataPrefix);

  /// Set the prefix for output files
  /**
   */
  void plotInterval(const int& a_plotInterval);


protected:
  // CFL number
  Real m_cfl;

  string m_dataPrefix;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // Verbosity of the object
  int m_verbosity;

  // Refinement threshold for gradient
  Real m_refineThresh;

  // Tag buffer size
  int m_tagBufferSize;

  // Initial dt multiplier
  Real m_initialDtMultiplier;

  // Pointer to the object that supplies the physics needed by the integrator
  // (used as a factory here)
  const GodunovPhysics* m_godunovPhysics;

  // Order of the normal predictor (1 -> PLM, 2 -> PPM)
  int m_normalPredOrder;

  // Use 4th other slope computations (otherwise, use 2nd order)
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
  
  // Has this object been defined
  bool m_isDefined;

  // Save the fault stations we will need
  Vector<Real> m_xFaultStations;
  Vector<Real> m_zFaultStations;

  // Save the body stations we will need
  Vector<Real> m_xBodyStations;
  Vector<Real> m_yBodyStations;
  Vector<Real> m_zBodyStations;

  // Where is the domain center
  Vector<Real> m_domainCenter;

  // How often to make plots files
  int m_plotInterval;

private:
  bool m_usePlasticity;

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelLinElastFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelLinElastFactory(const AMRLevelLinElastFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
