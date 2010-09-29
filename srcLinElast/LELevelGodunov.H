#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LELEVELGODUNOV_H_
#define _LELEVELGODUNOV_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "PiecewiseLinearFillPatch.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "Copier.H"

#include "PatchGodunov.H"
#include "LinElastPhysics.H"
#include "NamespaceHeader.H"

/// Level Godunov - really a hyperbolic level integrator
/**
*/
class LELevelGodunov
{
public:
    /// Default constructor
    /**
      Object requires define() to be called before all other functions.
      */
    LELevelGodunov();

    /// Destructor
    /**
      Destroys all objects created by define(). Passed in data references
      of define() are left alone.
      */
    ~LELevelGodunov();

    /// Actual constructor.
    /**
      Inside the routine, we cast away const-ness on the data members
      for the assignment.  The arguments passed in are maintained const
      (coding standards).  a_refine is the refinement ratio between
      this level and the next coarser level.  For the coarsest level, an
      empty DisjointBoxLayout is passed in for coarserDisjointBoxLayout.

      This has is just the inherited version plues the boundary data.
      */
    void define(const DisjointBoxLayout&    a_thisDisjointBoxLayout,
        //BD const DisjointBoxLayout&    a_bdryDisjointBoxLayout,
        const DisjointBoxLayout&    a_coarserDisjointBoxLayout,
        const ProblemDomain&        a_domain,
        const int&                  a_refineCoarse,
        const Real&                 a_dx,
        const LinElastPhysics* const a_godunovFactory,
        const int&                  a_normalPredOrder,
        const bool&                 a_useFourthOrderSlopes,
        const bool&                 a_usePrimLimiting,
        const bool&                 a_useCharLimiting,
        const bool&                 a_useFlattening,
        const bool&                 a_useArtificialViscosity,
        const Real&                 a_artificialViscosity,
        const bool&                 a_hasCoarser,
        const bool&                 a_hasFiner);

    /// Take one timestep for this grid level.
    /**
      For the coarsest level, UcoarseOld and UCoarseNew are empty
      LevelData<FArrayBox> objects.  Also, either UCoarseOld or UCoarseNew
      might be empty to indicate that t(nf) = t(nc) the one grid is at the
      current time and interpolation is not required for boundary condition
      generation.  U must be defined on the same DisjointBoxLayouts as were
      used in define().  a_flux is is used more generally to pass
      face-centered data in and out of the function. Coarse flux register
      is flux register with the next coarser level.  Fine flux register is
      the flux register with the next finer level.  To the finer level FR,
      this level is the coarse level.  To the coarser level FR, this level
      is the fine level.  If source terms do not exist, "a_S" should be null
      constructed and not defined (i.e. its define() should not be called).

      This has is just the inherited version plues the boundary data.
      */
    Real step(LevelData<FArrayBox>&       a_U,
        //BD LevelData<FArrayBox>&       a_bdryPsi,
        LevelData<FArrayBox>        a_flux[CH_SPACEDIM],
        LevelFluxRegister&          a_finerFluxRegister,
        LevelFluxRegister&          a_coarserFluxRegister,
        const LevelData<FArrayBox>& a_S,
        const LevelData<FArrayBox>& a_UCoarseOld,
        const Real&                 a_TCoarseOld,
        const LevelData<FArrayBox>& a_UCoarseNew,
        const Real&                 a_TCoarseNew,
        const Real&                 a_time,
        const Real&                 a_dt);

    /// To be added...
    /**
    */
    void computeWHalf(LayoutData<FluxBox>&        a_WHalf,
        LevelData<FArrayBox>&       a_U,
        const LevelData<FArrayBox>& a_S,
        const LevelData<FArrayBox>& a_UCoarseOld,
        const Real&                 a_TCoarseOld,
        const LevelData<FArrayBox>& a_UCoarseNew,
        const Real&                 a_TCoarseNew,
        const Real&                 a_time,
        const Real&                 a_dt);

    /// To be added...
    /**
    */
    Real computeUpdate(LevelData<FArrayBox>&       a_dU,
        LevelFluxRegister&          a_finerFluxRegister,
        LevelFluxRegister&          a_coarserFluxRegister,
        const LevelData<FArrayBox>& a_U,
        const LayoutData<FluxBox>&  a_WHalf,
        const Real&                 a_time,
        const Real&                 a_dt);

    /// Get maximum wave speed
    /**
    */
    Real getMaxWaveSpeed(const LevelData<FArrayBox>& a_U);

    /// get the godunovPhysics pointer from the PatchGodunov member
    /**
      Returns the pointer to the LinElastPhysics object which is
      used in the PatchGodunov class.
      */
    LinElastPhysics* getGodunovPhysicsPtr();


    const LinElastPhysics* getGodunovPhysicsPtrConst() const ;

    void highOrderLimiter(bool a_highOrderLimiter);

protected:
    // Box layout for this level
    DisjointBoxLayout m_grids;

    // Patch integrator
    PatchGodunov    m_LEpatchGodunov;

    // Number of ghost cells need locally for this level
    int m_numGhost;

    // Exchange copier
    Copier m_exchangeCopier;

    // Temporary storage space for conserved variables
    LevelData<FArrayBox> m_U;

    // Interpolator for filling in ghost cells from the next coarser level
    PiecewiseLinearFillPatch m_patcher;

    // Grid spacing
    Real m_dx;

    // Problem domain - index space for this level
    ProblemDomain m_domain;

    // Refinement ratio between this level and the next coarser
    int m_refineCoarse;

    // Number of conserved variables
    int m_numCons;

    // Number of conservative fluxes and other face centered quantities
    int m_numFluxes;

    // Flags for the existence of coarser and finer levels, respectively
    bool m_hasCoarser;
    bool m_hasFiner;

    // Order of the normal predictor (1 -> PLM, 2-> PPM)
    int m_normalPredOrder;

    // Use 4th order slope computations (otherwise, use 2nd order)
    bool m_useFourthOrderSlopes;

    // Do slope limiting in the primative or characteristic variables, respect.
    bool m_usePrimLimiting;
    bool m_useCharLimiting;

    // Do slope flattening - MUST BE USING 4th order slopes
    bool m_useFlattening;

    // Apply artificial viscosity of a set value
    bool m_useArtificialViscosity;
    Real m_artificialViscosity;

    // Has this object been defined
    bool m_isDefined;

    // Box layout for this level
    //BD DisjointBoxLayout m_bdryGrids;

    // Temporary storage space for boundary data variables
    //BD LevelData<FArrayBox> m_bdryPsi;

    // Exchange copier
    //BD Copier m_bdryExchangeCopier;


private:
    // Disallowed for all the usual reasons
    void operator=(const LELevelGodunov& a_input)
    {
        MayDay::Error("invalid operator");
    }

    // Disallowed for all the usual reasons
    LELevelGodunov(const LELevelGodunov& a_input)
    {
        MayDay::Error("invalid operator");
    }
};

#include "NamespaceFooter.H"
#endif
