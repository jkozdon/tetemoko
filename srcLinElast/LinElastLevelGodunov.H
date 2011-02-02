#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LINELASTLEVELGODUNOV_H_
#define _LINELASTLEVELGODUNOV_H_

#include "FArrayBox.H"
#include "GridData.H"
#include "FluxBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "PiecewiseLinearFillPatch.H"
#include "D1PiecewiseLinearFillPatch.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "Copier.H"
#include "LevelGodunov.H"

#include "LinElastPatchGodunov.H"
#include "NamespaceHeader.H"


/// Linear Elastic Level Godunov - really a hyperbolic level integrator
/**
*/
class LinElastLevelGodunov
{
public:
    /// Default constructor
    /**
      Object requires define() to be called before all other functions.
      */
    LinElastLevelGodunov();

    /// Destructor
    /**
      Destroys all objects created by define(). Passed in data references
      of define() are left alone.
      */
    ~LinElastLevelGodunov();

    /// Actual constructor.
    /**
      Inside the routine, we cast away const-ness on the data members
      for the assignment.  The arguments passed in are maintained const
      (coding standards).  a_refine is the refinement ratio between
      this level and the next coarser level.  For the coarsest level, an
      empty DisjointBoxLayout is passed in for coarserDisjointBoxLayout.
      */
    void define(const DisjointBoxLayout&    a_thisDisjointBoxLayout,
        const DisjointBoxLayout&    a_bdryDisjointBoxLayout,
        const DisjointBoxLayout&    a_coarserDisjointBoxLayout,
        const DisjointBoxLayout&    a_bdryCoarserDisjointBoxLayout,
        const ProblemDomain&        a_domain,
        const Box&                  a_bdryBox,
        const int&                  a_refineCoarse,
        const Real&                 a_dx,
        const GodunovPhysics* const a_godunovFactory,
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
      */
    Real step(LevelData<FArrayBox>&       a_U,
        LevelData<FArrayBox>&       a_B,
        LayoutData<GridData>&       a_relateUB,
        LevelData<FArrayBox>        a_flux[CH_SPACEDIM],
        LevelFluxRegister&          a_finerFluxRegister,
        LevelFluxRegister&          a_coarserFluxRegister,
        const LevelData<FArrayBox>& a_S,
        const LevelData<FArrayBox>& a_UCoarseOld,
        const LevelData<FArrayBox>& a_BCoarseOld,
        const Real&                 a_TCoarseOld,
        const LevelData<FArrayBox>& a_UCoarseNew,
        const LevelData<FArrayBox>& a_BCoarseNew,
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

    /// get the godunovPhysics pointer from the LinElastPatchGodunov member
    /**
      Returns the pointer to the GodunovPhysics object which is
      used in the LinElastPatchGodunov class.
      */
    GodunovPhysics* getGodunovPhysicsPtr();


    const GodunovPhysics* getGodunovPhysicsPtrConst() const ;

    void highOrderLimiter(bool a_highOrderLimiter);

protected:
    // Box layout for this level
    DisjointBoxLayout m_grids;

    // Boundary box layout for this level
    DisjointBoxLayout m_bdry;

    // Patch integrator
    LinElastPatchGodunov    m_patchGodunov;

    // Number of ghost cells need locally for this level
    int m_numGhost;

    // Exchange copier
    Copier m_exchangeCopier;
    Copier m_bdryExchangeCopier;

    // Temporary storage space for conserved variables
    LevelData<FArrayBox> m_U;

    // Temporary storage space for boundary variables
    LevelData<FArrayBox> m_B;

    // interpolator for filling in ghost cells from the next coarser level
    PiecewiseLinearFillPatch m_patcher;

    // interpolator for filling in ghost cells on the boundary from the next coarser level
    D1PiecewiseLinearFillPatch m_bndPatcher;

    // Grid spacing
    Real m_dx;

    // Problem domain - index space for this level
    ProblemDomain m_domain;

    // Boundary Box, face centered
    Box m_bdryFaceBox;

    // Refinement ratio between this level and the next coarser
    int m_refineCoarse;

    // Number of boundary variables
    int m_numBdryVars;

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


private:
    // Disallowed for all the usual reasons
    void operator=(const LinElastLevelGodunov& a_input)
    {
        MayDay::Error("LinElastLevelGodunov invalid operator");
    }

    // Disallowed for all the usual reasons
    LinElastLevelGodunov(const LinElastLevelGodunov& a_input)
    {
        MayDay::Error("LinElastLevelGodunov invalid operator");
    }
};

#include "NamespaceFooter.H"
#endif
