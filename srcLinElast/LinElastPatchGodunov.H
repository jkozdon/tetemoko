#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LINELASTPATCHGODUNOV_H_
#define _LINELASTPATCHGODUNOV_H_

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "FluxBox.H"

#include "LinElastPhysics.H"
#include "GodunovUtilities.H"
#include "NamespaceHeader.H"

class HDF5HeaderData;

///
/**
  The base class LinElastPatchGodunov provides an implementation of a
  second-order, unsplit Godunov method acting on a single grid/patch with
  extension to boundary data.  LinElastPatchGodunov provides an interface to
  the level integrator, LevelGodunov, which manages the entire level and flux
  corrections (via flux registers) to the coarser and finer levels.
  */
class LinElastPatchGodunov
{
public:
    /// Constructor
    /**
    */
    LinElastPatchGodunov();

    /// Destructor
    /**
    */
    virtual ~LinElastPatchGodunov();

    /// Define the object
    /**
    */
    virtual void define(const ProblemDomain&        a_domain,
        const Real&                 a_dx,
        const LinElastPhysics* const a_gdnvPhysicsPtr,
        const int&                  a_normalPredOrder,
        const bool&                 a_useFourthOrderSlopes,
        const bool&                 a_usePrimLimiting,
        const bool&                 a_useCharLimiting,
        const bool&                 a_useFlattening,
        const bool&                 a_useArtificialViscosity,
        const Real&                 a_artificialViscosity);

    /// Set the current time before calling updateState, computeWHalf, computeUpdate.
    /**
    */
    virtual void setCurrentTime(const Real& a_currentTime);

    /// Set the current box before calling updateState, computeWHalf, computeUpdate.
    /**
    */
    virtual void setCurrentBox(const Box& a_currentBox);

    /// Update the conserved variables and return the fluxes used to do this.
    /**
      Update the conserved variables and return the final fluxes that were used
      for this.  Compute the fluxes using a second-order, unsplit Godunov method
      based on the input conserved variables, a_U, and source terms, a_S.  If
      there are no source terms then a_S should be null constructed.  Also
      return the maximum wave speed.
      */
    virtual void updateState(FArrayBox&       a_U,
        FluxBox&         a_F,
        Real&            a_maxWaveSpeed,
        const FArrayBox& a_S,
        const Real&      a_dt,
        const Box&       a_box);

    virtual void plasticUpdate(FArrayBox& curU,
        const Real& a_dt,
        const Box&  curBox);

    virtual void updateState(FArrayBox&       a_U,
        FluxBox&         a_F,
        FluxBox&         a_wHalf,
        Real&            a_maxWaveSpeed,
        const FArrayBox& a_S,
        const Real&      a_dt,
        const Box&       a_box);

    /// Compute the time-centered values of the primitive variables on cell faces.
    /**
      This API is used in cases where some operation over the whole
      level must be performed on the face-centered variables prior to the
      final difference update. Examples include incompressible flow and MHD,
      in which it is necessary to compute the projection of a face-centered
      vector field on its divergence-free part. To complete the differencing,
      it is necessary to call the member function computeUpdate.
      */
    virtual void computeWHalf(FluxBox&         a_WHalf,
        const FArrayBox& a_U,
        const FArrayBox& a_S,
        const Real&      a_dt,
        const Box&       a_box);

    /// Compute the increment in the conserved variables from face variables.
    /**
      Compute dU = dt*dUdt, the change in the conserved variables over
      the time step. The fluxes are returned are suitable for use in refluxing.
      */
    virtual void computeUpdate(FArrayBox&       a_dU,
        FluxBox&         a_F,
        const FArrayBox& a_U,
        const FluxBox&   a_WHalf,
        const Real&      a_dt,
        const Box&       a_box);

    /// Compute fluxes from primitive variables on faces
    /**
    */
    void computeFluxes(FluxBox&        a_F,
        const FluxBox&  a_WHalf,
        const Box&      a_box);

    /// CTU normal predictor
    /**
      Compute the increments in the characteristic amplitudes using CTU
      (for CTU, increments are zero)
      */
    void CTUNormalPred(FArrayBox&       a_WMinus,
        FArrayBox&       a_WPlus,
        const Real&      a_dt,
        const Real&      a_dx,
        const FArrayBox& a_W,
        const FArrayBox& a_flat,
        const int&       a_dir,
        const Box&       a_box);

    /// PLM normal predictor
    /**
      Compute the increments in the characteristic amplitudes using PLM
      */
    void PLMNormalPred(FArrayBox&       a_WMinus,
        FArrayBox&       a_WPlus,
        const Real&      a_dt,
        const Real&      a_dx,
        const FArrayBox& a_W,
        const FArrayBox& a_flat,
        const int&       a_dir,
        const Box&       a_box);

    /// PPM normal predictor
    /**
      Compute the increments in the characteristic amplitudes using PPM
      */
    void PPMNormalPred(FArrayBox&       a_WMinus,
        FArrayBox&       a_WPlus,
        const Real&      a_dt,
        const Real&      a_dx,
        const FArrayBox& a_W,
        const FArrayBox& a_flat,
        const int&       a_dir,
        const Box&       a_box);

    /// Accessor to get at the GodunovPhysics object
    /**
    */
    GodunovPhysics* getGodunovPhysicsPtr();

    /// Accessor to get at the GodunovPhysics object
    /**
    */
    LinElastPhysics* getLinElastPhysicsPtr();

    /// Is the object completely defined
    /**
      Return true if the object is completely defined.
      */
    virtual bool isDefined() const;

    /// Set whether to use high-order limiter.
    void highOrderLimiter(bool a_highOrderLimiter);

#ifdef CH_USE_HDF5
    virtual void expressions(HDF5HeaderData& a_holder) const {;}
#endif

protected:
    // Problem domain and grid spacing
    ProblemDomain m_domain;
    Real          m_dx;

    // Object providing physics information
    LinElastPhysics*  m_linElastPhysics;

    // Object containing various methods for the Godunov computation
    GodunovUtilities m_util;

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

    // Current time and has it been set
    Real m_currentTime;
    bool m_isCurrentTimeSet;

    // Current box and has it been set
    Box  m_currentBox;
    bool m_isCurrentBoxSet;

    // Has this object been defined
    bool m_isDefined;

    // Use a high-order limiter?  (default false)
    bool m_highOrderLimiter;

private:
    // Disallowed for all the usual reasons
    void operator=(const LinElastPatchGodunov& a_input)
    {
        MayDay::Error("invalid operator");
    }

    // Disallowed for all the usual reasons
    LinElastPatchGodunov(const LinElastPatchGodunov& a_input)
    {
        MayDay::Error("invalid operator");
    }
};

#include "NamespaceFooter.H"
#endif
