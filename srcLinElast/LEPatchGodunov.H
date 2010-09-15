#ifndef _LEPATCHGODUNOV_H_
#define _LEPATCHGODUNOV_H_

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "FluxBox.H"
#include "PatchGodunov.H"

#include "GodunovPhysics.H"
#include "LinElastPhysics.H"
#include "GodunovUtilities.H"
#include "NamespaceHeader.H"

class HDF5HeaderData;

///
/**
  The base class PatchGodunov provides an implementation of a second-order,
  unsplit Godunov method acting on a single grid/patch.  PatchGodunov
  provides an interface to the level integrator, LevelGodunov, which manages
  the entire level and flux corrections (via flux registers) to the coarser
  and finer levels.
  */
class LEPatchGodunov: public PatchGodunov
{
public:
    /// Define the object
    /**
     * Error because user must use the new constructor with LinElastPhysics.
     */
    virtual void define(const ProblemDomain&        a_domain,
        const Real&                 a_dx,
        const GodunovPhysics* const a_gdnvPhysicsPtr,
        const int&                  a_normalPredOrder,
        const bool&                 a_useFourthOrderSlopes,
        const bool&                 a_usePrimLimiting,
        const bool&                 a_useCharLimiting,
        const bool&                 a_useFlattening,
        const bool&                 a_useArtificialViscosity,
        const Real&                 a_artificialViscosity);

    /// Define the object
    /**
     * Error because user must use the new constructor with LinElastPhysics.
     */
    virtual void define(const ProblemDomain&        a_domain,
        const Real&                 a_dx,
        const LinElastPhysics* const a_linElastPhysicsPtr,
        const int&                  a_normalPredOrder,
        const bool&                 a_useFourthOrderSlopes,
        const bool&                 a_usePrimLimiting,
        const bool&                 a_useCharLimiting,
        const bool&                 a_useFlattening,
        const bool&                 a_useArtificialViscosity,
        const Real&                 a_artificialViscosity);

    void updateStateBdry(FArrayBox&       a_U,
        FArrayBox&       a_Psi,
        FluxBox&         a_F,
        Real&            a_maxWaveSpeed,
        const FArrayBox& a_S,
        const Real&      a_dt,
        const Box&       a_box,
        const Box&       a_bdryBox);

    /// Compute the time-centered values of the primitive variables on cell faces.
    /**
      This API is used in cases where some operation over the whole
      level must be performed on the face-centered variables prior to the
      final difference update. Examples include incompressible flow and MHD,
      in which it is necessary to compute the projection of a face-centered
      vector field on its divergence-free part. To complete the differencing,
      it is necessary to call the member function computeUpdate.
      */
    void computeWHalfBdry(FluxBox&         a_WHalf,
        const FArrayBox& a_U,
        const FArrayBox& a_Psi,
        const FArrayBox& a_S,
        const Real&      a_dt,
        const Box&       a_box,
        const Box&       a_bdryBox);
protected:

private:
};

#include "NamespaceFooter.H"
#endif

