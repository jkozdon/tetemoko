#ifndef _LINELASTPHYSICS_H_
#define _LINELASTPHYSICS_H_

#include <string>
using std::string;

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "CH_HDF5.H"

#include "GodunovPhysics.H"

#include "UsingNamespace.H"

///
/**
  The base class LinElastPhysics provides the physics-dependent components
  for a higher-order Godunov method for a single patch: characteristic
  analysis, Riemann solver, quasilinear update, conservative update,
  and transformations between conserved, primitive, and flux variables.
  This class is essentially pure virtual, i.e. all of its member functions
  are virtual, and the ones that have default implementations are ones
  which are optionally defined, i.e. the default definition is to send
  an error message. Physics-dependent versions of this class that are
  required in real applications are derived from this class by inheritance.
  */
class LinElastPhysics: public GodunovPhysics
{
public:
    /// Constructor
    /**
    */
    LinElastPhysics(const int junk);

    /// Destructor
    /**
    */
    ~LinElastPhysics();

    /// Compute the maximum wave speed
    /**
    */
    Real getMaxWaveSpeed(const FArrayBox& a_U,
        const Box&       a_box);

    /// Object factory for this class
    /**
    */
    virtual LinElastPhysics* new_godunovPhysics() const;

    /// Number of conserved variables
    /**
      Return the number of conserved variables.
      */
    int numConserved();

    /// Names of the conserved variables
    /**
      Return the names of the conserved variables.
      */
    Vector<string> stateNames();

    /// Number of flux variables
    /**
      Return the  number of flux variables.  This can be greater than the number
      of conserved variables if addition fluxes/face-centered quantities are
      computed.
      */
    int numFluxes();

    /// Compute the increment in the conserved variables from face variables.
    /**
      Compute dU = dt*dUdt, the change in the conserved variables over
      the time step. The fluxes are returned are suitable for use in refluxing.
      This has a default implementation but can be redefined as needed.
      */
    void computeUpdate(FArrayBox&       a_dU,
        FluxBox&         a_F,
        const FArrayBox& a_U,
        const FluxBox&   a_WHalf,
        const bool&      a_useArtificialViscosity,
        const Real&      a_artificialViscosity,
        const Real&      a_currentTime,
        const Real&      a_dx,
        const Real&      a_dt,
        const Box&       a_box);

    /// Compute a flux from primitive variable values on a face
    /**
    */
    void getFlux(FArrayBox&       a_flux,
        const FArrayBox& a_whalf,
        const int&       a_dir,
        const Box&       a_box);

    /// Number of primitive variables
    /**
      Return the number of primitive variables.  This may be greater than the
      number of conserved variables if derived/redundant quantities are also
      stored for convenience.
      */
    int numPrimitives();

    /// Transform a_dW from primitive to characteristic variables
    /**
      On input, a_dW contains the increments of the primitive variables. On
      output, it contains the increments in the characteristic variables.

      IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
      smallest eigenvalue first, the largest eigenvalue last, and orders the
      characteristic variables accordingly.
      */
    void charAnalysis(FArrayBox&       a_dW,
        const FArrayBox& a_W,
        const int&       a_dir,
        const Box&       a_box);

    /// Transform a_dW from characteristic to primitive variables
    /**
      On input, a_dW contains the increments of the characteristic variables.
      On output, it contains the increments in the primitive variables.

      IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
      smallest eigenvalue first, the largest eigenvalue last, and orders the
      characteristic variables accordingly.
      */
    void charSynthesis(FArrayBox&       a_dW,
        const FArrayBox& a_W,
        const int&       a_dir,
        const Box&       a_box);

    /// Compute the characteristic values (eigenvalues)
    /**
      Compute the characteristic values (eigenvalues)

      IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
      smallest eigenvalue first, the largest eigenvalue last, and orders the
      characteristic variables accordingly.
      */
    void charValues(FArrayBox&       a_lambda,
        const FArrayBox& a_W,
        const int&       a_dir,
        const Box&       a_box);

    /// Add to (increment) the source terms given the current state
    /**
      On input, a_S contains the current source terms.  On output, a_S has
      had any additional source terms (based on the current state, a_W)
      added to it.  This should all be done on the region defined by a_box.
      */
    void incrementSource(FArrayBox&       a_S,
        const FArrayBox& a_W,
        const Box&       a_box);

    /// Compute the solution to the Riemann problem.
    /**
      Given input left and right states in a direction, a_dir, compute a
      Riemann problem and generate fluxes at the faces within a_box.
      */
    void riemann(/// face-centered solution to Riemann problem
        FArrayBox&       a_WStar,
        /// left state, on cells to left of each face
        const FArrayBox& a_WLeft,
        /// right state, on cells to right of each face
        const FArrayBox& a_WRight,
        /// state on cells, used to set boundary conditions
        const FArrayBox& a_W,
        /// current time
        const Real&      a_time,
        /// direction of faces
        const int&       a_dir,
        /// face-centered box on which to set a_WStar
        const Box&       a_box);

    /// Post-normal predictor calculation.
    /**
      Add increment to normal predictor, e.g. to account for source terms due to
      spatially-varying coefficients, to bound primitive variable ranges.
      */
    virtual void postNormalPred(FArrayBox&       a_dWMinus,
        FArrayBox&       a_dWPlus,
        const FArrayBox& a_W,
        const Real&      a_dt,
        const Real&      a_dx,
        const int&       a_dir,
        const Box&       a_box);

    /// Compute a temporary update to the boundary data
    /**
    */
    void quasilinearBoundaryUpdate(
        const FArrayBox& a_wHalf,
        const Real&      a_dt,
        const int&       a_dir,
        const Box&       a_box,
        const Real&      a_currentTime);

    /// Compute the quasilinear update A*dW/dx.
    /**
    */
    void quasilinearUpdate(FArrayBox&       a_AdWdx,
        const FArrayBox& a_wHalf,
        const FArrayBox& a_W,
        const Real&      a_scale,
        const int&       a_dir,
        const Box&       a_box);

    /// Compute the plastic correction
    /**
    */
    void plasticUpdate(FArrayBox& a_U,
        const Real&      a_dt,
        const Box&       a_box);

    /// Compute primitive variables from conserved variables.
    /**
    */
    void consToPrim(FArrayBox&       a_W,
        const FArrayBox& a_U,
        const Box&       a_box);

    /// Interval within the primitive variables corresponding to the velocities
    /**
      Return the interval of component indices within the primitive variable
      of the velocities.  Used for slope flattening (slope computation) and
      computing the divergence of the velocity (artificial viscosity).
      */
    virtual Interval velocityInterval();

    /// Component index within the primitive variables of the pressure
    /**
      Return the component index withn the primitive variables for the
      pressure.  Used for slope flattening (slope computation).
      */
    virtual int pressureIndex();

    /// Used to limit the absolute value of a "pressure" difference
    /**
      Return a value that is used by slope flattening to limit (away from
      zero) the absolute value of a slope in the pressureIndex() component
      (slope computation).
      */
    virtual Real smallPressure();

    /// Component index within the primitive variables of the bulk modulus
    /**
      Return the component index withn the primitive variables for the
      bulk modulus.  Used for slope flattening (slope computation) used
      as a normalization to measure shock strength.
      */
    virtual int bulkModulusIndex();

    /// Compute the artificial viscosity contribution to the flux
    /**
      Compute the artificial viscosity contribution to the flux.  This has
      a default implementation but this can be overridded as needed.
      */
    virtual void artVisc(FArrayBox&       a_F,
        const FArrayBox& a_U,
        const Real&      a_artificialViscosity,
        const Real&      a_currentTime,
        const int&       a_dir,
        const Box&       a_box);

#ifdef CH_USE_HDF5
    virtual void expressions(HDF5HeaderData& a_holder) const;
#endif


protected:

private:
    // Disallowed for all the usual reasons
    void operator=(const LinElastPhysics& a_input)
    {
        MayDay::Error("invalid operator");
    }

    // Disallowed for all the usual reasons
    LinElastPhysics(const LinElastPhysics& a_input)
    {
        MayDay::Error("invalid operator");
    }
};

#endif
