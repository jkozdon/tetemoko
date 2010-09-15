#ifndef _RSIBC_H_
#define _RSIBC_H_

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "PhysIBC.H"
#include "LEPhysIBC.H"

#include "UsingNamespace.H"

class RSIBC: public LEPhysIBC
{
public:
    /// Null Constructor
    /**
    */
    RSIBC();

    /// Constructor which defines parameters used by Fortran routines
    /**
    */
    RSIBC(const Real& cs,
        const Real& cp,
        const Real& mu,
        const vector<Real> background,
        const Real& r0,
        const Real& sigma,
        const Real& ntime,
        const Real& psi,
        const Real& a,
        const Real& b,
        const Real& V0,
        const Real& f0,
        const Real& L,
        const Real& fw,
        const Real& Vw);

    /// Destructor
    /**
    */
    virtual ~RSIBC();

    /// Define ?

    /// Factory method - this object is its own factory
    /**
      Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
      its define() must be called before it is used) and m_isFortranCommonSet
      set to value of m_isFortranCommonset in the current (factory) object.
      */
    PhysIBC *new_physIBC();

    void setFortranCommonSet();

    /// Set up initial conditions
    /**
    */
    virtual void initialize(LevelData<FArrayBox>& a_U);

    /// Set up initial boundary data
    /**
    */
    virtual void initializeBdry(LevelData<FArrayBox>& a_U);

    /// Set up initial boundary data
    /**
    */
    virtual bool hasBdryData();

    /// Set boundary primitive values.
    /**
    */
    virtual void primBC(FArrayBox& a_WGdnv,
        const FArrayBox&      a_WshiftInside,
        const FArrayBox&      a_W,
        const int&            a_dir,
        const Side::LoHiSide& a_side,
        const Real&           a_time);

    /// Set boundary primitive values.
    /**
    */
    virtual void primBC(FArrayBox& a_WGdnv,
        const FArrayBox&      a_WshiftInside,
        const FArrayBox&      a_W,
        const FArrayBox&      a_Psi,
        const int&            a_dir,
        const Side::LoHiSide& a_side,
        const Real&           a_time);

    /// Set boundary slopes
    /**
      The boundary slopes in a_dW are already set to one sided difference
      approximations.  If this function doesn't change them they will be
      used for the slopes at the boundaries.
      */
    virtual
        void setBdrySlopes(FArrayBox&       a_dW,
            const FArrayBox& a_W,
            const int&       a_dir,
            const Real&      a_time);

    /// Adjust boundary fluxes to account for artificial viscosity
    /**
    */
    virtual
        void artViscBC(FArrayBox&       a_F,
            const FArrayBox& a_U,
            const FArrayBox& a_divVel,
            const int&       a_dir,
            const Real&      a_time);


protected:
    bool m_isFortranCommonSet;
    Real m_psi;
    // Vector<int> m_boundaryType;

private:
};

#endif
