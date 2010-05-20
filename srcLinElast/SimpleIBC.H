#ifndef _SIMPLEIBC_H_
#define _SIMPLEIBC_H_

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "PhysIBC.H"

#include "UsingNamespace.H"

class SimpleIBC: public PhysIBC
{
public:
    /// Null Constructor
    /**
    */
    SimpleIBC();

    /// Constructor which defines parameters used by Fortran routines
    /**
    */
    SimpleIBC(const Real& cs,
        const Real& cp,
        const Real& mu,
        const Real& r0,
        const Real& mag,
        const Real& sig);

    /// Destructor
    /**
    */
    virtual ~SimpleIBC();

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

    /// Set boundary primitive values.
    /**
    */
    virtual void primBC(FArrayBox&            a_WGdnv,
        const FArrayBox&      a_Wextrap,
        const FArrayBox&      a_W,
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
};

#endif
