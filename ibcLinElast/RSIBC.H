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
        const Real& x0,
        const Real& y0,
        const Real& sigma,
        const Real& ntime,
        const Real& psi,
        const Real& a,
        const Real& b,
        const Real& V0,
        const Real& f0,
        const Real& L,
        const Real& fw,
        const Real& Vw,
        const Real& fExp,
        const Real ruptureVelocityThreshold,
        const Vector<int>& boundaryType);

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

    /// Store pointer to boundary data
    /**
    */
    virtual void setBdryData(FArrayBox* a_bdryData);

    /// Set boundary primitive values.
    /**
    */
    virtual void primBC(FArrayBox& a_WGdnv,
        const FArrayBox&      a_WshiftInside,
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

    /// Store pointer to boundary data
    /**
    */
    virtual void updateBoundary(const FArrayBox& a_WHalf,int dir,const Real& a_dt,const Real& a_time,const bool a_final);

    /// Do the initial tagging of cells
    /**
    */
    virtual bool tagCellsInit(FArrayBox& markFAB);

    /// Do the physics dependent tagging of cells
    /**
    */
    virtual bool tagCells(FArrayBox& markFAB, Real& time)
    {
        return tagCellsInit(markFAB);
    }

    virtual int numBdryVars()
    {
        return m_numBdryVars;
    }

    virtual Vector<string> bdryNames()
    {
        Vector<string> bdryNames;
        bdryNames.push_back("Vx");
        bdryNames.push_back("Vz");
        bdryNames.push_back("sxy");
        bdryNames.push_back("syz");
        bdryNames.push_back("V");
        bdryNames.push_back("slip_x");
        bdryNames.push_back("slip_z");
        bdryNames.push_back("slip_total");
        bdryNames.push_back("rupture");
        bdryNames.push_back("syy");
        bdryNames.push_back("psi");
        return bdryNames;
    }

protected:
    bool m_isFortranCommonSet;
    Real m_psi;
    Vector<int> m_boundaryType;

    bool m_tmpBdryDataSet;
    FArrayBox* m_tmpBdryData;
private:
};

#endif
