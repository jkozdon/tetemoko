#ifndef _SWIBC_H_
#define _SWIBC_H_

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "PhysIBC.H"
#include "LEPhysIBC.H"

#include "UsingNamespace.H"

class SWIBC: public LEPhysIBC
{
public:
    /// Null Constructor
    /**
    */
    SWIBC();

    /// Constructor which defines parameters used by Fortran routines
    /**
    */
    SWIBC(const Real& fricS,
        const Real& fricD,
        const Real& weakD,
        const Real& smoothValue,
        const int numPatches,
        const vector<Real> xloPatches,
        const vector<Real> xwPatches,
        const vector<Real> zloPatches,
        const vector<Real> zwPatches,
        const vector<Real> tauPatches,
        const vector<Real> fricBoxCenter,
        const vector<Real> fricBoxWidth,
        const Real outsideFriction,
        const Real ruptureVelocityThreshold,
        const Vector<int>& boundaryType);

    /// Destructor
    /**
    */
    virtual ~SWIBC();

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

    virtual void initializeBdry(LevelData<FArrayBox>& a_B,int a_side);

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
    virtual void updateBoundary(const FArrayBox& a_WHalf,int dir,const Real& a_dt,const Real& a_dx,const Real& a_time,const bool a_final);

    /// Do the initial tagging of cells
    /**
    */
    virtual bool tagCellsInit(FArrayBox& markFAB,const Real& threshold);

    /// Do the physics dependent tagging of cells
    /**
    */
    virtual bool tagCells(FArrayBox& markFAB, Real& time,const Real& threshold)
    {
        return tagCellsInit(markFAB,threshold);
    }

    virtual Vector<string> bdryNames(int a_side)
    {
      Vector<string> bdryNames;
      if(a_side == 2)
      {
        // MUST MATCH SWINDEX.H ORDER!!!
        bdryNames.push_back("Vx");
        bdryNames.push_back("Vz");
        bdryNames.push_back("sxy");
        bdryNames.push_back("syz");
        bdryNames.push_back("V");
        bdryNames.push_back("slip_x");
        bdryNames.push_back("slip_z");
        bdryNames.push_back("slip_total");
        bdryNames.push_back("rupture");
      }
      else
      {
        // WOULD BE BETTER TO CALL ORIGINAL FUNCTION!

        // MUST MATCH LEPhysIBC.H ORDER!!!
        bdryNames.push_back("vx_1");
        bdryNames.push_back("vy_1");
        bdryNames.push_back("vz_1");
        bdryNames.push_back("Tx_1");
        bdryNames.push_back("Ty_1");
        bdryNames.push_back("Tz_1");
        bdryNames.push_back("un_1");
        bdryNames.push_back("um_1");
        bdryNames.push_back("uz_1");

        bdryNames.push_back("vx_2");
        bdryNames.push_back("vy_2");
        bdryNames.push_back("vz_2");
        bdryNames.push_back("Tx_2");
        bdryNames.push_back("Ty_2");
        bdryNames.push_back("Tz_2");
        bdryNames.push_back("un_2");
        bdryNames.push_back("um_2");
        bdryNames.push_back("uz_2");
      }
      return bdryNames;
    }

    virtual void dumpBdryData(FILE * a_boundaryDataFile);

protected:
    bool m_isFortranCommonSet;
    Vector<int> m_boundaryType;
    bool m_isPatchBoxSet;

    int m_numPatches;
    Vector<Real> m_xcPatches;
    Vector<Real> m_xwPatches;
    Vector<Real> m_zcPatches;
    Vector<Real> m_zwPatches;
    Vector<Real> m_tauPatches;

    Real m_outsideFriction;
    Real m_smoothValue;
    int  m_smoothWidthNumCells;
    Vector<Real> m_fricBoxCenter;
    Vector<Real> m_fricBoxWidth;
    Vector<Box>  m_patchBoxes;
private:
};

#endif
