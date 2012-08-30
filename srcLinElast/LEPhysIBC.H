//
//  LEPhysIBC.H
// ============
//
// Extension of the PhysIBC class to include boundary data
//

#ifndef _LEPHYSIBC_H_
#define _LEPHYSIBC_H_

#define NUM_BDRY_HAT_VARS (9)
#define BX_VN             (0)
#define BX_VM             (1)
#define BX_VX             (2)
#define BX_SN             (3)
#define BX_SM             (4)
#define BX_SZ             (5)
#define BX_UN             (6)
#define BX_UM             (7)
#define BX_UZ             (8)

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "PhysIBC.H"
#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**
  Add boundary data to the PhysIBC
  */
class LEPhysIBC: public PhysIBC
{
public:

    /// Null Constructor
    /**
    */
    LEPhysIBC();

    /// Set up initial conditions for the boundary data
    /**
    */
    virtual void initializeBdry(LevelData<FArrayBox>& a_Psi,int a_side)
    {
        // MayDay::Error("No boundary data for this problem");
    }

    /// Allow the user to query whether there is boundary data
    /**
    */
    virtual bool hasBdryData()
    {
        return (0 == 1);
    }

    /// Store pointer to boundary data
    /**
    */
    virtual void setBdryData(FArrayBox* a_bdryData, int a_side)
    {
      m_bdryData[a_side]    = a_bdryData;
      m_tmpBdryData[a_side] = new FArrayBox(m_bdryData[a_side]->box(), m_numBdryVars[a_side]);
      m_tmpBdryDataSet[a_side] = false;
    }

    /// Store pointer to boundary data
    /**
    */
    virtual void updateBoundary(const FArrayBox& a_WHalf,int dir,const Real& a_dt,const Real& a_dx,const Real& a_time,const bool a_final)
    {
        // Default implementation is to do nothing
    }

    /// Do the initial tagging of cells
    /**
    */
    virtual bool tagCellsInit(FArrayBox& markFAB,const Real& threshold)
    {
        // Default implementation is to do nothing
        return false;
    }

    /// Do the physics dependent tagging of cells
    /**
    */
    virtual bool tagCells(FArrayBox& markFAB, Real& time,const Real& threshold)
    {
        return false;
    }

    virtual Vector<int> numBdryVars()
    {
        return m_numBdryVars;
    }

    virtual Vector<string> bdryNames(int a_side)
    {
      Vector<string> bdryNames;
      // MUST MATCH LEPhysIBC.H ORDER!!!
      bdryNames.push_back("vn");
      bdryNames.push_back("vm");
      bdryNames.push_back("vz");
      bdryNames.push_back("sn");
      bdryNames.push_back("sm");
      bdryNames.push_back("sz");
      bdryNames.push_back("un");
      bdryNames.push_back("um");
      bdryNames.push_back("uz");
      return bdryNames;
    }

    virtual void dumpBdryData(FILE * a_boundaryDataFile)
    {
    }

    void setFortranCommonLE(const Real& a_cs, const Real& a_cp, const Real& a_mu,
        const Real& a_sxx0,const Real& a_syy0,const Real& a_szz0,
        const Real& a_sxy0,const Real& a_sxz0,const Real& a_syz0);

    void setFortranCommonPlastic(const Real& a_mu, const Real& a_beta, const Real& a_eta);

protected:
    FArrayBox* m_bdryData[2*CH_SPACEDIM];
    FArrayBox* m_tmpBdryData[2*CH_SPACEDIM];
    Vector<int> m_numBdryVars;
    bool       m_tmpBdryDataSet[2*CH_SPACEDIM];
    bool       m_isFortranCommonLESet;
    bool       m_isFortranCommonPlasticSet;
private:
};

#include "NamespaceFooter.H"
#endif
