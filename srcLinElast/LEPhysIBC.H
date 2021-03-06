//
//  LEPhysIBC.H
// ============
//
// Extension of the PhysIBC class to include boundary data
//

#ifndef _LEPHYSIBC_H_
#define _LEPHYSIBC_H_

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
    virtual void initializeBdry(LevelData<FArrayBox>& a_Psi)
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
    virtual void setBdryData(FArrayBox* a_bdryData)
    {
        m_bdryData = a_bdryData;
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

    virtual int numBdryVars()
    {
        return 1;
    }

    virtual Vector<string> bdryNames()
    {
        Vector<string> junk;
        junk.push_back("junk");
        return junk;
    }

    virtual void dumpBdryData(FILE * a_boundaryDataFile)
    {
    }

    void setFortranCommonLE(const Real& a_cs, const Real& a_cp, const Real& a_mu,
        const Real& a_sxx0,const Real& a_syy0,const Real& a_szz0,
        const Real& a_sxy0,const Real& a_sxz0,const Real& a_syz0);

    void setFortranCommonPlastic(const Real& a_mu, const Real& a_beta, const Real& a_eta);

protected:
    FArrayBox* m_bdryData;
    int        m_numBdryVars;
    bool       m_isFortranCommonLESet;
    bool       m_isFortranCommonPlasticSet;
private:
};

#include "NamespaceFooter.H"
#endif
