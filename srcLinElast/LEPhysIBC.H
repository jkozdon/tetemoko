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
    virtual void updateBoundary(const FArrayBox& a_WHalf,int dir,const Real& a_dt)
    {
        // Default implementation is to do nothing
    }

    /// Do the initial tagging of cells
    /**
    */
    virtual bool tagCellsInit(FArrayBox& markFAB)
    {
        // Default implementation is to do nothing
        return false;
    }

    /// Do the physics dependent tagging of cells
    /**
    */
    virtual bool tagCells(FArrayBox& markFAB, Real& time)
    {
        return false;
    }

protected:
    FArrayBox* m_bdryData;
    int        m_numBdryVars;
private:
};

#include "NamespaceFooter.H"
#endif
