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

    /// Set boundary primitive values.
    /**
    */
    virtual void primBC(FArrayBox& a_WGdnv,
        const FArrayBox&      a_WshiftInside,
        const FArrayBox&      a_W,
        const FArrayBox&      a_Psi,
        const int&            a_dir,
        const Side::LoHiSide& a_side,
        const Real&           a_time)
    {
        MayDay::Error("No boundary data for this problem");
    }


    /// Set up initial conditions for the boundary data
    /**
    */
    virtual void initializeBdry(LevelData<FArrayBox>& a_Psi)
    {
        MayDay::Error("No boundary data for this problem");
    }

    /// Allow the user to query whether there is boundary data
    /**
    */
    virtual bool hasBdryData()
    {
        return (0 == 1);
    }

protected:
private:
};

#include "NamespaceFooter.H"
#endif
