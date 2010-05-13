#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "SimpleIBC.H"
#include "SimpleIBCF_F.H"

/// Null Constructor
SimpleIBC::SimpleIBC()
{
}

/// Constructor which defines parameters used by Fortran routines
SimpleIBC::SimpleIBC(const Real& a_cs,
    const Real& a_cp,
    const Real& a_mu)
{
    FORT_SIMPLESETF(CHF_CONST_REAL(a_cs),CHF_CONST_REAL(a_cp),CHF_CONST_REAL(a_mu));
    m_isFortranCommonSet = true;
}

/// Destructor
SimpleIBC::~SimpleIBC()
{
}

/// Define ?

/// Factory method - this object is its own factory
PhysIBC* SimpleIBC::new_physIBC()
{
    return NULL;
}

/// Set up initial conditions
void SimpleIBC::initialize(LevelData<FArrayBox>& a_U)
{
}

/// Set boundary primitive values.
void SimpleIBC::primBC(FArrayBox&            a_WGdnv,
    const FArrayBox&      a_Wextrap,
    const FArrayBox&      a_W,
    const int&            a_dir,
    const Side::LoHiSide& a_side,
    const Real&           a_time)
{
}

/// Set boundary slopes

void SimpleIBC::setBdrySlopes(FArrayBox&       a_dW,
    const FArrayBox& a_W,
    const int&       a_dir,
    const Real&      a_time)
{
}

void SimpleIBC::artViscBC(FArrayBox&       a_F,
    const FArrayBox& a_U,
    const FArrayBox& a_divVel,
    const int&       a_dir,
    const Real&      a_time)
{
}
