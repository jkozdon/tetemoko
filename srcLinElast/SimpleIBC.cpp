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

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_physIBC() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void SimpleIBC::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

/// Factory method - this object is its own factory
PhysIBC* SimpleIBC::new_physIBC()
{
    SimpleIBC* retval = new SimpleIBC();
    if(m_isFortranCommonSet == true)
    {
        retval->setFortranCommonSet();
    }
    return static_cast<PhysIBC*>(retval);
}

/// Set up initial conditions
void SimpleIBC::initialize(LevelData<FArrayBox>& a_U)
{
    pout() << "NOT SETUP :: SimpleIBC::initialize" << endl;
    CH_assert(m_isFortranCommonSet == true);
    CH_assert(m_isDefined == true);

    // Iterator of all grids in this level
    for (DataIterator dit = a_U.dataIterator();
        dit.ok(); ++dit)
    {
        // Storage for current grid
        FArrayBox& U = a_U[dit()];

        // Box of current grid
        Box uBox = U.box();
        uBox &= m_domain;

        // Set up initial condition in this grid
        //JK FORT_RAMPINITF(CHF_FRA(U),
        //JK     CHF_CONST_REAL(m_dx),
        //JK     CHF_BOX(uBox));
    }
}

/// Set boundary primitive values.
void SimpleIBC::primBC(FArrayBox&            a_WGdnv,
    const FArrayBox&      a_Wextrap,
    const FArrayBox&      a_W,
    const int&            a_dir,
    const Side::LoHiSide& a_side,
    const Real&           a_time)
{
    pout() << "NOT SETUP :: SimpleIBC::primBC" << endl;
}

/// Set boundary slopes

void SimpleIBC::setBdrySlopes(FArrayBox&       a_dW,
    const FArrayBox& a_W,
    const int&       a_dir,
    const Real&      a_time)
{
    pout() << "NOT SETUP :: SimpleIBC::setBdrySlopes" << endl;
}

void SimpleIBC::artViscBC(FArrayBox&       a_F,
    const FArrayBox& a_U,
    const FArrayBox& a_divVel,
    const int&       a_dir,
    const Real&      a_time)
{
    pout() << "NOT SETUP :: SimpleIBC::artViscBC" << endl;
}
