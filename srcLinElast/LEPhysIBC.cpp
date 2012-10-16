#include "LEPhysIBC.H"
#include "LinElastPhysicsF_F.H"
#include "PlasticPhysicsF_F.H"

/// Null Constructor
LEPhysIBC::LEPhysIBC()
{
    m_isFortranCommonLESet = false;
    m_isFortranCommonPlasticSet = false;
}

/// Set the fortran parameters for linear elasticity
void LEPhysIBC::setFortranCommonLE(const Real& a_cs,
    const Real& a_cp,
    const Real& a_mu,
    const Real& a_sxx0,
    const Real& a_syy0,
    const Real& a_szz0,
    const Real& a_sxy0,
    const Real& a_sxz0,
    const Real& a_syz0)
{
    FORT_LINELASTSETF(CHF_CONST_REAL(a_cs),CHF_CONST_REAL(a_cp),CHF_CONST_REAL(a_mu),
	   CHF_CONST_REAL(a_sxx0),CHF_CONST_REAL(a_syy0),CHF_CONST_REAL(a_szz0),
	   CHF_CONST_REAL(a_sxy0),CHF_CONST_REAL(a_sxz0),CHF_CONST_REAL(a_syz0));
    m_isFortranCommonLESet = true;
}

/// Set the fortran parameters for plasticity
void LEPhysIBC::setFortranCommonPlastic(const Real& a_mu,
    const Real& a_beta,
    const Real& a_eta)
{
    FORT_PLASTICSETF(CHF_CONST_REAL(a_mu),CHF_CONST_REAL(a_beta),CHF_CONST_REAL(a_eta));
    m_isFortranCommonPlasticSet = true;
}
