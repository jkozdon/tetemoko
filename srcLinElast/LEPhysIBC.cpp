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
    const vector<Real> a_back)
{
    FORT_LINELASTSETF(CHF_CONST_REAL(a_cs),CHF_CONST_REAL(a_cp),CHF_CONST_REAL(a_mu),CHF_CONST_VR(a_back));
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
