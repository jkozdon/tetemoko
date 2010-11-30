#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "SWIBC.H"
#include "SWIBCF_F.H"
#include "LinElastPhysicsF_F.H"

/// Null Constructor
SWIBC::SWIBC()
{
}

/// Constructor which defines parameters used by Fortran routines
SWIBC::SWIBC(const Real& a_cs,
    const Real& a_cp,
    const Real& a_mu,
    const vector<Real> a_back,
    const Real& a_fricS,
    const Real& a_fricD,
    const Real& a_weakD,
    const Real& a_tau_nuc,
    const Real& a_width,
    const vector<Real> a_nucPatch,
    const int a_numPatches,
    const vector<Real> a_xcPatches,
    const vector<Real> a_xwPatches,
    const vector<Real> a_zcPatches,
    const vector<Real> a_zwPatches,
    const vector<Real> a_tauPatches,
    const Vector<int>& a_boundaryType)
{
    FORT_LINELASTSETF(CHF_CONST_REAL(a_cs),CHF_CONST_REAL(a_cp),CHF_CONST_REAL(a_mu),CHF_CONST_VR(a_back));
    FORT_SWSETF(CHF_CONST_REAL(a_fricS),CHF_CONST_REAL(a_fricD),CHF_CONST_REAL(a_weakD),CHF_CONST_REAL(a_tau_nuc),
        CHF_CONST_VR(a_nucPatch),CHF_CONST_REAL(a_width));
    m_nucPatch           = a_nucPatch;
    m_boundaryType       = a_boundaryType;
    m_isFortranCommonSet = true;
    m_isNucBoxSet = false;

    m_numPatches = a_numPatches;
    m_xcPatches = a_xcPatches;
    m_xwPatches = a_xwPatches;
    m_zcPatches = a_zcPatches;
    m_zwPatches = a_zwPatches;
    m_tauPatches = a_tauPatches;
}

/// Destructor
SWIBC::~SWIBC()
{
}

/// Define ?

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_physIBC() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void SWIBC::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

/// Factory method - this object is its own factory
PhysIBC* SWIBC::new_physIBC()
{
    SWIBC* retval = new SWIBC();
    if(m_isFortranCommonSet == true)
    {
        retval->setFortranCommonSet();
    }
    retval->m_boundaryType = m_boundaryType;
    retval->m_nucPatch = m_nucPatch;
    retval->m_isNucBoxSet = m_isNucBoxSet;
    retval->m_nucBox = m_nucBox;
    retval->m_numPatches = m_numPatches;
    retval->m_xcPatches  = m_xcPatches;
    retval->m_xwPatches  = m_xwPatches;
    retval->m_zcPatches  = m_zcPatches;
    retval->m_zwPatches  = m_zwPatches;
    retval->m_tauPatches = m_tauPatches;
    return static_cast<PhysIBC*>(retval);
}

/// Set up initial conditions
void SWIBC::initialize(LevelData<FArrayBox>& a_U)
{
    // pout() << "SWIBC::initialize" << endl;
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
        FORT_SWINITF(CHF_FRA(U),
            CHF_CONST_REAL(m_dx),
            CHF_BOX(uBox));
    }
}

void SWIBC::initializeBdry(LevelData<FArrayBox>& a_B)
{
    const Real tmpVal = 0.0;
    for (DataIterator dit = a_B.dataIterator(); dit.ok(); ++dit)
    {
        // Storage for current grid
        FArrayBox& B = a_B[dit()];

        // Box of current grid
        Box bBox = B.box();
        bBox &= m_domain;

        // Set up initial condition in this grid
        FORT_LINELASTSETFAB(CHF_FRA1(B,2),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
    }
}

bool SWIBC::hasBdryData()
{
    return false;
}

/// Set boundary primitive values.
void SWIBC::primBC(FArrayBox& a_WGdnv,
    const FArrayBox&      a_WShiftInside,
    const FArrayBox&      a_W,
    const int&            a_dir,
    const Side::LoHiSide& a_side,
    const Real&           a_time)
{
    CH_assert(m_isFortranCommonSet == true);
    CH_assert(m_isDefined == true);
    Box boundaryBox;
    getBoundaryFaces(boundaryBox, a_WGdnv.box(), a_dir, a_side);

    // In periodic case, this doesn't do anything
    if (!m_domain.isPeriodic(a_dir))
    {
        int lohisign;
        Box tmp = a_WGdnv.box();

        // Determine which side and thus shifting directions
        lohisign = sign(a_side);
        tmp.shiftHalf(a_dir,lohisign);

        // Is there a domain boundary next to this grid
        if (!m_domain.contains(tmp))
        {
            tmp &= m_domain;

            Box boundaryBox;

            // Find the strip of cells next to the domain boundary
            if (a_side == Side::Lo)
            {
                boundaryBox = bdryLo(tmp,a_dir);
            }
            else
            {
                boundaryBox = bdryHi(tmp,a_dir);
            }

            if(lohisign == -1 && a_dir == 1)
            {
                FORT_SWFAULTBCF(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_WShiftInside),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_INT(lohisign),
                    CHF_CONST_REAL(m_dx),
                    CHF_CONST_REAL(a_time),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_FRA((*m_bdryData)),
                    CHF_CONST_INT(m_numPatches),
                    CHF_CONST_VR(m_xcPatches),
                    CHF_CONST_VR(m_xwPatches),
                    CHF_CONST_VR(m_zcPatches),
                    CHF_CONST_VR(m_zwPatches),
                    CHF_CONST_VR(m_tauPatches),
                    CHF_BOX(boundaryBox));
            }
            else if(m_boundaryType[a_dir*2 + (lohisign+1)/2] == 0)
            {
                FORT_LINELASTOUTBCF(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_WShiftInside),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_INT(lohisign),
                    CHF_CONST_REAL(m_dx),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox));
            }
            else if(m_boundaryType[a_dir*2 + (lohisign+1)/2] == 1)
            {
                FORT_LINELASTFREEBCF(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_WShiftInside),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_INT(lohisign),
                    CHF_CONST_REAL(m_dx),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox));
            }
            else
            {
                MayDay::Error("Invalid Boundary Type");
            }
        }
    }
}

/// Set boundary slopes
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void SWIBC::setBdrySlopes(FArrayBox&       a_dW,
    const FArrayBox& a_W,
    const int&       a_dir,
    const Real&      a_time)
{
}

void SWIBC::artViscBC(FArrayBox&       a_F,
    const FArrayBox& a_U,
    const FArrayBox& a_divVel,
    const int&       a_dir,
    const Real&      a_time)
{
    pout() << "NOT SETUP :: SWIBC::artViscBC" << endl;
}

void SWIBC::updateBoundary(const FArrayBox& a_WHalf,int a_dir,const Real& a_dt)
{
    if(a_dir == 1 && bdryLo(m_domain,1).contains(bdryLo(a_WHalf.box(),1)))
    {
        FORT_SWSETBND(CHF_FRA((*m_bdryData)),
            CHF_BOX(bdryLo(a_WHalf.box(),1)),
            CHF_CONST_FRA(a_WHalf),
            CHF_CONST_REAL(a_dt));
    }
}

/// Do the initial tagging of cells
/**
*/
bool SWIBC::tagCellsInit(FArrayBox& markFAB)
{
    if(!m_isNucBoxSet)
    {
        IntVect nucSm;
        IntVect nucBg;
        int offSet = 0;
        for(int itor = 0; itor < (SpaceDim); itor++)
        {
            if(itor == 1)
            {
                nucSm.setVal(itor,0);
                nucBg.setVal(itor,0);
            }
            else
            {
                nucSm.setVal(itor,floor(m_nucPatch[offSet  ]/m_dx));
                nucBg.setVal(itor, ceil(m_nucPatch[offSet+1]/m_dx));
                offSet+= 2;
            }
        }
        m_nucBox.define(Box(nucSm,nucBg));
        m_isNucBoxSet = true;
    }
    markFAB.setVal(1,m_nucBox & markFAB.box(),0);

    // FORT_BOUNDREFINE(
    //     CHF_FRA1(markFAB,0),
    //     CHF_CONST_REAL(refLocation),
    //     CHF_CONST_REAL(m_dx),
    //     CHF_BOX(b));
    return true;
}
