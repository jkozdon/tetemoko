#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "RSIBC.H"
#include "RSIBCF_F.H"
#include "LinElastPhysicsF_F.H"

/// Null Constructor
RSIBC::RSIBC()
{
}

/// Constructor which defines parameters used by Fortran routines
RSIBC::RSIBC(const Real& a_cs,
    const Real& a_cp,
    const Real& a_mu,
    const vector<Real> a_back,
    const Real& a_r0,
    const Real& a_sigma,
    const Real& a_ntime,
    const Real& a_psi,
    const Real& a_a,
    const Real& a_b,
    const Real& a_V0,
    const Real& a_f0,
    const Real& a_L,
    const Real& a_fw,
    const Real& a_Vw)
{
    FORT_LINELASTSETF(CHF_CONST_REAL(a_cs),CHF_CONST_REAL(a_cp),CHF_CONST_REAL(a_mu),CHF_CONST_VR(a_back));
    FORT_RSSETF(CHF_CONST_REAL(a_r0),CHF_CONST_REAL(a_sigma),CHF_CONST_REAL(a_ntime),
        CHF_CONST_REAL(a_a),CHF_CONST_REAL(a_b),CHF_CONST_REAL(a_V0),CHF_CONST_REAL(a_f0),
        CHF_CONST_REAL(a_L),CHF_CONST_REAL(a_fw),CHF_CONST_REAL(a_Vw));
    m_isFortranCommonSet = true;
    m_psi = a_psi;
}

/// Destructor
RSIBC::~RSIBC()
{
}

/// Define ?

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_physIBC() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void RSIBC::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

/// Factory method - this object is its own factory
PhysIBC* RSIBC::new_physIBC()
{
    RSIBC* retval = new RSIBC();
    if(m_isFortranCommonSet == true)
    {
        retval->setFortranCommonSet();
    }
    retval->m_psi = m_psi;

    return static_cast<PhysIBC*>(retval);
}

/// Set up initial conditions
void RSIBC::initialize(LevelData<FArrayBox>& a_U)
{
    // pout() << "RSIBC::initialize" << endl;
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
        FORT_RSINITF(CHF_FRA(U),
            CHF_CONST_REAL(m_dx),
            CHF_BOX(uBox));
    }
}

bool RSIBC::hasBndryData()
{
    return true;
}

/// Set up initial conditions
void RSIBC::initializeBndry(LevelData<FArrayBox>& a_Psi)
{
    CH_assert(m_isFortranCommonSet == true);
    CH_assert(m_isDefined == true);

    // Iterator of all grids in this level
    for (DataIterator dit = a_Psi.dataIterator();
        dit.ok(); ++dit)
    {
        // Storage for current grid
        FArrayBox& Psi = a_Psi[dit()];

        // Box of current grid
        Box uBox = Psi.box();
        uBox &= m_domain;

        // Set up initial condition in this grid
        FORT_RSSETFAB(CHF_FRA1(Psi,0),CHF_BOX(uBox),CHF_REAL(m_psi));
    }
}

/// Set boundary primitive values.
void RSIBC::primBC(FArrayBox&            a_WGdnv,
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
                MayDay::Error("RSIBC:: Fault boundary called with no boundary data");
            }
            else
            {
                FORT_LINELASTOUTBCF(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_WShiftInside),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_INT(lohisign),
                    CHF_CONST_REAL(m_dx),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox));
            }
        }
    }
}

/// Set boundary primitive values.
void RSIBC::primBC(FArrayBox&            a_WGdnv,
    const FArrayBox&      a_WShiftInside,
    const FArrayBox&      a_W,
    const FArrayBox&      a_Psi,
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
                FORT_RSFAULTBCF(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_WShiftInside),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_FRA1(a_Psi,0),
                    CHF_CONST_INT(lohisign),
                    CHF_CONST_REAL(m_dx),
                    CHF_CONST_REAL(a_time),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox));
            }
            else
            {
                FORT_LINELASTOUTBCF(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_WShiftInside),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_INT(lohisign),
                    CHF_CONST_REAL(m_dx),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox));
            }
        }
    }
}

/// Set boundary slopes
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void RSIBC::setBdrySlopes(FArrayBox&       a_dW,
    const FArrayBox& a_W,
    const int&       a_dir,
    const Real&      a_time)
{
    //JK if (! m_domain.isPeriodic(a_dir))
    //JK {
    //JK     Box loBox,hiBox,centerBox,domain;
    //JK     int hasLo,hasHi;
    //JK     Box slopeBox = a_dW.box();
    //JK     slopeBox.grow(a_dir,1);

    //JK     // Generate the domain boundary boxes, loBox and hiBox, if there are
    //JK     // domain boundarys there
    //JK     loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
    //JK         slopeBox,m_domain,a_dir);

    //JK     // Set the boundary slopes if necessary
    //JK     if ((hasLo != 0) || (hasHi != 0))
    //JK     {
    //JK         //JK FORT_RAMPSLOPEBCSF(CHF_FRA(a_dW),
    //JK         //JK     CHF_CONST_FRA(a_W),
    //JK         //JK     CHF_CONST_REAL(m_dx),
    //JK         //JK     CHF_CONST_INT(a_dir),
    //JK         //JK     CHF_BOX(loBox),
    //JK         //JK     CHF_CONST_INT(hasLo),
    //JK         //JK     CHF_BOX(hiBox),
    //JK         //JK     CHF_CONST_INT(hasHi));
    //JK     }
    //JK     MayDay::Error("RSIBC::setBdrySlopes: All directions must be periodic");
    //JK }
}

void RSIBC::artViscBC(FArrayBox&       a_F,
    const FArrayBox& a_U,
    const FArrayBox& a_divVel,
    const int&       a_dir,
    const Real&      a_time)
{
    pout() << "NOT SETUP :: RSIBC::artViscBC" << endl;
}
