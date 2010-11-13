#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "PseudoPulseIBC.H"
#include "PseudoPulseIBCF_F.H"
#include "LinElastPhysicsF_F.H"

/// Null Constructor
PseudoPulseIBC::PseudoPulseIBC()
{
}

/// Constructor which defines parameters used by Fortran routines
PseudoPulseIBC::PseudoPulseIBC(const Real& a_cs,
    const Real& a_cp,
    const Real& a_mu,
    const vector<Real> a_back,
    const Real& a_r0,
    const Real& a_slope,
    const Real& a_dTau,
    const Real& a_Vh,
    const Real& a_Vr,
    const Real& a_stopTime,
    const Real& a_jumpTime)
{
    FORT_LINELASTSETF(CHF_CONST_REAL(a_cs),CHF_CONST_REAL(a_cp),CHF_CONST_REAL(a_mu),CHF_CONST_VR(a_back));
    FORT_PSEUDOPULSESETF(CHF_CONST_REAL(a_r0),CHF_CONST_REAL(a_slope),CHF_CONST_REAL(a_dTau),CHF_CONST_REAL(a_Vh),CHF_CONST_REAL(a_Vr),CHF_CONST_REAL(a_jumpTime),CHF_CONST_REAL(a_stopTime));
    m_isFortranCommonSet = true;
}

/// Destructor
PseudoPulseIBC::~PseudoPulseIBC()
{
}

/// Define ?

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_physIBC() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void PseudoPulseIBC::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

/// Factory method - this object is its own factory
PhysIBC* PseudoPulseIBC::new_physIBC()
{
    PseudoPulseIBC* retval = new PseudoPulseIBC();
    if(m_isFortranCommonSet == true)
    {
        retval->setFortranCommonSet();
    }
    return static_cast<PhysIBC*>(retval);
}

/// Set up initial conditions
void PseudoPulseIBC::initialize(LevelData<FArrayBox>& a_U)
{
    // pout() << "PseudoPulseIBC::initialize" << endl;
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
        FORT_PSEUDOPULSEINITF(CHF_FRA(U),
            CHF_CONST_REAL(m_dx),
            CHF_BOX(uBox));
    }
}

bool PseudoPulseIBC::hasBdryData()
{
    return false;
}

/// Set boundary primitive values.
void PseudoPulseIBC::primBC(FArrayBox& a_WGdnv,
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
                FORT_PSEUDOPULSEFAULTBCF(CHF_FRA(a_WGdnv),
                    CHF_CONST_FRA(a_WShiftInside),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_INT(lohisign),
                    CHF_CONST_REAL(m_dx),
                    CHF_CONST_REAL(a_time),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox));
                // pout() << "Boundary Box :: " << m_bdryData->box() << " :: for box :: " << a_WGdnv.box() << endl;
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
void PseudoPulseIBC::setBdrySlopes(FArrayBox&       a_dW,
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
    //JK     MayDay::Error("PseudoPulseIBC::setBdrySlopes: All directions must be periodic");
    //JK }
}

void PseudoPulseIBC::artViscBC(FArrayBox&       a_F,
    const FArrayBox& a_U,
    const FArrayBox& a_divVel,
    const int&       a_dir,
    const Real&      a_time)
{
    pout() << "NOT SETUP :: PseudoPulseIBC::artViscBC" << endl;
}

void PseudoPulseIBC::updateBoundary(const FArrayBox& a_WHalf,int a_dir)
{
    if(a_dir == 1 && bdryLo(m_domain,1).contains(bdryLo(a_WHalf.box(),1)))
    {
        FORT_PSEUDOPULSESETBND(CHF_FRA((*m_bdryData)),
            CHF_BOX(bdryLo(a_WHalf.box(),1)),
            CHF_CONST_FRA(a_WHalf));
    }
}
