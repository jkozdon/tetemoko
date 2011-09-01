#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "RSIBC.H"
#include "RSIBCF_F.H"
#include "LinElastPhysicsF_F.H"
#include "BoxIterator.H"

#define IX_VX          (0)
#define IX_VZ          (1)
#define IX_SXY         (2)
#define IX_SYZ         (3)
#define IX_V           (4)
#define IX_SLIP_X      (5)
#define IX_SLIP_Z      (6)
#define IX_SLIP_T      (7)
#define IX_RT          (8)
#define IX_SYY         (9)
#define IX_PSI         (10)

/// Null Constructor
RSIBC::RSIBC()
{
}

/// Constructor which defines parameters used by Fortran routines
RSIBC::RSIBC(const Real& a_r0,
    const Real& a_x0,
    const Real& a_y0,
    const Real& a_nsig,
    const Real& a_ntime,
    const Real& a_psi,
    const Real& a_a,
    const Real& a_b,
    const Real& a_V0,
    const Real& a_f0,
    const Real& a_L,
    const Real& a_fw,
    const Real& a_Vw,
    const Real& a_fExp,
    const Real a_ruptureVelocityThreshold,
    const Vector<int>& a_boundaryType)
{
    FORT_RSSETF(CHF_CONST_REAL(a_r0),CHF_CONST_REAL(a_x0),CHF_CONST_REAL(a_y0),CHF_CONST_REAL(a_nsig),
        CHF_CONST_REAL(a_ntime),CHF_CONST_REAL(a_a),CHF_CONST_REAL(a_b),CHF_CONST_REAL(a_V0),CHF_CONST_REAL(a_f0),
        CHF_CONST_REAL(a_L),CHF_CONST_REAL(a_fw),CHF_CONST_REAL(a_Vw),CHF_CONST_REAL(a_fExp),CHF_CONST_REAL(a_ruptureVelocityThreshold));
    m_boundaryType       = a_boundaryType;
    m_isFortranCommonSet = true;
    m_psi = a_psi;
    m_numBdryVars = 11;
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
    retval->m_boundaryType = m_boundaryType;
    retval->m_psi          = m_psi;
    retval->m_numBdryVars  = m_numBdryVars;
    return static_cast<PhysIBC*>(retval);
}

/// Set up initial conditions
void RSIBC::initialize(LevelData<FArrayBox>& a_U)
{
    // pout() << "RSIBC::initialize" << endl;
    CH_assert(m_isFortranCommonSet == true);
    CH_assert(m_isFortranCommonLESet == true);
    CH_assert(m_isFortranCommonPlasticSet == true);
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
        FORT_LEINITF(CHF_FRA(U),
            CHF_CONST_REAL(m_dx),
            CHF_BOX(uBox));
    }
}

/// Set up initial conditions
void RSIBC::initializeBdry(LevelData<FArrayBox>& a_B)
{
    const Real tmpVal  =  0.0;
    const Real tmpVal2 = 100;
    for (DataIterator dit = a_B.dataIterator(); dit.ok(); ++dit)
    {
        // Storage for current grid
        FArrayBox& B = a_B[dit()];

        // Box of current grid
        Box bBox = B.box();
        bBox &= m_domain;

        // Set up initial condition in this grid
        FORT_LINELASTSETFAB(CHF_FRA1(B,IX_VX),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
        FORT_LINELASTSETFAB(CHF_FRA1(B,IX_VZ),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
        FORT_LINELASTSETFAB(CHF_FRA1(B,IX_SXY),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
        FORT_LINELASTSETFAB(CHF_FRA1(B,IX_SYZ),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
        FORT_LINELASTSETFAB(CHF_FRA1(B,IX_V),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
        FORT_LINELASTSETFAB(CHF_FRA1(B,IX_SLIP_X),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
        FORT_LINELASTSETFAB(CHF_FRA1(B,IX_SLIP_Z),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
        FORT_LINELASTSETFAB(CHF_FRA1(B,IX_SLIP_T),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
        FORT_LINELASTSETFAB(CHF_FRA1(B,IX_RT),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal2));
        FORT_LINELASTSETFAB(CHF_FRA1(B,IX_SYY),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
        FORT_LINELASTSETFAB(CHF_FRA1(B,IX_PSI),
            CHF_BOX(bBox),
            CHF_CONST_REAL(m_psi));
    }
}

// SET pointer and make a temporary storage array
void RSIBC::setBdryData(FArrayBox* a_bdryData)
{
    m_bdryData    = a_bdryData;
    m_tmpBdryData = new FArrayBox(m_bdryData->box(), m_numBdryVars);
    m_tmpBdryDataSet = false;
}

bool RSIBC::hasBdryData()
{
    return true;
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
    CH_assert(m_isFortranCommonLESet == true);
    CH_assert(m_isFortranCommonPlasticSet == true);
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
                if(m_tmpBdryDataSet)
                {
                    FORT_RSFAULTBCF(CHF_FRA(a_WGdnv),
                        CHF_CONST_FRA(a_WShiftInside),
                        CHF_CONST_FRA(a_W),
                        CHF_CONST_FRA1((*m_tmpBdryData),IX_PSI),
                        CHF_CONST_INT(lohisign),
                        CHF_CONST_REAL(m_dx),
                        CHF_CONST_REAL(a_time),
                        CHF_CONST_INT(a_dir),
                        CHF_BOX(boundaryBox));
                }
                else
                {
                    FORT_RSFAULTBCF(CHF_FRA(a_WGdnv),
                        CHF_CONST_FRA(a_WShiftInside),
                        CHF_CONST_FRA(a_W),
                        CHF_CONST_FRA1((*m_bdryData),IX_PSI),
                        CHF_CONST_INT(lohisign),
                        CHF_CONST_REAL(m_dx),
                        CHF_CONST_REAL(a_time),
                        CHF_CONST_INT(a_dir),
                        CHF_BOX(boundaryBox));
                }
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
void RSIBC::setBdrySlopes(FArrayBox&       a_dW,
    const FArrayBox& a_W,
    const int&       a_dir,
    const Real&      a_time)
{
}

void RSIBC::artViscBC(FArrayBox&       a_F,
    const FArrayBox& a_U,
    const FArrayBox& a_divVel,
    const int&       a_dir,
    const Real&      a_time)
{
    pout() << "NOT SETUP :: RSIBC::artViscBC" << endl;
}

void RSIBC::updateBoundary(const FArrayBox& a_WHalf,int a_dir,const Real& a_dt,const Real& a_time,const bool a_final)
{
    if(a_dir == 1 && bdryLo(m_domain,a_dir).contains(bdryLo(a_WHalf.box(),a_dir)))
    {
        if(a_final)
        {
            FORT_RSSETBND(
                CHF_FRA((*m_bdryData)),
                CHF_BOX(bdryLo(a_WHalf.box(),a_dir)),
                CHF_CONST_FRA((*m_bdryData)),
                CHF_CONST_FRA(a_WHalf),
                CHF_CONST_REAL(a_dt),
                CHF_CONST_REAL(a_time));
        }
        // THIS MAY BE NECESSARY IN 3-D, NOT SURE!!!
        // else if(m_tmpBdryDataSet)
        // {
        //     FORT_RSSETBND(
        //         CHF_FRA((*m_tmpBdryData)),
        //         CHF_BOX(bdryLo(a_WHalf.box(),a_dir)),
        //         CHF_CONST_FRA((*m_tmpBdryData)),
        //         CHF_CONST_FRA(a_WHalf),
        //         CHF_CONST_REAL(a_dt),
        //         CHF_CONST_REAL(a_time));
        // }
        else
        {
            // m_tmpBdryData = new FArrayBox(m_bdryData->box(), m_numBdryVars);
            FORT_RSSETBND(
                CHF_FRA((*m_tmpBdryData)),
                CHF_BOX(bdryLo(a_WHalf.box(),a_dir)),
                CHF_CONST_FRA((*m_bdryData)),
                CHF_CONST_FRA(a_WHalf),
                CHF_CONST_REAL(a_dt),
                CHF_CONST_REAL(a_time));
            m_tmpBdryDataSet = true;
        }
    }
}

/// Do the initial tagging of cells
/**
*/
bool RSIBC::tagCellsInit(FArrayBox& markFAB)
{
    // markFAB.setVal(1,m_patchBoxes[itor] & markFAB.box(),0);
    return true;
}
