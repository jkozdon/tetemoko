#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "SWIBC.H"
#include "SWIBCF_F.H"
#include "LinElastPhysicsF_F.H"
#include "BoxIterator.H"

/// Null Constructor
SWIBC::SWIBC()
{
}

/// Constructor which defines parameters used by Fortran routines
SWIBC::SWIBC(const Real& a_fricS,
    const Real& a_fricD,
    const Real& a_weakD,
    const Real& a_smoothValue,
    const int a_numPatches,
    const vector<Real> a_xcPatches,
    const vector<Real> a_xwPatches,
    const vector<Real> a_zcPatches,
    const vector<Real> a_zwPatches,
    const vector<Real> a_tauPatches,
    const vector<Real> a_fricBoxCenter,
    const vector<Real> a_fricBoxWidth,
    const Real a_outsideFriction,
    const Real a_ruptureVelocityThreshold,
    const Vector<int>& a_boundaryType)
{
    FORT_SWSETF(CHF_CONST_REAL(a_fricS),CHF_CONST_REAL(a_fricD),CHF_CONST_REAL(a_weakD),
        CHF_CONST_REAL(a_smoothValue), CHF_CONST_REAL(a_ruptureVelocityThreshold));
    m_boundaryType       = a_boundaryType;
    m_isFortranCommonSet = true;
    m_isPatchBoxSet = false;

    m_numPatches = a_numPatches;
    m_xcPatches = a_xcPatches;
    m_xwPatches = a_xwPatches;
    m_zcPatches = a_zcPatches;
    m_zwPatches = a_zwPatches;
    m_tauPatches = a_tauPatches;

    m_smoothValue     = a_smoothValue;
    m_fricBoxCenter   = a_fricBoxCenter;
    m_fricBoxWidth    = a_fricBoxWidth; 
    m_outsideFriction = a_outsideFriction;

    m_numBdryVars = 9;
}

/// Destructor
SWIBC::~SWIBC()
{
}


/// Define ?
//
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
    retval->m_boundaryType        = m_boundaryType;
    retval->m_isPatchBoxSet       = m_isPatchBoxSet;
    retval->m_numPatches          = m_numPatches;
    retval->m_xcPatches           = m_xcPatches;
    retval->m_xwPatches           = m_xwPatches;
    retval->m_zcPatches           = m_zcPatches;
    retval->m_zwPatches           = m_zwPatches;
    retval->m_tauPatches          = m_tauPatches;
    retval->m_fricBoxCenter       = m_fricBoxCenter;
    retval->m_fricBoxWidth        = m_fricBoxWidth;
    retval->m_outsideFriction     = m_outsideFriction;
    retval->m_smoothValue         = m_smoothValue;
    retval->m_patchBoxes          = m_patchBoxes;
    retval->m_smoothValue         = m_smoothValue;
    retval->m_smoothWidthNumCells = m_smoothWidthNumCells;
    retval->m_numBdryVars         = m_numBdryVars;
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
        FORT_LEINITF(CHF_FRA(U),
            CHF_CONST_REAL(m_dx),
            CHF_BOX(uBox));
    }
}

/// Set up initial conditions
void SWIBC::initializeBdry(LevelData<FArrayBox>& a_B)
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
        FORT_LINELASTSETFAB(CHF_FRA1(B,4),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
        FORT_LINELASTSETFAB(CHF_FRA1(B,5),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
        FORT_LINELASTSETFAB(CHF_FRA1(B,6),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal2));
    }
}

// SET pointer and make a temporary storage array
void SWIBC::setBdryData(FArrayBox* a_bdryData)
{
    m_bdryData    = a_bdryData;
    m_tmpBdryData = new FArrayBox(m_bdryData->box(), m_numBdryVars);
    m_tmpBdryDataSet = false;
}

bool SWIBC::hasBdryData()
{
    return true;
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
                if(m_tmpBdryDataSet)
                {
                    FORT_SWFAULTBCF(CHF_FRA(a_WGdnv),
                        CHF_CONST_FRA(a_WShiftInside),
                        CHF_CONST_FRA(a_W),
                        CHF_CONST_INT(lohisign),
                        CHF_CONST_REAL(m_dx),
                        CHF_CONST_REAL(a_time),
                        CHF_CONST_INT(a_dir),
                        CHF_CONST_FRA((*m_tmpBdryData)),
                        CHF_CONST_INT(m_numPatches),
                        CHF_CONST_VR(m_xcPatches),
                        CHF_CONST_VR(m_xwPatches),
                        CHF_CONST_VR(m_zcPatches),
                        CHF_CONST_VR(m_zwPatches),
                        CHF_CONST_VR(m_tauPatches),
                        CHF_CONST_VR(m_fricBoxCenter),
                        CHF_CONST_VR(m_fricBoxWidth),
                        CHF_CONST_REAL(m_outsideFriction),
                        CHF_BOX(boundaryBox));
                }
                else
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
                        CHF_CONST_VR(m_fricBoxCenter),
                        CHF_CONST_VR(m_fricBoxWidth),
                        CHF_CONST_REAL(m_outsideFriction),
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

void SWIBC::updateBoundary(const FArrayBox& a_WHalf,int a_dir,const Real& a_dt,const Real& a_time,const bool a_final)
{
    if(a_dir == 1 && bdryLo(m_domain,a_dir).contains(bdryLo(a_WHalf.box(),a_dir)))
    {
        if(a_final)
        {
            FORT_SWSETBND(
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
        //     FORT_SWSETBND(
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
            FORT_SWSETBND(
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
bool SWIBC::tagCellsInit(FArrayBox& markFAB)
{
    // We do this here becauase we need m_dx to be set, and it isn't set when
    // the object is defined
    if(!m_isPatchBoxSet)
    {
        m_patchBoxes.resize(m_numPatches);
        // Loop over all the patches and figure out the boxes
        for(int itor = 0; itor < m_numPatches; itor ++)
        {
            IntVect nucSm;
            IntVect nucBg;
            int offSet = 0;
            if(SpaceDim > 0)
            {
                nucSm.setVal(0,floor((m_fricBoxCenter[0]+m_xcPatches[itor]-m_xwPatches[itor])/m_dx));
                nucBg.setVal(0, ceil((m_fricBoxCenter[0]+m_xcPatches[itor]+m_xwPatches[itor])/m_dx));
            }
            if(SpaceDim > 1)
            {
                nucSm.setVal(1,0);
                nucBg.setVal(1,0);
            }
            if(SpaceDim > 2)
            {
                nucSm.setVal(2,floor((m_fricBoxCenter[1]+m_zcPatches[itor]-m_zwPatches[itor])/m_dx));
                nucBg.setVal(2, ceil((m_fricBoxCenter[1]+m_zcPatches[itor]+m_zwPatches[itor])/m_dx));
            }
            m_patchBoxes[itor] = Box(nucSm,nucBg);
            m_smoothWidthNumCells = ceil(m_smoothValue / m_dx / 2);
            m_isPatchBoxSet = true;
        }
    }

    for(int itor = 0; itor < m_numPatches; itor ++)
    {
        markFAB.setVal(1,m_patchBoxes[itor] & markFAB.box(),0);
    }
    // markFAB.setVal(1,markFAB.box(),0);
    // for(int itor = 0; itor < m_patchBoxes.capacity(); itor++)
    // {
    //     markFAB.setVal(1,adjCellLo(m_patchBoxes[itor],0, m_smoothWidthNumCells) & markFAB.box(),0);
    //     markFAB.setVal(1,adjCellHi(m_patchBoxes[itor],0, m_smoothWidthNumCells) & markFAB.box(),0);
    //     markFAB.setVal(1,adjCellLo(m_patchBoxes[itor],0,-m_smoothWidthNumCells) & markFAB.box(),0);
    //     markFAB.setVal(1,adjCellHi(m_patchBoxes[itor],0,-m_smoothWidthNumCells) & markFAB.box(),0);
    //     if(SpaceDim > 2)
    //     {
    //         markFAB.setVal(1,adjCellLo(m_patchBoxes[itor],2, m_smoothWidthNumCells) & markFAB.box(),0);
    //         markFAB.setVal(1,adjCellHi(m_patchBoxes[itor],2, m_smoothWidthNumCells) & markFAB.box(),0);
    //         markFAB.setVal(1,adjCellLo(m_patchBoxes[itor],2,-m_smoothWidthNumCells) & markFAB.box(),0);
    //         markFAB.setVal(1,adjCellHi(m_patchBoxes[itor],2,-m_smoothWidthNumCells) & markFAB.box(),0);
    //     }
    // }

    // FORT_BOUNDREFINE(
    //     CHF_FRA1(markFAB,0),
    //     CHF_CONST_REAL(refLocation),
    //     CHF_CONST_REAL(m_dx),
    //     CHF_BOX(b));
    return true;
}

void SWIBC::dumpBdryData(FILE * a_boundaryDataFile)
{
    // Get the box that defines this layer
    Box b = m_bdryData->box();

    // Remove the ghost cells
    b.grow(-(IntVect::Unit - BASISV(1)));


    // loop over the box saving the values
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit)
    {
        const IntVect& iv = bit();
        Real x = (iv[0]+0.5)*m_dx-m_fricBoxCenter[0];
        Real z = (iv[2]+0.5)*m_dx-m_fricBoxCenter[1];
        if(abs(x) <= m_fricBoxWidth[0] & abs(z) <= m_fricBoxWidth[1])
        {
            fprintf(a_boundaryDataFile,"%E %E %E\n", x, z, m_bdryData->get(iv,6));
        }
    }
}
