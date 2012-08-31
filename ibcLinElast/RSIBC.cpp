#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "RSIBC.H"
#include "RSINDEX.H"
#include "RSIBCF_F.H"
#include "LinElastPhysicsF_F.H"
#include "BoxIterator.H"

/// Null Constructor
RSIBC::RSIBC()
{
  m_isFortranCommonSet = false;
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
    const Real& a_ramp_x,
    const Real& a_ramp_w,
    const Real& a_ramp_a,
    const Real& a_ramp_Vw,
    const Real a_ruptureVelocityThreshold,
    const Vector<int>& a_boundaryType)
{
  FORT_RSSETF(CHF_CONST_REAL(a_r0),CHF_CONST_REAL(a_x0),CHF_CONST_REAL(a_y0),CHF_CONST_REAL(a_nsig),
      CHF_CONST_REAL(a_ntime),CHF_CONST_REAL(a_a),CHF_CONST_REAL(a_b),CHF_CONST_REAL(a_V0),CHF_CONST_REAL(a_f0),
      CHF_CONST_REAL(a_L),CHF_CONST_REAL(a_fw),CHF_CONST_REAL(a_Vw),CHF_CONST_REAL(a_fExp),
      CHF_CONST_REAL(a_ramp_x),CHF_CONST_REAL(a_ramp_w),CHF_CONST_REAL(a_ramp_a),CHF_CONST_REAL(a_ramp_Vw),
      CHF_CONST_REAL(a_ruptureVelocityThreshold));
  m_boundaryType       = a_boundaryType;
  m_isFortranCommonSet = true;
  m_psi = a_psi;
  m_R0 = 1e10;
  m_x0 = a_x0;
  m_numBdryVars.resize(2*CH_SPACEDIM,NUM_BDRY_HAT_VARS);
  for(int ix = 0; ix < 2*CH_SPACEDIM; ix++) m_tmpBdryDataSet[ix] = false;

  m_numBdryVars[2] = 11;

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


void RSIBC::setR0(const Real& a_R0)
{
  m_R0 = a_R0;
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
  retval->m_R0           = m_R0;
  retval->m_x0           = m_x0;
  for(int ix = 0; ix < 2*CH_SPACEDIM; ix++)
    retval->m_numBdryVars[ix]  = m_numBdryVars[ix];
  return static_cast<PhysIBC*>(retval);
}

/// Set up initial conditions
void RSIBC::initialize(LevelData<FArrayBox>& a_U)
{
  // pout() << "RSIBC::initialize" << endl;
  CH_assert(m_isFortranCommonSet == true);
  // CH_assert(m_isFortranCommonLESet == true);
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
void RSIBC::initializeBdry(LevelData<FArrayBox>& a_B,int a_side)
{
  const Real tmpVal  =  0.0;
  const Real tmpVal2 = 1e10;
  if(a_side!=2)
  {
    for (DataIterator dit = a_B.dataIterator(); dit.ok(); ++dit)
    {
      // Storage for current grid
      FArrayBox& B = a_B[dit()];

      // Box of current grid
      Box bBox = B.box();
      bBox &= m_domain;

      // Set up initial condition in this grid
      for(int ix = 0; ix < NUM_BDRY_HAT_VARS; ix++)
        FORT_LINELASTSETFAB(CHF_FRA1(B,ix),
            CHF_BOX(bBox),
            CHF_CONST_REAL(tmpVal));
    }
  }
  else
  {
    for (DataIterator dit = a_B.dataIterator(); dit.ok(); ++dit)
    {
      // Storage for current grid
      FArrayBox& B = a_B[dit()];

      // Box of current grid
      Box bBox = B.box();
      bBox &= m_domain;

      // Set up initial condition in this grid
      FORT_LINELASTSETFAB(CHF_FRA1(B,RX_VX),
          CHF_BOX(bBox),
          CHF_CONST_REAL(tmpVal));
      FORT_LINELASTSETFAB(CHF_FRA1(B,RX_VZ),
          CHF_BOX(bBox),
          CHF_CONST_REAL(tmpVal));
      FORT_LINELASTSETFAB(CHF_FRA1(B,RX_SXY),
          CHF_BOX(bBox),
          CHF_CONST_REAL(tmpVal));
      FORT_LINELASTSETFAB(CHF_FRA1(B,RX_SYZ),
          CHF_BOX(bBox),
          CHF_CONST_REAL(tmpVal));
      FORT_LINELASTSETFAB(CHF_FRA1(B,RX_V),
          CHF_BOX(bBox),
          CHF_CONST_REAL(tmpVal));
      FORT_LINELASTSETFAB(CHF_FRA1(B,RX_DX),
          CHF_BOX(bBox),
          CHF_CONST_REAL(tmpVal));
      FORT_LINELASTSETFAB(CHF_FRA1(B,RX_DZ),
          CHF_BOX(bBox),
          CHF_CONST_REAL(tmpVal));
      FORT_LINELASTSETFAB(CHF_FRA1(B,RX_DT),
          CHF_BOX(bBox),
          CHF_CONST_REAL(tmpVal));
      FORT_LINELASTSETFAB(CHF_FRA1(B,RX_RT),
          CHF_BOX(bBox),
          CHF_CONST_REAL(tmpVal2));
      FORT_LINELASTSETFAB(CHF_FRA1(B,RX_SYY),
          CHF_BOX(bBox),
          CHF_CONST_REAL(tmpVal));
      FORT_LINELASTSETFAB(CHF_FRA1(B,RX_PSI),
          CHF_BOX(bBox),
          CHF_CONST_REAL(m_psi));
    }
  }
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
  // CH_assert(m_isFortranCommonLESet == true);
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
        if(m_tmpBdryDataSet[2])
        {
          FORT_RSFAULTBCF(CHF_FRA(a_WGdnv),
              CHF_CONST_FRA(a_WShiftInside),
              CHF_CONST_FRA(a_W),
              CHF_CONST_FRA1((*m_tmpBdryData[2]),RX_PSI),
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
              CHF_CONST_FRA1((*m_bdryData[2]),RX_PSI),
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

void RSIBC::updateBoundary(const FArrayBox& a_WHalf,int a_dir,const Real& a_dt,const Real& a_dx,const Real& a_time,const bool a_final)
{
  if(bdryLo(m_domain,a_dir).contains(bdryLo(a_WHalf.box(),a_dir)))
  {
    if(a_dir == 1)
    {
      if(a_final)
      {
        FORT_RSSETBND(
            CHF_FRA((*m_bdryData[2*a_dir])),
            CHF_BOX(bdryLo(a_WHalf.box(),a_dir)),
            CHF_CONST_FRA((*m_bdryData[2*a_dir])),
            CHF_CONST_FRA(a_WHalf),
            CHF_CONST_REAL(a_dt),
            CHF_CONST_REAL(a_dx),
            CHF_CONST_REAL(a_time));
      }
      // THIS MAY BE NECESSARY IN 3-D, NOT SURE!!!
      // else if(m_tmpBdryDataSet[2*a_dir])
      // {
      //     FORT_RSSETBND(
      //         CHF_FRA((*m_tmpBdryData[2*a_dir])),
      //         CHF_BOX(bdryLo(a_WHalf.box(),a_dir)),
      //         CHF_CONST_FRA((*m_tmpBdryData[2*a_dir])),
      //         CHF_CONST_FRA(a_WHalf),
      //         CHF_CONST_REAL(a_dt),
      //         CHF_CONST_REAL(a_dx),
      //         CHF_CONST_REAL(a_time));
      // }
      else
      {
        // m_tmpBdryData[2*a_dir] = new FArrayBox(m_bdryData[2*a_dir]->box(), m_numBdryVars[2*a_dir]);
        FORT_RSSETBND(
            CHF_FRA((*m_tmpBdryData[2*a_dir])),
            CHF_BOX(bdryLo(a_WHalf.box(),a_dir)),
            CHF_CONST_FRA((*m_bdryData[2*a_dir])),
            CHF_CONST_FRA(a_WHalf),
            CHF_CONST_REAL(a_dt),
            CHF_CONST_REAL(a_dx),
            CHF_CONST_REAL(a_time));
        m_tmpBdryDataSet[2*a_dir] = true;
      }
    }
    else
    {
      int side = 2*a_dir;
      if(a_final)
      {
        FORT_LESETBND(
            CHF_FRA((*m_bdryData[side])),
            CHF_BOX(bdryLo(a_WHalf.box(),a_dir)),
            CHF_CONST_FRA((*m_bdryData[side])),
            CHF_CONST_FRA(a_WHalf),
            CHF_CONST_REAL(a_dt),
            CHF_CONST_REAL(a_dx),
            CHF_CONST_INT(side),
            CHF_CONST_REAL(a_time));
      }
      // THIS MAY BE NECESSARY IN 3-D, NOT SURE!!!
      // else if(m_tmpBdryDataSet[side])
      // {
      //     FORT_LESETBND(
      //         CHF_FRA((*m_tmpBdryData[side])),
      //         CHF_BOX(bdryLo(a_WHalf.box(),a_dir)),
      //         CHF_CONST_FRA((*m_tmpBdryData[side])),
      //         CHF_CONST_FRA(a_WHalf),
      //         CHF_CONST_REAL(a_dt),
      //         CHF_CONST_REAL(a_dx),
      //         CHF_CONST_REAL(a_time));
      // }
      else
      {
        // m_tmpBdryData[side] = new FArrayBox(m_bdryData[side]->box(), m_numBdryVars[side]);
        FORT_LESETBND(
            CHF_FRA((*m_tmpBdryData[side])),
            CHF_BOX(bdryLo(a_WHalf.box(),a_dir)),
            CHF_CONST_FRA((*m_bdryData[side])),
            CHF_CONST_FRA(a_WHalf),
            CHF_CONST_REAL(a_dt),
            CHF_CONST_REAL(a_dx),
            CHF_CONST_INT(side),
            CHF_CONST_REAL(a_time));
        m_tmpBdryDataSet[side] = true;
      }
    }
  }

  if(bdryHi(m_domain,a_dir).contains(bdryHi(a_WHalf.box(),a_dir)))
  {
    int side = 2*a_dir+1;
    if(a_final)
    {
      FORT_LESETBND(
          CHF_FRA((*m_bdryData[side])),
          CHF_BOX(bdryHi(a_WHalf.box(),a_dir)),
          CHF_CONST_FRA((*m_bdryData[side])),
          CHF_CONST_FRA(a_WHalf),
          CHF_CONST_REAL(a_dt),
          CHF_CONST_REAL(a_dx),
          CHF_CONST_INT(side),
          CHF_CONST_REAL(a_time));
    }
    // THIS MAY BE NECESSARY IN 3-D, NOT SURE!!!
    // else if(m_tmpBdryDataSet[side])
    // {
    //     FORT_LESETBND(
    //         CHF_FRA((*m_tmpBdryData[side])),
    //         CHF_BOX(bdryHi(a_WHalf.box(),a_dir)),
    //         CHF_CONST_FRA((*m_tmpBdryData[side])),
    //         CHF_CONST_FRA(a_WHalf),
    //         CHF_CONST_REAL(a_dt),
    //         CHF_CONST_REAL(a_dx),
    //         CHF_CONST_REAL(a_time));
    // }
    else
    {
      // m_tmpBdryData[side] = new FArrayBox(m_bdryData[side]->box(), m_numBdryVars[side]);
      FORT_LESETBND(
          CHF_FRA((*m_tmpBdryData[side])),
          CHF_BOX(bdryHi(a_WHalf.box(),a_dir)),
          CHF_CONST_FRA((*m_bdryData[side])),
          CHF_CONST_FRA(a_WHalf),
          CHF_CONST_REAL(a_dt),
          CHF_CONST_REAL(a_dx),
          CHF_CONST_INT(side),
          CHF_CONST_REAL(a_time));
      m_tmpBdryDataSet[side] = true;
    }
  }
}

/// Do the initial tagging of cells
/**
*/
bool RSIBC::tagCellsInit(FArrayBox& markFAB,const Real& threshold)
{
  return true;
  // If grid spacing > R0 refine otherwise only refine the nucleation patch
  if(m_dx > m_R0)
  {
    // pout() << m_dx << "  " << m_R0 << endl;
    markFAB.setVal(1,bdryLo(m_domain.domainBox(),1,1) & markFAB.box(),0);
  }
  else
  {
    IntVect nucSm;
    IntVect nucBg;
    if(SpaceDim > 0)
    {
      nucSm.setVal(0,floor(m_x0/m_dx));
      nucBg.setVal(0, ceil(m_x0/m_dx));
    }
    if(SpaceDim > 1)
    {
      nucSm.setVal(1,0);
      nucBg.setVal(1,0);
    }
    if(SpaceDim > 2)
    {
      nucSm.setVal(2,floor(m_x0/m_dx));
      nucBg.setVal(2, ceil(m_x0/m_dx));
    }
    markFAB.setVal(1,Box(nucSm,nucBg)& markFAB.box(),0);
  }

  // pout() << m_domain << endl;
  // FORT_ALLBOUNDREFINE(
  //     CHF_FRA1(markFAB,0),
  //     CHF_CONST_REAL(threshold),
  //     CHF_BOX(markFAB.box()));
  return true;
}

// Names for the boundary data
Vector<string> RSIBC::bdryNames(int a_side)
{
  Vector<string> bdryNames;
  if(a_side == 2)
  {
    // MUST MATCH RSINDEX.H ORDER!!!
    bdryNames.push_back("Vx");
    bdryNames.push_back("Vz");
    bdryNames.push_back("V");
    bdryNames.push_back("sxy");
    bdryNames.push_back("syz");
    bdryNames.push_back("syy");
    bdryNames.push_back("slip x");
    bdryNames.push_back("slip z");
    bdryNames.push_back("slip total");
    bdryNames.push_back("rupture time");
    bdryNames.push_back("psi");
  }
  else
  {
    // WOULD BE BETTER TO CALL ORIGINAL FUNCTION!

    // MUST MATCH LEPhysIBC.H ORDER!!!
    bdryNames.push_back("vx_1");
    bdryNames.push_back("vy_1");
    bdryNames.push_back("vz_1");
    bdryNames.push_back("Tx_1");
    bdryNames.push_back("Ty_1");
    bdryNames.push_back("Tz_1");
    bdryNames.push_back("un_1");
    bdryNames.push_back("um_1");
    bdryNames.push_back("uz_1");

    bdryNames.push_back("vx_2");
    bdryNames.push_back("vy_2");
    bdryNames.push_back("vz_2");
    bdryNames.push_back("Tx_2");
    bdryNames.push_back("Ty_2");
    bdryNames.push_back("Tz_2");
    bdryNames.push_back("un_2");
    bdryNames.push_back("um_2");
    bdryNames.push_back("uz_2");
  }
  return bdryNames;
}
