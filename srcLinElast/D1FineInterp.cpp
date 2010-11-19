#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "REAL.H"
#include "DataIterator.H"
#include "Tuple.H"
#include "InterpF_F.H"

#include "D1FineInterp.H"
#include "NamespaceHeader.H"

D1FineInterp::D1FineInterp()
    :
        is_defined(false)
{
}

D1FineInterp::~D1FineInterp()
{
}

void
D1FineInterp::define(const DisjointBoxLayout& a_fine_domain,
    const int& a_numcomps,
    const int& a_ref_ratio,
    const int& a_dir,
    const Box& a_fine_box)
{
    CH_TIME("D1FineInterp::define");

    m_dir = a_dir;
    m_ref_ratio = a_ref_ratio;
    m_numcomps = a_numcomps;
    m_fine_box = a_fine_box;
    is_defined = true;

    // Figure out the fine domain
    m_coarse_box.define(m_fine_box.smallEnd()/m_ref_ratio,
        (m_fine_box.bigEnd()+IntVect::Unit-BASISV(m_dir))/m_ref_ratio-(IntVect::Unit-BASISV(m_dir)),
        BASISV(m_dir));

    // Loop through the fine data boxes and define a new set of coarse data
    // boxes
    const DisjointBoxLayout& constGrids = a_fine_domain;
    DisjointBoxLayout coarsened_fine_domain;
    coarsened_fine_domain.deepCopy(constGrids);
    for(LayoutIterator lit = a_fine_domain.layoutIterator(); lit.ok(); ++lit)
    {
        const Box tmpFineBox = a_fine_domain[lit()];
        const Box tmpCoarseBox(tmpFineBox.smallEnd()/m_ref_ratio,
            (tmpFineBox.bigEnd()+IntVect::Unit-BASISV(m_dir))/m_ref_ratio-(IntVect::Unit-BASISV(m_dir)),
            BASISV(m_dir));
        coarsened_fine_domain.ref(lit()) = tmpCoarseBox;
    }
    coarsened_fine_domain.close();
    m_coarsened_fine_data.define ( coarsened_fine_domain,
        m_numcomps,
        IntVect::Unit-BASISV(m_dir));
}

bool
D1FineInterp::isDefined() const
{
    return ( is_defined );
}

// interpolate from coarse level to fine level
void
D1FineInterp::interpToFine(LevelData<FArrayBox>& a_fine_data,
    const LevelData<FArrayBox>& a_coarse_data)
{
    CH_TIME("D1FineInterp::interpToFine");
    CH_assert(is_defined);

    // this should handle all the periodic BCs as well,
    // by filling in the ghost cells in an appropriate way
    a_coarse_data.copyTo(a_coarse_data.interval(),
        m_coarsened_fine_data,
        m_coarsened_fine_data.interval() );

    const BoxLayout fine_domain = a_fine_data.boxLayout();
    DataIterator dit = fine_domain.dataIterator();

    // Loop over each fine box we need to fill
    for (dit.begin(); dit.ok(); ++dit)
    {
        const BaseFab<Real>& coarsened_fine = m_coarsened_fine_data[dit()];
        const Box& coarsened_fine_box = m_coarsened_fine_data.getBoxes()[dit()];
        BaseFab<Real>& fine = a_fine_data[dit()];
        // interpGridData interpolates from an entire coarse grid onto an
        // entire fine grid.
        interpGridData(fine,
            coarsened_fine,
            coarsened_fine_box,
            m_ref_ratio);
    }
    // CH_assert(false);
}

//JK void
//JK D1FineInterp::pwcinterpToFine(LevelData<FArrayBox>& a_fine_data,
//JK     const LevelData<FArrayBox>& a_coarse_data)
//JK {
//JK    CH_TIME("D1FineInterp::pwcinterpToFine");
//JK    CH_assert(is_defined);
//JK    // this should handle all the periodic BCs as well,
//JK    // by filling in the ghost cells in an appropriate way
//JK    a_coarse_data.copyTo(a_coarse_data.interval(),
//JK        m_coarsened_fine_data,
//JK        m_coarsened_fine_data.interval() );
//JK
//JK    const BoxLayout fine_domain = a_fine_data.boxLayout();
//JK    DataIterator dit = fine_domain.dataIterator();
//JK
//JK    for (dit.begin(); dit.ok(); ++dit)
//JK    {
//JK        const BaseFab<Real>& coarsened_fine = m_coarsened_fine_data[dit()];
//JK        const Box& coarsened_fine_box = m_coarsened_fine_data.getBoxes()[dit()];
//JK        BaseFab<Real>& fine = a_fine_data[dit()];
//JK        // interpGridData interpolates from an entire coarse grid onto an
//JK        // entire fine grid.
//JK        pwcinterpGridData(fine,
//JK            coarsened_fine,
//JK            coarsened_fine_box,
//JK            m_ref_ratio);
//JK    }
//JK }

//JK void
//JK D1FineInterp::pwcinterpGridData(BaseFab<Real>& a_fine,
//JK     const BaseFab<Real>& a_coarse,
//JK     const Box& a_coarsened_fine_box,
//JK     int a_ref_ratio) const
//JK {
//JK     CH_TIME("D1FineInterp::pwcinterpGridData");
//JK     // fill fine data with piecewise constant coarse data
//JK     const Box& b = a_coarsened_fine_box;
//JK     Box refbox(IntVect::Zero,
//JK         (a_ref_ratio-1)*IntVect::Unit);
//JK 
//JK     FORT_INTERPCONSTANT ( CHF_FRA(a_fine),
//JK         CHF_CONST_FRA(a_coarse),
//JK         CHF_BOX(b),
//JK         CHF_CONST_INT(a_ref_ratio),
//JK         CHF_BOX(refbox)
//JK         );
//JK }

// interpolate from fine grid to coarse grid.  prerequisite:
// coarsened.box contains coarsen(fine.box).
//
// uses piecewise bilinear interpolation with multidimensional-limited
// slopes.  see design document for details.
void
D1FineInterp::interpGridData(BaseFab<Real>& a_fine,
    const BaseFab<Real>& a_coarse,
    const Box& a_coarsened_fine_box,
    int a_ref_ratio)
const
{
    CH_TIME("D1FineInterp::interpGridData");

    // fill fine data with piecewise constant coarse data
    const Box& b = a_coarsened_fine_box;
    const int num_comp = a_fine.nComp ();
    Box refbox(IntVect::Zero,
        (a_ref_ratio-1)*(IntVect::Unit-BASISV(m_dir)));


    FORT_INTERPCONSTANT ( CHF_FRA(a_fine),
        CHF_CONST_FRA(a_coarse),
        CHF_BOX(b),
        CHF_CONST_INT(a_ref_ratio),
        CHF_BOX(refbox)
        );


    //  Tuple<BaseFab<Real>, SpaceDim> slopes;
    //  for (int dir = 0; dir < SpaceDim; ++dir)
    // hardwired to 3 due to lack of variable number of arguments in chfpp
    BaseFab<Real> slopes[3];
    for (int dir = 0; dir < 3; ++dir)
    {
        BaseFab<Real>& dir_slope = slopes[dir];
        dir_slope.resize(b, num_comp);
    }
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
        if(dir == m_dir)
        {
            slopes[m_dir].setVal(0.0);
            continue;
        }

        BaseFab<Real>& dir_slope = slopes[dir];

        const Box bcenter = grow(m_coarse_box,-BASISV(dir)) & b;
        if (!bcenter.isEmpty())
        {
            FORT_INTERPCENTRALSLOPE ( CHF_FRA ( dir_slope ),
                CHF_CONST_FRA ( a_coarse ),
                CHF_BOX ( bcenter ),
                CHF_CONST_INT ( dir )
                );
        }
        const Box blo = b & adjCellLo(grow(m_coarse_box,-BASISV(dir)),dir);
        if (!blo.isEmpty())
        {
            FORT_INTERPHISIDESLOPE ( CHF_FRA ( dir_slope ),
                CHF_CONST_FRA ( a_coarse ),
                CHF_BOX ( blo ),
                CHF_CONST_INT ( dir )

                );
        }
        const Box bhi = b & adjCellHi(grow(m_coarse_box,-BASISV(dir)),dir);
        if (!bhi.isEmpty())
        {
            FORT_INTERPLOSIDESLOPE ( CHF_FRA ( dir_slope ),
                CHF_CONST_FRA ( a_coarse ),
                CHF_BOX ( bhi ),
                CHF_CONST_INT ( dir )
                );
        }
    }
    
    // to do limits, we need to have a box which includes
    // the neighbors of a given point (to check for the
    // local maximum...

    //JK Modify so that the box excludes dim were not interpolating in
    Box neighborBox(-1*(IntVect::Unit-BASISV(m_dir)),
        (IntVect::Unit-BASISV(m_dir)),
        BASISV(m_dir));
    
    // GHM 7/12/01
    // interplimit iterates over box b_mod (was b), but cells within
    // 1 of the physical boundary never enter result (and this
    // wasted calculation may call upon uninitialized memory).
    // DFM 10/8/01
    // note that this turns off slope limiting for cells adjacent to the
    // boundary -- may want to revisit this in the future
    Box b_mod(b);
    b_mod.grow(IntVect::Unit-BASISV(m_dir));
    b_mod = m_coarse_box & b_mod;
    b_mod.grow(-(IntVect::Unit-BASISV(m_dir)));
    
    // create a box grown big enough to remove periodic BCs from domain
    Box domBox = grow(b, 2*(IntVect::Unit-BASISV(m_dir)));
    domBox = m_coarse_box & domBox;

    FORT_INTERPLIMIT ( CHF_FRA ( slopes[0] ),
        CHF_FRA ( slopes[1] ),
        CHF_FRA ( slopes[2] ),
        CHF_CONST_FRA ( a_coarse ),
        CHF_BOX ( b_mod ),
        CHF_BOX ( neighborBox ),
        CHF_BOX (domBox)
        );

    for (int dir = 0; dir < SpaceDim; ++dir)
    {
        if(dir == m_dir) continue;

        BaseFab<Real>& dir_slope = slopes[dir];
    
        FORT_INTERPLINEAR ( CHF_FRA ( a_fine ),
            CHF_CONST_FRA ( dir_slope ),
            CHF_BOX ( b ),
            CHF_CONST_INT ( dir ),
            CHF_CONST_INT ( a_ref_ratio ),
            CHF_BOX ( refbox )
            );
    }
}
#include "NamespaceFooter.H"
