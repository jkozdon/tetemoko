#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "InterpF_F.H"
#include "CH_Timer.H"
#include "MayDay.H"
using std::cout;
using std::endl;

#include "D1PiecewiseLinearFillPatch.H"
#include "NamespaceHeader.H"

#ifndef copysign
template <class T>
inline
T
copysign (const T& a,
    const T& b)
{
    return ( b >= 0 ) ? ( ( a >= 0 ) ? a : -a) : ( (a >= 0 ) ? -a : a);
}

#endif

const int D1PiecewiseLinearFillPatch::s_stencil_radius = 1;

D1PiecewiseLinearFillPatch::D1PiecewiseLinearFillPatch()
    :
        m_is_defined(false)
{
}

D1PiecewiseLinearFillPatch::~D1PiecewiseLinearFillPatch()
{
}

D1PiecewiseLinearFillPatch::D1PiecewiseLinearFillPatch(
    const DisjointBoxLayout& a_fine_domain,
    const DisjointBoxLayout& a_coarse_domain,
    int a_num_comps,
    const Box& a_fine_box,
    int a_dir,
    int a_ref_ratio,
    int a_interp_radius)
:
    m_is_defined(false)
{
    define(a_fine_domain,
        a_coarse_domain,
        a_num_comps,
        a_fine_box,
        a_dir,
        a_ref_ratio,
        a_interp_radius);
}

void
D1PiecewiseLinearFillPatch::define(
    const DisjointBoxLayout& a_fine_domain,
    const DisjointBoxLayout& a_coarse_domain,
    int a_num_comps,
    const Box& a_fine_box,
    int a_dir,
    int a_ref_ratio,
    int a_interp_radius
    )
{
    CH_TIME("D1PiecewiseLinearFillPatch::define");
    m_ref_ratio = a_ref_ratio;
    m_interp_radius = a_interp_radius;
    m_fine_box = a_fine_box;
    m_dir = a_dir;
    m_crse_box.define(m_fine_box.smallEnd()/a_ref_ratio,
        (m_fine_box.bigEnd()+IntVect::Unit-BASISV(m_dir))/a_ref_ratio-(IntVect::Unit-BASISV(m_dir)),
        BASISV(m_dir));

    // const ProblemDomain  fine_problem_domain = refine(m_crse_problem_domain,
    //     m_ref_ratio);

    // quick sanity checks
    //JK CH_assert (a_fine_domain.checkPeriodic(fine_problem_domain));
    if (a_coarse_domain.isClosed())
    {
        //JK CH_assert (a_coarse_domain.checkPeriodic(a_crse_problem_domain));

        //
        // create the work array
        DisjointBoxLayout coarsened_fine_domain;
        coarsen ( coarsened_fine_domain,
            a_fine_domain,
            m_ref_ratio );
        //JK pout() << a_fine_domain << endl;
        //JK pout() << ":::::" << endl;
        //JK pout() << coarsened_fine_domain << endl;


        const int coarse_slope_radius =
            (m_interp_radius + m_ref_ratio - 1) / m_ref_ratio;
        const int coarse_ghost_radius = coarse_slope_radius + s_stencil_radius;

        m_coarsenCopier.define(a_coarse_domain,
            coarsened_fine_domain,
            coarse_ghost_radius * (IntVect::Unit-BASISV(m_dir)));

        const IntVect coarse_slope = coarse_slope_radius * (IntVect::Unit-BASISV(m_dir));
        // (wasteful) extra storage here, but who cares?
        for(int dir=0; dir<3; dir++)
            m_slopes[dir].define(coarsened_fine_domain,
                a_num_comps,
                coarse_slope);
        const IntVect coarse_ghost = coarse_ghost_radius * (IntVect::Unit-BASISV(m_dir));
        m_coarsened_fine_data.define(coarsened_fine_domain,
            a_num_comps,
            coarse_ghost);

        // allocate intvectsets
        m_fine_interp.define(a_fine_domain);
        for (int dir = 0; dir < SpaceDim; ++dir)
        {
            m_coarse_centered_interp[dir].define(coarsened_fine_domain);
            m_coarse_lo_interp[dir].define(coarsened_fine_domain);
            m_coarse_hi_interp[dir].define(coarsened_fine_domain);
        }
        // Compute intvectsets. We do this mostly in the coarsened domain, with
        // only an intersection with the fine domain at the end.

        // first, create a box which will determine whether a given box
        // adjoins a periodic boundary
        //JK CANNOT HANDLE PERIODIC!!!
        //JK Box periodicTestBox(m_crse_problem_domain.domainBox());
        //JK if (m_crse_problem_domain.isPeriodic())
        //JK {
        //JK     for (int idir=0; idir<SpaceDim; idir++)
        //JK     {
        //JK         if (m_crse_problem_domain.isPeriodic(idir))
        //JK             periodicTestBox.grow(idir,-1);
        //JK     }
        //JK }

        // create regions in which to interpolate, then intersect with borders
        // to form regions for centered, one-sided low, and one-sided high
        // differences, per coordinate direction.
        
        DataIterator dit = coarsened_fine_domain.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
            const Box& fine_box = a_fine_domain[dit()];
            Box coarsened_fine_box = m_crse_box
                & coarsen(grow(fine_box, m_interp_radius*(IntVect::Unit-BASISV(m_dir))),m_ref_ratio);
            IntVectSet coarsened_fine_interp(coarsened_fine_box);
        
            // Iterate over boxes in coarsened fine domain, and subtract off
            // from the set of coarse cells from which the fine ghost cells
            // will be interpolated.
        
            LayoutIterator other_lit = coarsened_fine_domain.layoutIterator();
            for (other_lit.begin(); other_lit.ok(); ++other_lit)
            {
                const Box& other_coarsened_box
                    = coarsened_fine_domain.get(other_lit());
                coarsened_fine_interp -= other_coarsened_box;
                // also need to remove periodic images from list of cells
                // to be filled, since they will be filled through exchange
                // as well
                //JK CANNOT HANDLE PERIODIC!!!
                //JK if (m_crse_problem_domain.isPeriodic()
                //JK     && !periodicTestBox.contains(other_coarsened_box)
                //JK     && !periodicTestBox.contains(coarsened_fine_box))
                //JK {
                //JK     ShiftIterator shiftIt = m_crse_problem_domain.shiftIterator();
                //JK     IntVect shiftMult(m_crse_problem_domain.domainBox().size());
                //JK     Box shiftedBox(other_coarsened_box);
                //JK     for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                //JK     {
                //JK         IntVect shiftVect = shiftMult*shiftIt();
                //JK         shiftedBox.shift(shiftVect);
                //JK         coarsened_fine_interp -= shiftedBox;
                //JK         shiftedBox.shift(-shiftVect);
                //JK     }
                //JK }
            }
        
            // Now that we have the coarsened cells required for interpolation,
            // construct IntvectSets specifying the one-sided and centered
            // stencil locations.
        
            for (int dir = 0; dir < SpaceDim; ++dir)
            {
                if(dir != m_dir)
                {
                    IntVectSet& coarse_centered_interp
                        = m_coarse_centered_interp[dir][dit()];
                    coarse_centered_interp = coarsened_fine_interp;

                    IntVectSet& coarse_lo_interp = m_coarse_lo_interp[dir][dit()];
                    coarse_lo_interp = coarse_centered_interp;
                    coarse_lo_interp.shift(BASISV(dir));

                    IntVectSet& coarse_hi_interp = m_coarse_hi_interp[dir][dit()];
                    coarse_hi_interp = coarse_centered_interp;
                    coarse_hi_interp.shift(-BASISV(dir));

                    // We iterate over the coarse grids and subtract them off of the
                    // one-sided stencils.
                    LayoutIterator coarse_lit = a_coarse_domain.layoutIterator();
                    for (coarse_lit.begin();coarse_lit.ok();++coarse_lit)
                    {
                        Box bx = a_coarse_domain.get(coarse_lit());
                        coarse_lo_interp -= bx;
                        coarse_hi_interp -= bx;
                        //JK // once again, need to do periodic images, too
                        //JK CANNOT HANDLE PERIODIC!!!
                        //JK if (m_crse_problem_domain.isPeriodic()
                        //JK     && !periodicTestBox.contains(bx)
                        //JK     && !periodicTestBox.contains(coarsened_fine_box))
                        //JK {
                        //JK     ShiftIterator shiftIt = m_crse_problem_domain.shiftIterator();
                        //JK     IntVect shiftMult(m_crse_problem_domain.domainBox().size());
                        //JK     Box shiftedBox(bx);
                        //JK     for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                        //JK     {
                        //JK         IntVect shiftVect = shiftMult*shiftIt();
                        //JK         shiftedBox.shift(shiftVect);
                        //JK         coarse_lo_interp -= shiftedBox;
                        //JK         coarse_hi_interp -= shiftedBox;
                        //JK         shiftedBox.shift(-shiftVect);
                        //JK     }
                        //JK }
                    }

                    coarse_lo_interp.shift(-BASISV(dir));
                    coarse_hi_interp.shift(BASISV(dir));
                    coarse_centered_interp -= coarse_lo_interp;
                    coarse_centered_interp -= coarse_hi_interp;
                }
        
            }
            // Finally, we construct the fine cells that are going to be
            // interpolated by intersecting them with the refined version of the
            // coarse IntVectSet.
        
            IntVectSet& fine_interp = m_fine_interp[dit()];
            fine_interp = refine(coarsened_fine_interp,m_ref_ratio);
            fine_interp &= m_fine_box & grow(fine_box, m_interp_radius*(IntVect::Unit - BASISV(m_dir)));
        
        }
        m_is_defined = true;
    } // end if coarser level is well-defined
}

bool
D1PiecewiseLinearFillPatch::isDefined() const
{
    return ( m_is_defined );
}

// fill the interpolation region of the fine level ghost cells
void
D1PiecewiseLinearFillPatch::fillInterp(
    LevelData<FArrayBox>& a_fine_data,
    const LevelData<FArrayBox>& a_old_coarse_data,
    const LevelData<FArrayBox>& a_new_coarse_data,
    Real a_time_interp_coef,
    int a_src_comp,
    int a_dest_comp,
    int a_num_comp
    )
{
    CH_TIME("D1PiecewiseLinearFillPatch::fillInterp");
    // sanity checks
    CH_assert (m_is_defined);
    CH_assert (a_time_interp_coef >= 0.);
    CH_assert (a_time_interp_coef <= 1.);
    const DisjointBoxLayout oldCrseGrids = a_old_coarse_data.getBoxes();
    const DisjointBoxLayout newCrseGrids = a_new_coarse_data.getBoxes();
    const DisjointBoxLayout fineGrids = a_fine_data.getBoxes();

    //JK CH_assert (oldCrseGrids.checkPeriodic(m_crse_problem_domain));
    //JK CH_assert (newCrseGrids.checkPeriodic(m_crse_problem_domain));
    //JK CH_assert (fineGrids.checkPeriodic(refine(m_crse_problem_domain,
    //JK             m_ref_ratio)));


    // time interpolation of coarse level data, to coarsened fine level work array
    timeInterp(a_old_coarse_data,
        a_new_coarse_data,
        a_time_interp_coef,
        a_src_comp,
        a_dest_comp,
        a_num_comp);

    // piecewise contant interpolation, from coarsened fine level work
    // array, to fine level
    fillConstantInterp(a_fine_data,
        a_src_comp,
        a_dest_comp,
        a_num_comp);
    //
    // increment fine level data with per-direction linear terms
    computeSlopes(a_src_comp,
        a_num_comp);

    incrementLinearInterp(a_fine_data,
        a_src_comp,
        a_dest_comp,
        a_num_comp);
}

//
// time interpolation of coarse level data, to coarsened fine level work array
void
D1PiecewiseLinearFillPatch::timeInterp(
    const LevelData<FArrayBox>& a_old_coarse_data,
    const LevelData<FArrayBox>& a_new_coarse_data,
    Real a_time_interp_coef,
    int a_src_comp,
    int a_dest_comp,
    int a_num_comp
    )
{
    CH_TIME("D1PiecewiseLinearFillPatch::timeInterp");
    Interval src_interval (a_src_comp,  a_src_comp  + a_num_comp - 1);
    Interval dest_interval(a_dest_comp, a_dest_comp + a_num_comp - 1);

    if ( (a_old_coarse_data.boxLayout().size() == 0)  &&
        (a_new_coarse_data.boxLayout().size() == 0) )
    {
        MayDay::Error ( "D1PiecewiseLinearFillPatch::fillInterp: no old coarse data and no new coarse data" );
    }
    else if ( (a_time_interp_coef == 1.) ||
        (a_old_coarse_data.boxLayout().size() == 0) )
    {
        // old coarse data is absent, or fine time level is the new coarse time level
        a_new_coarse_data.copyTo( src_interval,
            m_coarsened_fine_data,
            dest_interval,
            m_coarsenCopier
            );
    }
    else if ( (a_time_interp_coef == 0.) ||
        (a_new_coarse_data.boxLayout().size() == 0) )
    {
        // new coarse data is absent, or fine time level is the old coarse time level
        a_old_coarse_data.copyTo( src_interval,
            m_coarsened_fine_data,
            dest_interval,
            m_coarsenCopier
            );
    }
    else
    {
        // linearly interpolate between old and new time levels
        a_new_coarse_data.copyTo( src_interval,
            m_coarsened_fine_data,
            dest_interval,
            m_coarsenCopier
            );
        const DisjointBoxLayout&
            coarsened_fine_layout = m_coarsened_fine_data.disjointBoxLayout();
        LevelData<FArrayBox>
            tmp_coarsened_fine_data(coarsened_fine_layout,
                m_coarsened_fine_data.nComp(),
                m_coarsened_fine_data.ghostVect());

        a_old_coarse_data.copyTo( src_interval,
            tmp_coarsened_fine_data,
            dest_interval,
            m_coarsenCopier
            );
        DataIterator dit = coarsened_fine_layout.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
            FArrayBox& coarsened_fine_fab = m_coarsened_fine_data[dit()];
            FArrayBox& tmp_coarsened_fine_fab = tmp_coarsened_fine_data[dit()];
            coarsened_fine_fab.mult(a_time_interp_coef,a_src_comp,a_num_comp);
            tmp_coarsened_fine_fab.mult(1.0 - a_time_interp_coef,a_src_comp,a_num_comp); 
            coarsened_fine_fab.plus(tmp_coarsened_fine_fab,a_src_comp,a_dest_comp,a_num_comp);
        }
    }
}

// fill the fine interpolation region piecewise-constantly
void
D1PiecewiseLinearFillPatch::fillConstantInterp(
    LevelData<FArrayBox>& a_fine_data,
    int a_src_comp,
    int a_dest_comp,
    int a_num_comp
    )
const
{
    CH_TIME("D1PiecewiseLinearFillPatch::fillConstantInterp");
    DataIterator dit = a_fine_data.boxLayout().dataIterator();

    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox& fine_fab = a_fine_data[dit()];
        const FArrayBox& coarse_fab = m_coarsened_fine_data[dit()];
        const IntVectSet& local_fine_interp = m_fine_interp[dit()];
        IVSIterator ivsit(local_fine_interp);
        for (ivsit.begin(); ivsit.ok(); ++ivsit)
        {
            const IntVect& fine_iv = ivsit();
            IntVect coarse_iv = coarsen(fine_iv, m_ref_ratio);
            int coarse_comp = a_src_comp;
            int fine_comp   = a_dest_comp;
            for(; coarse_comp < a_src_comp + a_num_comp; ++fine_comp, ++coarse_comp)
                fine_fab(fine_iv, fine_comp) = coarse_fab(coarse_iv, coarse_comp);
        }
    }
}

// compute slopes at the coarse interpolation sites in the specified direction
void
D1PiecewiseLinearFillPatch::computeSlopes(int a_src_comp,
    int a_num_comp)
{
    CH_TIME("D1PiecewiseLinearFillPatch::computeSlopes");

    DataIterator dit = m_coarsened_fine_data.boxLayout().dataIterator();

    for (int dir=0; dir<SpaceDim; dir++)
    {
        for (dit.begin(); dit.ok(); ++dit)
        {
            if(dir == m_dir)
            {
                m_slopes[m_dir][dit()].setVal(0.0);
                continue;
            }
            FArrayBox& slope_fab = m_slopes[dir][dit()];
            const FArrayBox& dataFab = m_coarsened_fine_data[dit()];
            const IntVectSet& local_centered_interp = m_coarse_centered_interp[dir][dit()];
            const IntVectSet& local_lo_interp = m_coarse_lo_interp[dir][dit()];
            const IntVectSet& local_hi_interp = m_coarse_hi_interp[dir][dit()];

            computeSimpleSlopesFab(slope_fab,
                a_src_comp,
                a_num_comp,
                dir,
                dataFab,
                local_centered_interp,
                local_lo_interp,
                local_hi_interp);
        }
    } 

    for (dit.begin(); dit.ok(); ++dit)
    {
        const Box& slopeBox = m_slopes[0][dit()].box();
        const FArrayBox& dataFab = m_coarsened_fine_data[dit()];
        computeMultiDimSlopes(m_slopes[0][dit()],
            m_slopes[1][dit()] ,
            m_slopes[2][dit()],
            dataFab,
            a_src_comp,
            a_num_comp,
            slopeBox);
    }
}

void D1PiecewiseLinearFillPatch::computeSimpleSlopesFab(FArrayBox & a_slopeFab,
    const int       & a_src_comp,
    const int       & a_num_comp,
    const int       & a_dir,
    const FArrayBox & a_dataFab,
    const IntVectSet& a_local_centered_interp,
    const IntVectSet& a_local_lo_interp,
    const IntVectSet& a_local_hi_interp)
{
    CH_TIME("D1PiecewiseLinearFillPatch::computeSimpleSlopes");
    CH_assert(a_dir != m_dir);
    // van leer limited central difference
    IVSIterator centered_ivsit(a_local_centered_interp);
    for (centered_ivsit.begin(); centered_ivsit.ok(); ++centered_ivsit)
    {
        const IntVect& iv = centered_ivsit();
        const IntVect ivlo = iv - BASISV(a_dir);
        const IntVect ivhi = iv + BASISV(a_dir);
        for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
        {
            Real dcenter = 0.5 * (a_dataFab(ivhi,comp) - a_dataFab(ivlo,comp));
            a_slopeFab(iv,comp) = dcenter;
        }
    }

    // one-sided difference (low)
    IVSIterator lo_ivsit(a_local_lo_interp);
    for (lo_ivsit.begin(); lo_ivsit.ok(); ++lo_ivsit)
    {
        const IntVect& iv = lo_ivsit();
        const IntVect ivlo = iv - BASISV(a_dir);
        for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
        {
            Real dlo = a_dataFab(iv,comp) - a_dataFab(ivlo,comp);
            a_slopeFab(iv,comp) = dlo;
        }
    }

    // one-sided difference (high)
    IVSIterator hi_ivsit(a_local_hi_interp);
    for (hi_ivsit.begin(); hi_ivsit.ok(); ++hi_ivsit)
    {
        const IntVect& iv = hi_ivsit();
        const IntVect ivhi = iv + BASISV(a_dir);
        for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
        {
            Real dhi = a_dataFab(ivhi,comp) - a_dataFab(iv,comp);
            a_slopeFab(iv,comp) = dhi;
        }
    }
}

void  D1PiecewiseLinearFillPatch::computeMultiDimSlopes(FArrayBox      & a_slopes0,
    FArrayBox      & a_slopes1,
    FArrayBox      & a_slopes2,
    const FArrayBox& a_dataFab,
    const int      & a_src_comp,
    const int      & a_num_comp,
    const Box      & a_slopeBox)
{
    // this is the same stuff that is in FineInterp.cpp
    Box b_mod(a_slopeBox);
    b_mod.grow(1);
    b_mod = m_crse_box & b_mod;
    b_mod.grow(-1);

    // create a box big enough to remove periodic BCs from domain
    Box domBox = grow(a_slopeBox,2);
    domBox = m_crse_box & domBox;

    // to do limits, we need to have a box which includes
    // the neighbors of a given point (to check for the
    // local maximum...
    //JK Remove the direction we are not interpolating in
    Box neighborBox(-1*(IntVect::Unit-BASISV(m_dir)),
        (IntVect::Unit-BASISV(m_dir)));


    FORT_INTERPLIMIT( CHF_FRA(a_slopes0),
        CHF_FRA(a_slopes1),
        CHF_FRA(a_slopes2),
        CHF_CONST_FRA(a_dataFab),
        CHF_BOX(b_mod),
        CHF_BOX(neighborBox),
        CHF_BOX(domBox));
}
// increment the fine interpolation sites with linear term for the
// specified coordinate direction
void
D1PiecewiseLinearFillPatch::incrementLinearInterp(
    LevelData<FArrayBox>& a_fine_data,
    int a_src_comp,
    int a_dest_comp,
    int a_num_comp
    )
const
{
    CH_TIME("D1PiecewiseLinearFillPatch::incrementLinearInterp");
    for (int dir=0; dir<SpaceDim; dir++)
    {
        DataIterator dit = a_fine_data.boxLayout().dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
            const FArrayBox& slope_fab = m_slopes[dir][dit()];
            FArrayBox& fine_data_fab = a_fine_data[dit()];
            const IntVectSet& fine_interp = m_fine_interp[dit()];
            IVSIterator ivsit(fine_interp);
            for (ivsit.begin(); ivsit.ok(); ++ivsit)
            {
                const IntVect& fine_iv = ivsit();
                const IntVect coarse_iv = coarsen(fine_iv,m_ref_ratio);
                const int offset = fine_iv[dir] - m_ref_ratio * coarse_iv[dir];
                Real interp_coef = -.5 + (offset +.5) / m_ref_ratio;
                int coarse_comp = a_src_comp;
                int fine_comp   = a_dest_comp;
                for (; coarse_comp < a_src_comp + a_num_comp; ++coarse_comp, ++fine_comp)
                {
                    fine_data_fab(fine_iv,fine_comp)
                        += interp_coef * slope_fab(coarse_iv,coarse_comp);
                }
            }
        }
    }
}

void
D1PiecewiseLinearFillPatch::printIntVectSets() const
{
    DataIterator lit = m_fine_interp.boxLayout().dataIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
        cout << "grid " << lit().intCode() << ": " << endl;
        cout << "fine ivs" << endl;
        cout << m_fine_interp[lit()] << endl;

        for (int dir = 0; dir < SpaceDim; ++dir)
        {
            cout << "coarse centered ivs [" << dir << "]: " << endl;
            cout << m_coarse_centered_interp[dir][lit()] << endl;
            cout << "coarse lo ivs [" << dir << "]: " << endl;
            cout << m_coarse_lo_interp[dir][lit()] << endl;
            cout << "coarse hi ivs [" << dir << "]: " << endl;
            cout << m_coarse_hi_interp[dir][lit()] << endl;
        }

    }
}
#include "NamespaceFooter.H"
