#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "REAL.H"
#include "FArrayBox.H"
#include "DataIterator.H"
#include "LayoutIterator.H"
#include "parstream.H"
#include "AverageF_F.H"

#include "D1CoarseAverage.H"
#include "NamespaceHeader.H"

D1CoarseAverage::D1CoarseAverage()
  :
    is_defined(false)
{
}

D1CoarseAverage::~D1CoarseAverage()
{
}

  void
D1CoarseAverage::define(const DisjointBoxLayout& a_fine_domain,
    int a_numcomps,
    int a_ref_ratio,
    int a_dir)
{
  CH_TIME("D1CoarseAverage::define1");
  m_ref_ratio = a_ref_ratio;
  m_dir = a_dir;
  m_numcomps = a_numcomps;

  const DisjointBoxLayout& constGrids = a_fine_domain;
  DisjointBoxLayout coarsened_fine_domain;

  // manually coarsen the domain all all dims except m_dir
  coarsened_fine_domain.deepCopy(constGrids);
  for(LayoutIterator lit = a_fine_domain.layoutIterator(); lit.ok(); ++lit)
  {
    const Box tmpFineBox = a_fine_domain[lit()];
    const Box tmpCoarseBox(tmpFineBox.smallEnd()/m_ref_ratio,
        (tmpFineBox.bigEnd()+IntVect::Unit-BASISV(m_dir))/m_ref_ratio-(IntVect::Unit-BASISV(m_dir)),
        BASISV(m_dir));
    // coarsened_fine_domain[lit()] = tmpCoarseBox; //3.1
    coarsened_fine_domain.ref(lit()) = tmpCoarseBox;
  }
  coarsened_fine_domain.close();

  // define
  m_coarsened_fine_data.define ( coarsened_fine_domain,
      m_numcomps);

  m_is_copier_defined = false;
  is_defined = true;
}

bool
D1CoarseAverage::isDefined() const
{
  return ( is_defined );
}

  void
D1CoarseAverage::averageToCoarse(LevelData<FArrayBox>& a_coarse_data,
    const LevelData<FArrayBox>& a_fine_data)
{
  computeAverages(a_coarse_data, a_fine_data, arithmetic);
}

  void
D1CoarseAverage::computeAverages(LevelData<FArrayBox>& a_coarse_data,
    const LevelData<FArrayBox>& a_fine_data,
    int a_averageType)
{
  CH_TIME("D1CoarseAverage::computeAverages");
  CH_assert(is_defined);
  // it would be nice if this could check for validity of a_fine_data.
  // this could be done with a redundant DisjointBoxLayout
  DataIterator dit = a_fine_data.boxLayout().dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    BaseFab<Real>& coarsened_fine = m_coarsened_fine_data[dit()];
    const BaseFab<Real>& fine = a_fine_data[dit()];
    // coarsenGridData coarsens from the entire fine grid onto the entire
    // coarse grid.
    averageGridData(coarsened_fine,
        fine,
        m_ref_ratio,
        a_averageType);
  }

  if (m_is_copier_defined)
  {
    // we can use the pre-defined copier to make things faster
    m_coarsened_fine_data.copyTo(m_coarsened_fine_data.interval(),
        a_coarse_data,
        a_coarse_data.interval(),
        m_copier);
  }
  else
  {
    m_coarsened_fine_data.copyTo(m_coarsened_fine_data.interval(),
        a_coarse_data,
        a_coarse_data.interval() );
  }
}

  void
D1CoarseAverage::averageGridData(BaseFab<Real>& a_coarse,
    const BaseFab<Real>& a_fine,
    int a_ref_ratio,
    int a_averageType)
const
{
  const Box& b = a_coarse.box();
  Box refbox(IntVect::Zero,
      (a_ref_ratio-1)*(IntVect::Unit-BASISV(m_dir)),
      BASISV(m_dir));
  if (a_averageType == arithmetic)
  {
    FORT_D1AVERAGE( CHF_FRA(a_coarse),
        CHF_CONST_FRA(a_fine),
        CHF_BOX(b),
        CHF_CONST_INT(a_ref_ratio),
        CHF_BOX(refbox)
        );
  }
  //JK else if (a_averageType == harmonic)
  //JK {
  //JK     FORT_AVERAGEHARMONIC( CHF_FRA(a_coarse),
  //JK         CHF_CONST_FRA(a_fine),
  //JK         CHF_BOX(b),
  //JK         CHF_CONST_INT(a_ref_ratio),
  //JK         CHF_BOX(refbox)
  //JK         );
  //JK }
  else
  {
    MayDay::Error("D1CoarseAverage::averageGridData -- bad averageType");
  }
}
#include "NamespaceFooter.H"
