#ifndef _D1FINEINTERP_H_
#define _D1FINEINTERP_H_

#include "REAL.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "FArrayBox.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

class DisjointBoxLayout;

/// replaces fine level data with interpolation of coarse level data.
/**
  This class replaces data at a fine level of refinement with data
  interpolated from a coarser level of refinement.  Interpolation is
  piecewise bi(tri)linear, with van Leer slopes if there is room for
  the stencil, with lower-order slopes if there isn't.  See the
  design document.

*/
class D1FineInterp
{
public:
    ///
    /**
      Default constructor.  User must subsequently call define().
      */
    D1FineInterp();

    ///
    /**
      Destructor.
      */
    ~D1FineInterp();

    //JK ///
    //JK /**
    //JK   Defining constructor.  Constructs a valid object.
    //JK   Equivalent to default construction followed by define().

    //JK   {\bf Arguments:}\\
    //JK   a_fine_domain (not modified): the fine level domain.\\
    //JK   a_numcomps (not modified): the number of components.\\
    //JK   a_ref_ratio (not modified): the refinement ratio.\\
    //JK   a_fine_problem_domain (not modified): problem domain at the fine level.\\
    //JK   */
    //JK FineInterp(const DisjointBoxLayout& a_fine_domain,
    //JK     const int& a_numcomps,
    //JK     const int& a_ref_ratio,
    //JK     const Box& a_fine_problem_domain);

    //JK ///
    //JK /**
    //JK   Defining constructor.  Constructs a valid object.
    //JK   Equivalent to default construction followed by define().

    //JK   {\bf Arguments:}\\
    //JK   a_fine_domain (not modified): the fine level domain.\\
    //JK   a_numcomps (not modified): the number of components.\\
    //JK   a_ref_ratio (not modified): the refinement ratio.\\
    //JK   a_fine_problem_domain (not modified): problem domain at the fine level.\\
    //JK   */
    //JK FineInterp(const DisjointBoxLayout& a_fine_domain,
    //JK     const int& a_numcomps,
    //JK     const int& a_ref_ratio,
    //JK     const ProblemDomain& a_fine_problem_domain);

    ///
    /**
      Defines this object.  Existing information is overriden.

      {\bf Arguments:}\\
      a_fine_domain (not modified): the fine level domain.\\
      a_numcomps (not modified): the number of components.\\
      a_ref_ratio (not modified): the refinement ratio.\\
      a_fine_problem_domain (not modified): problem domain at the fine level.\\

      {\bf This:}\\
      ---This object is modified.---
      */
    void
        define(const DisjointBoxLayout& a_fine_domain, // the fine level domain
            const int& a_numcomps,                   // the number of components
            const int& a_ref_ratio,                  // the refinement ratio
            const int& a_dir,                       // the direction we are not interpolating in
            const Box& a_fine_problem_domain);      // problem domain

    //JK ///
    //JK /**
    //JK   Defines this object.  Existing information is overriden.

    //JK   {\bf Arguments:}\\
    //JK   a_fine_domain (not modified): the fine level domain.\\
    //JK   a_numcomps (not modified): the number of components.\\
    //JK   a_ref_ratio (not modified): the refinement ratio.\\
    //JK   a_fine_problem_domain (not modified): problem domain at the fine level.\\

    //JK   {\bf This:}\\
    //JK   ---This object is modified.---
    //JK   */
    //JK void
    //JK     define(const DisjointBoxLayout& a_fine_domain, // the fine level domain
    //JK         const int& a_numcomps,                   // the number of components
    //JK         const int& a_ref_ratio,                  // the refinement ratio
    //JK         const ProblemDomain& a_fine_problem_domain);

    ///
    /**
      Returns true if this object was created with the defining
      constructor or if define() has been called.

      {\bf This:}\\
      This object is not modified.
      */
    bool
        isDefined() const;

    ///
    /**
      Replaces a_fine_data with data interpolated from a_coarse_data. It
      is an error to call if not this->isDefined().  The domain of
      a_fine_data should be the same as the fine domain specified in the
      most recent call to define().  It is expected that the coarse and
      fine level's domains are properly nested.  Both a_coarse_data and
      a_fine_data should have the same number of components specified in
      the most recent call to define().

      {\bf Arguments:}\\
      a_fine_data (modified): fine data. \\
      a_coarse_data (not modified): coarse data. \\

      {\bf This:}\\
      Well, it's complicated.  As far as the user is concerned, this object
      is not modified.  See the design document if you care for details.
      */
    void
        interpToFine(LevelData<FArrayBox>& a_fine_data,
            const LevelData<FArrayBox>& a_coarse_data);

    //JK /// Just do piecewise-constant interpolation.
    //JK void
    //JK     pwcinterpToFine(LevelData<FArrayBox>& a_fine_data,
    //JK         const LevelData<FArrayBox>& a_coarse_data);

protected:
    void
        interpGridData(BaseFab<Real>& a_fine,
            const BaseFab<Real>& a_coarse,
            const Box& a_coarsened_fine_box,
            int a_ref_ratio)
        const;
    //JK void
    //JK     pwcinterpGridData(BaseFab<Real>& a_fine,
    //JK         const BaseFab<Real>& a_coarse,
    //JK         const Box& a_coarsened_fine_box,
    //JK         int a_ref_ratio) const;

protected:
    bool is_defined;
    // the refinement ratio
    int m_ref_ratio;
    int m_numcomps;
    int m_dir;
    Box m_fine_box;
    Box m_coarse_box;
    // work array for the coarse level data in a domain that is the
    // outline of the fine level domain on the coarse level
    LevelData<FArrayBox> m_coarsened_fine_data;
    // coarse level problem domain
    ProblemDomain m_coarse_problem_domain;
};

#include "NamespaceFooter.H"
#endif
