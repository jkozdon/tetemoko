#ifdef CH_LANG_CC
/*
 * Indexing for the linear elasticity and plasticity code
 */
#endif


#ifndef _LEINDEX_H_
#define _LEINDEX_H_
#include "BaseNamespaceHeader.H"

#define IX_VX      ( 0)
#define IX_VY      ( 1)
#define IX_VZ      ( 2)
#define IX_SXX     ( 3)
#define IX_SYY     ( 4)
#define IX_SZZ     ( 5)
#define IX_SXY     ( 6)
#define IX_SXZ     ( 7)
#define IX_SYZ     ( 8)
#define IX_LAM     ( 9)
#define IX_GAM     (10)
#define IX_OFFSET  (11)

#define IX1_VX     (IX_VX )
#define IX1_VY     (IX_VY )
#define IX1_VZ     (IX_VZ )
#define IX1_SXX    (IX_SXX)
#define IX1_SYY    (IX_SYY)
#define IX1_SZZ    (IX_SZZ)
#define IX1_SXY    (IX_SXY)
#define IX1_SXZ    (IX_SXZ)
#define IX1_SYZ    (IX_SYZ)
#define IX1_LAM    (IX_LAM)
#define IX1_GAM    (IX_GAM)

#define IX2_VX     (IX_OFFSET + IX_VX )
#define IX2_VY     (IX_OFFSET + IX_VY )
#define IX2_VZ     (IX_OFFSET + IX_VZ )
#define IX2_SXX    (IX_OFFSET + IX_SXX)
#define IX2_SYY    (IX_OFFSET + IX_SYY)
#define IX2_SZZ    (IX_OFFSET + IX_SZZ)
#define IX2_SXY    (IX_OFFSET + IX_SXY)
#define IX2_SXZ    (IX_OFFSET + IX_SXZ)
#define IX2_SYZ    (IX_OFFSET + IX_SYZ)
#define IX2_LAM    (IX_OFFSET + IX_LAM)
#define IX2_GAM    (IX_OFFSET + IX_GAM)

#include "BaseNamespaceFooter.H"
#endif /*_LEINDEX_H_*/
