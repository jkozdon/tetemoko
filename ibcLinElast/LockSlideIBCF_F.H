#ifndef _LOCKSLIDEIBCF_F_H_
#define _LOCKSLIDEIBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDLOCKSLIDESETF 
#define GUARDLOCKSLIDESETF 
// Prototype for Fortran procedure LOCKSLIDESETF ...
//
void FORTRAN_NAME( LOCKSLIDESETF ,lockslidesetf )(
      CHFp_CONST_REAL(ar0)
      ,CHFp_CONST_REAL(aedge) );

#define FORT_LOCKSLIDESETF FORTRAN_NAME( inlineLOCKSLIDESETF, inlineLOCKSLIDESETF)
#define FORTNT_LOCKSLIDESETF FORTRAN_NAME( LOCKSLIDESETF, lockslidesetf)

inline void FORTRAN_NAME(inlineLOCKSLIDESETF, inlineLOCKSLIDESETF)(
      CHFp_CONST_REAL(ar0)
      ,CHFp_CONST_REAL(aedge) )
{
 CH_TIMELEAF("FORT_LOCKSLIDESETF");
 FORTRAN_NAME( LOCKSLIDESETF ,lockslidesetf )(
      CHFt_CONST_REAL(ar0)
      ,CHFt_CONST_REAL(aedge) );
}
#endif  // GUARDLOCKSLIDESETF 

#ifndef GUARDLOCKSLIDEINITF 
#define GUARDLOCKSLIDEINITF 
// Prototype for Fortran procedure LOCKSLIDEINITF ...
//
void FORTRAN_NAME( LOCKSLIDEINITF ,lockslideinitf )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_LOCKSLIDEINITF FORTRAN_NAME( inlineLOCKSLIDEINITF, inlineLOCKSLIDEINITF)
#define FORTNT_LOCKSLIDEINITF FORTRAN_NAME( LOCKSLIDEINITF, lockslideinitf)

inline void FORTRAN_NAME(inlineLOCKSLIDEINITF, inlineLOCKSLIDEINITF)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LOCKSLIDEINITF");
 FORTRAN_NAME( LOCKSLIDEINITF ,lockslideinitf )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLOCKSLIDEINITF 

#ifndef GUARDLOCKSLIDEFAULTBCF 
#define GUARDLOCKSLIDEFAULTBCF 
// Prototype for Fortran procedure LOCKSLIDEFAULTBCF ...
//
void FORTRAN_NAME( LOCKSLIDEFAULTBCF ,lockslidefaultbcf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WInside)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_LOCKSLIDEFAULTBCF FORTRAN_NAME( inlineLOCKSLIDEFAULTBCF, inlineLOCKSLIDEFAULTBCF)
#define FORTNT_LOCKSLIDEFAULTBCF FORTRAN_NAME( LOCKSLIDEFAULTBCF, lockslidefaultbcf)

inline void FORTRAN_NAME(inlineLOCKSLIDEFAULTBCF, inlineLOCKSLIDEFAULTBCF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WInside)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LOCKSLIDEFAULTBCF");
 FORTRAN_NAME( LOCKSLIDEFAULTBCF ,lockslidefaultbcf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WInside)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(t)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLOCKSLIDEFAULTBCF 

}

#endif
