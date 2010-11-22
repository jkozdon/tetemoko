#ifndef _PSEUDOPULSEIBCF_F_H_
#define _PSEUDOPULSEIBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDPSEUDOPULSESETF 
#define GUARDPSEUDOPULSESETF 
// Prototype for Fortran procedure PSEUDOPULSESETF ...
//
void FORTRAN_NAME( PSEUDOPULSESETF ,pseudopulsesetf )(
      CHFp_CONST_REAL(ar0)
      ,CHFp_CONST_REAL(aslope)
      ,CHFp_CONST_REAL(adTau)
      ,CHFp_CONST_REAL(aVh)
      ,CHFp_CONST_REAL(aVr)
      ,CHFp_CONST_REAL(aJumpTime)
      ,CHFp_CONST_REAL(aStopTime) );

#define FORT_PSEUDOPULSESETF FORTRAN_NAME( inlinePSEUDOPULSESETF, inlinePSEUDOPULSESETF)
#define FORTNT_PSEUDOPULSESETF FORTRAN_NAME( PSEUDOPULSESETF, pseudopulsesetf)

inline void FORTRAN_NAME(inlinePSEUDOPULSESETF, inlinePSEUDOPULSESETF)(
      CHFp_CONST_REAL(ar0)
      ,CHFp_CONST_REAL(aslope)
      ,CHFp_CONST_REAL(adTau)
      ,CHFp_CONST_REAL(aVh)
      ,CHFp_CONST_REAL(aVr)
      ,CHFp_CONST_REAL(aJumpTime)
      ,CHFp_CONST_REAL(aStopTime) )
{
 CH_TIMELEAF("FORT_PSEUDOPULSESETF");
 FORTRAN_NAME( PSEUDOPULSESETF ,pseudopulsesetf )(
      CHFt_CONST_REAL(ar0)
      ,CHFt_CONST_REAL(aslope)
      ,CHFt_CONST_REAL(adTau)
      ,CHFt_CONST_REAL(aVh)
      ,CHFt_CONST_REAL(aVr)
      ,CHFt_CONST_REAL(aJumpTime)
      ,CHFt_CONST_REAL(aStopTime) );
}
#endif  // GUARDPSEUDOPULSESETF 

#ifndef GUARDPSEUDOPULSEINITF 
#define GUARDPSEUDOPULSEINITF 
// Prototype for Fortran procedure PSEUDOPULSEINITF ...
//
void FORTRAN_NAME( PSEUDOPULSEINITF ,pseudopulseinitf )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_PSEUDOPULSEINITF FORTRAN_NAME( inlinePSEUDOPULSEINITF, inlinePSEUDOPULSEINITF)
#define FORTNT_PSEUDOPULSEINITF FORTRAN_NAME( PSEUDOPULSEINITF, pseudopulseinitf)

inline void FORTRAN_NAME(inlinePSEUDOPULSEINITF, inlinePSEUDOPULSEINITF)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PSEUDOPULSEINITF");
 FORTRAN_NAME( PSEUDOPULSEINITF ,pseudopulseinitf )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPSEUDOPULSEINITF 

#ifndef GUARDPSEUDOPULSESETBND 
#define GUARDPSEUDOPULSESETBND 
// Prototype for Fortran procedure PSEUDOPULSESETBND ...
//
void FORTRAN_NAME( PSEUDOPULSESETBND ,pseudopulsesetbnd )(
      CHFp_FRA(B)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_REAL(dt) );

#define FORT_PSEUDOPULSESETBND FORTRAN_NAME( inlinePSEUDOPULSESETBND, inlinePSEUDOPULSESETBND)
#define FORTNT_PSEUDOPULSESETBND FORTRAN_NAME( PSEUDOPULSESETBND, pseudopulsesetbnd)

inline void FORTRAN_NAME(inlinePSEUDOPULSESETBND, inlinePSEUDOPULSESETBND)(
      CHFp_FRA(B)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_REAL(dt) )
{
 CH_TIMELEAF("FORT_PSEUDOPULSESETBND");
 FORTRAN_NAME( PSEUDOPULSESETBND ,pseudopulsesetbnd )(
      CHFt_FRA(B)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_REAL(dt) );
}
#endif  // GUARDPSEUDOPULSESETBND 

#ifndef GUARDPSEUDOPULSEFAULTBCF 
#define GUARDPSEUDOPULSEFAULTBCF 
// Prototype for Fortran procedure PSEUDOPULSEFAULTBCF ...
//
void FORTRAN_NAME( PSEUDOPULSEFAULTBCF ,pseudopulsefaultbcf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WInside)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_PSEUDOPULSEFAULTBCF FORTRAN_NAME( inlinePSEUDOPULSEFAULTBCF, inlinePSEUDOPULSEFAULTBCF)
#define FORTNT_PSEUDOPULSEFAULTBCF FORTRAN_NAME( PSEUDOPULSEFAULTBCF, pseudopulsefaultbcf)

inline void FORTRAN_NAME(inlinePSEUDOPULSEFAULTBCF, inlinePSEUDOPULSEFAULTBCF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WInside)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PSEUDOPULSEFAULTBCF");
 FORTRAN_NAME( PSEUDOPULSEFAULTBCF ,pseudopulsefaultbcf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WInside)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(t)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPSEUDOPULSEFAULTBCF 

}

#endif
