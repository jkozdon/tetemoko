#ifndef _RSIBCF_F_H_
#define _RSIBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDRSSETF 
#define GUARDRSSETF 
// Prototype for Fortran procedure RSSETF ...
//
void FORTRAN_NAME( RSSETF ,rssetf )(
      CHFp_CONST_REAL(ar0)
      ,CHFp_CONST_REAL(ax0)
      ,CHFp_CONST_REAL(ay0)
      ,CHFp_CONST_REAL(asigma)
      ,CHFp_CONST_REAL(antime)
      ,CHFp_CONST_REAL(aRSa)
      ,CHFp_CONST_REAL(aRSb)
      ,CHFp_CONST_REAL(aRSV0)
      ,CHFp_CONST_REAL(aRSf0)
      ,CHFp_CONST_REAL(aRSL)
      ,CHFp_CONST_REAL(aRSfw)
      ,CHFp_CONST_REAL(aRSVw)
      ,CHFp_CONST_REAL(afExp)
      ,CHFp_CONST_REAL(aruptureThreshold) );

#define FORT_RSSETF FORTRAN_NAME( inlineRSSETF, inlineRSSETF)
#define FORTNT_RSSETF FORTRAN_NAME( RSSETF, rssetf)

inline void FORTRAN_NAME(inlineRSSETF, inlineRSSETF)(
      CHFp_CONST_REAL(ar0)
      ,CHFp_CONST_REAL(ax0)
      ,CHFp_CONST_REAL(ay0)
      ,CHFp_CONST_REAL(asigma)
      ,CHFp_CONST_REAL(antime)
      ,CHFp_CONST_REAL(aRSa)
      ,CHFp_CONST_REAL(aRSb)
      ,CHFp_CONST_REAL(aRSV0)
      ,CHFp_CONST_REAL(aRSf0)
      ,CHFp_CONST_REAL(aRSL)
      ,CHFp_CONST_REAL(aRSfw)
      ,CHFp_CONST_REAL(aRSVw)
      ,CHFp_CONST_REAL(afExp)
      ,CHFp_CONST_REAL(aruptureThreshold) )
{
 CH_TIMELEAF("FORT_RSSETF");
 FORTRAN_NAME( RSSETF ,rssetf )(
      CHFt_CONST_REAL(ar0)
      ,CHFt_CONST_REAL(ax0)
      ,CHFt_CONST_REAL(ay0)
      ,CHFt_CONST_REAL(asigma)
      ,CHFt_CONST_REAL(antime)
      ,CHFt_CONST_REAL(aRSa)
      ,CHFt_CONST_REAL(aRSb)
      ,CHFt_CONST_REAL(aRSV0)
      ,CHFt_CONST_REAL(aRSf0)
      ,CHFt_CONST_REAL(aRSL)
      ,CHFt_CONST_REAL(aRSfw)
      ,CHFt_CONST_REAL(aRSVw)
      ,CHFt_CONST_REAL(afExp)
      ,CHFt_CONST_REAL(aruptureThreshold) );
}
#endif  // GUARDRSSETF 

#ifndef GUARDRSSETBND 
#define GUARDRSSETBND 
// Prototype for Fortran procedure RSSETBND ...
//
void FORTRAN_NAME( RSSETBND ,rssetbnd )(
      CHFp_FRA(Bout)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(Bin)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(time) );

#define FORT_RSSETBND FORTRAN_NAME( inlineRSSETBND, inlineRSSETBND)
#define FORTNT_RSSETBND FORTRAN_NAME( RSSETBND, rssetbnd)

inline void FORTRAN_NAME(inlineRSSETBND, inlineRSSETBND)(
      CHFp_FRA(Bout)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FRA(Bin)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_RSSETBND");
 FORTRAN_NAME( RSSETBND ,rssetbnd )(
      CHFt_FRA(Bout)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FRA(Bin)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDRSSETBND 

#ifndef GUARDRSFAULTBCF 
#define GUARDRSFAULTBCF 
// Prototype for Fortran procedure RSFAULTBCF ...
//
void FORTRAN_NAME( RSFAULTBCF ,rsfaultbcf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WInside)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(GPsi)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_RSFAULTBCF FORTRAN_NAME( inlineRSFAULTBCF, inlineRSFAULTBCF)
#define FORTNT_RSFAULTBCF FORTRAN_NAME( RSFAULTBCF, rsfaultbcf)

inline void FORTRAN_NAME(inlineRSFAULTBCF, inlineRSFAULTBCF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WInside)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(GPsi)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RSFAULTBCF");
 FORTRAN_NAME( RSFAULTBCF ,rsfaultbcf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WInside)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA1(GPsi)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(t)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRSFAULTBCF 

#ifndef GUARDRSNEWTON 
#define GUARDRSNEWTON 
// Prototype for Fortran procedure RSNEWTON ...
//
void FORTRAN_NAME( RSNEWTON ,rsnewton )(
      CHFp_CONST_REAL(phi)
      ,CHFp_CONST_REAL(psi)
      ,CHFp_CONST_REAL(sn)
      ,CHFp_REAL(V)
      ,CHFp_REAL(tau) );

#define FORT_RSNEWTON FORTRAN_NAME( inlineRSNEWTON, inlineRSNEWTON)
#define FORTNT_RSNEWTON FORTRAN_NAME( RSNEWTON, rsnewton)

inline void FORTRAN_NAME(inlineRSNEWTON, inlineRSNEWTON)(
      CHFp_CONST_REAL(phi)
      ,CHFp_CONST_REAL(psi)
      ,CHFp_CONST_REAL(sn)
      ,CHFp_REAL(V)
      ,CHFp_REAL(tau) )
{
 CH_TIMELEAF("FORT_RSNEWTON");
 FORTRAN_NAME( RSNEWTON ,rsnewton )(
      CHFt_CONST_REAL(phi)
      ,CHFt_CONST_REAL(psi)
      ,CHFt_CONST_REAL(sn)
      ,CHFt_REAL(V)
      ,CHFt_REAL(tau) );
}
#endif  // GUARDRSNEWTON 

}

#endif
