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
      ,CHFp_CONST_REAL(asigma)
      ,CHFp_CONST_REAL(antime)
      ,CHFp_CONST_REAL(aRSa)
      ,CHFp_CONST_REAL(aRSb)
      ,CHFp_CONST_REAL(aRSV0)
      ,CHFp_CONST_REAL(aRSf0)
      ,CHFp_CONST_REAL(aRSL)
      ,CHFp_CONST_REAL(aRSfw)
      ,CHFp_CONST_REAL(aRSVw) );

#define FORT_RSSETF FORTRAN_NAME( inlineRSSETF, inlineRSSETF)
#define FORTNT_RSSETF FORTRAN_NAME( RSSETF, rssetf)

inline void FORTRAN_NAME(inlineRSSETF, inlineRSSETF)(
      CHFp_CONST_REAL(ar0)
      ,CHFp_CONST_REAL(asigma)
      ,CHFp_CONST_REAL(antime)
      ,CHFp_CONST_REAL(aRSa)
      ,CHFp_CONST_REAL(aRSb)
      ,CHFp_CONST_REAL(aRSV0)
      ,CHFp_CONST_REAL(aRSf0)
      ,CHFp_CONST_REAL(aRSL)
      ,CHFp_CONST_REAL(aRSfw)
      ,CHFp_CONST_REAL(aRSVw) )
{
 CH_TIMELEAF("FORT_RSSETF");
 FORTRAN_NAME( RSSETF ,rssetf )(
      CHFt_CONST_REAL(ar0)
      ,CHFt_CONST_REAL(asigma)
      ,CHFt_CONST_REAL(antime)
      ,CHFt_CONST_REAL(aRSa)
      ,CHFt_CONST_REAL(aRSb)
      ,CHFt_CONST_REAL(aRSV0)
      ,CHFt_CONST_REAL(aRSf0)
      ,CHFt_CONST_REAL(aRSL)
      ,CHFt_CONST_REAL(aRSfw)
      ,CHFt_CONST_REAL(aRSVw) );
}
#endif  // GUARDRSSETF 

#ifndef GUARDRSINITF 
#define GUARDRSINITF 
// Prototype for Fortran procedure RSINITF ...
//
void FORTRAN_NAME( RSINITF ,rsinitf )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_RSINITF FORTRAN_NAME( inlineRSINITF, inlineRSINITF)
#define FORTNT_RSINITF FORTRAN_NAME( RSINITF, rsinitf)

inline void FORTRAN_NAME(inlineRSINITF, inlineRSINITF)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RSINITF");
 FORTRAN_NAME( RSINITF ,rsinitf )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRSINITF 

#ifndef GUARDRSSETFAB 
#define GUARDRSSETFAB 
// Prototype for Fortran procedure RSSETFAB ...
//
void FORTRAN_NAME( RSSETFAB ,rssetfab )(
      CHFp_FRA1(U)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(val) );

#define FORT_RSSETFAB FORTRAN_NAME( inlineRSSETFAB, inlineRSSETFAB)
#define FORTNT_RSSETFAB FORTRAN_NAME( RSSETFAB, rssetfab)

inline void FORTRAN_NAME(inlineRSSETFAB, inlineRSSETFAB)(
      CHFp_FRA1(U)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(val) )
{
 CH_TIMELEAF("FORT_RSSETFAB");
 FORTRAN_NAME( RSSETFAB ,rssetfab )(
      CHFt_FRA1(U)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(val) );
}
#endif  // GUARDRSSETFAB 

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
      CHFp_CONST_REAL(psi)
      ,CHFp_REAL(vx)
      ,CHFp_CONST_REAL(vy)
      ,CHFp_REAL(vz)
      ,CHFp_REAL(sxy)
      ,CHFp_CONST_REAL(syy)
      ,CHFp_REAL(syz) );

#define FORT_RSNEWTON FORTRAN_NAME( inlineRSNEWTON, inlineRSNEWTON)
#define FORTNT_RSNEWTON FORTRAN_NAME( RSNEWTON, rsnewton)

inline void FORTRAN_NAME(inlineRSNEWTON, inlineRSNEWTON)(
      CHFp_CONST_REAL(psi)
      ,CHFp_REAL(vx)
      ,CHFp_CONST_REAL(vy)
      ,CHFp_REAL(vz)
      ,CHFp_REAL(sxy)
      ,CHFp_CONST_REAL(syy)
      ,CHFp_REAL(syz) )
{
 CH_TIMELEAF("FORT_RSNEWTON");
 FORTRAN_NAME( RSNEWTON ,rsnewton )(
      CHFt_CONST_REAL(psi)
      ,CHFt_REAL(vx)
      ,CHFt_CONST_REAL(vy)
      ,CHFt_REAL(vz)
      ,CHFt_REAL(sxy)
      ,CHFt_CONST_REAL(syy)
      ,CHFt_REAL(syz) );
}
#endif  // GUARDRSNEWTON 

}

#endif
