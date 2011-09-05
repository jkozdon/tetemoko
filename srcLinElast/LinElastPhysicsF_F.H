#ifndef _LINELASTPHYSICSF_F_H_
#define _LINELASTPHYSICSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMAXWAVESPEEDF 
#define GUARDMAXWAVESPEEDF 
// Prototype for Fortran procedure MAXWAVESPEEDF ...
//
void FORTRAN_NAME( MAXWAVESPEEDF ,maxwavespeedf )(
      CHFp_REAL(cmax)
      ,CHFp_CONST_FRA(U)
      ,CHFp_BOX(box) );

#define FORT_MAXWAVESPEEDF FORTRAN_NAME( inlineMAXWAVESPEEDF, inlineMAXWAVESPEEDF)
#define FORTNT_MAXWAVESPEEDF FORTRAN_NAME( MAXWAVESPEEDF, maxwavespeedf)

inline void FORTRAN_NAME(inlineMAXWAVESPEEDF, inlineMAXWAVESPEEDF)(
      CHFp_REAL(cmax)
      ,CHFp_CONST_FRA(U)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAXWAVESPEEDF");
 FORTRAN_NAME( MAXWAVESPEEDF ,maxwavespeedf )(
      CHFt_REAL(cmax)
      ,CHFt_CONST_FRA(U)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAXWAVESPEEDF 

#ifndef GUARDMAGNITUDEF 
#define GUARDMAGNITUDEF 
// Prototype for Fortran procedure MAGNITUDEF ...
//
void FORTRAN_NAME( MAGNITUDEF ,magnitudef )(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) );

#define FORT_MAGNITUDEF FORTRAN_NAME( inlineMAGNITUDEF, inlineMAGNITUDEF)
#define FORTNT_MAGNITUDEF FORTRAN_NAME( MAGNITUDEF, magnitudef)

inline void FORTRAN_NAME(inlineMAGNITUDEF, inlineMAGNITUDEF)(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAGNITUDEF");
 FORTRAN_NAME( MAGNITUDEF ,magnitudef )(
      CHFt_FRA1(magdata)
      ,CHFt_CONST_FRA(data)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAGNITUDEF 

#ifndef GUARDMAXF 
#define GUARDMAXF 
// Prototype for Fortran procedure MAXF ...
//
void FORTRAN_NAME( MAXF ,maxf )(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) );

#define FORT_MAXF FORTRAN_NAME( inlineMAXF, inlineMAXF)
#define FORTNT_MAXF FORTRAN_NAME( MAXF, maxf)

inline void FORTRAN_NAME(inlineMAXF, inlineMAXF)(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAXF");
 FORTRAN_NAME( MAXF ,maxf )(
      CHFt_FRA1(magdata)
      ,CHFt_CONST_FRA(data)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAXF 

#ifndef GUARDGETRELGRADF 
#define GUARDGETRELGRADF 
// Prototype for Fortran procedure GETRELGRADF ...
//
void FORTRAN_NAME( GETRELGRADF ,getrelgradf )(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_GETRELGRADF FORTRAN_NAME( inlineGETRELGRADF, inlineGETRELGRADF)
#define FORTNT_GETRELGRADF FORTRAN_NAME( GETRELGRADF, getrelgradf)

inline void FORTRAN_NAME(inlineGETRELGRADF, inlineGETRELGRADF)(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_GETRELGRADF");
 FORTRAN_NAME( GETRELGRADF ,getrelgradf )(
      CHFt_FRA1(du)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDGETRELGRADF 

#ifndef GUARDBOUNDREFINE 
#define GUARDBOUNDREFINE 
// Prototype for Fortran procedure BOUNDREFINE ...
//
void FORTRAN_NAME( BOUNDREFINE ,boundrefine )(
      CHFp_FRA1(markCells)
      ,CHFp_CONST_REAL(refLocation)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_BOUNDREFINE FORTRAN_NAME( inlineBOUNDREFINE, inlineBOUNDREFINE)
#define FORTNT_BOUNDREFINE FORTRAN_NAME( BOUNDREFINE, boundrefine)

inline void FORTRAN_NAME(inlineBOUNDREFINE, inlineBOUNDREFINE)(
      CHFp_FRA1(markCells)
      ,CHFp_CONST_REAL(refLocation)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_BOUNDREFINE");
 FORTRAN_NAME( BOUNDREFINE ,boundrefine )(
      CHFt_FRA1(markCells)
      ,CHFt_CONST_REAL(refLocation)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDBOUNDREFINE 

#ifndef GUARDGETFULLGRADF 
#define GUARDGETFULLGRADF 
// Prototype for Fortran procedure GETFULLGRADF ...
//
void FORTRAN_NAME( GETFULLGRADF ,getfullgradf )(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_GETFULLGRADF FORTRAN_NAME( inlineGETFULLGRADF, inlineGETFULLGRADF)
#define FORTNT_GETFULLGRADF FORTRAN_NAME( GETFULLGRADF, getfullgradf)

inline void FORTRAN_NAME(inlineGETFULLGRADF, inlineGETFULLGRADF)(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_GETFULLGRADF");
 FORTRAN_NAME( GETFULLGRADF ,getfullgradf )(
      CHFt_FRA1(du)
      ,CHFt_CONST_FRA(u)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDGETFULLGRADF 

#ifndef GUARDGETFLUXF 
#define GUARDGETFLUXF 
// Prototype for Fortran procedure GETFLUXF ...
//
void FORTRAN_NAME( GETFLUXF ,getfluxf )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wgdnv)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_GETFLUXF FORTRAN_NAME( inlineGETFLUXF, inlineGETFLUXF)
#define FORTNT_GETFLUXF FORTRAN_NAME( GETFLUXF, getfluxf)

inline void FORTRAN_NAME(inlineGETFLUXF, inlineGETFLUXF)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wgdnv)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETFLUXF");
 FORTRAN_NAME( GETFLUXF ,getfluxf )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Wgdnv)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETFLUXF 

#ifndef GUARDCHARANALYSISF 
#define GUARDCHARANALYSISF 
// Prototype for Fortran procedure CHARANALYSISF ...
//
void FORTRAN_NAME( CHARANALYSISF ,charanalysisf )(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_CHARANALYSISF FORTRAN_NAME( inlineCHARANALYSISF, inlineCHARANALYSISF)
#define FORTNT_CHARANALYSISF FORTRAN_NAME( CHARANALYSISF, charanalysisf)

inline void FORTRAN_NAME(inlineCHARANALYSISF, inlineCHARANALYSISF)(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARANALYSISF");
 FORTRAN_NAME( CHARANALYSISF ,charanalysisf )(
      CHFt_FRA(dW)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARANALYSISF 

#ifndef GUARDCHARSYNTHESISF 
#define GUARDCHARSYNTHESISF 
// Prototype for Fortran procedure CHARSYNTHESISF ...
//
void FORTRAN_NAME( CHARSYNTHESISF ,charsynthesisf )(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_CHARSYNTHESISF FORTRAN_NAME( inlineCHARSYNTHESISF, inlineCHARSYNTHESISF)
#define FORTNT_CHARSYNTHESISF FORTRAN_NAME( CHARSYNTHESISF, charsynthesisf)

inline void FORTRAN_NAME(inlineCHARSYNTHESISF, inlineCHARSYNTHESISF)(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARSYNTHESISF");
 FORTRAN_NAME( CHARSYNTHESISF ,charsynthesisf )(
      CHFt_FRA(dW)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARSYNTHESISF 

#ifndef GUARDCHARVALUESF 
#define GUARDCHARVALUESF 
// Prototype for Fortran procedure CHARVALUESF ...
//
void FORTRAN_NAME( CHARVALUESF ,charvaluesf )(
      CHFp_FRA(lambda)
      ,CHFp_BOX(box) );

#define FORT_CHARVALUESF FORTRAN_NAME( inlineCHARVALUESF, inlineCHARVALUESF)
#define FORTNT_CHARVALUESF FORTRAN_NAME( CHARVALUESF, charvaluesf)

inline void FORTRAN_NAME(inlineCHARVALUESF, inlineCHARVALUESF)(
      CHFp_FRA(lambda)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARVALUESF");
 FORTRAN_NAME( CHARVALUESF ,charvaluesf )(
      CHFt_FRA(lambda)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARVALUESF 

#ifndef GUARDRIEMANNF 
#define GUARDRIEMANNF 
// Prototype for Fortran procedure RIEMANNF ...
//
void FORTRAN_NAME( RIEMANNF ,riemannf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_RIEMANNF FORTRAN_NAME( inlineRIEMANNF, inlineRIEMANNF)
#define FORTNT_RIEMANNF FORTRAN_NAME( RIEMANNF, riemannf)

inline void FORTRAN_NAME(inlineRIEMANNF, inlineRIEMANNF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANNF");
 FORTRAN_NAME( RIEMANNF ,riemannf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANNF 

#ifndef GUARDCOPYF 
#define GUARDCOPYF 
// Prototype for Fortran procedure COPYF ...
//
void FORTRAN_NAME( COPYF ,copyf )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_BOX(box) );

#define FORT_COPYF FORTRAN_NAME( inlineCOPYF, inlineCOPYF)
#define FORTNT_COPYF FORTRAN_NAME( COPYF, copyf)

inline void FORTRAN_NAME(inlineCOPYF, inlineCOPYF)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_COPYF");
 FORTRAN_NAME( COPYF ,copyf )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCOPYF 

#ifndef GUARDGETADWDXF 
#define GUARDGETADWDXF 
// Prototype for Fortran procedure GETADWDXF ...
//
void FORTRAN_NAME( GETADWDXF ,getadwdxf )(
      CHFp_FRA(AdWdx)
      ,CHFp_CONST_FRA(WHalf)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_GETADWDXF FORTRAN_NAME( inlineGETADWDXF, inlineGETADWDXF)
#define FORTNT_GETADWDXF FORTRAN_NAME( GETADWDXF, getadwdxf)

inline void FORTRAN_NAME(inlineGETADWDXF, inlineGETADWDXF)(
      CHFp_FRA(AdWdx)
      ,CHFp_CONST_FRA(WHalf)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETADWDXF");
 FORTRAN_NAME( GETADWDXF ,getadwdxf )(
      CHFt_FRA(AdWdx)
      ,CHFt_CONST_FRA(WHalf)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETADWDXF 

#ifndef GUARDLINELASTFREEBCF 
#define GUARDLINELASTFREEBCF 
// Prototype for Fortran procedure LINELASTFREEBCF ...
//
void FORTRAN_NAME( LINELASTFREEBCF ,linelastfreebcf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WInside)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_LINELASTFREEBCF FORTRAN_NAME( inlineLINELASTFREEBCF, inlineLINELASTFREEBCF)
#define FORTNT_LINELASTFREEBCF FORTRAN_NAME( LINELASTFREEBCF, linelastfreebcf)

inline void FORTRAN_NAME(inlineLINELASTFREEBCF, inlineLINELASTFREEBCF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WInside)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LINELASTFREEBCF");
 FORTRAN_NAME( LINELASTFREEBCF ,linelastfreebcf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WInside)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLINELASTFREEBCF 

#ifndef GUARDLINELASTOUTBCF 
#define GUARDLINELASTOUTBCF 
// Prototype for Fortran procedure LINELASTOUTBCF ...
//
void FORTRAN_NAME( LINELASTOUTBCF ,linelastoutbcf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WInside)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_LINELASTOUTBCF FORTRAN_NAME( inlineLINELASTOUTBCF, inlineLINELASTOUTBCF)
#define FORTNT_LINELASTOUTBCF FORTRAN_NAME( LINELASTOUTBCF, linelastoutbcf)

inline void FORTRAN_NAME(inlineLINELASTOUTBCF, inlineLINELASTOUTBCF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WInside)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LINELASTOUTBCF");
 FORTRAN_NAME( LINELASTOUTBCF ,linelastoutbcf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WInside)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLINELASTOUTBCF 

#ifndef GUARDLINELASTSETF 
#define GUARDLINELASTSETF 
// Prototype for Fortran procedure LINELASTSETF ...
//
void FORTRAN_NAME( LINELASTSETF ,linelastsetf )(
      CHFp_CONST_REAL(aucs)
      ,CHFp_CONST_REAL(aucp)
      ,CHFp_CONST_REAL(aumu)
      ,CHFp_CONST_REAL(asxx0)
      ,CHFp_CONST_REAL(asyy0)
      ,CHFp_CONST_REAL(aszz0)
      ,CHFp_CONST_REAL(asxy0)
      ,CHFp_CONST_REAL(asxz0)
      ,CHFp_CONST_REAL(asyz0) );

#define FORT_LINELASTSETF FORTRAN_NAME( inlineLINELASTSETF, inlineLINELASTSETF)
#define FORTNT_LINELASTSETF FORTRAN_NAME( LINELASTSETF, linelastsetf)

inline void FORTRAN_NAME(inlineLINELASTSETF, inlineLINELASTSETF)(
      CHFp_CONST_REAL(aucs)
      ,CHFp_CONST_REAL(aucp)
      ,CHFp_CONST_REAL(aumu)
      ,CHFp_CONST_REAL(asxx0)
      ,CHFp_CONST_REAL(asyy0)
      ,CHFp_CONST_REAL(aszz0)
      ,CHFp_CONST_REAL(asxy0)
      ,CHFp_CONST_REAL(asxz0)
      ,CHFp_CONST_REAL(asyz0) )
{
 CH_TIMELEAF("FORT_LINELASTSETF");
 FORTRAN_NAME( LINELASTSETF ,linelastsetf )(
      CHFt_CONST_REAL(aucs)
      ,CHFt_CONST_REAL(aucp)
      ,CHFt_CONST_REAL(aumu)
      ,CHFt_CONST_REAL(asxx0)
      ,CHFt_CONST_REAL(asyy0)
      ,CHFt_CONST_REAL(aszz0)
      ,CHFt_CONST_REAL(asxy0)
      ,CHFt_CONST_REAL(asxz0)
      ,CHFt_CONST_REAL(asyz0) );
}
#endif  // GUARDLINELASTSETF 

#ifndef GUARDLEINITF 
#define GUARDLEINITF 
// Prototype for Fortran procedure LEINITF ...
//
void FORTRAN_NAME( LEINITF ,leinitf )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_LEINITF FORTRAN_NAME( inlineLEINITF, inlineLEINITF)
#define FORTNT_LEINITF FORTRAN_NAME( LEINITF, leinitf)

inline void FORTRAN_NAME(inlineLEINITF, inlineLEINITF)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LEINITF");
 FORTRAN_NAME( LEINITF ,leinitf )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLEINITF 

#ifndef GUARDLINELASTSETFAB 
#define GUARDLINELASTSETFAB 
// Prototype for Fortran procedure LINELASTSETFAB ...
//
void FORTRAN_NAME( LINELASTSETFAB ,linelastsetfab )(
      CHFp_FRA1(U)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(val) );

#define FORT_LINELASTSETFAB FORTRAN_NAME( inlineLINELASTSETFAB, inlineLINELASTSETFAB)
#define FORTNT_LINELASTSETFAB FORTRAN_NAME( LINELASTSETFAB, linelastsetfab)

inline void FORTRAN_NAME(inlineLINELASTSETFAB, inlineLINELASTSETFAB)(
      CHFp_FRA1(U)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(val) )
{
 CH_TIMELEAF("FORT_LINELASTSETFAB");
 FORTRAN_NAME( LINELASTSETFAB ,linelastsetfab )(
      CHFt_FRA1(U)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(val) );
}
#endif  // GUARDLINELASTSETFAB 

#ifndef GUARDGAUSSIANFUNC 
#define GUARDGAUSSIANFUNC 
// Prototype for Fortran procedure GAUSSIANFUNC ...
//
void FORTRAN_NAME( GAUSSIANFUNC ,gaussianfunc )(
      CHFp_REAL(val)
      ,CHFp_CONST_REAL(w)
      ,CHFp_CONST_REAL(r) );

#define FORT_GAUSSIANFUNC FORTRAN_NAME( inlineGAUSSIANFUNC, inlineGAUSSIANFUNC)
#define FORTNT_GAUSSIANFUNC FORTRAN_NAME( GAUSSIANFUNC, gaussianfunc)

inline void FORTRAN_NAME(inlineGAUSSIANFUNC, inlineGAUSSIANFUNC)(
      CHFp_REAL(val)
      ,CHFp_CONST_REAL(w)
      ,CHFp_CONST_REAL(r) )
{
 CH_TIMELEAF("FORT_GAUSSIANFUNC");
 FORTRAN_NAME( GAUSSIANFUNC ,gaussianfunc )(
      CHFt_REAL(val)
      ,CHFt_CONST_REAL(w)
      ,CHFt_CONST_REAL(r) );
}
#endif  // GUARDGAUSSIANFUNC 

#ifndef GUARDSMOOTHFUNC 
#define GUARDSMOOTHFUNC 
// Prototype for Fortran procedure SMOOTHFUNC ...
//
void FORTRAN_NAME( SMOOTHFUNC ,smoothfunc )(
      CHFp_REAL(val)
      ,CHFp_CONST_REAL(w)
      ,CHFp_CONST_REAL(r) );

#define FORT_SMOOTHFUNC FORTRAN_NAME( inlineSMOOTHFUNC, inlineSMOOTHFUNC)
#define FORTNT_SMOOTHFUNC FORTRAN_NAME( SMOOTHFUNC, smoothfunc)

inline void FORTRAN_NAME(inlineSMOOTHFUNC, inlineSMOOTHFUNC)(
      CHFp_REAL(val)
      ,CHFp_CONST_REAL(w)
      ,CHFp_CONST_REAL(r) )
{
 CH_TIMELEAF("FORT_SMOOTHFUNC");
 FORTRAN_NAME( SMOOTHFUNC ,smoothfunc )(
      CHFt_REAL(val)
      ,CHFt_CONST_REAL(w)
      ,CHFt_CONST_REAL(r) );
}
#endif  // GUARDSMOOTHFUNC 

#ifndef GUARDSMOOTHBOXCAR 
#define GUARDSMOOTHBOXCAR 
// Prototype for Fortran procedure SMOOTHBOXCAR ...
//
void FORTRAN_NAME( SMOOTHBOXCAR ,smoothboxcar )(
      CHFp_REAL(val)
      ,CHFp_CONST_REAL(w)
      ,CHFp_CONST_REAL(r)
      ,CHFp_CONST_REAL(p) );

#define FORT_SMOOTHBOXCAR FORTRAN_NAME( inlineSMOOTHBOXCAR, inlineSMOOTHBOXCAR)
#define FORTNT_SMOOTHBOXCAR FORTRAN_NAME( SMOOTHBOXCAR, smoothboxcar)

inline void FORTRAN_NAME(inlineSMOOTHBOXCAR, inlineSMOOTHBOXCAR)(
      CHFp_REAL(val)
      ,CHFp_CONST_REAL(w)
      ,CHFp_CONST_REAL(r)
      ,CHFp_CONST_REAL(p) )
{
 CH_TIMELEAF("FORT_SMOOTHBOXCAR");
 FORTRAN_NAME( SMOOTHBOXCAR ,smoothboxcar )(
      CHFt_REAL(val)
      ,CHFt_CONST_REAL(w)
      ,CHFt_CONST_REAL(r)
      ,CHFt_CONST_REAL(p) );
}
#endif  // GUARDSMOOTHBOXCAR 

#ifndef GUARDBOXCARFUNC 
#define GUARDBOXCARFUNC 
// Prototype for Fortran procedure BOXCARFUNC ...
//
void FORTRAN_NAME( BOXCARFUNC ,boxcarfunc )(
      CHFp_REAL(val)
      ,CHFp_CONST_REAL(w)
      ,CHFp_CONST_REAL(r) );

#define FORT_BOXCARFUNC FORTRAN_NAME( inlineBOXCARFUNC, inlineBOXCARFUNC)
#define FORTNT_BOXCARFUNC FORTRAN_NAME( BOXCARFUNC, boxcarfunc)

inline void FORTRAN_NAME(inlineBOXCARFUNC, inlineBOXCARFUNC)(
      CHFp_REAL(val)
      ,CHFp_CONST_REAL(w)
      ,CHFp_CONST_REAL(r) )
{
 CH_TIMELEAF("FORT_BOXCARFUNC");
 FORTRAN_NAME( BOXCARFUNC ,boxcarfunc )(
      CHFt_REAL(val)
      ,CHFt_CONST_REAL(w)
      ,CHFt_CONST_REAL(r) );
}
#endif  // GUARDBOXCARFUNC 

#ifndef GUARDSMOOTHFUNC2 
#define GUARDSMOOTHFUNC2 
// Prototype for Fortran procedure SMOOTHFUNC2 ...
//
void FORTRAN_NAME( SMOOTHFUNC2 ,smoothfunc2 )(
      CHFp_REAL(val)
      ,CHFp_CONST_REAL(w)
      ,CHFp_CONST_REAL(r) );

#define FORT_SMOOTHFUNC2 FORTRAN_NAME( inlineSMOOTHFUNC2, inlineSMOOTHFUNC2)
#define FORTNT_SMOOTHFUNC2 FORTRAN_NAME( SMOOTHFUNC2, smoothfunc2)

inline void FORTRAN_NAME(inlineSMOOTHFUNC2, inlineSMOOTHFUNC2)(
      CHFp_REAL(val)
      ,CHFp_CONST_REAL(w)
      ,CHFp_CONST_REAL(r) )
{
 CH_TIMELEAF("FORT_SMOOTHFUNC2");
 FORTRAN_NAME( SMOOTHFUNC2 ,smoothfunc2 )(
      CHFt_REAL(val)
      ,CHFt_CONST_REAL(w)
      ,CHFt_CONST_REAL(r) );
}
#endif  // GUARDSMOOTHFUNC2 

#ifndef GUARDSMOOTHRAMP 
#define GUARDSMOOTHRAMP 
// Prototype for Fortran procedure SMOOTHRAMP ...
//
void FORTRAN_NAME( SMOOTHRAMP ,smoothramp )(
      CHFp_REAL(val)
      ,CHFp_REAL(x)
      ,CHFp_REAL(k) );

#define FORT_SMOOTHRAMP FORTRAN_NAME( inlineSMOOTHRAMP, inlineSMOOTHRAMP)
#define FORTNT_SMOOTHRAMP FORTRAN_NAME( SMOOTHRAMP, smoothramp)

inline void FORTRAN_NAME(inlineSMOOTHRAMP, inlineSMOOTHRAMP)(
      CHFp_REAL(val)
      ,CHFp_REAL(x)
      ,CHFp_REAL(k) )
{
 CH_TIMELEAF("FORT_SMOOTHRAMP");
 FORTRAN_NAME( SMOOTHRAMP ,smoothramp )(
      CHFt_REAL(val)
      ,CHFt_REAL(x)
      ,CHFt_REAL(k) );
}
#endif  // GUARDSMOOTHRAMP 

#ifndef GUARDARTDISP 
#define GUARDARTDISP 
// Prototype for Fortran procedure ARTDISP ...
//
void FORTRAN_NAME( ARTDISP ,artdisp )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_ARTDISP FORTRAN_NAME( inlineARTDISP, inlineARTDISP)
#define FORTNT_ARTDISP FORTRAN_NAME( ARTDISP, artdisp)

inline void FORTRAN_NAME(inlineARTDISP, inlineARTDISP)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ARTDISP");
 FORTRAN_NAME( ARTDISP ,artdisp )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_REAL(coeff)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDARTDISP 

}

#endif
