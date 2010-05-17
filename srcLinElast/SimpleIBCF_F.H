#ifndef _SIMPLEIBCF_F_H_
#define _SIMPLEIBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSIMPLESETF 
#define GUARDSIMPLESETF 
// Prototype for Fortran procedure SIMPLESETF ...
//
void FORTRAN_NAME( SIMPLESETF ,simplesetf )(
      CHFp_CONST_REAL(aucs)
      ,CHFp_CONST_REAL(aucp)
      ,CHFp_CONST_REAL(aumu) );

#define FORT_SIMPLESETF FORTRAN_NAME( inlineSIMPLESETF, inlineSIMPLESETF)
#define FORTNT_SIMPLESETF FORTRAN_NAME( SIMPLESETF, simplesetf)

inline void FORTRAN_NAME(inlineSIMPLESETF, inlineSIMPLESETF)(
      CHFp_CONST_REAL(aucs)
      ,CHFp_CONST_REAL(aucp)
      ,CHFp_CONST_REAL(aumu) )
{
 CH_TIMELEAF("FORT_SIMPLESETF");
 FORTRAN_NAME( SIMPLESETF ,simplesetf )(
      CHFt_CONST_REAL(aucs)
      ,CHFt_CONST_REAL(aucp)
      ,CHFt_CONST_REAL(aumu) );
}
#endif  // GUARDSIMPLESETF 

}

#endif
