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
      CHFp_CONST_REAL(ar0)
      ,CHFp_CONST_VR(amag)
      ,CHFp_CONST_REAL(asig) );

#define FORT_SIMPLESETF FORTRAN_NAME( inlineSIMPLESETF, inlineSIMPLESETF)
#define FORTNT_SIMPLESETF FORTRAN_NAME( SIMPLESETF, simplesetf)

inline void FORTRAN_NAME(inlineSIMPLESETF, inlineSIMPLESETF)(
      CHFp_CONST_REAL(ar0)
      ,CHFp_CONST_VR(amag)
      ,CHFp_CONST_REAL(asig) )
{
 CH_TIMELEAF("FORT_SIMPLESETF");
 FORTRAN_NAME( SIMPLESETF ,simplesetf )(
      CHFt_CONST_REAL(ar0)
      ,CHFt_CONST_VR(amag)
      ,CHFt_CONST_REAL(asig) );
}
#endif  // GUARDSIMPLESETF 

#ifndef GUARDSIMPLEINITF 
#define GUARDSIMPLEINITF 
// Prototype for Fortran procedure SIMPLEINITF ...
//
void FORTRAN_NAME( SIMPLEINITF ,simpleinitf )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SIMPLEINITF FORTRAN_NAME( inlineSIMPLEINITF, inlineSIMPLEINITF)
#define FORTNT_SIMPLEINITF FORTRAN_NAME( SIMPLEINITF, simpleinitf)

inline void FORTRAN_NAME(inlineSIMPLEINITF, inlineSIMPLEINITF)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SIMPLEINITF");
 FORTRAN_NAME( SIMPLEINITF ,simpleinitf )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSIMPLEINITF 

}

#endif
