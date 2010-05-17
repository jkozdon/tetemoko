#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "LoHiSide.H"

#include "LinElastPhysics.H"

//JK #include "PolytropicPhysicsF_F.H"


/// Constructor
LinElastPhysics::LinElastPhysics(const int junk)
{
}

/// Destructor
LinElastPhysics::~LinElastPhysics()
{
}

/// Compute the maximum wave speed
Real LinElastPhysics::getMaxWaveSpeed(const FArrayBox& a_U,
    const Box&       a_box)
{
    pout() << "NOT DEFINED :: LinElastPhysics::getMaxWaveSpeed" << endl;
        CH_assert(isDefined());
    CH_assert(a_U.contains(a_box));

    Real speed = 0.0;

    //JK FORT_MAXWAVESPEEDF(CHF_REAL(speed),
    //JK                    CHF_CONST_FRA(a_U),
    //JK                    CHF_BOX(a_box));

    return speed;
}

/// Object factory for this class
GodunovPhysics* LinElastPhysics::new_godunovPhysics() const
{
    CH_assert(m_isBCSet);

    GodunovPhysics* retval = static_cast<GodunovPhysics*>
        (new LinElastPhysics(0 /*JUNK*/));

    retval->setPhysIBC(m_bc);
    return retval;
}


/// Number of conserved variables
/**
  Return the number of conserved variables.
  */
int LinElastPhysics::numConserved()
{
  CH_assert(isDefined());

  return 9;
}

/// Names of the conserved variables
/**
  Return the names of the conserved variables.
  */
Vector<string> LinElastPhysics::stateNames()
{
    CH_assert(isDefined());

    Vector<string> retval;

    retval.push_back("v_x");
    retval.push_back("v_y");
    retval.push_back("v_z");
    retval.push_back("sig_xx");
    retval.push_back("sig_yy");
    retval.push_back("sig_zz");
    retval.push_back("sig_xy");
    retval.push_back("sig_xz");
    retval.push_back("sig_yz");

    return retval;
}

/// Number of flux variables
/**
  Return the  number of flux variables.  This can be greater than the number
  of conserved variables if addition fluxes/face-centered quantities are
  computed.
  */
int LinElastPhysics::numFluxes()
{
    CH_assert(isDefined());

    // In some computations there may be more fluxes than conserved variables
    return numConserved();
}

/// Compute a flux from primitive variable values on a face
/**
*/
void LinElastPhysics::getFlux(FArrayBox&       a_flux,
    const FArrayBox& a_whalf,
    const int&       a_dir,
    const Box&       a_box)
{
    pout() << "NOT DEFINED :: LinElastPhysics::getFlux" << endl;
    CH_assert(isDefined());

    //JK FORT_GETFLUXF(CHF_FRA(a_flux),
    //JK     CHF_CONST_FRA(a_whalf),
    //JK     CHF_CONST_INT(a_dir),
    //JK     CHF_BOX(a_box));
}

/// Number of primitive variables
/**
  Return the number of primitive variables.  This may be greater than the
  number of conserved variables if derived/redundant quantities are also
  stored for convenience.
  */
int LinElastPhysics::numPrimitives()
{
    CH_assert(isDefined());
    pout() << "MAY NOT BE RIGHT :: LinElastPhysics::getFlux" << endl;

    //JK NOT SURE IF THIS SHOULD BE THE SAME
    return numConserved();
}

/// Transform a_dW from primitive to characteristic variables
/**
  On input, a_dW contains the increments of the primitive variables. On
  output, it contains the increments in the characteristic variables.

  IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
  smallest eigenvalue first, the largest eigenvalue last, and orders the
  characteristic variables accordingly.
  */
void LinElastPhysics::charAnalysis(FArrayBox&       a_dW,
    const FArrayBox& a_W,
    const int&       a_dir,
    const Box&       a_box)
{
    pout() << "NOT DEFINED :: LinElastPhysics::charAnalysis" << endl;
    CH_assert(isDefined());

    //JK FORT_CHARANALYSISF(CHF_FRA(a_dW),
    //JK     CHF_CONST_FRA(a_W),
    //JK     CHF_CONST_INT(a_dir),
    //JK     CHF_BOX(a_box));
}

/// Transform a_dW from characteristic to primitive variables
/**
  On input, a_dW contains the increments of the characteristic variables.
  On output, it contains the increments in the primitive variables.

  IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
  smallest eigenvalue first, the largest eigenvalue last, and orders the
  characteristic variables accordingly.
  */
void LinElastPhysics::charSynthesis(FArrayBox&       a_dW,
    const FArrayBox& a_W,
    const int&       a_dir,
    const Box&       a_box)
{
    pout() << "NOT DEFINED :: LinElastPhysics::charSynthesis" << endl;
    CH_assert(isDefined());

    //JK FORT_CHARSYNTHESISF(CHF_FRA(a_dW),
    //JK     CHF_CONST_FRA(a_W),
    //JK     CHF_CONST_INT(a_dir),
    //JK     CHF_BOX(a_box));
}

/// Compute the characteristic values (eigenvalues)
/**
  Compute the characteristic values (eigenvalues)

  IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
  smallest eigenvalue first, the largest eigenvalue last, and orders the
  characteristic variables accordingly.
  */
void LinElastPhysics::charValues(FArrayBox&       a_lambda,
    const FArrayBox& a_W,
    const int&       a_dir,
    const Box&       a_box)
{
    pout() << "NOT DEFINED :: LinElastPhysics::charValues" << endl;
    CH_assert(isDefined());

    //JK FORT_CHARVALUESF(CHF_FRA(a_lambda),
    //JK     CHF_CONST_FRA(a_W),
    //JK     CHF_CONST_INT(a_dir),
    //JK     CHF_BOX(a_box));
}

/// Add to (increment) the source terms given the current state
/**
  On input, a_S contains the current source terms.  On output, a_S has
  had any additional source terms (based on the current state, a_W)
  added to it.  This should all be done on the region defined by a_box.
  */
void LinElastPhysics::incrementSource(FArrayBox&       a_S,
    const FArrayBox& a_W,
    const Box&       a_box)
{
    pout() << "MAY NOT BE RIGHT :: LinElastPhysics::incrementSource" << endl;
}

/// Compute the solution to the Riemann problem.
/**
  Given input left and right states in a direction, a_dir, compute a
  Riemann problem and generate fluxes at the faces within a_box.
  */
void LinElastPhysics::riemann(/// face-centered solution to Riemann problem
    FArrayBox&       a_WStar,
    /// left state, on cells to left of each face
    const FArrayBox& a_WLeft,
    /// right state, on cells to right of each face
    const FArrayBox& a_WRight,
    /// state on cells, used to set boundary conditions
    const FArrayBox& a_W,
    /// current time
    const Real&      a_time,
    /// direction of faces
    const int&       a_dir,
    /// face-centered box on which to set a_WStar
    const Box&       a_box)
{
    pout() << "NOT DEFINED :: LinElastPhysics::riemann" << endl;
}

/// Post-normal predictor calculation.
/**
  Add increment to normal predictor, e.g. to account for source terms due to
  spatially-varying coefficients, to bound primitive variable ranges.
  */
void LinElastPhysics::postNormalPred(FArrayBox&       a_dWMinus,
    FArrayBox&       a_dWPlus,
    const FArrayBox& a_W,
    const Real&      a_dt,
    const Real&      a_dx,
    const int&       a_dir,
    const Box&       a_box)
{
    pout() << "NOT DEFINED :: LinElastPhysics::postNormalPred" << endl;
    // Bound "a_dWMinus" and "a_dWPlus" so that density and pressure will never
    // be less than "smallr" and "smallp", respectively
    //JK FORT_POSTNORMALPREDF(CHF_FRA(a_dWMinus),
    //JK     CHF_FRA(a_dWPlus),
    //JK     CHF_CONST_FRA(a_W),
    //JK     CHF_BOX(a_box));
}

/// Compute the quasilinear update A*dW/dx.
/**
*/
void LinElastPhysics::quasilinearUpdate(FArrayBox&       a_dWdx,
    const FArrayBox& a_wHalf,
    const FArrayBox& a_W,
    const Real&      a_scale,
    const int&       a_dir,
    const Box&       a_box)
{
    pout() << "NOT DEFINED :: LinElastPhysics::quasilinear" << endl;
    CH_assert(isDefined());
    CH_assert(a_dWdx.box().contains(a_box));

    //JK FORT_GETADWDXF(CHF_FRA(a_dWdx),
    //JK     CHF_CONST_FRA(a_WHalf),
    //JK     CHF_CONST_FRA(a_W),
    //JK     CHF_CONST_REAL(a_scale),
    //JK     CHF_CONST_INT(a_dir),
    //JK     CHF_BOX(a_box));
}

/// Compute primitive variables from conserved variables.
/**
*/
void LinElastPhysics::consToPrim(FArrayBox&       a_W,
    const FArrayBox& a_U,
    const Box&       a_box)
{
    pout() << "NOT DEFINED :: LinElastPhysics::consToPrim" << endl;
    CH_assert(isDefined());
    CH_assert(a_U.box().contains(a_box));
    CH_assert(a_W.box().contains(a_box));

    //JK FORT_CONSTOPRIMF(CHF_FRA(a_W),
    //JK     CHF_CONST_FRA(a_U),
    //JK     CHF_BOX(a_box));
}

/// Interval within the primitive variables corresponding to the velocities
/**
  Return the interval of component indices within the primitive variable
  of the velocities.  Used for slope flattening (slope computation) and
  computing the divergence of the velocity (artificial viscosity).
  */
Interval LinElastPhysics::velocityInterval()
{
    MayDay::Error("LinElastPhysics::velocityInterval - not defined");

    Interval retval(-1,-1);
    return retval;
}

/// Component index within the primitive variables of the pressure
/**
  Return the component index withn the primitive variables for the
  pressure.  Used for slope flattening (slope computation).
  */
int LinElastPhysics::pressureIndex()
{
    MayDay::Error("LinElastPhysics::pressureIndex - not defined");

    return -1;
}

/// Used to limit the absolute value of a "pressure" difference
/**
  Return a value that is used by slope flattening to limit (away from
  zero) the absolute value of a slope in the pressureIndex() component
  (slope computation).
  */
Real LinElastPhysics::smallPressure()
{
    MayDay::Error("LinElastPhysics::smallPressure - not defined");

    return -1.0;
}

/// Component index within the primitive variables of the bulk modulus
/**
  Return the component index withn the primitive variables for the
  bulk modulus.  Used for slope flattening (slope computation) used
  as a normalization to measure shock strength.
  */
int LinElastPhysics::bulkModulusIndex()
{
    MayDay::Error("LinElastPhysics::bulkModulusIndex - not defined");

    return -1;
}

#ifdef CH_USE_HDF5
void LinElastPhysics::expressions(HDF5HeaderData& a_expressions) const
{
  //JK a_expressions.m_string["scalar gamma"] = "1.4";

  //JK a_expressions.m_string["vector velocity"] = "momentum/density";
  //JK a_expressions.m_string["scalar kinetic_energy"] = "dot(velocity,velocity)/2";
  //JK a_expressions.m_string["scalar pressure"] = "(gamma-1)*(<energy-density>-kinetic_energy*density)";
  //JK a_expressions.m_string["scalar soundspeed"] = "sqrt(gamma*(pressure/density))";
  //JK a_expressions.m_string["scalar log10entropy"] = "log10(pressure) - gamma*log10(<density>)";
}

#endif
