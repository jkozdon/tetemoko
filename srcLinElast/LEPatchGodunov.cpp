#include "LoHiSide.H"

#include "PatchGodunov.H"
#include "LEPatchGodunov.H"
#include "LinElastPhysics.H"
#include "NamespaceHeader.H"


// Define this object and the boundary condition object
void LEPatchGodunov::define(const ProblemDomain&           a_domain,
    const Real&                    a_dx,
    const GodunovPhysics* const    a_physPtr,
    const int&                     a_normalPredOrder,
    const bool&                    a_useFourthOrderSlopes,
    const bool&                    a_usePrimLimiting,
    const bool&                    a_useCharLimiting,
    const bool&                    a_useFlattening,
    const bool&                    a_useArtificialViscosity,
    const Real&                    a_artificialViscosity)
{
    MayDay::Error("Must use new define with LinElastPhysics not GodunovPhysics");
}

// Define this object and the boundary condition object
// We just want to force the user to pass the right physics object
void LEPatchGodunov::define(const ProblemDomain&           a_domain,
    const Real&                    a_dx,
    const LinElastPhysics* const   a_physPtr,
    const int&                     a_normalPredOrder,
    const bool&                    a_useFourthOrderSlopes,
    const bool&                    a_usePrimLimiting,
    const bool&                    a_useCharLimiting,
    const bool&                    a_useFlattening,
    const bool&                    a_useArtificialViscosity,
    const Real&                    a_artificialViscosity)
{
    PatchGodunov::define(a_domain,
        a_dx,
        (GodunovPhysics*) a_physPtr,
        a_normalPredOrder,
        a_useFourthOrderSlopes,
        a_usePrimLimiting,
        a_useCharLimiting,
        a_useFlattening,
        a_useArtificialViscosity,
        a_artificialViscosity);
}

void LEPatchGodunov::updateStateBND(FArrayBox&       a_U,
    FArrayBox&       a_Psi,
    FluxBox&         a_F,
    Real&            a_maxWaveSpeed,
    const FArrayBox& a_S,
    const Real&      a_dt,
    const Box&       a_box,
    const Box&       a_bdryBox)
{
    //JK pout() << "updateStateBND" << endl;
    CH_assert(isDefined());
    CH_assert(a_box == m_currentBox);

    int numPrim = m_gdnvPhysics->numPrimitives();
    int numFlux = m_gdnvPhysics->numFluxes();

    FluxBox whalf(a_box,numPrim);
    whalf.setVal(0.0);

    a_F.resize(a_box,numFlux);
    a_F.setVal(0.0);

    computeWHalfBND(whalf, a_U, a_Psi, a_S, a_dt, a_box, a_bdryBox);

    FArrayBox dU(a_U.box(),a_U.nComp());
    computeUpdate(dU, a_F, a_U, whalf, a_dt, a_box);

    a_U += dU;

    // Get and return the maximum wave speed on this patch/grid
    a_maxWaveSpeed = m_gdnvPhysics->getMaxWaveSpeed(a_U, m_currentBox);
}

// Compute the time-centered values of the primitive variable on the
// interior cell faces of the cell-centered box a_box.
void LEPatchGodunov::computeWHalfBND(FluxBox&         a_WHalf,
    const FArrayBox& a_U,
    const FArrayBox& a_Psi,
    const FArrayBox& a_S,
    const Real&      a_dt,
    const Box&       a_box,
    const Box&       a_bdryBox)
{
    //JK pout() << "computeWHalfBND" << endl;
    CH_assert(isDefined());
    CH_assert(a_box == m_currentBox);

    // Get the number of various variables
    int numPrim = m_gdnvPhysics->numPrimitives();

    // The current box of valid cells
    Box curBox = m_currentBox;

    // Boxes for face-centered state - used for the riemann() and
    // artificialViscosity() calls
    Box faceBox[SpaceDim];

    // Boxes for face-centered fluxes
    Box fluxBox[SpaceDim];

    // Boxes for cell-centered state - used for the updatePrim() calls
    Box ccBox[SpaceDim];

    for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
    {
        // faceBox[dir1] is face-centered in direction "dir1",
        // is valid "curBox" grown by 1 in all directions except "dir1",
        // and stays one cell away, in "dir1", from the domain boundary.
        faceBox[dir1] = curBox; // valid cell-centered box
        faceBox[dir1].grow(1);
        faceBox[dir1] &= m_domain;
        faceBox[dir1].grow(dir1, -1);
        faceBox[dir1].surroundingNodes(dir1);

        // ccBox[dir1] is cell-centered,
        // is valid "curBox" grown by 1 in all directions except "dir1",
        // but only within the domain, "m_domain".
        ccBox[dir1] = curBox;
        ccBox[dir1].grow(1);
        ccBox[dir1].grow(dir1, -1);
        ccBox[dir1] &= m_domain;

        // fluxBox[dir1] is face-centered in direction "dir1",
        // consisting of all those faces of cells of "ccBox[dir1]".
        fluxBox[dir1] = ccBox[dir1];
        fluxBox[dir1].surroundingNodes(dir1);

        // The difference between faceBox[dir1] and fluxBox[dir1] is:
        // if curBox abuts the boundary of the domain in direction "dir1",
        // then fluxBox[dir1] contains faces along that boundary,
        // but faceBox[dir1] does not.
    }

    // cell-centered box:  but restrict it to within domain
    Box WBox = a_U.box();
    WBox &= m_domain;

    // Primitive variables
    FArrayBox W(WBox,numPrim);

    // Calculate the primitive variables W from the conserved variables a_U,
    // on cell-centered WBox.
    m_gdnvPhysics->consToPrim(W, a_U, WBox);

    // slopeBox is the cell-centered box where slopes will be needed:
    // it is one larger than the final update box "curBox" of valid cells,
    // but within domain.
    // On slopeBox we define the FABs flattening, WMinus, and WPlus.
    Box slopeBox = curBox;
    slopeBox.grow(1);
    slopeBox &= m_domain;

    // Compute flattening once for all slopes if needed
    FArrayBox flattening(slopeBox, 1); // cell-centered
    if (m_useFlattening)
    {
        Interval velInt    = m_gdnvPhysics->velocityInterval();
        int pressureIndex  = m_gdnvPhysics->pressureIndex();
        Real smallPressure = m_gdnvPhysics->smallPressure();
        int bulkIndex      = m_gdnvPhysics->bulkModulusIndex();

        m_util.computeFlattening(flattening,
            W,
            velInt,
            pressureIndex,
            smallPressure,
            bulkIndex,
            slopeBox);
    }

    // Intermediate, extrapolated primitive variables
    FArrayBox WMinus[SpaceDim];
    FArrayBox WPlus [SpaceDim];

    // Initial fluxes
    FArrayBox WHalf1[SpaceDim];

    // Source term computed from current state
    FArrayBox localSource;

    // If the source term is valid, make a local copy, increment it, and scale
    // it by 1/2 the timestep
    if (!a_S.box().isEmpty())
    {
        localSource.define(a_S.box(), a_S.nComp());
        localSource.copy(a_S);

        m_gdnvPhysics->incrementSource(localSource,W,slopeBox);

        localSource *= 0.5 * a_dt;
    }

    // Compute initial fluxes
    for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
        // Size the intermediate, extrapolated primitive variables
        WMinus[dir1].resize(slopeBox,numPrim); // cell-centered
        WPlus [dir1].resize(slopeBox,numPrim); // cell-centered

        // Compute predictor step to obtain extrapolated primitive variables
        if (m_normalPredOrder == 0)
        {
            CTUNormalPred(WMinus[dir1],WPlus[dir1],a_dt,m_dx,W,
                flattening,dir1,slopeBox);
        }
        else if (m_normalPredOrder == 1)
        {
            PLMNormalPred(WMinus[dir1],WPlus[dir1],a_dt,m_dx,W,
                flattening,dir1,slopeBox);
        }
        else if (m_normalPredOrder == 2)
        {
            PPMNormalPred(WMinus[dir1],WPlus[dir1],a_dt,m_dx,W,
                flattening,dir1,slopeBox);
        }
        else
        {
            MayDay::Error("PatchGodunov::computeWHalfBND: Normal predictor order must be 1 (PLM) or 2 (PPM)");
        }

        // If the source term is valid add it to the primitive quantities
        if (!localSource.box().isEmpty())
        {
            WMinus[dir1] += localSource;
            WPlus [dir1] += localSource;
        }

        // Solve the Riemann problem
        WHalf1[dir1].resize(fluxBox[dir1],numPrim); // face-centered

        // If the box contains the boundary with the fault call riemannBND
        if(a_bdryBox.sameType(faceBox[dir1]))
        {
            // we know from the define that it is really of this type
            ((LinElastPhysics*)m_gdnvPhysics)->riemannBND(WHalf1[dir1],WPlus[dir1],
                WMinus[dir1],W,a_Psi,m_currentTime,dir1,faceBox[dir1],a_bdryBox);
        }
        else
        {
            // we know from the define that it is really of this type
            m_gdnvPhysics->riemann(WHalf1[dir1],WPlus[dir1],WMinus[dir1],
                W,m_currentTime,dir1,faceBox[dir1]);
        }

    }

#if (CH_SPACEDIM == 3)
    // In 3D, compute some additional intermediate fluxes
    //
    // NOTE:  The diagonal entries of this array of fluxes are not
    // used and will not be defined.
    FArrayBox WHalf2[SpaceDim][SpaceDim];

    // Compute the intermediate, corrected fluxes in each direction
    for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
        // Correct fluxes using fluxes from a different direction
        for (int dir2 = 0; dir2 < SpaceDim; dir2++)
        {
            // A different direction has been found
            if (dir2 != dir1)
            {
                // Temporary primitive variables
                FArrayBox WTempMinus(WMinus[dir1].box(),numPrim);
                FArrayBox WTempPlus (WPlus [dir1].box(),numPrim);
                FArrayBox AdWdx(WPlus[dir1].box(),numPrim);

                // preventing uninitialized memory reads which cause
                // FPE's on some machines
                // AdWdx.setVal(666.666);

                // Copy data for in place modification
                WTempMinus.copy(WMinus[dir1]);
                WTempPlus .copy(WPlus [dir1]);

                // Update the current, extrapolated primitive variable using a flux
                // in a different direction

                m_gdnvPhysics->quasilinearUpdate(AdWdx,
                    WHalf1[dir2],WTempMinus,-(1.0/3.0) * a_dt / m_dx,dir2,ccBox[dir2]);
                WTempMinus += AdWdx;

                m_gdnvPhysics->quasilinearUpdate(AdWdx,
                    WHalf1[dir2],WTempPlus ,-(1.0/3.0) * a_dt / m_dx,dir2,ccBox[dir2]);
                WTempPlus  += AdWdx;

                // Solve the Riemann problem.

                WHalf2[dir1][dir2].resize(fluxBox[dir1],numPrim);
                // we know from the define that it is really of this type

                // If the box contains the boundary with the fault call riemannBND
                if(a_bdryBox.sameType(faceBox[dir1]))
                {
                    // we know from the define that it is really of this type
                    ((LinElastPhysics*)m_gdnvPhysics)->riemannBND(WHalf2[dir1][dir2],WTempPlus,WTempMinus,W,
                        a_Psi,m_currentTime,dir1,faceBox[dir1],a_bdryBox);
                }
                else
                {
                    // we know from the define that it is really of this type
                    ((LinElastPhysics*)m_gdnvPhysics)->riemann(WHalf2[dir1][dir2],WTempPlus,WTempMinus,W,
                        m_currentTime,dir1,faceBox[dir1]);
                }
            }
        }
    }
#endif

    // faceBox and fluxBox are now a bit smaller for the final corrections
    for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
    {
        faceBox[dir1] = curBox;
        faceBox[dir1].grow(dir1,1);
        faceBox[dir1] &= m_domain;
        faceBox[dir1].grow(dir1,-1);
        faceBox[dir1].surroundingNodes(dir1);

        fluxBox[dir1] = curBox;
        fluxBox[dir1].surroundingNodes(dir1);
    }

    // Do the final corrections to the fluxes
    for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
        // Correct the flux using fluxes in the remaining direction(s)
        for (int dir2 = 0; dir2 < SpaceDim; dir2++)
        {
            // A different direction has been found
            if (dir2 != dir1)
            {
#if (CH_SPACEDIM == 2)
                // In 2D, the current primitive state is updated by a flux in
                // the other direction
                FArrayBox AdWdx(WPlus[dir1].box(),numPrim);

                m_gdnvPhysics->quasilinearUpdate(AdWdx,WHalf1[dir2],WMinus[dir1],
                    -(1.0/2.0) * a_dt / m_dx,dir2,ccBox[dir2]);
                WMinus[dir1] += AdWdx;

                m_gdnvPhysics->quasilinearUpdate(AdWdx,WHalf1[dir2],WPlus [dir1],
                    -(1.0/2.0) * a_dt / m_dx,dir2,ccBox[dir2]);
                WPlus[dir1] += AdWdx;

#elif (CH_SPACEDIM == 3)
                // In 3D, find a direction different from the two above
                int dir3 = 3 - dir1 - dir2;

                // Update the conservative state using both corrected fluxes in
                // the other two directions
                FArrayBox AdWdx(WPlus[dir1].box(),numPrim);

                m_gdnvPhysics->quasilinearUpdate(AdWdx,WHalf2[dir2][dir3],WMinus[dir1],
                    -(1.0/2.0) * a_dt / m_dx,dir2,ccBox[dir2]);
                WMinus[dir1].plus(AdWdx,0,0,numPrim);

                m_gdnvPhysics->quasilinearUpdate(AdWdx,WHalf2[dir2][dir3],WPlus [dir1],
                    -(1.0/2.0) * a_dt / m_dx,dir2,ccBox[dir2]);
                WPlus [dir1].plus(AdWdx,0,0,numPrim);
#else
                // Only 2D and 3D should be possible
                MayDay::Error("PatchGodunov::computeWHalfBND: CH_SPACEDIM not 2 or 3!");
#endif
            }
        }

        // Solve the Riemann problem to obtain time-centered face values.
        // be returned
        FArrayBox& WHalf = a_WHalf[dir1];
        // we know from the define that it is really of this type

        // If the box contains the boundary with the fault call riemannBND
        if(a_bdryBox.sameType(faceBox[dir1]))
        {
            // we know from the define that it is really of this type
            ((LinElastPhysics*)m_gdnvPhysics)->riemannBND(WHalf,WPlus[dir1],WMinus[dir1],
                W,a_Psi,m_currentTime,dir1,faceBox[dir1],a_bdryBox);
        }
        else
        {
            // we know from the define that it is really of this type
            ((LinElastPhysics*)m_gdnvPhysics)->riemann(WHalf,WPlus[dir1],WMinus[dir1],W,
                 m_currentTime,dir1,faceBox[dir1]);
        }
    }
}

#include "NamespaceFooter.H"
