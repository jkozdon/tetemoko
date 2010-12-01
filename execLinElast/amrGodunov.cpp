#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::ifstream;
using std::ios;

// #define HALEM_PROC_SPEED
#ifdef HALEM_PROC_SPEED
#include <cstdio>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#endif

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "FABView.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "memusage.H"

#include "AMR.H"
#include "AMRLevel.H"
#include "AMRLevelLinElastFactory.H"
#include "AMRLevelLinElast.H"

#include "LinElastPhysics.H"

#include "SimpleIBC.H"
#include "LockSlideIBC.H"
#include "VelSlideAsinh1sIBC.H"
#include "SWIBC.H"
#include "PseudoPulseIBC.H"

#include "UsingNamespace.H"

#ifdef CH_Linux
// Should be undefined by default
#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

OldTimer Everything    ("gov Everything", 0);
OldTimer TimeReadInput ("gov Read Input",   Everything);
OldTimer TimeSetupAMR  ("gov Setup AMR",    Everything);
OldTimer TimeRun       ("gov Run",          Everything);
OldTimer TimeConclude  ("gov Conclude",     Everything);

// Possible pressure relationships for the initial condition
#define PRESSURE_ISENTROPIC 0
#define PRESSURE_CONSTANT   1

// amrGodunov is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void amrGodunov();

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
    const ProblemDomain&  a_domain,
    int                   a_maxLevel,
    int                   a_maxGridSize,
    int                   a_blockFactor,
    int                   a_verbosity,
    std::string           a_gridFile);

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
    // Start MPI
    MPI_Init(&a_argc,&a_argv);
#ifdef CH_AIX
    H5dont_atexit();
#endif
    // setChomboMPIErrorHandler();
#endif

    int rank, number_procs;

#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank = 0;
    number_procs = 1;
#endif

    if (rank == 0)
    {
        pout() << " number_procs = " << number_procs << endl;
    }

    OldTimer::TimerInit(rank);

    Everything.start();

    // Check for an input file
    char* inFile = NULL;

    if (a_argc > 1)
    {
        inFile = a_argv[1];
    }
    else
    {
        pout() << "Usage:  amrGodunov...ex <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
    }

    // Parse the command line and the input file (if any)
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

#ifdef TRAP_FPE
    enableFpExceptions();
#endif

    // Run amrGodunov, i.e., do the computation
    amrGodunov();

    Everything.stop();

#ifndef CH_NTIMER
    Real end_memory = get_memory_usage_from_OS();

    pout() << endl
        << "Everything completed --- "
        << "mem: "
        << setw(8) << setprecision(3)
        << setiosflags(ios::fixed)
        << end_memory
        << " MB, time: "
        << setw(8) << setprecision(3)
        << setiosflags(ios::fixed)
        << Everything.wc_time()
        << " sec (wall-clock)" << endl << endl;
#endif

#if !defined(CH_NTIMER) && defined(CH_MPI)
    Real avg_memory, min_memory, max_memory;
    gather_memory_from_procs(end_memory, avg_memory, min_memory, max_memory);
#endif

    // OldTimer::TimerSummary();
    CH_TIMER_REPORT();

#ifdef CH_MPI
    // Exit MPI
    dumpmemoryatexit();
    MPI_Finalize();
#endif
}

void amrGodunov()
{
    // Start timing the reading of the input file
    TimeReadInput.start();

    // Read inputs that are prefixed with "physics." These are the physical
    // parameters of the simulation
    ParmParse ppphysics("physics");

    // Read inputs that are prefixed with "comp." These are the computational
    // parameters of the simulation
    ParmParse ppcomp("comp");

    // Set the verbosity, i.e. how much diagnostic output
    int verbosity = 0;
    ppcomp.query("verbosity",verbosity);
    CH_assert(verbosity >= 0);

    // Get the physical parameters of the simulation
    // S-Wave speed
    Real cs = 3.0;
    ppphysics.query("cs",cs);
    CH_assert(cs >= 0);

    // P-Wave speed
    Real cp = cs*sqrt(3);
    ppphysics.query("cp",cp);
    CH_assert(cp >= 0);

    // Shear modulus
    Real mu = 30;
    ppphysics.query("mu",mu);
    CH_assert(mu >= 0);

    if(verbosity >= 1)
    {
           pout() << endl;
           pout() << ":::::::::::::::::::::::::::" << endl;
           pout() << ":: Simulation Parameters ::" << endl;
           pout() << ":::::::::::::::::::::::::::" << endl << endl;
           pout() << "S-wave Speed    = " << cs << endl;
           pout() << "P-wave Speed    = " << cp << endl;
           pout() << "Shear Modulus   = " << mu << endl << endl;
    }

    // Define IBC for ramp problem
    std::string problemString;
    bool isPeriodic[SpaceDim];

    // Set the physical size of the longest dimension of the domain
    Real domainLength = 1.0;
    ppcomp.get("domain_length",domainLength);

    // Stop when the simulation time get here
    Real maxTime = 0.0;
    ppcomp.get("max_time",maxTime);

    // Stop after this number of steps
    int nstop = 0;
    ppcomp.get("max_step",nstop);

    // Set the resolution of the coarsest level
    vector<int> numCells(SpaceDim);
    for (int i = 0; i < SpaceDim; ++i)
    {
        numCells[i] = 0;
    }
    ppcomp.getarr("num_cells",numCells,0,SpaceDim);

    CH_assert(D_TERM(   (numCells[0] > 0),
            && (numCells[1] > 0),
            && (numCells[2] > 0)));
    CH_assert(D_TERM(   (numCells[0] % 2 == 0),
            && (numCells[1] % 2 == 0),
            && (numCells[2] % 2 == 0)));

    // Maximum AMR level limit
    int maxLevel = 0;
    ppcomp.get("max_level",maxLevel);
    int numReadLevels = Max(maxLevel,1);

    // Refinement ratios between levels
    std::vector<int> refRatios;
    // Note: this requires a refRatio to be defined for the finest level
    // (even though it will never be used)
    ppcomp.getarr("ref_ratio",refRatios,0,numReadLevels+1);

    // Number of coarse time steps from one regridding to the next
    std::vector<int> regridIntervals;
    ppcomp.getarr("regrid_interval",regridIntervals,0,numReadLevels);

    // How far to extend refinement from cells newly tagged for refinement
    int tagBufferSize = 3;
    ppcomp.get("tag_buffer_size",tagBufferSize);

    // Threshold that triggers refinement
    Real refineThresh = 0.3;
    ppcomp.get ("refine_thresh",refineThresh);

    // Minimum dimension of a grid
    int blockFactor = 1;
    ppcomp.get("block_factor",blockFactor);

    // Maximum dimension of a grid
    int maxGridSize = 32;
    ppcomp.get("max_grid_size",maxGridSize);

    Real fillRatio = 0.75;
    ppcomp.get("fill_ratio",fillRatio);

    // Order of the normal predictor (CTU -> 0, PLM -> 1, PPM -> 2)
    std::string normalPred;
    int normalPredOrder;
    ppcomp.get("normal_predictor",normalPred);
    if (normalPred == "CTU" || normalPred == "ctu")
    {
        normalPredOrder = 0;
    }
    else if (normalPred == "PLM" || normalPred == "plm")
    {
        normalPredOrder = 1;
    }
    else if (normalPred == "PPM" || normalPred == "ppm")
    {
        normalPredOrder = 2;
    }
    else
    {
        MayDay::Error("Normal predictor must by PLM or PPM");
    }

    // Use fourth order slopes
    int inFourthOrderSlopes = 1;
    bool useFourthOrderSlopes;
    ppcomp.get("use_fourth_order_slopes",inFourthOrderSlopes);
    useFourthOrderSlopes = (inFourthOrderSlopes == 1);

    // Do slope limiting
    int inPrimLimiting = 1;
    bool usePrimLimiting;
    ppcomp.get("use_prim_limiting",inPrimLimiting);
    usePrimLimiting = (inPrimLimiting == 1);

    // Do slope limiting using characteristics
    int inCharLimiting = 0;
    bool useCharLimiting;
    ppcomp.get("use_char_limiting",inCharLimiting);
    useCharLimiting = (inCharLimiting == 1);

    // Do slope flattening
    int inFlattening = 0;
    bool useFlattening;
    //JK NOT CURRENTLY IMPLEMENTED
    // ppcomp.get("use_flattening",inFlattening);
    useFlattening = (inFlattening == 1);

    // Apply artificial viscosity
    int inArtificialViscosity = 0;
    bool useArtificialViscosity;
    //JK NOT CURRENTLY IMPLEMENTED
    // ppcomp.get("use_artificial_viscosity",inArtificialViscosity);
    useArtificialViscosity = (inArtificialViscosity == 1);

    // Artificial viscosity coefficient/multiplier
    Real artificialViscosity = 0.1;
    ppcomp.get("artificial_viscosity",artificialViscosity);

    // Don't use high-order limiter by default
    int inHighOrderLimiter = 0;
    bool highOrderLimiter;
    ppcomp.query("high_order_limiter", inHighOrderLimiter);
    highOrderLimiter = (inHighOrderLimiter == 1);

    // Set up checkpointing
    int checkpointInterval = 0;
    ppcomp.query("checkpoint_interval",checkpointInterval);

    // Set up plot file writing
    int plotInterval = 0;
    ppcomp.query("plot_interval",plotInterval);

    // CFL multiplier
    Real cfl = 0.8;
    ppcomp.get("cfl",cfl);

    // Initial CFL multiplier
    Real initialCFL = 0.1;
    ppcomp.get("initial_cfl",initialCFL);

    // Determine if a fixed or variable time step will be used
    Real fixedDt = -1;
    ppcomp.query("fixed_dt",fixedDt);

    // Limit the time step growth
    Real maxDtGrowth = 1.1;
    ppcomp.get("max_dt_growth",maxDtGrowth);

    // Let the time step grow by this factor above the "maximum" before
    // reducing it
    Real dtToleranceFactor = 1.1;
    ppcomp.get("dt_tolerance_factor",dtToleranceFactor);

    // End timing the reading of the input file
    TimeReadInput.stop();

    // Create and define IBC (initial and boundary condition) object
    PhysIBC* ibc;

    // Background Values
    vector<Real> backgroundVals(9,0);
    ppphysics.queryarr("background",backgroundVals,0,9);

    if (ppphysics.contains("problem"))
    {
        ppphysics.query("problem",problemString);
        if (problemString == "simple")
        {
            // Pulse Location
            Real r0 = 0.5;
            ppphysics.query("r0",r0);
            CH_assert(r0 >= 0);

            // Pulse Magnitude
            vector<Real> mag(9,0.0);
            ppphysics.getarr("mag",mag,0,9);

            // Pulse Width
            Real sig = 10;
            ppphysics.query("sig",sig);
            CH_assert(sig >= 0);

            // Determine which spatial directions are periodic
            vector<int> isPeriodica(SpaceDim,0);

            ppcomp.getarr("is_periodic",isPeriodica,0,SpaceDim);

            vector<int> xBoudaryT(6,0);
            ppcomp.queryarr("x_boundary",xBoudaryT,0,2);

            vector<int> yBoudaryT(2,0);
            ppcomp.queryarr("y_boundary",yBoudaryT,0,2);

            vector<int> zBoudaryT(2,0);
            ppcomp.queryarr("z_boundary",zBoudaryT,0,2);

            std::vector<int> boundaryType(6,0);
            boundaryType[0] = xBoudaryT[0];
            boundaryType[1] = xBoudaryT[1];
            boundaryType[2] = yBoudaryT[0];
            boundaryType[3] = yBoudaryT[1];
            boundaryType[4] = zBoudaryT[0];
            boundaryType[5] = zBoudaryT[1];

            // convert periodic from int->bool
            for (int dim = 0; dim < SpaceDim; dim++)
            {
                isPeriodic[dim] = (isPeriodica[dim] == 1);
                if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
                {
                    pout() << "Using Periodic BCs in direction: " << dim << endl;
                }
                else if(verbosity >= 2 && procID() == 0)
                {
                    pout() << "Using BCs " << boundaryType[dim*2] << " and " << boundaryType[dim*2+1] << " in direction: " << dim << endl;
                }
            }


            SimpleIBC* simpleibc = new SimpleIBC(cs,cp,mu,backgroundVals,r0,mag,sig,boundaryType);
            ibc = simpleibc;
            if(verbosity >= 1)
            {
                pout() << "Pulse Location  = " << r0 << endl;
                pout() << "Pulse Magnitude = " << 
                    mag[0] << " " <<
                    mag[1] << " " <<
                    mag[2] << " " <<
                    mag[3] << " " <<
                    mag[4] << " " <<
                    mag[5] << " " <<
                    mag[6] << " " <<
                    mag[7] << " " <<
                    mag[8] << " " <<
                    endl;
                pout() << "Pulse Width     = " << sig << endl << endl;
            }
        }
        else if (problemString == "lockSlide")
        {

            // where is the center of the slip region
            Real center = 0.5;
            ppphysics.query("center_scale",center);
            // CH_assert(center >= 0);
            if(center < 0)
            {
                center = -center;
            }
            else
            {
                center = center * domainLength;
            }

            // how big is the slip region
            Real edge = 0.25;
            ppphysics.query("edge_scale",edge);
            // CH_assert(edge >= 0);
            if(edge < 0)
            {
                edge = -edge;
            }
            else
            {
                edge = edge * domainLength;
            }



            // no dimensions are periodic with this problem
            for (int dim = 0; dim < SpaceDim; dim++)
            {
                isPeriodic[dim] = false;
            }


            LockSlideIBC* lockSlideibc = 
                new LockSlideIBC(cs,cp,mu,backgroundVals,center,edge);
            ibc = lockSlideibc;
            if(verbosity >= 1)
            {
                pout() << "Background Values = " << 
                    backgroundVals[0] << " " <<
                    backgroundVals[1] << " " <<
                    backgroundVals[2] << " " <<
                    backgroundVals[3] << " " <<
                    backgroundVals[4] << " " <<
                    backgroundVals[5] << " " <<
                    backgroundVals[6] << " " <<
                    backgroundVals[7] << " " <<
                    backgroundVals[8] << " " <<
                    endl << endl;
            }
        }
        else if (problemString == "velslideasinh1s")
        {

            // where is the center of the slip region
            Real center = 0.5;
            ppphysics.query("center_scale",center);
            // CH_assert(center >= 0);
            if(center < 0)
            {
                center = -center;
            }
            else
            {
                center = center * domainLength;
            }

            // how big is the slip region
            Real sigma = 10;
            ppphysics.query("sigma",sigma);
            CH_assert(sigma >= 0);

            // how big is the slip region
            Real ntime = 1;
            ppphysics.query("nucleation_time",ntime);
            CH_assert(ntime >= 0);

            // no dimensions are periodic with this problem
            for (int dim = 0; dim < SpaceDim; dim++)
            {
                isPeriodic[dim] = false;
            }


            VelSlideAsinh1sIBC* velSlideibc =
                new VelSlideAsinh1sIBC(cs,cp,mu,backgroundVals,center,sigma,ntime);
            ibc = velSlideibc;
            if(verbosity >= 1)
            {
                pout() << "Background Values = " << 
                    backgroundVals[0] << " " <<
                    backgroundVals[1] << " " <<
                    backgroundVals[2] << " " <<
                    backgroundVals[3] << " " <<
                    backgroundVals[4] << " " <<
                    backgroundVals[5] << " " <<
                    backgroundVals[6] << " " <<
                    backgroundVals[7] << " " <<
                    backgroundVals[8] << " " <<
                    endl << endl;
            }
        }
        else if (problemString == "slipweak")
        {
            // This is a very slimple IBC loosely based on the SCEC problems

            // Where is the nucleation patch
            int numPatches = 0;
            ppphysics.get("num_patches",numPatches);

            vector<Real> xcPatches(numPatches,0);
            vector<Real> xwPatches(numPatches,0);
            vector<Real> zcPatches(numPatches,0);
            vector<Real> zwPatches(numPatches,0);
            vector<Real> tauPatches(numPatches,0);

            ppphysics.getarr("xc_patches",xcPatches,0,numPatches);
            ppphysics.getarr("xw_patches",xwPatches,0,numPatches);
            ppphysics.getarr("zc_patches",zcPatches,0,numPatches);
            ppphysics.getarr("zw_patches",zwPatches,0,numPatches);
            ppphysics.getarr("tau_patches",tauPatches,0,numPatches);

            Real fricD = 0.525;
            ppphysics.query("f_dynamic",fricD);

            Real fricS = 0.677;
            ppphysics.query("f_static",fricS);

            Real weakDist = 0.4;
            ppphysics.query("d_weak",weakDist);


            // get the boundary conditiions
            // 0 : Outflow
            // 1 : Free surface
            // 2 : Fault (must be y_bound 2 X)
            vector<int> xBoudaryT(2,0);
            ppcomp.queryarr("x_boundary",xBoudaryT,0,2);

            vector<int> yBoudaryT(2,0);
            ppcomp.queryarr("y_boundary",yBoudaryT,0,2);

            vector<int> zBoudaryT(2,0);
            ppcomp.queryarr("z_boundary",zBoudaryT,0,2);

            std::vector<int> boundaryType(6,0);
            boundaryType[0] = xBoudaryT[0];
            boundaryType[1] = xBoudaryT[1];
            boundaryType[2] = yBoudaryT[0];
            boundaryType[3] = yBoudaryT[1];
            boundaryType[4] = zBoudaryT[0];
            boundaryType[5] = zBoudaryT[1];

            // no dimensions are periodic with this problem
            for (int dim = 0; dim < SpaceDim; dim++)
            {
                isPeriodic[dim] = false;
                if(verbosity >= 2 && procID() == 0)
                {
                    pout() << "Using BCs " << boundaryType[dim*2] << " and " << boundaryType[dim*2+1] << " in direction: " << dim << endl;
                }
            }

            Real dx = domainLength;
            for(int dim = 0; dim < SpaceDim; dim++)
            {
                dx = min(dx,domainLength/((Real) numCells[dim]));
            }
            // for(int lvl = 0; lvl < numReadLevels; lvl++)
            // {
            //     dx = dx / ((Real) refRatios[lvl]);
            // }

            Real width = 12*dx;

            SWIBC* swibc =
                new SWIBC(cs,cp,mu,backgroundVals,fricS,fricD,weakDist,width,
                    numPatches,xcPatches,xwPatches,zcPatches,zwPatches,tauPatches,
                    boundaryType);
            ibc = swibc;
            if(verbosity >= 1)
            {
                pout() << "Static Friction    = " << fricS << endl;
                pout() << "Dynamic Friction   = " << fricD << endl;
                pout() << "Weakening Distance = " << weakDist << endl;
                pout() << "num patches        = " << numPatches << endl;
                for(int itor = 0; itor < numPatches; itor++)
                {
                    pout() << "Patch " << itor << ": tau = " << tauPatches[itor];
                    pout() << " at " << "(" << xcPatches[itor] << ",0," << zcPatches[itor] << ")";
                    pout() << " +/- (" << xwPatches[itor] << ",0," << zwPatches[itor] << ")" << endl;
                }
                pout() << "Background Values = " << 
                    backgroundVals[0] << " " <<
                    backgroundVals[1] << " " <<
                    backgroundVals[2] << " " <<
                    backgroundVals[3] << " " <<
                    backgroundVals[4] << " " <<
                    backgroundVals[5] << " " <<
                    backgroundVals[6] << " " <<
                    backgroundVals[7] << " " <<
                    backgroundVals[8] << " " <<
                    endl << endl;
            }
        }
        else if (problemString == "pseudopulse")
        {

            // where is the center of the slip region
            Real center = 0.5;
            ppphysics.query("center_scale",center);
            // CH_assert(center >= 0);
            if(center < 0)
            {
                center = -center;
            }
            else
            {
                center = center * domainLength;
            }

            // how big is the slip region
            Real slope = 10;
            ppphysics.query("slope",slope);
            CH_assert(slope >= 0);

            Real dTau = 1;
            ppphysics.query("dTau",dTau);
            CH_assert(dTau >= 0);

            Real Vh = 1.2;
            ppphysics.query("Vh",Vh);
            CH_assert(Vh >= 0);

            Real Vr = 1.6;
            ppphysics.query("Vr",Vr);
            CH_assert(Vr >= 0);

            Real jumpTime = maxTime / 2.0;
            ppphysics.query("jump_time",jumpTime);
            CH_assert(jumpTime >= 0);

            Real stopTime = maxTime;
            ppphysics.query("stop_time",stopTime);
            CH_assert(stopTime >= 0);


            // no dimensions are periodic with this problem
            for (int dim = 0; dim < SpaceDim; dim++)
            {
                isPeriodic[dim] = false;
            }

            PseudoPulseIBC* pseudopulseibc =
                new PseudoPulseIBC(cs,cp,mu,backgroundVals,center,slope,dTau,Vh,Vr,stopTime,jumpTime);
            ibc = pseudopulseibc;
            if(verbosity >= 1)
            {
                pout() << "slope   = " << slope << endl;
                pout() << "dTau    = " << dTau << endl;
                pout() << "Vr / cs = " << Vr << endl;
                pout() << "Vh / cs = " << Vh << endl;
                pout() << "Background Values = " << 
                    backgroundVals[0] << " " <<
                    backgroundVals[1] << " " <<
                    backgroundVals[2] << " " <<
                    backgroundVals[3] << " " <<
                    backgroundVals[4] << " " <<
                    backgroundVals[5] << " " <<
                    backgroundVals[6] << " " <<
                    backgroundVals[7] << " " <<
                    backgroundVals[8] << " " <<
                    endl << endl;
            }
        }
        else
        {
            // The sample problem name given isn't valid
            pout() << "Invalid problem, \"" << problemString << "\", specified in input file" << endl << endl;
            return;
        }
    }
    else
    {
        // A sample problem must be specified
        pout() << "\"godunov.problem\" not specified in input file" << endl << endl;
        return;
    }


#ifndef CH_NTIMER
    pout() << "Input Read completed --- "
        << "mem: "
        << setw(8) << setprecision(3)
        << setiosflags(ios::fixed)
        << get_memory_usage_from_OS()
        << " MB, time: "
        << setw(8) << setprecision(3)
        << setiosflags(ios::fixed)
        << TimeReadInput.wc_time()
        << " sec (wall-clock)" << endl;
#endif
    
    // Start timing AMR solver setup
    TimeSetupAMR.start();

    //JK // A minimum pressure needed to construct PolytropicPhysics - used in slope
    //JK // flattening
    //JK Real smallPressure;

    //JK Don't think that we need this
    // Don't use source term by default
    bool useSourceTerm = false;
    
    // Source term multiplier
    Real sourceTermScaling = 0.0;


    if (verbosity >= 1)
    {
        pout() << endl;
        pout() << "::::::::::::::::::::::::::::::" << endl;
        pout() << ":: Computational Parameters ::" << endl;
        pout() << "::::::::::::::::::::::::::::::" << endl << endl;
        pout() << "verbosity            = " << verbosity << endl;

        pout() << "maximum_step         = " << nstop << endl;
        pout() << "maximum_time         = " << maxTime << endl;
        if (fixedDt > 0)
        {
            pout() << "fixed_dt             = " << fixedDt << endl;
        }

        pout() << "number_of_cells      = " << D_TERM(numCells[0] << "  " <<,
            numCells[1] << "  " <<,
            numCells[2] << ) endl;
        pout() << "is_period            = " << D_TERM(isPeriodic[0] << "  " <<,
            isPeriodic[1] << "  " <<,
            isPeriodic[2] << ) endl;

        pout() << "maximum_level        = " << maxLevel << endl;
        pout() << "refinement_ratio     = ";
        for (int i = 0; i < refRatios.size(); ++i)
        {
            pout() << refRatios[i] << " ";
        }
        pout() << endl;

        pout() << "regrid_interval      = ";
        for (int i = 0; i < regridIntervals.size(); ++i)
        {
            pout() << regridIntervals[i] << " ";
        }
        pout() << endl;
        pout() << "tag_buffer_size      = " << tagBufferSize << endl;

        pout() << "refinement_threshold = " << refineThresh << endl;

        pout() << "blocking_factor      = " << blockFactor << endl;
        pout() << "max_grid_size        = " << maxGridSize << endl;
        pout() << "fill_ratio           = " << fillRatio << endl;

        pout() << "normal_predictor = ";
        if (normalPredOrder == 0)
        {
            pout() << "CTU" << endl;
        }
        else if (normalPredOrder == 1)
        {
            pout() << "PLM" << endl;
        }
        else if (normalPredOrder == 2)
        {
            pout() << "PPM" << endl;
        }
        else
        {
            pout() << "Unknown (" << normalPredOrder << ")" << endl;
        }

        pout() << "slope_order     = "
            << (useFourthOrderSlopes ? "2nd" : "4th") << endl;
        pout() << "use_primitive_slope_limiting      = "
            << (usePrimLimiting ? "yes" : "no") << endl;
        pout() << "use_characteristic_slope_limiting = "
            << (useCharLimiting ? "yes" : "no") << endl;
        pout() << "use_slope_flattening              = "
            << (useFlattening ? "yes" : "no") << endl;

        //JK pout() << "use_artificial_viscosity = "
        //JK     << (useArtificialViscosity ? "yes" : "no") << endl;
        //JK if (useArtificialViscosity)
        //JK {
        //JK     pout() << "artificial_viscosity     = " << artificialViscosity << endl;
        //JK }

        //JK pout() << "use_source_term          = "
        //JK     << (useSourceTerm ? "yes" : "no") << endl;
        //JK if (useSourceTerm)
        //JK {
        //JK     pout() << "source_term_scaling  = " << sourceTermScaling << endl;
        //JK }

        pout() << "checkpoint_interval  = " << checkpointInterval << endl;
        pout() << "plot_interval        = " << plotInterval << endl;

        pout() << "CFL                  = " << cfl << endl;
        pout() << "initial_CFL          = " << initialCFL << endl;

        pout() << "maximum_dt_growth    = " << maxDtGrowth << endl;
        pout() << "dt_tolerance_factor  = " << dtToleranceFactor << endl;
    }

    ProblemDomain probDomain (IntVect::Zero,
        IntVect(D_DECL(numCells[0]-1,
                numCells[1]-1,
                numCells[2]-1)),
        isPeriodic);
    
    // Set up the physics for linear elasticitiy
    LinElastPhysics linElastPhysics(0 /*JUNK*/);
    linElastPhysics.setPhysIBC(ibc);

    // Cast to physics base class pointer for technical reasons
    GodunovPhysics* godunovPhysics = static_cast<GodunovPhysics*> (&linElastPhysics);

    // Set up the AMRLevel... factory
    AMRLevelLinElastFactory amrGodFact;

    amrGodFact.define(cfl,
        domainLength,
        verbosity,
        refineThresh,
        tagBufferSize,
        initialCFL,
        godunovPhysics,
        normalPredOrder,
        useFourthOrderSlopes,
        usePrimLimiting,
        useCharLimiting,
        useFlattening,
        useArtificialViscosity,
        artificialViscosity,
        useSourceTerm,
        sourceTermScaling,
        highOrderLimiter);

    AMR amr;
    
    // Set up the AMR object
    amr.define(maxLevel,refRatios,probDomain,&amrGodFact);
    
    if (fixedDt > 0)
    {
      amr.fixedDt(fixedDt);
    }
    
    // Set grid generation parameters
    amr.maxGridSize(maxGridSize);
    amr.blockFactor(blockFactor);
    amr.fillRatio(fillRatio);
    
    // The hyperbolic codes use a grid buffer of 1
    amr.gridBufferSize(1);
    
    // Set output parameters
    amr.checkpointInterval(checkpointInterval);
    amr.plotInterval(plotInterval);
    amr.regridIntervals(regridIntervals);
    amr.maxDtGrow(maxDtGrowth);
    amr.dtToleranceFactor(dtToleranceFactor);

    // Set up output files
    if (ppcomp.contains("plot_prefix"))
    {
        std::string prefix;
        ppcomp.query("plot_prefix",prefix);
        amr.plotPrefix(prefix);
    }

    if (ppcomp.contains("chk_prefix"))
    {
        std::string prefix;
        ppcomp.query("chk_prefix",prefix);
        amr.checkpointPrefix(prefix);
    }

    amr.verbosity(verbosity);

    // Set up input files
    if (!ppcomp.contains("restart_file"))
    {
        if (!ppcomp.contains("fixed_hierarchy"))
        {
            // initialize from scratch for AMR run
            // initialize hierarchy of levels
            amr.setupForNewAMRRun();
        }
        else
        {
            std::string gridFile;
            ppcomp.query("fixed_hierarchy",gridFile);

            // initialize from a list of grids in "gridFile"
            Vector<Vector<Box> > amrGrids(maxLevel+1);
            setupFixedGrids(amrGrids,
                probDomain,
                maxLevel,
                maxGridSize,
                blockFactor,
                verbosity,
                gridFile);
            amr.setupForFixedHierarchyRun(amrGrids,1);
        }
    }
    else
    {
        std::string restartFile;
        ppcomp.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
        HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
        // read from checkpoint file
        amr.setupForRestart(handle);
        handle.close();
#else
        MayDay::Error("amrGodunov restart only defined with hdf5");
#endif
    }

    // End timing AMR solver setup
    TimeSetupAMR.stop();

#ifndef CH_NTIMER
    pout() << "AMR Setup completed ---- "
        << "mem: "
        << setw(8) << setprecision(3)
        << setiosflags(ios::fixed)
        << get_memory_usage_from_OS()
        << " MB, time: "
        << setw(8) << setprecision(3)
        << setiosflags(ios::fixed)
        << TimeSetupAMR.wc_time()
        << " sec (wall-clock)" << endl;

    if (verbosity >= 1)
    {
        pout() << endl;
    }
#endif

    // Run and time the computation
    TimeRun.start();
    amr.run(maxTime,nstop);
    Vector<AMRLevel*> vecLevels = amr.getAMRLevels();
    for(int ator = 0; ator < vecLevels.size(); ator++)
    {
        if(((AMRLevelLinElast*)vecLevels[ator])->level() == 0)
        {
        }
    }
    TimeRun.stop();

#ifndef CN_NTIMER
    if (verbosity >= 1)
    {
        pout() << endl;
    }

    pout() << "AMR Run completed ------ "
        << "mem: "
        << setw(8) << setprecision(3)
        << setiosflags(ios::fixed)
        << get_memory_usage_from_OS()
        << " MB, time: "
        << setw(8) << setprecision(3)
        << setiosflags(ios::fixed)
        << TimeRun.wc_time()
        << " sec (wall-clock)" << endl;
#endif

    // Output the last plot file and statistics - time the process
    TimeConclude.start();
    amr.conclude();
    TimeConclude.stop();

#ifndef CH_NTIMER
    pout() << "AMR Conclude completed - "
        << "mem: "
        << setw(8) << setprecision(3)
        << setiosflags(ios::fixed)
        << get_memory_usage_from_OS()
        << " MB, time: "
        << setw(8) << setprecision(3)
        << setiosflags(ios::fixed)
        << TimeConclude.wc_time()
        << " sec (wall-clock)" << endl;
#endif
}

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
    const ProblemDomain&  a_domain,
    int                   a_maxLevel,
    int                   a_maxGridSize,
    int                   a_blockFactor,
    int                   a_verbosity,
    std::string           a_gridFile)
{
    MayDay::Error("setupFixedGrids not defined");
    //  // Run this task on one processor
    //  if (procID() == uniqueProc(SerialTask::compute))
    //  {
    //    a_amrGrids.push_back(Vector<Box>(1,a_domain.domainBox()));
    //
    //    // Read in predefined grids
    //    ifstream is(a_gridFile.c_str(), ios::in);
    //
    //    if (is.fail())
    //    {
    //      MayDay::Error("Cannot open grids file");
    //    }
    //
    //    // Format of file:
    //    //   number of levels, then for each level (starting with level 1):
    //    //   number of grids on level, list of boxes
    //
    //    int inNumLevels;
    //    is >> inNumLevels;
    //
    //    CH_assert (inNumLevels <= a_maxLevel+1);
    //
    //    if (a_verbosity >= 3)
    //    {
    //      pout() << "numLevels = " << inNumLevels << endl;
    //    }
    //
    //    while (is.get() != '\n');
    //
    //    a_amrGrids.resize(inNumLevels);
    //
    //    // Check to see if coarsest level needs to be broken up
    //    domainSplit(a_domain,a_amrGrids[0],a_maxGridSize,a_blockFactor);
    //
    //    if (a_verbosity >= 3)
    //    {
    //      pout() << "level 0: ";
    //      for (int n = 0; n < a_amrGrids[0].size(); n++)
    //      {
    //        pout() << a_amrGrids[0][0] << endl;
    //      }
    //    }
    //
    //    // Now loop over levels, starting with level 1
    //    int ngrid;
    //    for (int lev = 1; lev < inNumLevels; lev++)
    //    {
    //      is >> ngrid;
    //
    //      if (a_verbosity >= 3)
    //      {
    //        pout() << "level " << lev << " numGrids = " << ngrid << endl;
    //        pout() << "Grids: ";
    //      }
    //
    //      while (is.get() != '\n');
    //
    //      a_amrGrids[lev].resize(ngrid);
    //
    //      for (int i = 0; i < ngrid; i++)
    //      {
    //        Box bx;
    //        is >> bx;
    //
    //        while (is.get() != '\n');
    //
    //        // Quick check on box size
    //        Box bxRef(bx);
    //
    //        if (bxRef.longside() > a_maxGridSize)
    //        {
    //          pout() << "Grid " << bx << " too large" << endl;
    //          MayDay::Error();
    //        }
    //
    //        if (a_verbosity >= 3)
    //        {
    //          pout() << bx << endl;
    //        }
    //
    //        a_amrGrids[lev][i] = bx;
    //      } // End loop over boxes on this level
    //    } // End loop over levels
    //  }
    //
    //  // Broadcast results to all the processors
    //  broadcast(a_amrGrids,uniqueProc(SerialTask::compute));
}

#ifdef TRAP_FPE
#include <fenv.h>

// FE_INEXACT    inexact result
// FE_DIVBYZERO  division by zero
// FE_UNDERFLOW  result not representable due to underflow
// FE_OVERFLOW   result not representable due to overflow
// FE_INVALID    invalid operation

static void enableFpExceptions ()
{
    if (feclearexcept(FE_ALL_EXCEPT) != 0)
    {
        MayDay::Abort("feclearexcept failed");
    }

    int flags = FE_DIVBYZERO |
        FE_INVALID   |
        //              FE_UNDERFLOW |
        FE_OVERFLOW  ;

    if (feenableexcept(flags) == -1)
    {
        MayDay::Abort("feenableexcept failed");
    }
}
#endif
