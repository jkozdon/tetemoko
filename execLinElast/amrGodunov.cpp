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
// #include "AMRLevelPolytropicGasFactory.H"
// #include "AMRLevelPolytropicGas.H"

#include "LinElastPhysics.H"

#include "SimpleIBC.H"

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
    Real cp = 5.2;
    ppphysics.query("cp",cp);
    CH_assert(cp >= 0);

    // Shear modulus
    Real mu = 30;
    ppphysics.query("mu",mu);
    CH_assert(mu >= 0);

    // Output the physical parameters if doing output
    if(verbosity >= 1)
    {
           pout() << endl;
           pout() << ":::::::::::::::::::::::::::" << endl;
           pout() << ":: Simulation Parameters ::" << endl;
           pout() << ":::::::::::::::::::::::::::" << endl << endl;
           pout() << "S-wave Speed  = " << cs << endl;
           pout() << "P-wave Speed  = " << cp << endl;
           pout() << "Shear Modulus = " << mu << endl << endl;
    }


    // Stop after this number of steps
    int nstop = 0;
    ppcomp.get("max_step",nstop);

    // Stop when the simulation time get here
    Real stopTime = 0.0;
    ppcomp.get("max_time",stopTime);

    // Set the physical size of the longest dimension of the domain
    Real domainLength = 1.0;
    ppcomp.get("domain_length",domainLength);

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

    // Determine which spatial directions are periodic
    vector<int> isPeriodica(SpaceDim,0);
    bool isPeriodic[SpaceDim];

    ppcomp.getarr("is_periodic",isPeriodica,0,SpaceDim);
    // convert periodic from int->bool
    for (int dim = 0; dim < SpaceDim; dim++)
    {
        isPeriodic[dim] = (isPeriodica[dim] == 1);
        if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
        {
            pout() << "Using Periodic BCs in direction: " << dim << endl;
        }
    }

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

    //JK // Order of the normal predictor (CTU -> 0, PLM -> 1, PPM -> 2)
    //JK std::string normalPred;
    //JK int normalPredOrder;
    //JK ppcomp.get("normal_predictor",normalPred);
    //JK if (normalPred == "CTU" || normalPred == "ctu")
    //JK {
    //JK     normalPredOrder = 0;
    //JK }
    //JK else if (normalPred == "PLM" || normalPred == "plm")
    //JK {
    //JK     normalPredOrder = 1;
    //JK }
    //JK else if (normalPred == "PPM" || normalPred == "ppm")
    //JK {
    //JK     normalPredOrder = 2;
    //JK }
    //JK else
    //JK {
    //JK     MayDay::Error("Normal predictor must by PLM or PPM");
    //JK }

    //JK // Use fourth order slopes
    //JK int inFourthOrderSlopes = 1;
    //JK bool useFourthOrderSlopes;
    //JK ppcomp.get("use_fourth_order_slopes",inFourthOrderSlopes);
    //JK useFourthOrderSlopes = (inFourthOrderSlopes == 1);

    //JK // Do slope limiting
    //JK int inPrimLimiting = 1;
    //JK bool usePrimLimiting;
    //JK ppcomp.get("use_prim_limiting",inPrimLimiting);
    //JK usePrimLimiting = (inPrimLimiting == 1);

    //JK // Do slope limiting using characteristics
    //JK int inCharLimiting = 0;
    //JK bool useCharLimiting;
    //JK ppcomp.get("use_char_limiting",inCharLimiting);
    //JK useCharLimiting = (inCharLimiting == 1);

    //JK // Do slope flattening
    //JK int inFlattening = 1;
    //JK bool useFlattening;
    //JK ppcomp.get("use_flattening",inFlattening);
    //JK useFlattening = (inFlattening == 1);

    //JK // Apply artificial viscosity
    //JK int inArtificialViscosity = 1;
    //JK bool useArtificialViscosity;
    //JK ppcomp.get("use_artificial_viscosity",inArtificialViscosity);
    //JK useArtificialViscosity = (inArtificialViscosity == 1);

    //JK // Artificial viscosity coefficient/multiplier
    //JK Real artificialViscosity = 0.1;
    //JK ppcomp.get("artificial_viscosity",artificialViscosity);

    //JK // Don't use high-order limiter by default
    //JK int inHighOrderLimiter = 0;
    //JK bool highOrderLimiter;
    //JK ppcomp.query("high_order_limiter", inHighOrderLimiter);
    //JK highOrderLimiter = (inHighOrderLimiter == 1);

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
    
    // Create and define IBC (initial and boundary condition) object
    PhysIBC* ibc;

    //JK // A minimum pressure needed to construct PolytropicPhysics - used in slope
    //JK // flattening
    //JK Real smallPressure;
    //JK 
    //JK // Don't use source term by default
    //JK bool useSourceTerm = false;
    //JK 
    //JK // Source term multiplier
    //JK Real sourceTermScaling = 0.0;

    // Define IBC for ramp problem
    SimpleIBC* simpleibc = new SimpleIBC(cs,cp,mu);
    //simpleibc->setFortranCommon(smallPressure,
    //                          gamma,
    //                          alpha,
    //                          ms,
    //                          xcorner,
    //                          artificialViscosity);
    ibc = simpleibc;


    if (verbosity >= 1)
    {
        pout() << endl;
        pout() << "::::::::::::::::::::::::::::::" << endl;
        pout() << ":: Computational Parameters ::" << endl;
        pout() << "::::::::::::::::::::::::::::::" << endl << endl;
        pout() << "verbosity            = " << verbosity << endl;

        pout() << "maximum_step         = " << nstop << endl;
        pout() << "maximum_time         = " << stopTime << endl;
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

        //JK pout() << "normal_predictor = ";
        //JK if (normalPredOrder == 0)
        //JK {
        //JK     pout() << "CTU" << endl;
        //JK }
        //JK else if (normalPredOrder == 1)
        //JK {
        //JK     pout() << "PLM" << endl;
        //JK }
        //JK else if (normalPredOrder == 2)
        //JK {
        //JK     pout() << "PPM" << endl;
        //JK }
        //JK else
        //JK {
        //JK     pout() << "Unknown (" << normalPredOrder << ")" << endl;
        //JK }

        //JK pout() << "slope_order = "
        //JK     << (useFourthOrderSlopes ? "2nd" : "4th") << endl;
        //JK pout() << "use_primitive_slope_limiting = "
        //JK     << (usePrimLimiting ? "yes" : "no") << endl;
        //JK pout() << "use_characteristic_slope_limiting = "
        //JK     << (useCharLimiting ? "yes" : "no") << endl;
        //JK pout() << "use_slope_flattening = "
        //JK     << (useFlattening ? "yes" : "no") << endl;

        //JK pout() << "use_artificial_viscosity = "
        //JK     << (useArtificialViscosity ? "yes" : "no") << endl;
        //JK if (useArtificialViscosity)
        //JK {
        //JK     pout() << "artificial_viscosity = " << artificialViscosity << endl;
        //JK }

        //JK pout() << "use_source_term = "
        //JK     << (useSourceTerm ? "yes" : "no") << endl;
        //JK if (useSourceTerm)
        //JK {
        //JK     pout() << "source_term_scaling = " << sourceTermScaling << endl;
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

    //   // Set up the AMRLevel... factory
    //   AMRLevelPolytropicGasFactory amrGodFact;
    // 
    //   amrGodFact.define(cfl,
    //                     domainLength,
    //                     verbosity,
    //                     refineThresh,
    //                     tagBufferSize,
    //                     initialCFL,
    //                     godunovPhysics,
    //                     normalPredOrder,
    //                     useFourthOrderSlopes,
    //                     usePrimLimiting,
    //                     useCharLimiting,
    //                     useFlattening,
    //                     useArtificialViscosity,
    //                     artificialViscosity,
    //                     useSourceTerm,
    //                     sourceTermScaling,
    //                     highOrderLimiter);
    // 
    //   AMR amr;
    // 
    //   // Set up the AMR object
    //   amr.define(maxLevel,refRatios,probDomain,&amrGodFact);
    // 
    //   if (fixedDt > 0)
    //   {
    //     amr.fixedDt(fixedDt);
    //   }
    // 
    //   // Set grid generation parameters
    //   amr.maxGridSize(maxGridSize);
    //   amr.blockFactor(blockFactor);
    //   amr.fillRatio(fillRatio);
    // 
    //   // The hyperbolic codes use a grid buffer of 1
    //   amr.gridBufferSize(1);
    // 
    //   // Set output parameters
    //   amr.checkpointInterval(checkpointInterval);
    //   amr.plotInterval(plotInterval);
    //   amr.regridIntervals(regridIntervals);
    //   amr.maxDtGrow(maxDtGrowth);
    //   amr.dtToleranceFactor(dtToleranceFactor);
    // 
    //   // Set up output files
    //   if (ppgodunov.contains("plot_prefix"))
    //   {
    //     std::string prefix;
    //     ppgodunov.query("plot_prefix",prefix);
    //     amr.plotPrefix(prefix);
    //   }
    // 
    //   if (ppgodunov.contains("chk_prefix"))
    //   {
    //     std::string prefix;
    //     ppgodunov.query("chk_prefix",prefix);
    //     amr.checkpointPrefix(prefix);
    //   }
    // 
    //   amr.verbosity(verbosity);
    // 
    //   // Set up input files
    //   if (!ppgodunov.contains("restart_file"))
    //   {
    //     if (!ppgodunov.contains("fixed_hierarchy"))
    //     {
    //       // initialize from scratch for AMR run
    //       // initialize hierarchy of levels
    //       amr.setupForNewAMRRun();
    //     }
    //     else
    //     {
    //       std::string gridFile;
    //       ppgodunov.query("fixed_hierarchy",gridFile);
    // 
    //       // initialize from a list of grids in "gridFile"
    //       Vector<Vector<Box> > amrGrids(maxLevel+1);
    //       setupFixedGrids(amrGrids,
    //                       probDomain,
    //                       maxLevel,
    //                       maxGridSize,
    //                       blockFactor,
    //                       verbosity,
    //                       gridFile);
    //       amr.setupForFixedHierarchyRun(amrGrids,1);
    //     }
    //   }
    //   else
    //   {
    //     std::string restartFile;
    //     ppgodunov.query("restart_file",restartFile);
    // 
    // #ifdef CH_USE_HDF5
    //     HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
    //     // read from checkpoint file
    //     amr.setupForRestart(handle);
    //     handle.close();
    // #else
    //     MayDay::Error("amrGodunov restart only defined with hdf5");
    // #endif
    //   }
    // 
    //   // End timing AMR solver setup
    //   TimeSetupAMR.stop();
    // 
    // #ifndef CH_NTIMER
    //   pout() << "AMR Setup completed ---- "
    //          << "mem: "
    //          << setw(8) << setprecision(3)
    //          << setiosflags(ios::fixed)
    //          << get_memory_usage_from_OS()
    //          << " MB, time: "
    //          << setw(8) << setprecision(3)
    //          << setiosflags(ios::fixed)
    //          << TimeSetupAMR.wc_time()
    //          << " sec (wall-clock)" << endl;
    // 
    //   if (verbosity >= 1)
    //   {
    //     pout() << endl;
    //   }
    // #endif
    // 
    //   // Run and time the computation
    //   TimeRun.start();
    //   amr.run(stopTime,nstop);
    //   TimeRun.stop();
    // 
    // #ifndef CN_NTIMER
    //   if (verbosity >= 1)
    //   {
    //     pout() << endl;
    //   }
    // 
    //   pout() << "AMR Run completed ------ "
    //          << "mem: "
    //          << setw(8) << setprecision(3)
    //          << setiosflags(ios::fixed)
    //          << get_memory_usage_from_OS()
    //          << " MB, time: "
    //          << setw(8) << setprecision(3)
    //          << setiosflags(ios::fixed)
    //          << TimeRun.wc_time()
    //          << " sec (wall-clock)" << endl;
    // #endif
    // 
    //   // Output the last plot file and statistics - time the process
    //   TimeConclude.start();
    //   amr.conclude();
    //   TimeConclude.stop();
    // 
    // #ifndef CH_NTIMER
    //   pout() << "AMR Conclude completed - "
    //          << "mem: "
    //          << setw(8) << setprecision(3)
    //          << setiosflags(ios::fixed)
    //          << get_memory_usage_from_OS()
    //          << " MB, time: "
    //          << setw(8) << setprecision(3)
    //          << setiosflags(ios::fixed)
    //          << TimeConclude.wc_time()
    //          << " sec (wall-clock)" << endl;
    // #endif
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
