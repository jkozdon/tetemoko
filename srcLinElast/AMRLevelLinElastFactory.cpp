#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRLevel.H"

#include "AMRLevelLinElastFactory.H"
#include "AMRLevelLinElast.H"

AMRLevelLinElastFactory::AMRLevelLinElastFactory()
{
    m_linElastPhysics = NULL;
    m_isDefined = false;
}

AMRLevelLinElastFactory::~AMRLevelLinElastFactory()
{
    if (m_linElastPhysics != NULL)
    {
        delete m_linElastPhysics;
        m_linElastPhysics = NULL;
    }

    m_isDefined = false;
}

void AMRLevelLinElastFactory::define(const Real&                 a_cfl,
    const Real&                 a_domainLength,
    const int&                  a_verbosity,
    const Real&                 a_refineThresh,
    const Real&                 a_plasticThresh,
    const int&                  a_tagBufferSize,
    const Real&                 a_initialDtMultiplier,
    const LinElastPhysics* const a_linElastPhysics,
    const int&                  a_normalPredOrder,
    const bool&                 a_useFourthOrderSlopes,
    const bool&                 a_usePrimLimiting,
    const bool&                 a_useCharLimiting,
    const bool&                 a_useFlattening,
    const bool&                 a_useArtificialViscosity,
    const Real&                 a_artificialViscosity,
    const bool&                 a_useSourceTerm,
    const Real&                 a_sourceTermScaling,
    const bool&                 a_highOrderLimiter,
    const Vector<Real>&         a_xCoarsen,
    const Vector<Real>&         a_yCoarsen,
    const Vector<Real>&         a_zCoarsen,
    const Real&                 a_slopeCoarsen,
    const Real&                 a_widthCoarsen,
    const Vector<Real>&         a_xFaultStations,
    const Vector<Real>&         a_zFaultStations,
    const Vector<Real>&         a_xBodyStations,
    const Vector<Real>&         a_yBodyStations,
    const Vector<Real>&         a_zBodyStations,
    const Vector<Real>&         a_domainCenter,
    const bool&                 a_usePlasticity)
{
    // Store the CFL number
    m_cfl = a_cfl;

    // Store the physical dimension of the longest side of the domain
    m_domainLength = a_domainLength;

    // Store the verbosity of the object
    m_verbosity = a_verbosity;

    // Store the refinement threshold for gradient
    m_refineThresh = a_refineThresh;
    
    // Store the refinement threshold for plasticity
    m_plasticThresh = a_plasticThresh;

    // Store the tag buffer size
    m_tagBufferSize = a_tagBufferSize;

    // Store the initial dt multiplier
    m_initialDtMultiplier = a_initialDtMultiplier;

    // Delete any existing physics object
    if (m_linElastPhysics != NULL)
    {
        delete m_linElastPhysics;
        m_linElastPhysics = NULL;
    }

    // Store the object that supplies the physics needed by the integrator
    // (used as a factory)
    m_linElastPhysics = a_linElastPhysics->new_godunovPhysics();

    // Store the order of the normal predictor (1 -> PLM, 2 -> PPM)
    m_normalPredOrder = a_normalPredOrder;

    // Store the slope computation parameters
    m_useFourthOrderSlopes = a_useFourthOrderSlopes;
    m_usePrimLimiting      = a_usePrimLimiting;
    m_useCharLimiting      = a_useCharLimiting;
    m_useFlattening        = a_useFlattening;

    // Artificial viscosity coefficient must be greater than zero
    CH_assert(!a_useArtificialViscosity || (a_artificialViscosity >= 0.0));

    // Store the artificial viscosity flag and coefficient
    m_useArtificialViscosity = a_useArtificialViscosity;
    m_artificialViscosity    = a_artificialViscosity;

    // Supply a source term to the computation
    m_useSourceTerm     = a_useSourceTerm;
    m_sourceTermScaling = a_sourceTermScaling;

    // Use a high-order limiter?
    m_highOrderLimiter = a_highOrderLimiter;

    // define the fault stations
    m_xFaultStations   = a_xFaultStations;
    if(SpaceDim > 2)
    {
        m_zFaultStations   = a_zFaultStations;
    }
    else
    {
        m_zFaultStations   = a_xFaultStations;
    }

    // define the body stations
    m_xBodyStations   = a_xBodyStations;
    m_yBodyStations   = a_yBodyStations;
    if(SpaceDim > 2)
    {
        m_zBodyStations   = a_zBodyStations;
    }
    else
    {
        m_zBodyStations   = a_xBodyStations;
    }

    // Set the center of the domain
    m_domainCenter = a_domainCenter;

    m_dataPrefix = "";

    m_plotInterval = 0;

    m_usePlasticity = a_usePlasticity;

    // The object is defined
    m_isDefined = true;

    m_xCoarsen     = a_xCoarsen;
    m_yCoarsen     = a_yCoarsen;
    m_zCoarsen     = a_zCoarsen;
    m_slopeCoarsen = a_slopeCoarsen;
    m_widthCoarsen = a_widthCoarsen;
}

// Virtual constructor
AMRLevel* AMRLevelLinElastFactory::new_amrlevel() const
{
    // Make sure everything is defined
    CH_assert(isDefined());

    // Create a new AMRLevelLinElast
    AMRLevelLinElast* amrGodPtr = new AMRLevelLinElast();

    // Define the new object
    amrGodPtr->defineParams(m_cfl,
        m_domainLength,
        m_verbosity,
        m_refineThresh,
        m_plasticThresh,
        m_tagBufferSize,
        m_initialDtMultiplier,
        m_linElastPhysics,
        m_normalPredOrder,
        m_useFourthOrderSlopes,
        m_usePrimLimiting,
        m_useCharLimiting,
        m_useFlattening,
        m_useArtificialViscosity,
        m_artificialViscosity,
        m_useSourceTerm,
        m_sourceTermScaling,
        m_highOrderLimiter,
        m_xCoarsen,
        m_yCoarsen,
        m_zCoarsen,
        m_slopeCoarsen,
        m_widthCoarsen,
        m_xFaultStations,
        m_zFaultStations,
        m_xBodyStations,
        m_yBodyStations,
        m_zBodyStations,
        m_domainCenter,
        m_dataPrefix,
        m_plotInterval,
        m_usePlasticity);

    // Return it
    return (static_cast <AMRLevel*> (amrGodPtr));
}

// Check that everything is defined
bool AMRLevelLinElastFactory::isDefined() const
{
    return m_isDefined;
}

void AMRLevelLinElastFactory::dataPrefix(const string& a_dataPrefix)
{
    m_dataPrefix = a_dataPrefix;
}

void AMRLevelLinElastFactory::plotInterval(const int& a_plotInterval)
{
    m_plotInterval = a_plotInterval;
}
