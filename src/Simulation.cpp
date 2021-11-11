#include "Simulation.h" 
#include "WCSPHSolver.h"
#include "PCISPHSolver.h"
#include "DFSPHSolver.h"
#include "Poly6.h"
#include "Spiky.h"
#include "ViscoK.h"
#include "Cohesion.h"
#include "Adhesion.h"
#include "CubicSpline.h"
#include "Logger.h"
#include "SceneLoader.h"
#include <filesystem>

#include "AkinciBoundaryModel.h"
#include "PCISPHBoundaryModel.h"
#include "CubeBoundaryModel.h"

Simulation* Simulation::current = nullptr;

const int Simulation::WCSPH_METHOD = 0;
const int Simulation::PCISPH_METHOD = 1;
const int Simulation::DFSPH_METHOD = 2;

const int Simulation::CUBE_BOUNDARY_METHOD = 0;
const int Simulation::PCISPH_BOUNDARY_METHOD = 1;
const int Simulation::AKINCI_BOUNDARY_METHOD = 2;

const int Simulation::STANDARD_VISCOSITY_METHOD = 0;
const int Simulation::ARTIFICIAL_VISCOSITY_METHOD = 1;
const int Simulation::XSPH_VISCOSITY_METHOD = 2;
const int Simulation::AKINCI_SURFACE_TENSION_METHOD = 3;
const int Simulation::AKINCI_ADHESION_METHOD = 4;

const int Simulation::POLY6_KERNEL_METHOD = 0;
const int Simulation::SPIKY_KERNEL_METHOD = 1;
const int Simulation::CUBIC_KERNEL_METHOD = 2;

const int Simulation::POLY6_GRADIENT_METHOD = 0;
const int Simulation::SPIKY_GRADIENT_METHOD = 1;
const int Simulation::CUBIC_GRADIENT_METHOD = 2;

const int Simulation::SPLINE_LAPLACIAN_METHOD = 0;

Simulation::Simulation()
{
    // Inicializar  parametros por defecto
    activeSave = false;
    name = "";

    current_method = -1;
    current_bh_method = CUBE_BOUNDARY_METHOD;
    current_visco_method = -1;

    current_kernel_method = -1;
    current_gradient_method = -1;
    current_laplacian_method = -1;

    solver = nullptr;
    grid = nullptr;

    grid = new HashTable();
}

Simulation::~Simulation()
{
    // Eliminar punteros
    delete solver;
    delete grid;

    // Eliminar contenido de los punteros en fluidModels y limpiar el vector
    for (unsigned i = 0; i < fluidModels.size(); ++i)
        delete fluidModels[i];
    fluidModels.clear();

    // Eliminar contenido de los punteros en boundaryModels y limpiar el vector
    for (unsigned i = 0; i < boundaryModels.size(); ++i)
        delete boundaryModels[i];
    boundaryModels.clear();

    delete current;
}

Simulation* Simulation::getCurrent()
{
    if (current == nullptr)
        current = new Simulation();

    return current;
}

void Simulation::setCurrent(Simulation* sim)
{
    current = sim;
}

bool Simulation::hasCurrent()
{
    return current != nullptr;
}

void Simulation::init()
{
    initKernels();
    initGrid();
    initMasses();
    solver -> init();

    printInfo();
}

void Simulation::initGrid()
{
    unsigned int nTotalParticles = 0;
    unsigned int nFluidModels = numberFluidModels();
    for (unsigned int i = 0; i < nFluidModels; ++i)
    {
        FluidModel *fm = getFluidModel(i);
        nTotalParticles += fm -> getNumParticles();
    }
    grid -> reserve(2 * nTotalParticles);
}

void Simulation::initKernels()
{
    Poly6::setSupportRadius(supportRadius);
    Spiky::setSupportRadius(supportRadius);
    ViscoK::setSupportRadius(supportRadius);
    Cohesion::setSupportRadius(supportRadius);
    Adhesion::setSupportRadius(supportRadius);
    CubicSpline::setSupportRadius(supportRadius);
}

void Simulation::initMasses()
{
    Real diameter = 2.0 * getParticleRadius();
    Real particleVolume = pow(diameter, 3.0);

    for (unsigned int i = 0; i < numberFluidModels(); ++i)
    {
        FluidModel *fm = getFluidModel(i);
        Real density0 = fm -> getRefDensity();

        fm -> setMasses(density0 * particleVolume);
    }
}

bool Simulation::step()
{
    Real time = getTime();
    static Real compTime = 0;
    bool save = false;

    if (time >= getStartTime() && time <= getEndTime())
    {
        LOG("--TIME-------------- ", time, " s ------------------------------");
        LOG("  FRAME -> ", tm.getFrame());
        startCounting("Step computation ");
        solver -> step();
        stopCounting("Step computation ");

        compTime += getInterval("Step computation ");
 
        LOG("Avg time step     -> ", time / (solver -> getSteps() - 1), " s");
        LOG("Avg step comp     -> ", compTime / (solver -> getSteps() - 1), " s");

        LOG("----------------------------------------------------------------");
        LOG();

        save = tm.hasToSave();

        if (save)
        {
            if (activeSave)
                SceneLoader::writeFluid();

            tm.setFrame(tm.getFrame() + 1);
        }
    }

    return save;
}

void Simulation::run()
{
    Real time = getTime();

    while (time >= getStartTime() && time <= getEndTime())
        step();
}

void Simulation::setParticleRadius(Real particleRadius)
{
    this -> particleRadius = particleRadius;

    // Actualizar valor de la distancia de suavizado
    this -> supportRadius = 4.0 * particleRadius /** pow(kernelParticles, 1./3.)*/;
    grid -> setSmoothingLength(this -> supportRadius);
}

void Simulation::setSupportRadius(Real supportRadius)
{
    this -> supportRadius = supportRadius;
    grid -> setSmoothingLength(supportRadius);
}

void Simulation::setSimulationMethod(const int method)
{
    if (method == current_method || method < 0)
        return;

    current_method = method;

    if (solver != nullptr)
        delete solver;

    switch (method)
    {
        case WCSPH_METHOD:
            solver = new WCSPHSolver();
            break;

        case PCISPH_METHOD:
            solver = new PCISPHSolver();
            break;

        case DFSPH_METHOD:
            solver = new DFSPHSolver();
            break;
    }
}

void Simulation::setBoundaryMethod(const int method)
{
    current_bh_method = method;
}

void Simulation::setKernelMethod(const int method)
{
    if (method == current_kernel_method || method < 0)
        return;

    current_kernel_method = method;

    switch (method)
    {
        case POLY6_KERNEL_METHOD:
            value = Poly6::W;  
            break;

        case SPIKY_KERNEL_METHOD:
            value = Spiky::W;
            break;
    }
}

void Simulation::setGradientMethod(const int method)
{
    if (method == current_gradient_method || method < 0)
        return;

    current_gradient_method = method;
    
    switch (method)
    {
        case POLY6_GRADIENT_METHOD:
            grad = Poly6::gradW;   
            break;

        case SPIKY_GRADIENT_METHOD:
            grad = Spiky::gradW;
            break;
    }  
}

void Simulation::setLaplacianMethod(const int method)
{
    if (method == current_laplacian_method || method < 0)
        return;

    current_laplacian_method = method;

    switch (method)
    {
        case SPLINE_LAPLACIAN_METHOD:
            lapl = ViscoK::laplW;   
            break;
    }
}

FluidModel* Simulation::addFluidModel(std::vector<Vector3r>& fluidPoints, std::vector<Vector3r>& fluidVelocities)
{
    FluidModel *fm = new FluidModel();
    fm -> init(fluidPoints, fluidVelocities);
    fluidModels.push_back(fm);

    return fm;
}

FluidModel* Simulation::addFluidModel()
{
    FluidModel *fm = new FluidModel();
    fluidModels.push_back(fm);

    return fm;
}

void Simulation::addBoundaryModel(BoundaryModel *bm)
{
    boundaryModels.push_back(bm);
}

void Simulation::computeNonPressureForces()
{
    unsigned int nFluidModels = numberFluidModels();

    for (unsigned int i = 0; i < nFluidModels; ++i)
    {
        FluidModel *fm = getFluidModel(i);
        for (unsigned int i = 0; i < fm -> numberNonPressureForces(); ++i)
            fm -> getNonPressureForce(i) -> step();
    }   
}

void Simulation::emitParticles()
{
    for (unsigned int i = 0; i < numberFluidModels(); ++i)
        fluidModels[i] -> emitParticles();
}

bool Simulation::importScene(std::string path)
{
    LOG("Importing scene from ", path, "...");
    name = path;

    SimulationInfo simData;

    if (!SceneLoader::readConfiguration(simData, path))
        return false;

    const SceneInfo & sceneData = simData.sceneData;
    const FluidInfo & fluidData = simData.fluidData;
    const BoundaryInfo & boundaryData = simData.boundaryData;

    // Scene info
    setTimeStep(sceneData.timeStep); 
    setFPS(sceneData.fps);
    setMinTimeStep(sceneData.minTimeStep); 
    setMaxTimeStep(sceneData.maxTimeStep);  
    setStartTime(sceneData.startTime);
    setEndTime(sceneData.endTime); 
    setSimulationMethod(sceneData.simulationMethod);
    setBoundaryMethod(sceneData.boundaryMethod);
    setGravity(sceneData.gravity);
    setParticleRadius(sceneData.particleRadius);

    if (sceneData.simulationMethod == WCSPH_METHOD)
    {
        WCSPHSolver *wcsph = static_cast<WCSPHSolver*>(getSolver());
        wcsph -> setStiffness(sceneData.stiffness);
        wcsph -> setGamma(sceneData.gamma);
    }
    else if (sceneData.simulationMethod == PCISPH_METHOD)
    {
        PCISPHSolver *pcisph = static_cast<PCISPHSolver*>(getSolver());
        pcisph -> setMaxError(sceneData.eta);
        pcisph -> setMinIterations(sceneData.minIterations);
        pcisph -> setMaxIterations(sceneData.maxIterations);
    }
    else if (sceneData.simulationMethod == DFSPH_METHOD)
    {
        DFSPHSolver *dfsph = static_cast<DFSPHSolver*>(getSolver());
        dfsph -> setMaxError(sceneData.eta);
        dfsph -> setMaxErrorV(sceneData.etaV);
        dfsph -> setCFLFactor(sceneData.cflFactor);
        dfsph -> setMinIterations(sceneData.minIterations);
        dfsph -> setMaxIterations(sceneData.maxIterations);
        dfsph -> setMinIterationsV(sceneData.minIterationsV);
        dfsph -> setMaxIterationsV(sceneData.maxIterationsV);
    }

    // Fluid info
    unsigned int numFluids = fluidData.fluids.size();
    const std::vector<Fluid> & fluids = fluidData.fluids; // vector de fluidModels

    setViscosityMethod(fluidData.viscosityMethod);
    setSurfaceTensionMethod(fluidData.surfaceTensionMethod);
    setAdhesionMethod(fluidData.adhesionMethod);

    for (unsigned int i = 0; i < numFluids; ++i)
    {
        FluidModel *fm = nullptr;
        
        unsigned int numBlocks = fluids[i].fluidBlocks.size();
        unsigned int numEmitters = fluids[i].emitters.size();
        //unsigned int numGeometries = fluids[i].geometries.size();
        
        if (numBlocks > 0)
            fm = buildFluidBlock(fluids[i].fluidBlocks);

        // geomettry

        // emitters
        if (!fm) // si no hay fm porque no hay fluidBlocks se crea uno con posiciones y velocidades vacias
            fm = addFluidModel();

        for (unsigned int j = 0; j < numEmitters; ++j)
        {
            fm -> addEmitter(fluids[i].emitters[j].type,
                             fluids[i].emitters[j].numParticles,
                             fluids[i].emitters[j].r,
                             fluids[i].emitters[j].v,
                             fluids[i].emitters[j].rot,
                             fluids[i].emitters[j].startTime,
                             fluids[i].emitters[j].width,
                             fluids[i].emitters[j].height,
                             fluids[i].emitters[j].spacing);

        }

        fm -> setRefDensity(fluids[i].density0);

        if (fluidData.viscosityMethod > -1 && fluids[i].viscosity != 0.0)
            fm -> setViscosityForce(fluids[i].viscosity, fluids[i].boundaryViscosity);

        if (fluidData.surfaceTensionMethod > -1 && fluids[i].surfaceTension != 0.0)
            fm -> setSurfaceTensionForce(fluids[i].surfaceTension);

        if (fluidData.adhesionMethod > -1 && fluids[i].adhesion != 0.0)
            fm -> setAdhesionForce(fluids[i].adhesion);
    }

    // Boundary info
    unsigned int numBoundaries = boundaryData.boundaries.size();
    const std::vector<Boundary> & boundaries = boundaryData.boundaries;

    if (sceneData.boundaryMethod == AKINCI_BOUNDARY_METHOD)
    {
        for (unsigned int i = 0; i < numBoundaries; ++i)
        {
            AkinciBoundaryModel *bm = new AkinciBoundaryModel();

            unsigned int numBox = boundaries[i].box.size();
            unsigned int numSphere = boundaries[i].sphere.size();
            unsigned int numGeometry = boundaries[i].geometry.size();

            for (unsigned int j = 0; j < numBox; ++j)
                bm -> addCube(boundaries[i].box[j].first.min, boundaries[i].box[j].first.max);

            for (unsigned int j = 0; j < numSphere; ++j)
                bm -> addSphere(boundaries[i].sphere[j].first.pos, boundaries[i].sphere[j].first.radius);

            for (unsigned int j = 0; j < numGeometry; ++j)
                bm -> addGeometry(boundaries[i].geometry[j].path, boundaries[i].geometry[j].spacing * sceneData.particleRadius * 2.0);
            
            addBoundaryModel(bm);
        }
    }
    else if (sceneData.boundaryMethod == PCISPH_BOUNDARY_METHOD)
    {
        if (numberFluidModels() == 1)
        {
            for (unsigned int i = 0; i < numBoundaries; ++i)
            {
                PCISPHBoundaryModel *bm = new PCISPHBoundaryModel();
                Real density0 = getFluidModel(0) -> getRefDensity(); 
                Real particleVolume = pow(2.0 * getParticleRadius(), 3.0);
                bm -> setMass(density0 * particleVolume);
                bm -> setRefDensity(density0);
                bm -> setNormalFct(boundaries[i].normalFct);
                bm -> setTangentialFct(boundaries[i].tangFct);

                unsigned int numBox = boundaries[i].box.size();
                unsigned int numSphere = boundaries[i].sphere.size();

                for (unsigned int j = 0; j < numBox; ++j)
                    bm -> addCube(boundaries[i].box[j].first.min, boundaries[i].box[j].first.max, boundaries[i].box[j].second);

                for (unsigned int j = 0; j < numSphere; ++j)
                    bm -> addSphere(boundaries[i].sphere[j].first.pos, boundaries[i].sphere[j].first.radius, boundaries[i].sphere[j].second);

                addBoundaryModel(bm);
            }
        }
        else
            LOG("PCISPH BOUNDARY METHOD is not compatible with multiphase");
    }

    LOG("Scene imported!");

    init();

    return true;
}

std::vector<Real> & Simulation::getDivergenceError(const unsigned int fmIndex)
{
    if (current -> getSimulationMethod() == DFSPH_METHOD)
    {
        DFSPHSolver *solver = static_cast<DFSPHSolver*>(this -> solver);

        return solver -> getDivergenceError(fmIndex);
    }
    else
    {
        FluidModel *fm = getFluidModel(fmIndex);

        return fm -> getPressures();
    }
}

FluidModel* Simulation::buildFluidBlock(const std::vector<BlockInfo> & fluidBlocks)
{
    Real dist = 2.0 * getParticleRadius();

    std::vector<Vector3r> position;
    std::vector<Vector3r> velocity;

    for (const BlockInfo & block: fluidBlocks)
    {
        Vector3i ppedge = floor((block.max - block.min) / dist + 1e-5);

        for (int i = 0; i < ppedge.x; ++i)
                for (int j = 0; j < ppedge.y; ++j)
                    for (int k = 0; k < ppedge.z; ++k)
                    {
                        Vector3r pos(i, j, k);
                        pos *= dist;
                        pos += block.min + getParticleRadius();

                        position.push_back(pos);
                        velocity.push_back(Vector3r(0.0));
                    }
    }

    addFluidModel(position, velocity);

    return getFluidModel(numberFluidModels() - 1);
}

void Simulation::printInfo()
{
    std::string solverLabel = "none";
    std::string bhLabel = "none";
    std::string viscoLabel = "none";
    std::string surftenLabel = "none";
    std::string adhesionLabel = "none";

    if (current_method == WCSPH_METHOD)
        solverLabel = "Weakly Compressible SPH";
    else if (current_method == PCISPH_METHOD)
        solverLabel = "Predictive-Corrective Incompressible SPH";
    else if (current_method == DFSPH_METHOD)
        solverLabel = "Divergence-Free SPH";

    if (current_bh_method == CUBE_BOUNDARY_METHOD)
        bhLabel == "Simple box boundary handling method";
    else if (current_bh_method == PCISPH_BOUNDARY_METHOD)
        bhLabel = "PCISPH boundary handling method";
    else if (current_bh_method == AKINCI_BOUNDARY_METHOD)
        bhLabel = "Akinci boundary handling method";

    if (current_visco_method == STANDARD_VISCOSITY_METHOD)
        viscoLabel = "Standard SPH formulation of viscosity (Explicit)";
    else if (current_visco_method == ARTIFICIAL_VISCOSITY_METHOD)
        viscoLabel = "Artificial viscosity method (Explicit)";
    else if (current_visco_method == XSPH_VISCOSITY_METHOD)
        bhLabel = "XSPH viscosity method";

    if (current_surften_method == AKINCI_SURFACE_TENSION_METHOD)
        surftenLabel = "Akinci surface tension method";

    if (current_adhesion_method == AKINCI_ADHESION_METHOD)
        adhesionLabel = "Akinci adhesion method";    

    LOG("--------------------------------------------------------------------------------------");
    LOG("Launching simulation with the following configuration:");
    LOG("Number of fluid models:         ", numberFluidModels(), " fm");
    LOG("Number of boundary models:      ", numberBoundaryModels(), " bm");
    for (unsigned int i = 0; i < numberFluidModels(); ++i)
        LOG("Number of particles in fm ", i + 1,  ":    ", getFluidModel(i) -> getNumParticles(), " particles");
    for (unsigned int i = 0; i < numberBoundaryModels(); ++i)
        LOG("Number of particles in bm ", i + 1, ":    ", getBoundaryModel(i) -> getNumParticles(), " particles");
    LOG("Pressure solver:                " , solverLabel);
    if (current_method == WCSPH_METHOD)
    {
        WCSPHSolver *currentSolver = static_cast<WCSPHSolver*>(solver);
        LOG("Stiffness:                      ", currentSolver -> getStiffness());
        LOG("Gamma:                          ", currentSolver -> getGamma());
    }
    else if (current_method == PCISPH_METHOD)
    {
        PCISPHSolver *currentSolver = static_cast<PCISPHSolver*>(solver);
        LOG("Max. allowed density error:     ", currentSolver -> getMaxError());
        LOG("Min. iterations:                ", currentSolver -> getMinIterations());
        LOG("Max. iterations:                ", currentSolver -> getMaxIterations());
    }
    else if (current_method == DFSPH_METHOD)
    {
        DFSPHSolver *currentSolver = static_cast<DFSPHSolver*>(solver);
        LOG("Max. allowed density error:     ", currentSolver -> getMaxError());
        LOG("Max. allowed divergence error:  ", currentSolver -> getMaxErrorV());
        LOG("Min. density iterations:        ", currentSolver -> getMinIterations());
        LOG("Max. density iterations:        ", currentSolver -> getMaxIterations());
        LOG("Min divergence iterations:      ", currentSolver -> getMinIterationsV());
        LOG("Max divergence iterations:      ", currentSolver -> getMaxIterationsV());
    }
    LOG("Boundary handling method:       ", bhLabel);
    LOG("Viscosity method:               ", viscoLabel);
    LOG("Surface tension method:         ", surftenLabel);
    LOG("Adhesion method:                ", adhesionLabel);
    LOG("Particle radius:                ", particleRadius, " m");
    LOG("Support radius:                 ", supportRadius, " m");
    LOG("Gravity:                        ", "(" , gravity.x , "," , gravity.y , "," , gravity.z , ") m/s^2");
    LOG("--------------------------------------------------------------------------------------");
    LOG();
    LOG("Press any key to start the simulation");

    getchar();
}





