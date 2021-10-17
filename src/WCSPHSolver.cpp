#include "WCSPHSolver.h"
#include "Simulation.h"
#include "CubeBoundaryModel.h"
#include "PCISPHBoundaryModel.h"
#include "AkinciBoundaryModel.h"
#include "Spiky.h"
#include <iostream>

void WCSPHSolver::init()
{
    SPHSolver::init();
}

void WCSPHSolver::step()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();

    std::cout << "stiffness " << stiffness << std::endl;

    // Insertar las particulas de todos los fluidModels en el grid
    sim -> startCounting();
    insertFluidParticles();
    sim -> stopCounting();
    std::cout << "Fill grid                 -> " << sim -> getInterval() << std::endl;

    // Hacer la busqueda de vecinos de todos los fluidModels y boundaryModels
    sim -> startCounting();
    neighborhoodSearch();
    sim -> stopCounting();
    std::cout << "Search neighboorhoods     -> " << sim -> getInterval() << std::endl;

    // Calcular densidad 
    sim -> startCounting();
    computeDensities();
    sim -> stopCounting();
    std::cout << "Compute density           -> " << sim -> getInterval() << std::endl;

    // Calcular presion
    sim -> startCounting();
    computePressures();
    sim -> stopCounting();
    std::cout << "Compute pressure          -> " << sim -> getInterval() << std::endl;

    // Calcular fuerza de presion de cada fluidModel
    sim -> startCounting();
    for (unsigned int i = 0; i < nFluidModels; ++i)
    {
        computePressureForce(i);
    }
    sim -> stopCounting();
    std::cout << "Compute pressure force    -> " << sim -> getInterval() << std::endl;

    // Calcular fuerzas de no presion para todos los fluid models
    for (unsigned int i = 0; i < nFluidModels; ++i)
    {
        FluidModel *fm = sim -> getFluidModel(i);
        unsigned int nNonPressureForces = fm -> numberNonPressureForces();

        for (unsigned int j = 0; j < nNonPressureForces; ++j)
        {
            NonPressureForce *force = fm -> getNonPressureForce(j);
            sim -> startCounting();
            force -> step();
            sim -> stopCounting();
        }
    }
    
    std::cout << "Compute nonpressure force -> " << sim -> getInterval() << std::endl;

    updateTimeStep();

    sim -> startCounting();
    integrate();
    sim -> stopCounting();
    std::cout << "Compute integration       -> " << sim -> getInterval() << std::endl;

    sim -> startCounting();
    if (sim -> getBoundaryHandlingMethod() == Simulation::CUBE_BOUNDARY_METHOD)
        for (unsigned int nBoundary = 0; nBoundary < sim -> numberBoundaryModels(); ++nBoundary)
        {
            CubeBoundaryModel* bm = static_cast<CubeBoundaryModel*>(sim -> getBoundaryModel(nBoundary));

            bm -> correctPositionAndVelocity();
        }
    else if (sim -> getBoundaryHandlingMethod() == Simulation::PCISPH_BOUNDARY_METHOD)
    {
        for (unsigned int nBoundary = 0; nBoundary < sim -> numberBoundaryModels(); ++nBoundary)
        {
            PCISPHBoundaryModel* pcibm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nBoundary));

            pcibm -> correctPositions();
            pcibm -> correctVelocities();
        }
    }
    sim -> stopCounting();
    std::cout << "Boundary handling         -> " << sim -> getInterval() << std::endl;

    std::cout << "Pre emit" << std::endl;
    sim -> emitParticles();
    std::cout << "Post emit" << std::endl;

    sim -> setTime(sim -> getTime() + sim -> getTimeStep());
}

void WCSPHSolver::computeFluidPressures(const unsigned int fmIndex)
{
    Simulation *sim = Simulation::getCurrent();
    FluidModel *fm = sim -> getFluidModel(fmIndex);
    unsigned int numParticles = fm -> getNumActiveParticles();
    Real density0 = fm -> getRefDensity();

    #pragma omp parallel for
    for (unsigned int i = 0; i < numParticles; ++i)
    {
        Real & pressure = fm -> getPressure(i);
        Real & density = fm -> getDensity(i);

        pressure = glm::max(static_cast<Real>(0.0), stiffness * density0 / gamma * (glm::pow(density / density0, gamma) - 1));
    }
}

void WCSPHSolver::computeBoundaryPressures(const unsigned int bmIndex)
{
    Simulation *sim = Simulation::getCurrent();

    if (sim -> getBoundaryHandlingMethod() == Simulation::PCISPH_BOUNDARY_METHOD)
    {
        PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(bmIndex));
        unsigned int numParticles = bm -> getNumParticles();
        Real density0 = bm -> getRefDensity();

        #pragma omp parallel for
        for (unsigned int i = 0; i < numParticles; ++i)
        {
            Real & pressure = bm -> getPressure(i);
            Real & density = bm -> getDensity(i);

            pressure = glm::max(static_cast<Real>(0.0), stiffness * density0 / gamma * (glm::pow(density / density0, gamma) - 1));
        }
    }
}

void WCSPHSolver::computePressures()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    unsigned int nBoundaryModels = sim -> numberBoundaryModels();

    for (unsigned int i = 0; i < nFluidModels; ++i)
        computeFluidPressures(i);

    for (unsigned int i = 0; i < nBoundaryModels; ++i)
        computeBoundaryPressures(i);
}


void WCSPHSolver::computePressureForce(const unsigned int fmIndex)
{
    Simulation *sim = Simulation::getCurrent();
    FluidModel *fm = sim -> getFluidModel(fmIndex);
    HashTable *grid = sim -> getGrid();
    unsigned int numParticles = fm -> getNumActiveParticles();
    int boundaryMethod = sim -> getBoundaryHandlingMethod();
    Real density0 = fm -> getRefDensity();

    #pragma omp parallel for
    for (unsigned int i = 0; i < numParticles; ++i)
    {
        Vector3r & ri = fm -> getPosition(i);
        Vector3r & ai = fm -> getAcceleration(i);
        Real & press_i = fm -> getPressure(i);
        Real & dens_i = fm -> getDensity(i);

        Vector3i cellId = floor(ri / sim -> getSupportRadius());

        forall_fluid_neighbors_in_same_phase
        (
            Vector3r & rj = fm -> getPosition(j);
            Real & press_j = fm -> getPressure(j);
            Real & dens_j = fm -> getDensity(j);

            ai -= fm -> getMass(j) * (press_i / (dens_i * dens_i) + press_j / (dens_j * dens_j)) * Spiky::gradW(ri - rj);
        );

        if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
        {
            forall_boundary_neighbors
            (
                PCISPHBoundaryModel *nbm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                Vector3r & rb = nbm -> getPosition(b);
                Real & press_b = nbm -> getPressure(b);
                Real & dens_b = nbm -> getDensity(b);

                ai -= nbm -> getMass() * (press_i / (dens_i * dens_i) + press_b / (dens_b * dens_b)) * Spiky::gradW(ri - rb);     
            );
        }
        else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
        {
            forall_boundary_neighbors
            (
                AkinciBoundaryModel *nbm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                Vector3r & rb = nbm -> getPosition(b);

                ai -= density0 * nbm -> getVolume(b) * (press_i / (dens_i * dens_i) + press_i / (density0 * density0)) * Spiky::gradW(ri - rb); 
            );
        }
    }
}

void WCSPHSolver::updateTimeStep()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    Real ts = sim -> getTimeStep();

    maxVel = 0.1;
    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    {
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumActiveParticles();

        for (unsigned int i = 0; i < numParticles; ++i)
        {
            Vector3r & vi = fm -> getVelocity(i);
            Vector3r & ai = fm -> getAcceleration(i);
            Vector3r acc = ai + sim -> getGravity();

            Vector3r current_vi = vi + acc * ts;

            maxVel = glm::max(maxVel, length(current_vi));
        }
    }

    Real cflFactor = 0.5;
    Real newTs = cflFactor * 0.4 * 2.0 * sim -> getParticleRadius() / maxVel;

    std::cout << "Propuesta ts " << newTs << std::endl;

    if (newTs > ts)
        newTs = ts * 1.005;
    /*else if (newTs < ts)
        newTs = ;*/
        
    newTs = glm::max(newTs, sim -> getMinTimeStep()); 
    newTs = glm::min(newTs, sim -> getMaxTimeStep()); 

    sim -> setTimeStep(newTs);
    
    std::cout << "TimeStep " << sim -> getTimeStep() << std::endl;
}
