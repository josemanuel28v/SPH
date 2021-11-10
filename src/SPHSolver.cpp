#include "SPHSolver.h" 
#include "PCISPHBoundaryModel.h"
#include "AkinciBoundaryModel.h"
#include "Simulation.h"
#include "FluidModel.h"
#include "HashTable.h"
#include "Poly6.h"
#include "CubicSpline.h"
#include "Spiky.h"
#include <iostream>

SPHSolver::SPHSolver()
{
    steps = 0;
    iterations = 0;
    minIterations = 3;
    maxIterations = 100;
    sumIterations = 0;
    maxError = 0.01; // 1%
}

void SPHSolver::init()
{
    Simulation *sim = Simulation::getCurrent();

    insertBoundaryParticles();

    if (sim -> getBoundaryHandlingMethod() == Simulation::AKINCI_BOUNDARY_METHOD)
    {
        neighborhoodSearch();
        unsigned int nBoundaries = sim -> numberBoundaryModels();

        for (unsigned int bmIndex = 0; bmIndex< nBoundaries; ++bmIndex)
        {
            AkinciBoundaryModel *abm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(bmIndex));
            abm -> computeVolume();
        }
    }
}

SPHSolver::~SPHSolver()
{}

void SPHSolver::computeFluidDensities(const unsigned int fmIndex)
{
    Simulation *sim = Simulation::getCurrent();
    FluidModel *fm = sim -> getFluidModel(fmIndex);
    unsigned int numParticles = fm -> getNumActiveParticles();
    int boundaryMethod = sim -> getBoundaryHandlingMethod();
    Real density0 = fm -> getRefDensity();

    HashTable *grid = sim -> getGrid();

    // Calcular densidad para las particulas de fluido de un fluid model
    #pragma omp parallel for 
    for (unsigned int i = 0; i < numParticles; ++i)
    {
        Real & density = fm -> getDensity(i);
        Vector3r & ri = fm -> getPosition(i);

        Vector3i cellId = floor(ri / sim -> getSupportRadius());
        density = 0;

        // Contribucion de las partículas de fluido vecinas
        forall_fluid_neighbors_in_same_phase
        (
            Vector3r & rj = fm -> getPosition(j);
            density += fm -> getMass(j) * CubicSpline::W(ri - rj);
        );

        // Contribucion de todos los boundary models
        if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
        {
            forall_boundary_neighbors
            (
                PCISPHBoundaryModel *nbm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));
                Vector3r & rb = nbm -> getPosition(b);
                density += nbm -> getMass() * CubicSpline::W(ri - rb);
            );
        }
        else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
        {
            forall_boundary_neighbors
            (
                AkinciBoundaryModel *nbm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));
                Vector3r & rb = nbm -> getPosition(b);
                density += density0 * nbm -> getVolume(b) * CubicSpline::W(ri - rb);
            );
        }
    }
}

void SPHSolver::computeBoundaryDensities(const unsigned int bmIndex)
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();

    if (sim -> getBoundaryHandlingMethod() == Simulation::PCISPH_BOUNDARY_METHOD)
    {
        PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(bmIndex));
        unsigned int numParticles = bm -> getNumParticles();

        #pragma omp parallel for 
        for (unsigned int i = 0; i < numParticles; ++i)
        {
            Real & density = bm -> getDensity(i);
            Vector3r & ri = bm -> getPosition(i);

            Vector3i cellId = floor(ri / sim -> getSupportRadius());
            density = 0;

            forall_fluid_neighbors
            (
                Vector3r & rj = nfm -> getPosition(j);
                density += nfm -> getMass(j) * CubicSpline::W(ri - rj);
            );

            forall_boundary_neighbors
            (
                PCISPHBoundaryModel *nbm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));
                Vector3r & rb = nbm -> getPosition(b);
                density += nbm -> getMass() * CubicSpline::W(ri - rb);
            );
        }
    }
}

void SPHSolver::computeDensities()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    unsigned int nBoundaryModels = sim -> numberBoundaryModels();

    for (unsigned int i = 0; i < nFluidModels; ++i)
        computeFluidDensities(i);

    for (unsigned int i = 0; i < nBoundaryModels; ++i)
        computeBoundaryDensities(i);
}

void SPHSolver::integrate()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    Real ts = sim -> getTimeStep();

    maxVel = 0.0;
    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    {
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumActiveParticles();

        #pragma omp parallel for 
        for (unsigned int i = 0; i < numParticles; ++i)
        {
            Vector3r& r = fm -> getPosition(i);
            Vector3r& v = fm -> getVelocity(i);
            Vector3r& a = fm -> getAcceleration(i);

            a += sim -> getGravity();

            v += a * ts;
            r += v * ts;

            a = Vector3r(0, 0, 0);

            maxVel = glm::max(maxVel, length(v));
        }
    }
}

void SPHSolver::insertFluidParticles()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int numFluidModels = sim -> numberFluidModels();
    FluidModel *fm;

    // Vaciar el grid antes de insertar las particulas
    grid -> clear();

    // Insertar las particulas de fluido en el grid
    for (unsigned int fluidModelIndex = 0; fluidModelIndex < numFluidModels; ++fluidModelIndex)
    {
        fm = sim->getFluidModel(fluidModelIndex);
        unsigned int numParticles = fm -> getNumActiveParticles();

        for (unsigned int i = 0; i < numParticles; ++i)
        {
            grid -> insertFluidParticle(fm -> getPosition(i), i, fluidModelIndex); 
        }
    }
}

void SPHSolver::insertBoundaryParticles()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int numBoundaryModels = sim -> numberBoundaryModels();

    // Vaciar el grid antes de insertar las particulas
    grid -> clearB();

    // Insertar las particulas de fluido en el grid
    for (unsigned int bmIndex = 0; bmIndex < numBoundaryModels; ++bmIndex)
    {
        BoundaryModel *bm = sim->getBoundaryModel(bmIndex); 
        unsigned int numParticles = bm -> size();

        for (unsigned int i = 0; i < numParticles; ++i)
            grid -> insertBoundaryParticle(bm -> getPosition(i), i, bmIndex); 
    }
}

void SPHSolver::neighborhoodSearch()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();

    grid -> neighborhoodSearch();
}



