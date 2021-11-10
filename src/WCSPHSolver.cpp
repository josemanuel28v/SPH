#include "WCSPHSolver.h"
#include "Simulation.h"
#include "CubeBoundaryModel.h"
#include "PCISPHBoundaryModel.h"
#include "AkinciBoundaryModel.h"
#include "Spiky.h"
#include "CubicSpline.h"

void WCSPHSolver::init()
{
    SPHSolver::init();
}

void WCSPHSolver::step()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    ++steps;

    sim -> startCounting("Fill grid        ");
    insertFluidParticles();
    sim -> stopCounting("Fill grid        ");

    sim -> startCounting("Neigh search     ");
    neighborhoodSearch();
    sim -> stopCounting("Neigh search     ");

    sim -> startCounting("Densities        ");
    computeDensities();
    sim -> stopCounting("Densities        ");

    sim -> startCounting("Pressures        ");
    computePressures();
    sim -> stopCounting("Pressures        ");

    sim -> startCounting("Pressure force   ");
    for (unsigned int i = 0; i < nFluidModels; ++i)
        computePressureForce(i);
    sim -> stopCounting("Pressure force   ");

    sim -> startCounting("NPForces         ");
    sim -> computeNonPressureForces();
    sim -> stopCounting("NPForces         ");
    
    //updateTimeStep();

    sim -> startCounting("Integratio       ");
    integrate();
    sim -> stopCounting("Integration      ");

    if (sim -> getBoundaryHandlingMethod() == Simulation::CUBE_BOUNDARY_METHOD)
    {
        sim -> startCounting("Boundary Handling");
        for (unsigned int nBoundary = 0; nBoundary < sim -> numberBoundaryModels(); ++nBoundary)
        {
            CubeBoundaryModel* bm = static_cast<CubeBoundaryModel*>(sim -> getBoundaryModel(nBoundary));

            bm -> correctPositionAndVelocity();
        }
        sim -> startCounting("Boundary Handling");
    }
    else if (sim -> getBoundaryHandlingMethod() == Simulation::PCISPH_BOUNDARY_METHOD)
    {
        sim -> startCounting("Boundary Handling");
        for (unsigned int nBoundary = 0; nBoundary < sim -> numberBoundaryModels(); ++nBoundary)
        {
            PCISPHBoundaryModel* pcibm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nBoundary));

            pcibm -> correctPositions();
            pcibm -> correctVelocities();
        }
        sim -> stopCounting("Boundary Handling");
    }

    sim -> emitParticles();

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

        Real densPress_i = press_i / (dens_i * dens_i);

        forall_fluid_neighbors_in_same_phase
        (
            Vector3r & rj = fm -> getPosition(j);
            Real & press_j = fm -> getPressure(j);
            Real & dens_j = fm -> getDensity(j);

            ai -= fm -> getMass(j) * (densPress_i + press_j / (dens_j * dens_j)) * CubicSpline::gradW(ri - rj);
        );

        if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
        {
            forall_boundary_neighbors
            (
                PCISPHBoundaryModel *nbm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                Vector3r & rb = nbm -> getPosition(b);
                Real & press_b = nbm -> getPressure(b);
                Real & dens_b = nbm -> getDensity(b);

                ai -= nbm -> getMass() * (densPress_i + press_b / (dens_b * dens_b)) * CubicSpline::gradW(ri - rb);     
            );
        }
        else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
        {
            forall_boundary_neighbors
            (
                AkinciBoundaryModel *nbm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                Vector3r & rb = nbm -> getPosition(b);

                // Segun el paper es 1.0 * densPress_i, pero como este metodo permite compresion, ese valor no es suficiente evitar la penetracioin
                ai -= density0 * nbm -> getVolume(b) * (densPress_i * 1.0) * CubicSpline::gradW(ri - rb); 
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

    if (newTs > ts)
        newTs = ts * 1.005;
    /*else if (newTs < ts)
        newTs = ;*/
        
    newTs = glm::max(newTs, sim -> getMinTimeStep()); 
    newTs = glm::min(newTs, sim -> getMaxTimeStep()); 

    sim -> setTimeStep(newTs);
    
}
