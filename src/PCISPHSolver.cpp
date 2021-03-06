#include "PCISPHSolver.h"
#include "Simulation.h"
#include "CubeBoundaryModel.h"
#include "PCISPHBoundaryModel.h"
#include "AkinciBoundaryModel.h"
#include "Poly6.h"
#include "Spiky.h"
#include "CubicSpline.h"
#include "Logger.h"
#include <iostream>

void PCISPHSolver::init()
{
    //maxIterations = 10;
    resizeData();
    computeScalingFactor();
    
    SPHSolver::init();
}

void PCISPHSolver::step()
{
    Simulation *sim = Simulation::getCurrent();
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

    sim -> startCounting("npForces         ");
    sim -> computeNonPressureForces();
    sim -> stopCounting("npForces         ");

    // Inicializar presion y fuerza de presion
    initPressure();

    // Bucle predictivo correctivo para calcular la presion 
    sim -> startCounting("Pressure solver  ");
    pressureSolver();
    sim -> stopCounting("Pressure solver  ");

    sim -> startCounting("Integration      ");
    integrate();
    sim -> stopCounting("Integration      ");

    if (sim -> getBoundaryHandlingMethod() == Simulation::CUBE_BOUNDARY_METHOD)
    {
        sim -> startCounting("Boundary handling");
        for (unsigned int nBoundary = 0; nBoundary < sim -> numberBoundaryModels(); ++nBoundary)
        {
            CubeBoundaryModel* bm = static_cast<CubeBoundaryModel*>(sim -> getBoundaryModel(nBoundary));

            bm -> correctPositionAndVelocity();
        }
        sim -> stopCounting("Boundary handling");
    }
    else if (sim -> getBoundaryHandlingMethod() == Simulation::PCISPH_BOUNDARY_METHOD)
    {
        sim -> startCounting("Boundary handling");
        for (unsigned int nBoundary = 0; nBoundary < sim -> numberBoundaryModels(); ++nBoundary)
        {
            PCISPHBoundaryModel* pcibm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nBoundary));

            pcibm -> correctPositions();
            pcibm -> correctVelocities();
        }
        sim -> stopCounting("Boundary handling");
    }

    sim -> emitParticles();

    sim -> setTime(sim -> getTime() + sim -> getTimeStep());
}

void PCISPHSolver::initPressure()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    unsigned int nBoundaryModels = sim -> numberBoundaryModels();

    // Inicializar a 0 presion y aceleracion de presion de todos los fluidModels
    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    { 
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumActiveParticles();

        #pragma omp parallel for
        for (unsigned int i = 0; i < numParticles; ++i)
        {
            fm -> setPressure(i, 0.0);
            Vector3r zero(0.0, 0.0, 0.0);
            setPressureAcc(fmIndex, i, zero);
        }
    }

    // Inicializar a 0 la presion de las boundary particles en el metodo de bh de pcisph
    if (sim -> getBoundaryHandlingMethod() == Simulation::PCISPH_BOUNDARY_METHOD)
    {
        for (unsigned int bmIndex = 0; bmIndex < nBoundaryModels; ++bmIndex)
        { 
            PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(bmIndex));
            unsigned int numParticles = bm -> getNumParticles();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                Real & pressure = bm -> getPressure(i);
                pressure = 0.0;
            }
        }
    }
}

void PCISPHSolver::predictVelocityAndPosition()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    Real ts = sim -> getTimeStep();

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
            Vector3r& pacc = getPressureAcc(fmIndex, i);

            Vector3r & predR = getPredR(fmIndex, i);
            Vector3r & predV = getPredV(fmIndex, i);
            Vector3r acc = a + pacc + sim -> getGravity();

            predV = v + acc * ts;
            predR = r + predV * ts;
        }
    }
}

void PCISPHSolver::predictDensities()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    unsigned int nBoundaryModels = sim -> numberBoundaryModels();

    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        predictFluidDensities(fmIndex);

    for (unsigned int bmIndex = 0; bmIndex < nBoundaryModels; ++bmIndex)
        predictBoundaryDensities(bmIndex);
}

void PCISPHSolver::predictFluidDensities(const unsigned int fmIndex)
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
        Vector3r & ri = getPredR(fmIndex, i);

        Vector3i cellId = floor(fm -> getPosition(i) / sim -> getSupportRadius());

        density = 0;

        // contribuci??n de las part??culas de fluido vecinas
        forall_fluid_neighbors_in_same_phase
        (
            Vector3r & rj = getPredR(fmIndex, j);
            density += fm -> getMass(j) * CubicSpline::W(ri - rj);
        );

        // casos de los distintos boudarymodels  
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

void PCISPHSolver::predictBoundaryDensities(const unsigned int bmIndex)
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
                Vector3r & rj = getPredR(nfmIndex, j);
                density += nfm -> getMass(j) * Poly6::W(ri - rj);
            );

            forall_boundary_neighbors
            (
                PCISPHBoundaryModel *nbm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));
                Vector3r & rb = nbm -> getPosition(b);
                density += nbm -> getMass() * Poly6::W(ri - rb);
            );
        }
    }
}

void PCISPHSolver::updatePressure()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    unsigned int nBoundaryModels = sim -> numberBoundaryModels();

    maxDensityError = 0;

    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    { 
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumActiveParticles();
        Real density0 = fm -> getRefDensity();

        #pragma omp parallel for
        for (unsigned int i = 0; i < numParticles; ++i)
        {
            Real & pressure = fm -> getPressure(i);
            Real error = glm::max(static_cast<Real>(0.0), fm -> getDensity(i) - density0);

            pressure += getScalingFactor(fmIndex) * error;

            // maxError esta normalizado por lo que tambien se normaliza maxDensityError
            maxDensityError = glm::max(maxDensityError, error / density0);
        }
    }

    if (sim -> getBoundaryHandlingMethod() == Simulation::PCISPH_BOUNDARY_METHOD)
    {
        for (unsigned int bmIndex = 0; bmIndex < nBoundaryModels; ++bmIndex)
        {
            PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(bmIndex));
            unsigned int numParticles = bm -> getNumParticles();
            Real density0 = bm -> getRefDensity();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                Real & pressure = bm -> getPressure(i);
                Real error = glm::max(static_cast<Real>(0.0), bm -> getDensity(i) - density0);

                // Como el scaling factor depende de la density0 coger un scaling factor de un fluido que tenga la misma density0
                pressure += getScalingFactor(0) * error;
            }
        }        
    }
}

void PCISPHSolver::computePressureAcc()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int nFluidModels = sim -> numberFluidModels();
    int boundaryMethod = sim -> getBoundaryHandlingMethod();

    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    { 
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumActiveParticles();
        Real density0 = fm -> getRefDensity();

        #pragma omp parallel for
        for (unsigned int i = 0; i < numParticles; ++i)
        {
            Vector3r & ri = fm -> getPosition(i);
            Vector3r & pacc = getPressureAcc(fmIndex, i);
            Real & press_i = fm -> getPressure(i);
            Real & dens_i = fm -> getDensity(i);

            Vector3i cellId = floor(ri / sim -> getSupportRadius());

            pacc = Vector3r(0.0, 0.0, 0.0);

            Real densPress_i = press_i / (dens_i * dens_i);
            forall_fluid_neighbors_in_same_phase
            (
                Vector3r & rj = fm -> getPosition(j);
                Real & press_j = fm -> getPressure(j);
                Real & dens_j = fm -> getDensity(j);

                pacc -= fm -> getMass(j) * (densPress_i + press_j / (dens_j * dens_j)) * Spiky::gradW(ri - rj);
            );

            if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
            {
                forall_boundary_neighbors
                (
                    PCISPHBoundaryModel *nbm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                    Vector3r & rb = nbm -> getPosition(b);
                    Real & press_b = nbm -> getPressure(b);
                    Real & dens_b = nbm -> getDensity(b);

                    pacc -= nbm -> getMass() * (densPress_i + press_b / (dens_b * dens_b)) * Spiky::gradW(ri - rb);                
                );
            }
            else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
            {
                forall_boundary_neighbors
                (
                    AkinciBoundaryModel *nbm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                    Vector3r & rb = nbm -> getPosition(b);

                    pacc -= density0 * nbm -> getVolume(b) * (densPress_i) * Spiky::gradW(ri - rb);                
                );
            }
        }
    }
}

void PCISPHSolver::pressureSolver()
{
    Simulation *sim = Simulation::getCurrent();

    iterations = 0;
    while ((maxDensityError > maxError || iterations < minIterations) && iterations < maxIterations)
    {
        predictVelocityAndPosition();

        if (sim -> getBoundaryHandlingMethod() == Simulation::PCISPH_BOUNDARY_METHOD)
        {
            for (unsigned int nBoundary = 0; nBoundary < sim -> numberBoundaryModels(); ++nBoundary)
            {
                PCISPHBoundaryModel* pcibm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nBoundary));

                pcibm -> correctPredPositions();
            }
        }

        predictDensities();
        updatePressure();
        computePressureAcc();

        ++iterations;

        LOG("It -> ", iterations, "           -> ", maxDensityError * 100.0, "%");
    }
}

void PCISPHSolver::computeScalingFactor()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    Real ts = sim -> getTimeStep();

    Real dist = 2.0 * sim -> getParticleRadius(); 
    Vector3r sumGrad(0, 0, 0);      // Sumatorio del gradiente
    Real sumSqGrad = 0;             // Sumatorio del cuadrado del gradiente

    for (int i = -1; i <= 1; i++)
        for (int j = -1; j <= 1; j++)
            for (int k = -1; k <= 1; k++)
            {
                Vector3r rij(-i, -j, -k);
                rij *= dist;

                Vector3r grad = Poly6::gradW(rij); // Parece que funciona mejor con el standard kernel que con el pkernel

                sumGrad += grad;
                sumSqGrad += glm::dot(grad, grad);
            }

    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    {
        FluidModel *fm = sim -> getFluidModel(fmIndex);

        Real beta = 2 * pow(fm -> getMass(0) * ts / fm -> getRefDensity(), 2);
        Real dotSumGrad = glm::dot(sumGrad, sumGrad);

        setScalingFactor(fmIndex, - 1.0 / (beta * (- dotSumGrad - sumSqGrad)));
    }
}

void PCISPHSolver::resizeData()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();

    predR.resize(nFluidModels);
    predV.resize(nFluidModels);
    pressureAcc.resize(nFluidModels);
    scalingFactor.resize(nFluidModels);

    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    {
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumParticles();

        predR[fmIndex].resize(numParticles);
        predV[fmIndex].resize(numParticles);
        pressureAcc[fmIndex].resize(numParticles);
    }
}

void PCISPHSolver::integrate()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    Real ts = sim -> getTimeStep();

    maxVel = 0; // Antes de la integracion de todos los modelos o antes de la integracion de cada modelo

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

            a += sim -> getGravity() + getPressureAcc(fmIndex, i);

            v += a * ts;
            r += v * ts;

            a = Vector3r(0, 0, 0);

            maxVel = glm::max(maxVel, length(v));
        }
    }
}