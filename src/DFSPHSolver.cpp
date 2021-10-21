#include "DFSPHSolver.h"
#include "Simulation.h"
#include "Spiky.h"
#include "Poly6.h"
#include "CubicSpline.h"
#include "CubeBoundaryModel.h"
#include "PCISPHBoundaryModel.h"
#include "AkinciBoundaryModel.h"

void DFSPHSolver::init()
{
    minIterations = 2; // Iteraciones minimas del bucle de correccion de densidad
    minIterationsV = 1; // Iteraciones minimas del bucle de correccion de divergencia

    maxIterations = 100;
    sumIterationsV = 0;

    eps = 1e-5;

    resizeData();

    SPHSolver::init();

    insertFluidParticles();
    neighborhoodSearch();
    computeDensities();
    computeAlpha();
}

void DFSPHSolver::step()
{
    Simulation *sim = Simulation::getCurrent();
    steps ++;
    
    sim -> startCounting("npForces         ");
    sim -> computeNonPressureForces();
    sim -> stopCounting("npForces         ");

    updateTimeStep();

    predictVelocities(); 

    if (sim -> getBoundaryHandlingMethod() == Simulation::PCISPH_BOUNDARY_METHOD)
        for (unsigned int nBoundary = 0; nBoundary < sim -> numberBoundaryModels(); ++nBoundary)
        {
            PCISPHBoundaryModel* bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nBoundary));

            bm -> correctVelocities(); // cuidado: necesita las normales calculadas en predictPositions() y si se calcula fuerza de st se sobreescriben
                                       // en el primer step las normales estaran a 0 pero igualmente las velocidades tambien estaran a 0 entonces no pasaria nada
                                       // y a la siguiente vuelta ya se habran calculado las normales en predictPositions()
        }

    sim -> startCounting("Density solver   ");
    correctDensityError();
    sim -> stopCounting("Density solver   ");

    updatePositions();

    if (sim -> getBoundaryHandlingMethod() == Simulation::PCISPH_BOUNDARY_METHOD)
    {
        sim -> startCounting("Boundary handling");
        for (unsigned int nBoundary = 0; nBoundary < sim -> numberBoundaryModels(); ++nBoundary)
        {
            PCISPHBoundaryModel* bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nBoundary));

            bm -> correctPositions();
        }
        sim -> stopCounting("Boundary handling");
    }

    sim -> emitParticles();

    sim -> startCounting("Fill grid        ");
    insertFluidParticles();
    sim -> stopCounting("Fill grid        ");

    sim -> startCounting("Neigh search     ");
    neighborhoodSearch();
    sim -> stopCounting("Neigh search     ");

    sim -> startCounting("Densities        ");
    computeDensities();
    sim -> stopCounting("Densities        ");    

    sim -> startCounting("Alpha factor     ");
    computeAlpha();
    sim -> stopCounting("Alpha factor     ");

    sim -> startCounting("Divergence solver");
    correctDivergenceError();
    sim -> stopCounting("Divergence solver");

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

    sim -> setTime(sim -> getTime() + sim -> getTimeStep());
}

void DFSPHSolver::computeAlpha()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int nFluidModels = sim -> numberFluidModels();
    Real supportRadius = sim -> getSupportRadius();
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
            Real & alpha = getAlpha(fmIndex, i);

            Vector3i cellId = floor(ri / supportRadius);

            Vector3r sumGrad(0.0, 0.0, 0.0);
            Real sumSqGrad = 0.0;

            forall_fluid_neighbors_in_same_phase
            (
                Vector3r & rj = fm -> getPosition(j);
                Vector3r massGrad = fm -> getMass(j) * CubicSpline::gradW(ri - rj);
                Real massGradMag = length(massGrad);

                sumGrad += massGrad;
                sumSqGrad +=  massGradMag * massGradMag;
            );

            if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
            {
                forall_boundary_neighbors
                (
                    PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                    Vector3r & rb = bm -> getPosition(b);
                    Vector3r massGrad = bm -> getMass() * CubicSpline::gradW(ri - rb);

                    //Real massGradMag = length(massGrad);

                    sumGrad += massGrad;
                    //sumSqGrad += massGradMag * massGradMag; // aumenta el numero de iteraciones del bucle corrector de densidad
                );
            }
            else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
            {
                forall_boundary_neighbors
                (
                    AkinciBoundaryModel *bm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                    Vector3r & rb = bm -> getPosition(b);
                    Vector3r massGrad = density0 * bm -> getVolume(b) * CubicSpline::gradW(ri - rb);

                    //Real massGradMag = length(massGrad);

                    sumGrad += massGrad;
                    //sumSqGrad += massGradMag * massGradMag; // aumenta el numero de iteraciones del bucle corrector de densidad
                );
            }

            Real sumGradMag = length(sumGrad);
            Real denominator = sumGradMag * sumGradMag + sumSqGrad;

            if (fabs(denominator) > eps)
                alpha = fm -> getDensity(i) / denominator; // importante tiene signo negativo pero a empezado a funcionar al quitarselo
            else 
                alpha = 0.0;
        }
    }
}

void DFSPHSolver::predictVelocities()
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
            Vector3r & vi = fm -> getVelocity(i);
            Vector3r & ai = fm -> getAcceleration(i);
            Vector3r acc = ai + sim -> getGravity();

            vi += acc * ts;

            ai = Vector3r(0.0, 0.0, 0.0); // limpiar la aceleracion ya que es el unico sitio donde se usa
        }
    }
}

void DFSPHSolver::updatePositions()
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
            Vector3r & ri = fm -> getPosition(i);
            Vector3r & vi = fm -> getVelocity(i);

            ri += vi * ts;
        }
    }
}

void DFSPHSolver::predictDensities()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int nFluidModels = sim -> numberFluidModels();
    Real supportRadius = sim -> getSupportRadius();
    Real ts = sim -> getTimeStep();
    unsigned int totalNumParticles = 0;
    int boundaryMethod = sim -> getBoundaryHandlingMethod();

    avgDensError = 0.0;
    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    {
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumActiveParticles();
        Real density0 = fm -> getRefDensity();
        totalNumParticles += numParticles;

        #pragma omp parallel for reduction(+:avgDensError)
        for (unsigned int i = 0; i < numParticles; ++i)
        {
            Vector3r & ri = fm -> getPosition(i);
            Vector3r & vi = fm -> getVelocity(i);
            Real & pred_density = getPredDensity(fmIndex, i);

            Real densityChange = 0;

            Vector3i cellId = floor(ri / supportRadius);

            forall_fluid_neighbors_in_same_phase
            (
                Vector3r & rj = fm -> getPosition(j);
                Vector3r & vj = fm -> getVelocity(j);

                densityChange += fm -> getMass(j) * dot(vi - vj, CubicSpline::gradW(ri - rj));
            );

            if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
            {
                forall_boundary_neighbors
                (
                    PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                    Vector3r & rb = bm -> getPosition(b);
                    densityChange += bm -> getMass() * dot(vi, CubicSpline::gradW(ri - rb));
                );
            }
            else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
            {
                forall_boundary_neighbors
                (
                    AkinciBoundaryModel *bm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                    Vector3r & rb = bm -> getPosition(b);
                    densityChange += density0 * bm -> getVolume(b) * dot(vi, CubicSpline::gradW(ri - rb));
                );
            }

            pred_density = fm -> getDensity(i) + densityChange * ts;
            pred_density = glm::max(pred_density, density0); // importante que la densidad predicha nunca sea menor de 1000 ya que solo se corrige la compresión

            avgDensError += (pred_density - density0) / density0;
        }
    }

    avgDensError /= totalNumParticles;

    LOG("Avg error         -> ", avgDensError * 100.0, "% | ", maxError * 100.0, "%");
}

void DFSPHSolver::correctDensityError()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int nFluidModels = sim -> numberFluidModels();
    int boundaryMethod =sim -> getBoundaryHandlingMethod();
    Real supportRadius = sim -> getSupportRadius();
    Real ts = sim -> getTimeStep();
    Real ts2 = ts * ts;

    densityWarmStart();

    predictDensities();

    iterations = 0;
    avgDensError = 0;
    while (((avgDensError > maxError) || (iterations < minIterations)) && iterations < maxIterations)
    {
        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();
            Real density0 = fm -> getRefDensity();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                Vector3r & ri = fm -> getPosition(i);
                Vector3r & vi = fm -> getVelocity(i);
                Real dens_i = fm -> getDensity(i);
                Real predDens_i = getPredDensity(fmIndex, i);
                Real alpha_i = getAlpha(fmIndex, i);

                Vector3r sum(0.0, 0.0, 0.0);

                Real ki = (predDens_i - density0) * alpha_i / ts2;

                // Warm start
                k[fmIndex][i] += (predDens_i - density0) * alpha_i;
                
                Vector3i cellId = floor(ri / supportRadius);

                forall_fluid_neighbors_in_same_phase
                (
                    Vector3r & rj = fm -> getPosition(j);
                    Real dens_j = fm -> getDensity(j);
                    Real predDens_j = getPredDensity(fmIndex, j);
                    Real alpha_j = getAlpha(fmIndex, j);

                    Real kj = (predDens_j - density0) * alpha_j / ts2;

                    Real kSum = (ki / dens_i + kj / dens_j);

                    if (fabs(kSum) > eps) // ver si sirve de algo (probablemente solo para ahorrar tiempo)
                        sum += fm -> getMass(j) * kSum * CubicSpline::gradW(ri - rj);
                );

                if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
                {
                    forall_boundary_neighbors
                    (
                        PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                        Vector3r & rb = bm -> getPosition(b);

                        Real kSum = ki / dens_i;

                        if (fabs(kSum > eps))
                            sum += bm -> getMass() * kSum * CubicSpline::gradW(ri - rb);
                    );
                } 
                else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
                {
                    forall_boundary_neighbors
                    (
                        AkinciBoundaryModel *bm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                        Vector3r & rb = bm -> getPosition(b);

                        Real kSum = ki / dens_i; // densidad al cuadrado? (como en el paper)

                        if (fabs(kSum) > eps)
                            sum += density0 * bm -> getVolume(b) * kSum * CubicSpline::gradW(ri - rb);
                    );
                }
                            
                vi -= ts * sum;
            } 

            predictDensities();
        }

        ++iterations;
    }
    sumIterations += iterations;

    LOG("Solver its        -> ", iterations);
    LOG("Solver avg its    -> ", sumIterations / (Real) steps);
}

void DFSPHSolver::computeDivergenceError()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int nFluidModels = sim -> numberFluidModels();
    int boundaryMethod = sim -> getBoundaryHandlingMethod();
    Real supportRadius = sim -> getSupportRadius();
    unsigned int totalNumParticles = 0;

    avgDivError = 0;
    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    {
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumActiveParticles();
        Real density0 = fm -> getRefDensity();
        totalNumParticles += numParticles;

        #pragma omp parallel for reduction(+:avgDivError)
        for (unsigned int i = 0; i < numParticles; ++i)
        {
            Vector3r & ri = fm -> getPosition(i);
            Vector3r & vi = fm -> getVelocity(i);
            Real & divError = getDivError(fmIndex, i);

            divError = 0;
            unsigned int numNeighbors = 0;

            Vector3i cellId = floor(ri / supportRadius);

            forall_fluid_neighbors_in_same_phase
            (
                Vector3r & rj = fm -> getPosition(j);
                Vector3r & vj = fm -> getVelocity(j);

                Vector3r grad = CubicSpline::gradW(ri - rj);
                divError += fm -> getMass(j) * dot(vi - vj, grad);

                if (length(grad) > 0.0)
                    numNeighbors ++;
            );

            if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
            {
                forall_boundary_neighbors
                (
                    PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                    Vector3r & rb = bm -> getPosition(b);

                    Vector3r grad = CubicSpline::gradW(ri - rb);
                    divError += bm -> getMass() * dot(vi, grad);

                    if (length(grad) > 0.0)
                        numNeighbors ++;
                );
            } 
            else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
            {
                forall_boundary_neighbors
                (
                    AkinciBoundaryModel *bm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                    Vector3r & rb = bm -> getPosition(b);

                    Vector3r grad = CubicSpline::gradW(ri - rb);
                    divError += density0 * bm -> getVolume(b) * dot(vi, grad);

                    if (length(grad) > 0.0)
                        numNeighbors ++;
                );
            } 

            if (numNeighbors < 20) // No corregir la divergencia de las particulas que tienen deficiencia de vecinos (Importante contar solo como vecinos las particulas que estan a una distancia menor o igual que h ya que la tabla hash devuelve como vecinos mas particulas a parte de estas)
                divError = 0.0;
            else
                divError = glm::max(divError, static_cast<Real>(0.0));  // solo se toma el error de divergencia positivo ya que es el error de compresion al igual que en la densidad

            avgDivError += divError / density0;
        }
    }

    avgDivError /= totalNumParticles;

    LOG("Avg errorV        -> ", avgDivError * sim -> getTimeStep() * 100.0, "% | ", maxErrorV * 100.0, "%");
}

void DFSPHSolver::correctDivergenceError()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int nFluidModels = sim -> numberFluidModels();
    int boundaryMethod = sim -> getBoundaryHandlingMethod();
    Real supportRadius = sim -> getSupportRadius();
    Real ts = sim -> getTimeStep();

    //divergenceWarmStart(); // De momento inestabiliza la simulacion

    computeDivergenceError();

    iterationsV = 0;
    // en el codigo de dfsph se divide el maxErrorV entre el timeStep
    // puede ser porque el avgDivError es el drho/dt o sea la variacion de la densidad respecto al tiempo
    // pero por unidad de tiempo, sin embargo el maxErrorV es el error entre pasos de simulacion es decir ya esta implicitamente
    // multiplicado por el dt, entonces se tiene que dividir por el dt
    // en este caso no se divide porque ya se multiplica cada div error que se acumula en avgDivError por el dt
    avgDivError = 0;
    while (((avgDivError > maxErrorV / ts) || (iterationsV < minIterationsV)) && iterationsV < maxIterations)
    {
        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();
            Real density0 = fm -> getRefDensity();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                Vector3r & ri = fm -> getPosition(i);
                Vector3r & vi = fm -> getVelocity(i);
                Real dens_i = fm -> getDensity(i);
                Real divError_i = getDivError(fmIndex, i);
                Real alpha_i = getAlpha(fmIndex, i);

                Vector3r sum(0.0, 0.0, 0.0);

                Real kiv = divError_i * alpha_i / ts;

                // Warm start
                kv[fmIndex][i] += divError_i * alpha_i;

                Vector3i cellId = floor(ri / supportRadius);

                forall_fluid_neighbors_in_same_phase
                (
                    Vector3r & rj = fm -> getPosition(j);
                    Real dens_j = fm -> getDensity(j);
                    Real divError_j = getDivError(fmIndex, j);
                    Real alpha_j = getAlpha(fmIndex, j);

                    Real kjv = divError_j * alpha_j / ts;
                    Real kSumV = (kiv / dens_i + kjv / dens_j);

                    if (fabs(kSumV) > eps)
                        sum += fm -> getMass(j) * kSumV * CubicSpline::gradW(ri - rj);
                );

                if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
                {
                    forall_boundary_neighbors
                    (
                        PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                        Vector3r & rb = bm -> getPosition(b);

                        Real kSum = kiv / dens_i;

                        if (fabs(kSum) > eps)
                            sum += bm -> getMass() * kSum * CubicSpline::gradW(ri - rb);
                    );
                } 
                else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
                {
                    forall_boundary_neighbors
                    (
                        AkinciBoundaryModel *bm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                        Vector3r & rb = bm -> getPosition(b);

                        Real kSum = kiv / dens_i; // densidad al cuadrado? (como en el paper)

                        if (fabs(kSum) > eps)
                            sum += density0 * bm -> getVolume(b) * kSum * CubicSpline::gradW(ri - rb);
                    );
                } 

                vi -= ts * sum;
            } 

            computeDivergenceError();
        }

        ++iterationsV;
    }
    sumIterationsV += iterationsV;

    LOG("SolverV its       -> ", iterationsV);
    LOG("SolverV avg its   -> ", sumIterationsV / (Real) steps);
}

void DFSPHSolver::densityWarmStart()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int nFluidModels = sim -> numberFluidModels();
    int boundaryMethod =sim -> getBoundaryHandlingMethod();
    Real supportRadius = sim -> getSupportRadius();
    Real ts = sim -> getTimeStep();
    Real ts2 = ts * ts;

    ////// WARM START ////// 
    if (sim -> getTimeStep() > 0.0)
    {
        //predictDensities();

        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();
            Real density0 = fm -> getRefDensity();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                //if (getPredDensity(fmIndex, i) > density0)
                //{
                    Real dens_i = fm -> getDensity(i);
                    //k[fmIndex][i] *= 1.0 / (ts2 * dens_i);
                    //k[fmIndex][i] = - glm::max(-k[fmIndex][i] / dens_i, -0.000000001) / ts2; // asi vibra poco y el paso de tiempo aumenta respecto a no usar warm start y respecto a utilizar la linea de arriba
                    k[fmIndex][i] = - glm::max(-k[fmIndex][i] / dens_i, -0.0000000005) / ts2; // asi vibran poco pero el paso de tiempo no aumenta tanto como con la linea de arriba
                //}
                //else
                //{
                //    k[fmIndex][i] = 0.0;
                //}
            }
        }

        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();
            Real density0 = fm -> getRefDensity();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                //if (getPredDensity(fmIndex, i) > density0)
                {
                    Vector3r & ri = fm -> getPosition(i);
                    Vector3r & vi = fm -> getVelocity(i);

                    Vector3r sum(0.0, 0.0, 0.0);

                    Real ki = k[fmIndex][i];

                    Vector3i cellId = floor(ri / supportRadius);


                    forall_fluid_neighbors_in_same_phase
                    (
                        Vector3r & rj = fm -> getPosition(j);

                        Real kj = k[fmIndex][j];

                        Real kSum = (ki + kj);

                        if (fabs(kSum) > eps) // ver si sirve de algo (probablemente solo para ahorrar tiempo)
                            sum += fm -> getMass(j) * kSum * CubicSpline::gradW(ri - rj);;
                    );

                    if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
                    {
                        forall_boundary_neighbors
                        (
                            PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                            Vector3r & rb = bm -> getPosition(b);

                            Real kSum = ki;

                            if (fabs(kSum > eps))
                                sum += bm -> getMass() * kSum * CubicSpline::gradW(ri - rb);
                        );
                    } 
                    else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
                    {
                        forall_boundary_neighbors
                        (
                            AkinciBoundaryModel *bm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                            Vector3r & rb = bm -> getPosition(b);

                            Real kSum = 2.0 * ki; 

                            if (fabs(kSum) > eps)
                                sum += density0 * bm -> getVolume(b) * kSum * CubicSpline::gradW(ri - rb);
                        );
                    }

                    vi -= sum * ts;
                }
            }
        }

        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                k[fmIndex][i] = 0.0;
            }
        }
    }

    /*if (sim -> getTimeStep() > 0.0)
    {
        //predictDensities();

        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();
            Real density0 = fm -> getRefDensity();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                if (getPredDensity(fmIndex, i) > density0)
                {
                    Real dens_i = fm -> getDensity(i);
                    //k[fmIndex][i] *= 1.0 / (ts2 * dens_i);
                    //k[fmIndex][i] = - glm::max(-k[fmIndex][i] / dens_i, -0.000000001) / ts2; // asi vibra poco y el paso de tiempo aumenta respecto a no usar warm start y respecto a utilizar la linea de arriba
                    k[fmIndex][i] = - glm::max(-k[fmIndex][i] / dens_i, -0.0000000005) / ts2; // asi vibran poco pero el paso de tiempo no aumenta tanto como con la linea de arriba
                }
                else
                {
                    k[fmIndex][i] = 0.0;
                }
            }
        }

        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();
            Real density0 = fm -> getRefDensity();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                //if (getPredDensity(fmIndex, i) > density0)
                {
                    Vector3r & ri = fm -> getPosition(i);
                    Vector3r & vi = fm -> getVelocity(i);

                    Vector3r sum(0.0, 0.0, 0.0);

                    Real ki = k[fmIndex][i];

                    Vector3i cellId = floor(ri / supportRadius);


                    forall_fluid_neighbors_in_same_phase
                    (
                        Vector3r & rj = fm -> getPosition(j);

                        Real kj = k[fmIndex][j];

                        Real kSum = (ki + kj);

                        if (fabs(kSum) > eps) // ver si sirve de algo (probablemente solo para ahorrar tiempo)
                            sum += fm -> getMass(j) * kSum * CubicSpline::gradW(ri - rj);;
                    );

                    if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
                    {
                        forall_boundary_neighbors
                        (
                            PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                            Vector3r & rb = bm -> getPosition(b);

                            Real kSum = ki;

                            if (fabs(kSum > eps))
                                sum += bm -> getMass() * kSum * CubicSpline::gradW(ri - rb);
                        );
                    } 
                    else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
                    {
                        forall_boundary_neighbors
                        (
                            AkinciBoundaryModel *bm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                            Vector3r & rb = bm -> getPosition(b);

                            Real kSum = 2.0 * ki; 

                            if (fabs(kSum) > eps)
                                sum += density0 * bm -> getVolume(b) * kSum * CubicSpline::gradW(ri - rb);
                        );
                    }

                    vi -= sum * ts;
                }
            }
        }

        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                k[fmIndex][i] = 0.0;
            }
        }
    }*/
}

void DFSPHSolver::divergenceWarmStart()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int nFluidModels = sim -> numberFluidModels();
    int boundaryMethod = sim -> getBoundaryHandlingMethod();
    Real supportRadius = sim -> getSupportRadius();
    Real ts = sim -> getTimeStep();

    /////// WARM START /////////
    if (sim -> getTimeStep() > 0.0)
    {
        computeDivergenceError();

        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                if (getDivError(fmIndex, i) > 0.0)
                {
                    Real dens_i = fm -> getDensity(i);
                    kv[fmIndex][i] *= 1.0 / (ts * dens_i);
                }
                else
                    kv[fmIndex][i] = 0.0;
            }
        }

        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();
            Real density0 = fm -> getRefDensity();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                //if (getDivError(fmIndex, i) > 0.0) 
                {
                    Vector3r & ri = fm -> getPosition(i);
                    Vector3r & vi = fm -> getVelocity(i);

                    Vector3r sum(0.0, 0.0, 0.0);

                    Real kiv = kv[fmIndex][i];

                    Vector3i cellId = floor(ri / supportRadius);

                    forall_fluid_neighbors_in_same_phase
                    (
                        Vector3r & rj = fm -> getPosition(j);

                        Real kjv = kv[fmIndex][j];
                        Real kSumV = (kiv + kjv);

                        if (fabs(kSumV) > eps)
                            sum += fm -> getMass(j) * kSumV * CubicSpline::gradW(ri - rj);
                    );

                    if (boundaryMethod == Simulation::PCISPH_BOUNDARY_METHOD)
                    {
                        forall_boundary_neighbors
                        (
                            PCISPHBoundaryModel *bm = static_cast<PCISPHBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                            Vector3r & rb = bm -> getPosition(b);

                            Real kSum = kiv;

                            if (fabs(kSum) > eps)
                                sum += bm -> getMass() * kSum * CubicSpline::gradW(ri - rb);
                        );
                    } 
                    else if (boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
                    {
                        forall_boundary_neighbors
                        (
                            AkinciBoundaryModel *bm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                            Vector3r & rb = bm -> getPosition(b);

                            Real kSum = kiv; 

                            if (fabs(kSum) > eps)
                                sum += density0 * bm -> getVolume(b) * kSum * CubicSpline::gradW(ri - rb);
                        );
                    }

                    vi -= ts * sum;
                }
            } 
        }

        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();

            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                kv[fmIndex][i] = 0.0;
            }
        }
    }
}

void DFSPHSolver::resizeData()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();

    alpha.resize(nFluidModels);
    divError.resize(nFluidModels);
    predDensity.resize(nFluidModels);    

    k.resize(nFluidModels);
    kv.resize(nFluidModels);

    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    {
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumParticles();

        alpha[fmIndex].resize(numParticles);
        divError[fmIndex].resize(numParticles);
        predDensity[fmIndex].resize(numParticles);    
        k[fmIndex].resize(numParticles);
        kv[fmIndex].resize(numParticles);
    }
}


void DFSPHSolver::updateTimeStep()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();
    Real ts = sim -> getTimeStep();

    maxVel = 0.1;
    //maxAcc = 0.00;
    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    {
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumActiveParticles();

        #pragma omp parallel for reduction(max: maxVel)
        for (unsigned int i = 0; i < numParticles; ++i)
        {
            Vector3r & vi = fm -> getVelocity(i);
            Vector3r & ai = fm -> getAcceleration(i);
            Vector3r acc = ai + sim -> getGravity();

            Vector3r current_vi = vi + acc * ts;

            maxVel = glm::max(maxVel, sqrt(length(current_vi)));
            //maxAcc = glm::max(maxAcc, length(acc));
        }
    }

    Real newTs;
    Real cflTs = 1.0 * 0.4 * 2.0 * sim -> getParticleRadius() / maxVel; 

    if (cflTs < ts)
        newTs = cflTs;
    else if (iterations > 3)
        newTs = ts * 0.99;
    else 
        newTs = ts * 1.001;
        
    newTs = glm::max(newTs, sim -> getMinTimeStep()); 
    newTs = glm::min(newTs, sim -> getMaxTimeStep()); 

    sim -> setTimeStep(newTs);
    
    LOG("Current time step -> ", newTs, " s (", cflTs, " s)");
}

