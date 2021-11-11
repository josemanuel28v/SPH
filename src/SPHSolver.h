#ifndef _SPH_SOLVER_H_
#define _SPH_SOLVER_H_

#include "types.h"

class SPHSolver
{
    protected:

        unsigned int iterations;
        unsigned int minIterations;
        unsigned int maxIterations;

        unsigned int steps;
        unsigned int sumIterations;

        Real maxError;

        Real maxVel; // Ver si la mayoria de metodos lo necesitan, sino quitar de aqui y meter en la clase derivada que lo necesite        

    public:

        SPHSolver();
        virtual ~SPHSolver();

        // Calculo de densidad fijo para todos los solvers
        virtual void computeFluidDensities(const unsigned int);
        virtual void computeBoundaryDensities(const unsigned int);
        virtual void computeDensities();

        virtual void integrate();

        virtual void insertFluidParticles();
        virtual void insertBoundaryParticles();
        virtual void neighborhoodSearch();

        virtual void init() = 0;
        virtual void step() = 0;
        virtual void resizeData() = 0;

        void setMaxError(Real error) { maxError = error; }
        void setMinIterations(unsigned int minIt) { minIterations = minIt; }
        void setMaxIterations(unsigned int maxIt) { maxIterations = maxIt; }

        Real getMaxError() { return maxError; }
        Real getMaxVel() { return maxVel; }
        unsigned int getSteps() { return steps; }
        unsigned int getMinIterations() { return minIterations; }
        unsigned int getMaxIterations() { return maxIterations; }
};

#endif
