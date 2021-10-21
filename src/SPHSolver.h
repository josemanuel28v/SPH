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

        // Suponiendo que todos los metodos calculan la densidad de la misma forma
        virtual void computeFluidDensities(const unsigned int);
        virtual void computeBoundaryDensities(const unsigned int);
        virtual void computeDensities();

        virtual void integrate();

        // Suponiendo que todos los metodos obtienen los vecinos de la misma forma
        // De momento el solver inserta siempre en la tabla hash las particulas 
        // de todos los fluidModels y fluidBoundaries
        // Estos metodos cambiaran si se modifica el metodo de busqueda de vecinos
        virtual void insertFluidParticles();
        virtual void insertBoundaryParticles();
        virtual void neighborhoodSearch();

        virtual void init() = 0;
        virtual void step() = 0;
        virtual void resizeData() = 0;

        void setMaxError(Real error) { maxError = error; }

        Real getMaxError() { return maxError; }
        Real getMaxVel() { return maxVel; }
        unsigned int getSteps() { return steps; }
};

#endif
