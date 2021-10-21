#ifndef _DFSPH_SOLVER_H_
#define _DFSPH_SOLVER_H_

#include <vector>
#include "SPHSolver.h"
#include "types.h"

class DFSPHSolver: public SPHSolver
{
    protected:

        std::vector<std::vector<Real>> alpha;
        std::vector<std::vector<Real>> divError;
        std::vector<std::vector<Real>> predDensity;    

        std::vector<std::vector<Real>> k;    
        std::vector<std::vector<Real>> kv;    

        Real avgDensError;
        Real avgDivError;

        Real maxErrorV; // Maximo error de divergencia permitido

        unsigned int sumIterationsV;
        unsigned int minIterationsV;
        unsigned int iterationsV;

        Real maxAcc;
        
        Real eps;

    public:

        void init();
        void step();

        void computeAlpha();
        void predictVelocities();
        void updatePositions();
        void correctDensityError();
        void correctDivergenceError();

        void predictDensities();
        void computeDivergenceError();

        void divergenceWarmStart();
        void densityWarmStart();

        void setMaxErrorV(Real etaV) { maxErrorV = etaV; }

        Real & getAlpha(const unsigned int fmIndex, const unsigned int i) { return alpha[fmIndex][i]; }
        Real & getDivError(const unsigned int fmIndex, const unsigned int i) { return divError[fmIndex][i]; }
        Real & getPredDensity(const unsigned int fmIndex, const unsigned int i) { return predDensity[fmIndex][i]; }

        void updateTimeStep();

        void resizeData();

        // para pintar en opengl
        std::vector<Real> & getDivergenceError(const unsigned fmIndex) { return divError[fmIndex]; }
};

#endif  
