#ifndef _WCSPH_SOLVER_H_
#define _WCSPH_SOLVER_H_

#include <vector>
#include "SPHSolver.h"
#include "types.h"

class WCSPHSolver: public SPHSolver
{
    protected:

        Real stiffness;
        Real gamma;

    public:

        void init();
        void step();

        void computeFluidPressures(const unsigned int);
        void computeBoundaryPressures(const unsigned int);
        void computePressures();
        void computePressureForce(const unsigned int);

        void setStiffness(Real stiffness) { this -> stiffness = stiffness; }
        void setGamma(Real gamma) { this -> gamma = gamma; }

        Real getStiffness() { return stiffness; }
        Real getGamma() { return gamma; }

        void updateTimeStep();

        void resizeData() {}
};

#endif
