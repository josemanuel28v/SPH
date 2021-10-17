#ifndef _PCISPH_SOLVER_H_
#define _PCISPH_SOLVER_H_

#include <vector>
#include "SPHSolver.h"
#include "types.h"

class PCISPHSolver: public SPHSolver
{
    protected:

        std::vector<std::vector<Vector3r>> predR;
        std::vector<std::vector<Vector3r>> predV;
        std::vector<std::vector<Vector3r>> pressureAcc;
        std::vector<Real> scalingFactor;

        Real maxDensityError;

    public:

        void init();
        void step();

        void computeScalingFactor();
        void pressureSolver();
        void initPressure();
        void predictVelocityAndPosition();
        void predictDensities();
        void predictFluidDensities(const unsigned int);
        void predictBoundaryDensities(const unsigned int);
        void updatePressure();
        void computePressureAcc();
        void integrate();

        void setPressureAcc(const unsigned int fmIndex, const unsigned int i, Vector3r & a) { pressureAcc[fmIndex][i] = a; }
        void setPredR(const unsigned int fmIndex, const unsigned int i, Vector3r & r) { predR[fmIndex][i] = r; }
        void setPredV(const unsigned int fmIndex, const unsigned int i, Vector3r & v) { predV[fmIndex][i] = v; }
        void setScalingFactor(const unsigned int fmIndex, Real sf) { scalingFactor[fmIndex] = sf; }

        Vector3r & getPressureAcc(const unsigned int fmIndex, const unsigned int i) { return pressureAcc[fmIndex][i]; }
        Vector3r & getPredR(const unsigned int fmIndex, const unsigned int i) { return predR[fmIndex][i]; }
        Vector3r & getPredV(const unsigned int fmIndex, const unsigned int i) { return predV[fmIndex][i]; }
        Real & getScalingFactor(const unsigned int fmIndex) { return scalingFactor[fmIndex]; }

        std::vector<std::vector<Vector3r>>* getPredR() { return &predR; }
        std::vector<std::vector<Vector3r>>* getPredV() { return &predV; }

        void resizeData();
};

#endif  
