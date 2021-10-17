#ifndef _XSPH_VISCOSITY_H_
#define _XSPH_VISCOSITY_H_

#include "Simulation.h"
#include "NonPressureForce.h"
#include "AkinciBoundaryModel.h"
#include "Spiky.h"
#include "CubicSpline.h"
#include "types.h"

#include <iostream>

class XSPHViscosity: public NonPressureForce
{
    protected:

        Real viscosity;
        Real boundaryViscosity;

    public:

        XSPHViscosity(FluidModel *fm) :NonPressureForce(fm) {}

        void init(Real viscosity, Real boundaryViscosity)
        {
            this -> viscosity = viscosity;
            this -> boundaryViscosity = boundaryViscosity;
        }

        void step()
        {
            Simulation *sim = Simulation::getCurrent();
            HashTable *grid = sim -> getGrid();
            unsigned int numParticles = fm -> getNumActiveParticles();
            int boundaryMethod = sim -> getBoundaryHandlingMethod();
            Real supportRadius = sim -> getSupportRadius();
            Real density0 = fm -> getRefDensity();

            #pragma omp parallel for 
            for (uint i = 0; i < numParticles; i++)
            {
                Vector3r& ri = fm -> getPosition(i);
                Vector3r& vi = fm -> getVelocity(i);
                Vector3r& ai = fm -> getAcceleration(i);

                Vector3i cellId = floor(ri / supportRadius);

                forall_fluid_neighbors_in_same_phase
                (
                    const Vector3r &rj = fm -> getPosition(j);
                    const Vector3r &vj = fm -> getVelocity(j);

                    ai +=  viscosity * (fm -> getMass(j) / fm -> getDensity(j)) * (vj - vi) * CubicSpline::W(ri - rj) / sim -> getTimeStep();
                );

                if (boundaryViscosity != 0.0 && boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
                {
                    forall_boundary_neighbors
                    (
                        AkinciBoundaryModel *nbm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));
                        Vector3r & rb = nbm -> getPosition(b);

                        ai +=  boundaryViscosity * (density0 * nbm -> getVolume(b) / fm -> getDensity(i)) * (- vi) * CubicSpline::W(ri - rb) / sim -> getTimeStep();
                    );
                }
            }
        }

        void resize(const unsigned int ) {}
};

#endif 
