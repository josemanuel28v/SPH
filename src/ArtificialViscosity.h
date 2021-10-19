#ifndef _ARTIFICIAL_VISCOSITY_H_
#define _ARTIFICIAL_VISCOSITY_H_

#include "Simulation.h"
#include "NonPressureForce.h"
#include "AkinciBoundaryModel.h"
#include "Spiky.h"
#include "CubicSpline.h"
#include "types.h"

#include <iostream>

class ArtificialViscosity: public NonPressureForce
{
    protected:

        Real viscosity;
        Real boundaryViscosity;

    public:

        ArtificialViscosity(FluidModel *fm) :NonPressureForce(fm) {}

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
            Real h2 = supportRadius * supportRadius;

            #pragma omp parallel for 
            for (unsigned int i = 0; i < numParticles; i++)
            {
                Vector3r& ri = fm -> getPosition(i);
                Vector3r& vi = fm -> getVelocity(i);
                Vector3r& ai = fm -> getAcceleration(i);
                Real dens_i = fm -> getDensity(i);

                Vector3i cellId = floor(ri / supportRadius);

                forall_fluid_neighbors_in_same_phase
                (
                    const Vector3r &rj = fm -> getPosition(j);
                    const Vector3r &vj = fm -> getVelocity(j);

                    const Vector3r rirj = ri - rj;
                    Real dotProd = dot(vi - vj, rirj);
				    const Real dens_j = fm -> getDensity(j);
                    const Real sqDist = length(rirj) * length(rirj);

                    //Real kij = 2 * density0 / (dens_i + dens_j);

                    ai += /*kij **/ 10.0 * 0.1 * viscosity * (fm -> getMass(j) / dens_j) * dotProd / (sqDist + 0.01 * h2) * CubicSpline::gradW(rirj);
                );

                if (boundaryViscosity != 0.0 && boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
                {
                    forall_boundary_neighbors
                    (
                        AkinciBoundaryModel *nbm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));
                        Vector3r & rb = nbm -> getPosition(b);

                        const Vector3r rirb = ri - rb;
                        Real dotProd = dot(vi, rirb);
                        const Real sqDist = length(rirb) * length(rirb);

                        ai += 10.0 * 0.1 * boundaryViscosity * (density0 * nbm -> getVolume(b) / dens_i) * dotProd / (sqDist + 0.01 * h2) * CubicSpline::gradW(rirb);
                    );
                }


            }
        }

        void resize(const unsigned int ) {}
};

#endif 
