#ifndef _STANDARD_VISCOSITY_FORCE_H_
#define _STANDARD_VISCOSITY_FORCE_H_

#include "Simulation.h"
#include "NonPressureForce.h"
#include "AkinciBoundaryModel.h"
#include "ViscoK.h"
#include "types.h"

#include <iostream>

class StandardViscosityForce: public NonPressureForce
{
    protected:

        Real viscosity;
        Real boundaryViscosity;

    public:

        StandardViscosityForce(FluidModel *fm) :NonPressureForce(fm) {}

        // init?
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
            Real density0 = fm -> getRefDensity();

            // Calcular densidad para las particulas de fluido vecinas del mismo fluidModel
            #pragma omp parallel for num_threads(16)
            for (unsigned int i = 0; i < numParticles; ++i)
            {
                Vector3r& ri = fm -> getPosition(i);
                Vector3r& vi = fm -> getVelocity(i);
                Vector3r& ai = fm -> getAcceleration(i);

                Vector3r viscoForce(0.0);

                Vector3i cellId = floor(ri / sim -> getSupportRadius());

                // contribución de las partículas de fluido vecinas
                forall_fluid_neighbors_in_same_phase
                (
                    Vector3r & rj = fm -> getPosition(j);
                    Vector3r & vj = fm -> getVelocity(j);
                    viscoForce += fm -> getMass(j) * (vj - vi) / fm -> getDensity(j) * ViscoK::laplW(ri - rj);
                );

                Real kinematicVisco = viscosity / fm -> getRefDensity(); 
                viscoForce *= kinematicVisco;

                ai += viscoForce;

                viscoForce = Vector3r(0.0, 0.0, 0.0);

                if (boundaryViscosity != 0.0 && boundaryMethod == Simulation::AKINCI_BOUNDARY_METHOD)
                {
                    forall_boundary_neighbors
                    (
                        AkinciBoundaryModel *bm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                        Vector3r & rb = bm -> getPosition(b);

                        viscoForce += density0 * bm -> getVolume(b) * (- vi) / fm -> getDensity(i) * ViscoK::laplW(ri - rb);
                    );

                    Real kinematicVisco = boundaryViscosity / fm -> getRefDensity(); 
                    viscoForce *= kinematicVisco;

                    ai += viscoForce;
                }
                // casos de los distintos boudarymodels
            }
        }

        void resize(const unsigned int ) {}
};

#endif 

 
