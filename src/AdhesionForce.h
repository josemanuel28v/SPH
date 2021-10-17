 
#ifndef _ADHESION_FORCE_H_
#define _ADHESION_FORCE_H_

#include "Simulation.h"
#include "NonPressureForce.h"
#include "Adhesion.h"
#include "AkinciBoundaryModel.h"
#include "types.h"
#include <vector>

#include <iostream>

class AdhesionForce: public NonPressureForce
{
    protected:

        Real beta;


    public:

        AdhesionForce(FluidModel *fm) :NonPressureForce(fm) { resize(fm -> getNumParticles()); }

        void init(Real beta)
        {
            this -> beta = beta;
        }

        void step()
        {
            Simulation *sim = Simulation::getCurrent();
            HashTable *grid = sim -> getGrid();
            unsigned int numParticles = fm -> getNumActiveParticles();
            Real density0 = fm -> getRefDensity();
            Real supportRadius = sim -> getSupportRadius();

            #pragma omp parallel for 
            for (unsigned int i = 0; i < numParticles; i++)
            {
                Vector3r & ri = fm -> getPosition(i);
                Vector3r & ai = fm -> getAcceleration(i);

                Vector3i cellId = floor(ri / supportRadius);

                forall_boundary_neighbors
                (
                    AkinciBoundaryModel *bm = static_cast<AkinciBoundaryModel*>(sim -> getBoundaryModel(nbmIndex));

                    Vector3r & rb = bm -> getPosition(b);

                    Vector3r rib = ri - rb;
                    Real mod_rib = length(rib);

                    if (mod_rib > 0.0)
                    {
                        ai -= beta * density0 * bm -> getVolume(b) * Adhesion::W(rib) * rib / mod_rib;
                    }
                );

            }
        }

        void resize(const unsigned int size) { }
};

#endif 