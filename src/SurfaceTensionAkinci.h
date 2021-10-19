#ifndef _SURFACE_TENSION_AKINCI_H_
#define _SURFACE_TENSION_AKINCI_H_

#include "Simulation.h"
#include "NonPressureForce.h"
#include "CubicSpline.h"
#include "Cohesion.h"
#include "types.h"
#include <vector>

#include <iostream>

class SurfaceTensionAkinci: public NonPressureForce
{
    protected:

        Real stCoef;

        std::vector<Vector3r> normal;

    public:

        SurfaceTensionAkinci(FluidModel *fm) :NonPressureForce(fm) { resize(fm -> getNumParticles()); }

        // init?
        void init(Real stCoef)
        {
            this -> stCoef = stCoef;
        }

        void step()
        {
            Simulation *sim = Simulation::getCurrent();
            HashTable *grid = sim -> getGrid();
            unsigned int numParticles = fm -> getNumActiveParticles();
            Real density0 = fm -> getRefDensity();
            Real supportRadius = sim -> getSupportRadius();

            // Calcular normales
            #pragma omp parallel for
            for (unsigned int i = 0; i < numParticles; i++)
            {
                Vector3r & ri = fm -> getPosition(i);
                Vector3r & ni = getNormal(i);

                ni = Vector3r(0.0, 0.0, 0.0);

                Vector3r cellId = floor(ri / supportRadius);
                
                forall_fluid_neighbors_in_same_phase
                (
                    Vector3r & rj = fm -> getPosition(j);
                    Real vol = fm -> getMass(j) / fm -> getDensity(j);
                    ni += vol * CubicSpline::gradW(ri - rj);
                );

                ni *= supportRadius;
            }

            // Compute curvature and cohesion force
            #pragma omp parallel for 
            for (unsigned int i = 0; i < numParticles; i++)
            {
                Vector3r & ri = fm -> getPosition(i);
                Vector3r & ai = fm -> getAcceleration(i);
                Vector3r & ni = getNormal(i);
                Real dens_i = fm -> getDensity(i);
                
                Vector3r cellId = floor(ri / supportRadius);

                Vector3r stAcc(0.0, 0.0, 0.0);
                
                forall_fluid_neighbors_in_same_phase
                (
                    Vector3r cohesion(0.0, 0.0, 0.0);
                    Vector3r curvature(0.0, 0.0, 0.0);

                    Vector3r & rj = fm -> getPosition(j);
                    Vector3r & nj = getNormal(j);
                    Real dens_j = fm -> getDensity(j);

                    Vector3r rij = ri - rj;
                    Real rijMag = glm::length(rij);
                    Real kij = 2.0 * density0 / (dens_i + dens_j);

                    if (rijMag > 0.0 && rijMag <= supportRadius)
                    {
                        cohesion  = - fm -> getMass(j) * Cohesion::W(rij) * rij / rijMag;
                        curvature = - 1.0 * (ni - nj); // Se ha multiplicado por 2.0 porque da mejores resultados

                        ai += stCoef * kij * (cohesion + curvature);
                    }
                );

                //ai += stAcc;
            }
        }

        void resize(const unsigned int size) { normal.resize(size); }

        Vector3r & getNormal(const unsigned int i) { return normal[i]; }
};

#endif 