#include "CubeBoundaryModel.h"
#include "Simulation.h"
#include <iostream>

void CubeBoundaryModel::init(std::vector<Vector3r> & points)
{
    BoundaryModel::init(points);

    min = &getPosition(0);
    max = &getPosition(1);

    normalFct = 0.0;
    tangentialFct = 1.0;
}

void CubeBoundaryModel::correctPositionAndVelocity()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();

    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    {
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumActiveParticles();

        #pragma omp parallel for
        for (unsigned i = 0; i < numParticles; i++)
        {
            Vector3r & pos = fm -> getPosition(i);
            Vector3r & vel = fm -> getVelocity(i);

            // Coordenada X
            if (pos.x < min -> x)
            {
                pos.x = min -> x;
                if (vel.x < 0.f)
                    vel.x *= - normalFct;

                vel.y *= tangentialFct;
                vel.z *= tangentialFct;
            }
            else if (pos.x > max -> x)
            {
                pos.x = max -> x;
                if (vel.x > 0.f)
                    vel.x *= - normalFct;

                vel.y *= tangentialFct;
                vel.z *= tangentialFct;
            }

            // Coordenada Y
            if (pos.y < min -> y)
            {
                pos.y = min -> y;
                if (vel.y < 0.f)
                    vel.y *= - normalFct;

                vel.x *= tangentialFct;
                vel.z *= tangentialFct;
            }
            else if (pos.y > max -> y)
            {
                pos.y = max -> y;
                if (vel.y > 0.f)
                    vel.y *= - normalFct;

                vel.x *= tangentialFct;
                vel.z *= tangentialFct;
            }

            // Coordenada Z
            if (pos.z < min -> z)
            {
                pos.z = min -> z;
                if (vel.z < 0.f)
                    vel.z *= - normalFct;

                vel.x *= tangentialFct;
                vel.y *= tangentialFct;
            }
            else if (pos.z > max -> z)
            {
                pos.z = max -> z;
                if (vel.z > 0.f)
                    vel.z *= - normalFct;

                vel.x *= tangentialFct;
                vel.y *= tangentialFct;
                
            }
        }
    }
} 
