#include "AkinciBoundaryModel.h" 
#include "Simulation.h"
#include "PCISPHSolver.h"
#include "DFSPHSolver.h"
#include "Poly6.h"
#include "CubicSpline.h"
#include <iostream>

void AkinciBoundaryModel::init(std::vector<Vector3r> & points)
{
    BoundaryModel::init(points);

    normalFct = 0.0;
    tangentialFct = 1.0;

    sample();
}

void AkinciBoundaryModel::computeVolume()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    Real supportRadius = sim -> getSupportRadius();

    /**
     * Si se van a solapar distintos boundary models ya sean dinamicos
     * animados o estaticos se deberan recorrer los boundary neighbors
     * de todos los boundary models no solo los propios para calcular un 
     * volumen correcto en los solapamientos
     */
    #pragma omp parallel for
    for (unsigned int i = 0; i < r.size(); ++i)
    {
        Vector3r & ri = getPosition(i);
        Real & volume = getVolume(i);
        Vector3i cellId = floor (ri / supportRadius);

        Real delta = 0.0;

        forsame_boundary_neighbors
        (
            Vector3r & rb = getPosition(b);
            delta += CubicSpline::W(ri - rb);
        );

        volume = 1.0 / delta;
    }
}

void AkinciBoundaryModel::sample()
{
    if (r.size() == 1)
        sampleSphere();
    else if (r.size() == 2)
        sampleCube();
}

void AkinciBoundaryModel::sampleCube()
{
    Simulation *sim = Simulation::getCurrent();

    // Longitud del cubo en cada eje
    Vector3r min = getPosition(0);
    Vector3r max = getPosition(1);
    Real supportRadius = sim -> getSupportRadius();

    Real lx = glm::abs(min.x - max.x);
    Real ly = glm::abs(min.y - max.y); 
    Real lz = glm::abs(min.z - max.z); 
    
    unsigned nx = glm::ceil(lx / (0.5 * supportRadius));
    unsigned ny = glm::ceil(ly / (0.5 * supportRadius));
    unsigned nz = glm::ceil(lz / (0.5 * supportRadius));
    
    // Distancia entre particulas de cada eje
    Real dx = lx / nx; 
    Real dy = ly / ny; 
    Real dz = lz / nz;

    vector<Vector3r> boundary_position;
    unsigned count = 0;
    for (uint i = 0; i <= nx; i++)
        for (uint j = 0; j <= ny; j++)
            for (uint k = 0; k <= nz; k++)
                if ((i == 0 || i == nx) ||
                    (j == 0 || j == ny) ||
                    (k == 0 || k == nz))
                {
                    Vector3r position(i * dx, j * dy, k * dz);
                    position += min;
                    boundary_position.push_back(position);

                    count++;
                }

    resizeBoundary(boundary_position.size()); // Añadir las particulas fantasma;

    for (uint i = 0; i < boundary_position.size(); i++)
    {
        r[i] = boundary_position[i];
        v[i] = Vector3r(0, 0, 0);
    }
}

void AkinciBoundaryModel::sampleSphere()
{
    Simulation *sim = Simulation::getCurrent();

    Vector3r origen = getPosition(0);
    Real radius = getRadius(); // mitad de dimensions en blender

    Real separacion = 1.5; // 1 indica sin solaparse, cuanto mayor sea, mas se solaparan
    unsigned num_parts = floor(4.0 * M_PI * radius * radius / (M_PI * sim -> getParticleRadius() * sim -> getParticleRadius()) * separacion); // Superficie de la boundary sphere / superficie de un circulo con el radio de la particula del fluido

    vector<Vector3r> boundary_position;

    // Equidistant sphere
    Real a = 4.0 * M_PI  / num_parts;
    Real d = sqrt(a);
    int m_theta = int(round(M_PI / d));
    Real d_theta = M_PI / m_theta;
    Real d_phi = a / d_theta;

    for (int m = 0; m < m_theta; ++m)
    {
        Real theta = M_PI * (m + 0.5) / m_theta;
        int m_phi = int(round(2.0 * M_PI * sin(theta) / d_phi));
        for (int n = 0; n < m_phi; ++n)
        {
            Real phi = 2.0 * M_PI * n / m_phi;
            Real x = sin(theta) * cos(phi);
            Real y = sin(theta) * sin(phi);
            Real z = cos(theta);

            boundary_position.push_back(radius * Vector3r(x, y, z) + origen);
        }
    }

    resizeBoundary(boundary_position.size()); // Añadir las particulas fantasma;

    for (uint i = 0; i < boundary_position.size(); i++)
    {
        r[i] = boundary_position[i];
        v[i] = Vector3r(0, 0, 0);
    }
}

void AkinciBoundaryModel::resizeBoundary(const unsigned int size)
{
    BoundaryModel::resizeBoundary(size);
    v.resize(size);
    volume.resize(size);
}
 
