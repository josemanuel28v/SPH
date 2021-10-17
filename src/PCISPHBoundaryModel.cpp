#include "PCISPHBoundaryModel.h" 
#include "Simulation.h"
#include "PCISPHSolver.h"
#include "DFSPHSolver.h"
#include <iostream>

void PCISPHBoundaryModel::init(std::vector<Vector3r> & points)
{
    BoundaryModel::init(points);

    normalFct = 0.0;
    tangentialFct = 1.0;

    sample();
}

void PCISPHBoundaryModel::sample()
{
    if (r.size() == 1)
        sampleSphere();
    else if (r.size() == 2)
        sampleCube();
}

void PCISPHBoundaryModel::sampleCube()
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

    // Face normals
    Vector3r leftNormal(1, 0, 0);
    Vector3r rightNormal(-1, 0, 0);
    Vector3r topNormal(0, -1, 0);
    Vector3r bottomNormal(0, 1, 0);
    Vector3r frontNormal(0, 0, -1);
    Vector3r behindNormal(0, 0, 1);
    
    // Edge normals
    Vector3r leftBottomNormal = glm::normalize((leftNormal + bottomNormal) * static_cast<Real>(0.5));
    Vector3r behindBottomNormal = glm::normalize((behindNormal + bottomNormal) * static_cast<Real>(0.5));
    Vector3r frontBottomNormal = glm::normalize((frontNormal + bottomNormal) * static_cast<Real>(0.5));
    Vector3r rightBottomNormal = glm::normalize((rightNormal + bottomNormal) * static_cast<Real>(0.5));

    Vector3r leftTopNormal = glm::normalize((leftNormal + topNormal) * static_cast<Real>(0.5));
    Vector3r behindTopNormal = glm::normalize((behindNormal + topNormal) * static_cast<Real>(0.5));
    Vector3r frontTopNormal = glm::normalize((frontNormal + topNormal) * static_cast<Real>(0.5));
    Vector3r rightTopNormal = glm::normalize((rightNormal + topNormal) * static_cast<Real>(0.5));

    Vector3r leftBehindNormal = glm::normalize((leftNormal + behindNormal) * static_cast<Real>(0.5));
    Vector3r rightBehindNormal = glm::normalize((rightNormal + behindNormal) * static_cast<Real>(0.5));
    Vector3r rightFrontNormal = glm::normalize((rightNormal + frontNormal) * static_cast<Real>(0.5));
    Vector3r leftFrontNormal = glm::normalize((leftNormal + frontNormal) * static_cast<Real>(0.5));

    // Corner normals
    Vector3r leftBehindBottomNormal = glm::normalize((leftNormal + behindNormal + bottomNormal) / static_cast<Real>(3.0));
    Vector3r rightBehindBottomNormal = glm::normalize((rightNormal + behindNormal + bottomNormal) / static_cast<Real>(3.0));
    Vector3r leftFrontBottomNormal = glm::normalize((leftNormal + frontNormal + bottomNormal) / static_cast<Real>(3.0));
    Vector3r rightFrontBottomNormal = glm::normalize((rightNormal + frontNormal + bottomNormal) / static_cast<Real>(3.0));

    Vector3r leftBehindTopNormal = glm::normalize((leftNormal + behindNormal + topNormal) / static_cast<Real>(3.0));
    Vector3r rightBehindTopNormal = glm::normalize((rightNormal + behindNormal + topNormal) / static_cast<Real>(3.0));
    Vector3r leftFrontTopNormal = glm::normalize((leftNormal + frontNormal + topNormal) / static_cast<Real>(3.0));
    Vector3r rightFrontTopNormal = glm::normalize((rightNormal + frontNormal + topNormal) / static_cast<Real>(3.0));

    vector<Vector3r> boundary_position;
    vector<Vector3r> boundary_normal;
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

                    // Corners
                    if (i == 0 && j == 0 && k == 0)
                        boundary_normal.push_back(leftBehindBottomNormal);
                    else if (i == nx && j == 0 && k == 0)
                        boundary_normal.push_back(rightBehindBottomNormal);
                    else if (i == 0 && j == 0 && k == nz)
                        boundary_normal.push_back(leftFrontBottomNormal);
                    else if (i == nx && j == 0 && k == nz)
                        boundary_normal.push_back(rightFrontBottomNormal);
                    else if (i == 0 && j == ny && k == 0)
                        boundary_normal.push_back(leftBehindTopNormal);
                    else if (i == nx && j == ny && k == 0)
                        boundary_normal.push_back(rightBehindTopNormal);
                    else if (i == 0 && j == ny && k == nz)
                        boundary_normal.push_back(leftFrontTopNormal);
                    else if (i == nx && j == ny && k == nz)
                        boundary_normal.push_back(rightFrontTopNormal);

                    // Edges
                    else if (i == 0 && j == 0 && k > 0 && k < nz)
                        boundary_normal.push_back(leftBottomNormal);
                    else if (i > 0 && i < nx && j == 0 && k == 0)
                        boundary_normal.push_back(behindBottomNormal);
                    else if (i > 0 && i < nx && j == 0 && k == nz)
                        boundary_normal.push_back(frontBottomNormal);
                    else if (i == nx && j == 0 && k > 0 && k < nz)
                        boundary_normal.push_back(rightBottomNormal);

                    else if (i == 0 && j == ny && k > 0 && k < nz)
                        boundary_normal.push_back(leftTopNormal);
                    else if (i > 0 && i < nx && j == ny && k == 0)
                        boundary_normal.push_back(behindTopNormal);
                    else if (i > 0 && i < nx && j == ny && k == nz)
                        boundary_normal.push_back(frontTopNormal);
                    else if (i == nx && j == ny && k > 0 && k < nz)
                        boundary_normal.push_back(rightTopNormal);

                    else if (i == 0 && j > 0 && j < ny && k == 0)
                        boundary_normal.push_back(leftBehindNormal);
                    else if (i == nx && j > 0 && j < ny && k == 0)
                        boundary_normal.push_back(rightBehindNormal);
                    else if (i == nx && j > 0 && j < ny && k == nz)
                        boundary_normal.push_back(rightFrontNormal);
                    else if (i == 0 && j > 0 && j < ny && k == nz)
                        boundary_normal.push_back(leftFrontNormal);

                    // Faces
                    else if (i == 0 && j > 0 && j < ny && k > 0 && k < nz)
                        boundary_normal.push_back(leftNormal);
                    else if (i == nx && j > 0 && j < ny && k > 0 && k < nz)
                        boundary_normal.push_back(rightNormal);
                    else if (i > 0 && i < nx && j == ny && k > 0 && k < nz)
                        boundary_normal.push_back(topNormal);
                    else if (i > 0 && i < nx && j == 0 && k > 0 && k < nz)
                        boundary_normal.push_back(bottomNormal);
                    else if (i > 0 && i < nx && j > 0 && j < ny && k == nz)
                        boundary_normal.push_back(frontNormal);

                    else if (i > 0 && i < nx && j > 0 && j < ny && k == 0)
                        boundary_normal.push_back(behindNormal);                    

                    count++;
                }

    resizeBoundary(boundary_position.size()); // Añadir las particulas fantasma;

    for (uint i = 0; i < boundary_position.size(); i++)
    {
        r[i] = boundary_position[i];
        n[i] = boundary_normal[i];
        v[i] = Vector3r(0, 0, 0);
        density[i] = 0;
        pressure[i] = 0;
    }

    std::cout << "Cube particles " << boundary_position.size() << std::endl;
}

void PCISPHBoundaryModel::sampleSphere()
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
        n[i] = glm::normalize(origen - boundary_position[i]);
        density[i] = 0;
        pressure[i] = 0;
    }
}

void PCISPHBoundaryModel::correctPositions()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int nFluidModels = sim -> numberFluidModels();
    Real supportRadius = sim -> getSupportRadius();
    Real restRadius = 0.5 * supportRadius;

    if (sim -> getSimulationMethod() == Simulation::PCISPH_METHOD)
    {
        PCISPHSolver *solver = static_cast<PCISPHSolver*>(sim -> getSolver());
        std::vector<std::vector<Vector3r>> *pred_r = solver -> getPredR();;

        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();

            //#pragma omp parallel for 
            for (uint i = 0; i < numParticles; i++)
            {       
                Vector3r & ri = fm -> getPosition(i); 
                Vector3r & ni = fm -> getNormal(i);
                Vector3r avg_normal(0, 0, 0);
                Real w_sum = 0;
                Real w2_sum = 0;

                Vector3i cellId = floor(ri / supportRadius);

                //#pragma omp critical asi parece que se soluciona pero ralentiza la paralelización y parece que va mas rapido sin paralelizar
                forsame_boundary_neighbors
                (
                    Vector3r & rb = r[b];
                    Vector3r predict_rib = (*pred_r)[fmIndex][i] - rb;

                    if (glm::length(ri - rb) < restRadius) // Colision
                    {
                        Real w = glm::max(static_cast<Real>(0.0), (restRadius - glm::length(predict_rib)) / restRadius);
                        w_sum += w;
                        w2_sum += w * (restRadius - length(predict_rib));
                        avg_normal += w * n[b];
                    } 
                );

                //#pragma omp critical

                ni = glm::normalize(avg_normal);

                //#pragma omp critical
                if (w_sum > 0)
                    ri += static_cast<Real>(1.0) / w_sum * w2_sum * ni;
            }
        }
    }
    else if (sim -> getSimulationMethod() == Simulation::WCSPH_METHOD || sim -> getSimulationMethod() == Simulation::DFSPH_METHOD)  
    {
        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();

            //#pragma omp parallel for 
            for (uint i = 0; i < numParticles; i++)
            {       
                Vector3r & ri = fm -> getPosition(i); 
                Vector3r & ni = fm -> getNormal(i);
                Vector3r avg_normal(0, 0, 0);
                Real w_sum = 0;
                Real w2_sum = 0;

                Vector3i cellId = floor(ri / supportRadius);

                //#pragma omp critical asi parece que se soluciona pero ralentiza la paralelización y parece que va mas rapido sin paralelizar
                forsame_boundary_neighbors
                (
                    Vector3r & rb = r[b];
                    Vector3r rib = ri - rb;

                    if (glm::length(rib) < restRadius) // Colision
                    {
                        Real w = glm::max(static_cast<Real>(0.0), (restRadius - length(rib)) / restRadius);
                        w_sum += w;
                        w2_sum += w * (restRadius - length(rib));
                        avg_normal += w * n[b];
                    } 
                );

                //#pragma omp critical

                ni = glm::normalize(avg_normal);

                //#pragma omp critical
                if (w_sum > 0)
                    ri += static_cast<Real>(1.0) / w_sum * w2_sum * ni;
            }
        }
    }  
}

void PCISPHBoundaryModel::correctVelocities()
{
    Simulation *sim = Simulation::getCurrent();
    unsigned int nFluidModels = sim -> numberFluidModels();    

    if (sim -> getSimulationMethod() == Simulation::PCISPH_METHOD) // con las velocidades predichas se corrige la velocidad
    {
        PCISPHSolver *solver = static_cast<PCISPHSolver*>(sim -> getSolver());
        std::vector<std::vector<Vector3r>> *pred_v = solver -> getPredV();

        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();

            #pragma omp parallel for
            for (uint i = 0; i < numParticles; i++)
            {       
                Vector3r & ni = fm -> getNormal(i);
                Vector3r & vi = fm -> getVelocity(i);

                if (glm::length(ni) > 0)
                {
                    Vector3r nVel = dot((*pred_v)[fmIndex][i], ni) * ni;
                    Vector3r tVel = (*pred_v)[fmIndex][i] - nVel;
                    vi = tangentialFct * tVel - normalFct * nVel;
                }
            }
        }
    }
    else if (sim -> getSimulationMethod() == Simulation::WCSPH_METHOD)  // con la velocidad se corrige la velocidad
    {
        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();

            #pragma omp parallel for
            for (uint i = 0; i < numParticles; i++)
            {       
                Vector3r & ni = fm -> getNormal(i);
                Vector3r & vi = fm -> getVelocity(i);

                if (glm::length(ni) > 0)
                {
                    Vector3r nVel = dot(vi, ni) * ni;
                    Vector3r tVel = vi - nVel;
                    vi = tangentialFct * tVel - normalFct * nVel;
                }
            }
        }
    }
    else if (sim -> getSimulationMethod() == Simulation::DFSPH_METHOD) // Con la velocidad predicha se corrige la velocidad predicha
    {
        for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
        {
            FluidModel *fm = sim -> getFluidModel(fmIndex);
            unsigned int numParticles = fm -> getNumActiveParticles();

            #pragma omp parallel for
            for (uint i = 0; i < numParticles; i++)
            {       
                Vector3r & ni = fm -> getNormal(i);
                Vector3r & vi = fm -> getVelocity(i);

                if (glm::length(ni) > 0)
                {
                    Vector3r nVel = dot(vi, ni) * ni;
                    Vector3r tVel = vi - nVel;
                    vi = tangentialFct * tVel - normalFct * nVel;
                }
            }
        }
    }
}

void PCISPHBoundaryModel::correctPredPositions()
{
    Simulation *sim = Simulation::getCurrent();
    HashTable *grid = sim -> getGrid();
    unsigned int nFluidModels = sim -> numberFluidModels();
    Real supportRadius = sim -> getSupportRadius();
    Real restRadius = 0.5 * supportRadius;

    PCISPHSolver *solver = static_cast<PCISPHSolver*>(sim -> getSolver());
    std::vector<std::vector<Vector3r>> *pred_r = solver -> getPredR();

    for (unsigned int fmIndex = 0; fmIndex < nFluidModels; ++fmIndex)
    {
        FluidModel *fm = sim -> getFluidModel(fmIndex);
        unsigned int numParticles = fm -> getNumActiveParticles();

        //#pragma omp parallel for num_threads(16) 
        for (uint i = 0; i < numParticles; i++)
        {       
            Vector3r & ri = fm -> getPosition(i); 
            Vector3r & ni = fm -> getNormal(i);
            Vector3r avg_normal(0, 0, 0);
            Real w_sum = 0;
            Real w2_sum = 0;

            Vector3i cellId = floor(ri / supportRadius);

            //#pragma omp critical asi parece que se soluciona pero ralentiza la paralelización y parece que va mas rapido sin paralelizar

            forsame_boundary_neighbors
            (
                Vector3r & rb = r[b];
                Vector3r & nb = n[b];
                Vector3r predict_rij = (*pred_r)[fmIndex][i] - rb;

                if (glm::length(ri - rb) < restRadius) // Colision
                {
                    Real w = glm::max(static_cast<Real>(0.0), (restRadius - glm::length(predict_rij)) / restRadius);
                    w_sum += w;
                    w2_sum += w * (restRadius - length(predict_rij));
                    avg_normal += w * nb;
                } 
            );

            //#pragma omp critical

            ni = glm::normalize(avg_normal);

            //#pragma omp critical
            
            if (w_sum > 0)
                (*pred_r)[fmIndex][i] += static_cast<Real>(1.0) / w_sum * w2_sum * ni;
        }
    }
}

void PCISPHBoundaryModel::resizeBoundary(const unsigned int size)
{
    BoundaryModel::resizeBoundary(size);
    v.resize(size);
    n.resize(size);
    density.resize(size);
    pressure.resize(size);
}
