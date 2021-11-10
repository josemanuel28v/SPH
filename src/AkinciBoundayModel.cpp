#include "AkinciBoundaryModel.h" 
#include "Simulation.h"
#include "PCISPHSolver.h"
#include "DFSPHSolver.h"
#include "Poly6.h"
#include "CubicSpline.h"
#include "../extern/OBJLoader.h"
#include "../extern/RegularTriangleSampling.h"
#include <glm/gtc/matrix_transform.hpp>

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

void AkinciBoundaryModel::addCube(Vector3r min, Vector3r max)
{
    unsigned int currentBoundaryParticles = size();
    std::vector<Vector3r> points;

    sampleCube(min, max, points);

    resizeBoundary(currentBoundaryParticles + points.size());

    for (unsigned int i = 0; i < points.size(); ++i)
    {
        unsigned int id = currentBoundaryParticles + i;

        r[id] = points[i];
        v[id] = Vector3r(0.0, 0.0, 0.0);
    }
}

void AkinciBoundaryModel::addSphere(Vector3r pos, Real radius)
{
    unsigned int currentBoundaryParticles = size();
    std::vector<Vector3r> points;

    sampleSphere(pos, radius, points);

    resizeBoundary(currentBoundaryParticles + points.size());

    for (unsigned int i = 0; i < points.size(); ++i)
    {
        unsigned int id = currentBoundaryParticles + i;

        r[id] = points[i];
        v[id] = Vector3r(0.0, 0.0, 0.0);
    }
}

void AkinciBoundaryModel::addGeometry(std::string path, Real particleRadius)
{
    unsigned int currentBoundaryParticles = size();
    std::vector<Vector3r> points;

    Vector3r scale = Vector3r(1.0, 1.0, 1.0);
    Vector3r translate = Vector3r(0.0, 0.0, 0.0);
    Vector3r rotate = Vector3r(M_PI * 0.5, 0.0, 0.0);

    sampleGeometry(path, particleRadius, scale, translate, rotate, points);

    resizeBoundary(currentBoundaryParticles + points.size());

    for (unsigned int i = 0; i < points.size(); ++i)
    {
        unsigned int id = currentBoundaryParticles + i;

        r[id] = points[i];
        v[id] = Vector3r(0.0, 0.0, 0.0);
    }
}

void AkinciBoundaryModel::sampleCube(Vector3r min, Vector3r max, std::vector<Vector3r> & points)
{
    Simulation *sim = Simulation::getCurrent();

    Real supportRadius = sim -> getSupportRadius();

    // Longitud del cubo en cada eje
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
                    points.push_back(position);

                    count++;
                }
}

void AkinciBoundaryModel::sampleSphere(Vector3r origen, Real radius, std::vector<Vector3r> & points)
{
    Simulation *sim = Simulation::getCurrent();

    Real separacion = 1.5; // 1 sin solaparse, cuanto mayor sea, mas se solaparan
    unsigned num_parts = floor(4.0 * M_PI * radius * radius / (M_PI * sim -> getParticleRadius() * sim -> getParticleRadius()) * separacion); // Superficie de la boundary sphere / superficie de un circulo con el radio de la particula del fluido

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

            points.push_back(radius * Vector3r(x, y, z) + origen);
        }
    }
}

void AkinciBoundaryModel::sampleGeometry(std::string path, Real maxDistance, Vector3r scale, Vector3r translate, Vector3r rotate, std::vector<Vector3r> & points)
{
    // Regular triangle sampling
	Utilities::OBJLoader::Vec3f scale_ = {(float) scale.x, (float) scale.y, (float) scale.z};
	std::vector<Utilities::OBJLoader::Vec3f> x;
	std::vector<Utilities::MeshFaceIndices> f;
	std::vector<Utilities::OBJLoader::Vec3f> n;
	std::vector<Utilities::OBJLoader::Vec2f> tc;

	Utilities::OBJLoader::loadObj(path, &x, &f, &n, &tc, scale_);

	// Cambio de los tipos de loader (Vec3f) a los tipos del sampler (Vector3r de eigen)
	std::vector<SPH::Vector3r> x_(x.size());
	for (unsigned int i = 0; i < x.size(); ++i)
	{
		x_[i][0] = x[i][0];
		x_[i][1] = x[i][1];
		x_[i][2] = x[i][2];
	}

	std::vector<unsigned int> f_(f.size() * 3);
	for (unsigned int i = 0; i < f.size(); ++i)
	{
		f_[3 * i] = f[i].posIndices[0] - 1;
		f_[3 * i + 1] = f[i].posIndices[1] - 1;
		f_[3 * i + 2] = f[i].posIndices[2] - 1;
	}

	std::vector<SPH::Vector3r> samples;

	SPH::RegularTriangleSampling::sampleMesh(x.size(), &x_[0], f.size(), &f_[0], maxDistance, samples);

    points.resize(samples.size());

    Vector3r axisX(1.0, 0.0, 0.0);
    Vector3r axisY(0.0, 1.0, 0.0);
    Vector3r axisZ(0.0, 0.0, 1.0);
	for (unsigned int i = 0; i < points.size(); ++i)
	{
		points[i].x = samples[i][0];
		points[i].y = samples[i][1];
		points[i].z = samples[i][2];

        Vector4r tmp(points[i].x, points[i].y, points[i].z, 1.0);

        Matrix4r rotM = glm::rotate(Matrix4r(1.0), rotate.x, axisX);
        rotM = glm::rotate(rotM, rotate.y, axisY);
        rotM = glm::rotate(rotM, rotate.z, axisZ);
        
        tmp = rotM * tmp;

        points[i] = Vector3r(tmp.x, tmp.y, tmp.z);

        points[i] += translate;
	}
}


void AkinciBoundaryModel::resizeBoundary(const unsigned int size)
{
    BoundaryModel::resizeBoundary(size);
    v.resize(size);
    volume.resize(size);
}
 
