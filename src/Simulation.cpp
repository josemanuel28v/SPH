#include "Simulation.h" 
#include "WCSPHSolver.h"
#include "PCISPHSolver.h"
#include "DFSPHSolver.h"
#include "Poly6.h"
#include "Spiky.h"
#include "ViscoK.h"
#include "Cohesion.h"
#include "Adhesion.h"
#include "CubicSpline.h"
#include <iostream>

#include "../Extern/OBJLoader.h"
#include "../Extern/RegularTriangleSampling.h"
#include "AkinciBoundaryModel.h"
#include <glm/gtc/matrix_transform.hpp>

Simulation* Simulation::current = nullptr;

const int Simulation::WCSPH_METHOD = 0;
const int Simulation::PCISPH_METHOD = 1;
const int Simulation::DFSPH_METHOD = 2;

const int Simulation::CUBE_BOUNDARY_METHOD = 0;
const int Simulation::PCISPH_BOUNDARY_METHOD = 1;
const int Simulation::AKINCI_BOUNDARY_METHOD = 2;

const int Simulation::STANDARD_VISCOSITY_METHOD = 0;
const int Simulation::ARTIFICIAL_VISCOSITY_METHOD = 1;
const int Simulation::XSPH_VISCOSITY_METHOD = 2;
const int Simulation::AKINCI_SURFACE_TENSION_METHOD = 3;
const int Simulation::AKINCI_ADHESION_METHOD = 4;

const int Simulation::POLY6_KERNEL_METHOD = 0;
const int Simulation::SPIKY_KERNEL_METHOD = 1;

const int Simulation::POLY6_GRADIENT_METHOD = 0;
const int Simulation::SPIKY_GRADIENT_METHOD = 1;

const int Simulation::SPLINE_LAPLACIAN_METHOD = 0;

Simulation::Simulation()
{
    // inicializar  parametros por defecto
    current_method = -1;
    current_bh_method = CUBE_BOUNDARY_METHOD;
    current_visco_method = -1;

    current_kernel_method = -1;
    current_gradient_method = -1;
    current_laplacian_method = -1;

    solver = nullptr;
    grid = nullptr;

    grid = new HashTable();
}

Simulation::~Simulation()
{
    //eliminar punteros
    delete solver;
    delete grid;

    // Eliminar contenido de los punteros en fluidModels y limpiar el vector
    for (unsigned i = 0; i < fluidModels.size(); ++i)
        delete fluidModels[i];
    fluidModels.clear();

    // Eliminar contenido de los punteros en boundaryModels y limpiar el vector
    for (unsigned i = 0; i < boundaryModels.size(); ++i)
        delete boundaryModels[i];
    boundaryModels.clear();

    //current = nullptr; // delete current?
    delete current;
}

Simulation* Simulation::getCurrent()
{
    if (current == nullptr)
        current = new Simulation();

    return current;
}

void Simulation::setCurrent(Simulation* sim)
{
    current = sim;
}

bool Simulation::hasCurrent()
{
    return current != nullptr;
}

void Simulation::init()
{
    initKernels();
    initGrid();
    solver -> init();
}

void Simulation::initGrid()
{
    unsigned int nTotalParticles = 0;
    unsigned int nFluidModels = numberFluidModels();
    for (unsigned int i = 0; i < nFluidModels; ++i)
    {
        FluidModel *fm = getFluidModel(i);
        nTotalParticles += fm -> getNumParticles();
    }
    grid -> reserve(2 * nTotalParticles);
}

void Simulation::initKernels()
{
    Poly6::setSupportRadius(supportRadius);
    Spiky::setSupportRadius(supportRadius);
    ViscoK::setSupportRadius(supportRadius);
    Cohesion::setSupportRadius(supportRadius);
    Adhesion::setSupportRadius(supportRadius);
    CubicSpline::setSupportRadius(supportRadius);
}

bool Simulation::step()
{
    std::cout << "STEP" << " " << getTime() << std::endl;
    solver -> step();
    std::cout << std::endl;

    return tm.hasToSave();
}

void Simulation::run()
{
    while (getTime() < 5.0)
        step();
}

void Simulation::setParticleRadius(Real particleRadius)
{
    this -> particleRadius = particleRadius;

    // Actualizar valor de la distancia de suavizado
    this -> supportRadius = 4.0 * particleRadius /** pow(kernelParticles, 1./3.)*/;
    grid -> setSmoothingLength(this -> supportRadius);
}

void Simulation::setSupportRadius(Real supportRadius)
{
    this -> supportRadius = supportRadius;
    grid -> setSmoothingLength(supportRadius);
}

void Simulation::setSimulationMethod(const int method)
{
    if (method == current_method || method < 0)
        return;

    current_method = method;

    if (solver != nullptr)
        delete solver;

    switch (method)
    {
        case WCSPH_METHOD:
            solver = new WCSPHSolver();
            break;

        case PCISPH_METHOD:
            solver = new PCISPHSolver();
            break;

        case DFSPH_METHOD:
            solver = new DFSPHSolver();
            break;
    }
}

void Simulation::setBoundaryMethod(const int method)
{
    current_bh_method = method;
}

void Simulation::setKernelMethod(const int method)
{
    if (method == current_kernel_method || method < 0)
        return;

    current_kernel_method = method;

    switch (method)
    {
        case POLY6_KERNEL_METHOD:
            value = Poly6::W;  
            break;

        case SPIKY_KERNEL_METHOD:
            value = Spiky::W;
            break;
    }
}

void Simulation::setGradientMethod(const int method)
{
    if (method == current_gradient_method || method < 0)
        return;

    current_gradient_method = method;
    
    switch (method)
    {
        case POLY6_GRADIENT_METHOD:
            grad = Poly6::gradW;   
            break;

        case SPIKY_GRADIENT_METHOD:
            grad = Spiky::gradW;
            break;
    }  
}

void Simulation::setLaplacianMethod(const int method)
{
    if (method == current_laplacian_method || method < 0)
        return;

    current_laplacian_method = method;

    switch (method)
    {
        case SPLINE_LAPLACIAN_METHOD:
            lapl = ViscoK::laplW;   
            break;
    }
}

void Simulation::addFluidModel(std::vector<Vector3r>& fluidPoints, std::vector<Vector3r>& fluidVelocities)
{
    FluidModel *fm = new FluidModel();
    fm -> init(fluidPoints, fluidVelocities);
    fluidModels.push_back(fm);
}

void Simulation::addBoundaryModel(BoundaryModel *bm)
{
    boundaryModels.push_back(bm);
}

void Simulation::emitParticles()
{
    for (unsigned int i = 0; i < numberFluidModels(); ++i)
        fluidModels[i] -> emitParticles();
}

std::vector<Real> & Simulation::getDivergenceError(const unsigned int fmIndex)
{
    if (current -> getSimulationMethod() == DFSPH_METHOD)
    {
        DFSPHSolver *solver = static_cast<DFSPHSolver*>(this -> solver);

        return solver -> getDivergenceError(fmIndex);
    }
    else
    {
        FluidModel *fm = getFluidModel(fmIndex);

        return fm -> getPressures();
    }
}

/** 
 * Esta funcion utiliza el sampling de la biblioteca SPlisHSPlasH
 * cuando se implemente un metodo de sampleo de superficie se sustituira por el propio
 */
void Simulation::addBoundaryModelFromOBJ(std::string path, Real maxDistance, Vector3r scale, Vector3r translate, Vector3r rotate)
{
    if (current_bh_method != AKINCI_BOUNDARY_METHOD)
    {
        std::cout << "Boundary method must be AKINCI_BOUNDARY_METHOD " << std::endl;
        return; 
    }

    // Regular triangle sampling
	Utilities::OBJLoader::Vec3f scale_ = {scale.x, scale.y, scale.z};
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

	std::vector<Vector3r> points(samples.size());

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

	AkinciBoundaryModel *abm = new AkinciBoundaryModel();

	abm -> init(points);

	addBoundaryModel(abm);

	std::cout << "Boundary particles from obj " << points.size() << std::endl;
}

FluidModel* Simulation::buildFluidBlock(Real radius, std::vector<FluidBlockInfo> fluidBlockInfo)
{
    Real dist = 2.0 * radius;

    std::vector<Vector3r> position;
    std::vector<Vector3r> velocity;

    for (FluidBlockInfo block: fluidBlockInfo)
    {    
        block.origin.x -= dist * (block.nx - 0.5) * 0.5;
        block.origin.y -= dist * (block.ny - 0.5) * 0.5;
        block.origin.z -= dist * (block.nz - 0.5) * 0.5;

        for (unsigned i = 0; i < block.nx; ++i)
            for (unsigned j = 0; j < block.ny; ++j)
                for (unsigned k = 0; k < block.nz; ++k)
                {
                    Vector3r pos(i, j, k);
                    pos *= dist;
                    pos += block.origin;

                    position.push_back(pos);
                    velocity.push_back(Vector3r(0.0));
                }
    }

    addFluidModel(position, velocity);

    return getFluidModel(numberFluidModels() - 1);
}





