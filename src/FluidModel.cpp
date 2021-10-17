#include "FluidModel.h"
#include "NonPressureForce.h"
#include "Simulation.h"
#include "StandardViscosityForce.h"
#include "ArtificialViscosity.h"
#include "XSPHViscosity.h"
#include "SurfaceTensionAkinci.h"
#include "AdhesionForce.h"
#include <iostream>

FluidModel::FluidModel()
{
    //inicializar variables y punteros que no se usen al principio a nullptr
    density_0 = 1000;
}

FluidModel::~FluidModel()
{
    // borrar punteros

    cleanFluid();
}

void FluidModel::setMasses(const Real mass)
{
    for (unsigned int i = 0; i < getNumParticles(); ++i)
        setMass(i, mass);
}

void FluidModel::init(std::vector<Vector3r> & fluidPoints, std::vector<Vector3r> & fluidVelocities)
{
    resizeFluid(fluidPoints.size());

    for (unsigned int i = 0; i < fluidPoints.size(); ++i)
    {
        r[i] = fluidPoints[i];
        v[i] = fluidVelocities[i];
    }

    // establecer el número aqui de momento hasta que se añadan los emisores 
    numParticles = fluidPoints.size();
    numActiveParticles = numParticles;
}

void FluidModel::resizeFluid(const unsigned int size)
{
    mass.resize(size);
    density.resize(size);
    pressure.resize(size);
    r.resize(size);
    v.resize(size);
    a.resize(size);
    n.resize(size);
}

void FluidModel::cleanFluid()
{
    mass.clear();
    density.clear();
    pressure.clear();
    r.clear();
    v.clear();
    a.clear();
    n.clear();
}

void FluidModel::setViscosityForce(Real viscosity, Real bViscosity)
{
    Simulation *sim = Simulation::getCurrent();
    const int current_method = sim -> getViscosityMethod();

    if (current_method == Simulation::STANDARD_VISCOSITY_METHOD)
    {
        StandardViscosityForce *svf = new StandardViscosityForce(this);
        svf -> init(viscosity, bViscosity);
        npForces.push_back(svf);
    }
    else if (current_method == Simulation::ARTIFICIAL_VISCOSITY_METHOD)
    {
        ArtificialViscosity *av = new ArtificialViscosity(this);
        av -> init(viscosity, bViscosity);
        npForces.push_back(av);
    }
    else if (current_method == Simulation::XSPH_VISCOSITY_METHOD)
    {
        XSPHViscosity *xv = new XSPHViscosity(this);
        xv -> init(viscosity, bViscosity);
        npForces.push_back(xv);
    }
}

void FluidModel::setSurfaceTensionForce(Real stCoef)
{
    Simulation *sim = Simulation::getCurrent();
    const int current_method = sim -> getSurfaceTensionMethod();

    if (current_method == Simulation::AKINCI_SURFACE_TENSION_METHOD)
    {
        SurfaceTensionAkinci *sta = new SurfaceTensionAkinci(this);
        sta -> init(stCoef);
        npForces.push_back(sta);
    }
}

void FluidModel::setAdhesionForce(Real beta)
{
    Simulation *sim = Simulation::getCurrent();
    const int current_method = sim -> getAdhesionMethod();

    if (current_method == Simulation::AKINCI_ADHESION_METHOD)
    {
        AdhesionForce *ad = new AdhesionForce(this);
        ad -> init(beta);
        npForces.push_back(ad);
    }
}

void FluidModel::addEmitter(unsigned int type, unsigned int numParticles, Vector3r r, Real v, Matrix4r rot, Real startTime, Real w, Real h)
{
    Simulation *sim = Simulation::getCurrent();
    SPHSolver *solver = sim -> getSolver();
    HashTable *grid = sim -> getGrid();

    Emitter emitter(this, type, numParticles, r, v, rot, startTime, w, h);

    // Cada vez que se añade un emisor se modifica el numero total de particulas de fluidModel y se redimensiona el solver y las nonpresureforces
    this -> numParticles += numParticles;
    resizeFluid(this -> numParticles);

    solver -> resizeData();

    grid -> reserve(2 * this -> numParticles);

    setMasses(density_0 * 4.0 / 3.0 * M_PI * pow(sim -> getParticleRadius(), 3.0));

    for (unsigned int i = 0; i < npForces.size(); ++i)
        npForces[i] -> resize(this -> numParticles);

    emitters.push_back(emitter);
}

void FluidModel::emitParticles()
{
    for (unsigned int i = 0; i < emitters.size(); ++i)
        emitters[i].emitParticles();
}

std::vector<Vector3r> & FluidModel::getPositions()
{
    return r;
}

std::vector<Real> & FluidModel::getPressures()
{
    return pressure;
}

std::vector<Real> FluidModel::getVelocities()
{    
    std::vector<Real> v_mag(v.size());

    for (unsigned int i = 0; i < v.size(); ++i)
        v_mag[i] = length(v[i]);

    return v_mag;
}