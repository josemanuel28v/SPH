#ifndef _FLUID_MODEL_H_
#define _FLUID_MODEL_H_

#include <vector>
#include "types.h"
#include "Emitter.h"

class NonPressureForce;

class FluidModel
{
    protected:

        Real density_0;
        Real volume;
        std::vector<Real> mass;
        std::vector<Real> density;
        std::vector<Real> pressure;
        std::vector<Vector3r> r;
        std::vector<Vector3r> v;
        std::vector<Vector3r> a;
        std::vector<Vector3r> n; // Mover a PCISPHBoundaryModel

        unsigned int numParticles;
        unsigned int numActiveParticles;

        std::vector<NonPressureForce*> npForces;
        std::vector<Emitter> emitters;

    public:

        FluidModel();
        ~FluidModel();

        void init(std::vector<Vector3r> & fluidPoints, std::vector<Vector3r> & fluidVelocities);
        void resizeFluid(const unsigned int);
        void cleanFluid();

        void setMasses(const Real);
        void setRefDensity(Real density_0) { this -> density_0 = density_0; }
        void setVolume(Real volume) { this -> volume = volume; }
        void setMass(const unsigned int i, const Real mass) { this -> mass[i] = mass; }
        void setDensity(const unsigned int i, const Real density) { this -> density[i] = density; }
        void setPressure(const unsigned int i, const Real pressure) { this -> pressure[i] = pressure; }
        void setPosition(const unsigned int i, const Vector3r & r) { this -> r[i] = r; }
        void setVelocity(const unsigned int i, const Vector3r & v) { this -> v[i] = v; }
        void setAcceleration(const unsigned int i, const Vector3r & a) { this -> a[i] = a; }
        void setNormal(const unsigned int i, const Vector3r & n) { this -> n[i] = n; }

        Real& getRefDensity() { return density_0; }
        Real& getVolume() { return volume; }
        Real& getMass(const unsigned int i) { return mass[i]; }
        Real& getDensity(const unsigned int i) { return density[i]; }
        Real& getPressure(const unsigned int i) { return pressure[i]; }
        Vector3r& getPosition(const unsigned int i) { return r[i]; }
        Vector3r& getVelocity(const unsigned int i) { return v[i]; }
        Vector3r& getAcceleration(const unsigned int i) { return a[i]; }
        Vector3r& getNormal(const unsigned int i) { return n[i]; }

        unsigned int getNumParticles() { return mass.size(); }
        unsigned int getNumActiveParticles() { return numActiveParticles; }

        void setNumActiveParticles(const unsigned int n) { if (n <= getNumParticles()) { numActiveParticles = n; } }

        NonPressureForce* getNonPressureForce(const unsigned int i) { return npForces[i]; }
        unsigned int numberNonPressureForces() { return npForces.size(); }

        void setViscosityForce(Real, Real);
        void setSurfaceTensionForce(Real);
        void setAdhesionForce(Real);

        void addEmitter(unsigned int type, unsigned int numParticles, Vector3r r, Real v, Quat4r rot, Real startTime, Real w, Real h, Real s);
        void emitParticles();

        // Necesario para pintar colores en OpenGL
        std::vector<Vector3r> & getPositions();
        std::vector<Real> & getPressures();
        std::vector<Real>  getVelocities();

        // Pausar emisores desde la interfaz de OpenGL
        void setEmittersPause(bool pause) { for (auto & emitter: emitters) { emitter.setPause(pause); } }
};

#endif
