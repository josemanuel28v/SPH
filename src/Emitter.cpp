#include "Emitter.h"
#include "Simulation.h"

#include <iostream>
#include <glm/gtc/matrix_transform.hpp>

unsigned int Emitter::SQUARE_EMITTER = 0;
unsigned int Emitter::CIRCLE_EMITTER = 1;

Emitter::Emitter(FluidModel *fm, unsigned int type, unsigned int numParticles, Vector3r r, Real v, Quat4r rot, Real startTime, Real w, Real h, Real spacing)
{
    this -> fm = fm;
    this -> type = type;
    this -> numParticles = numParticles;
    this -> r = r;
    this -> v = v;
    this -> rot = rot;
    this -> startTime = startTime;
    this -> nextTimeEmit = startTime;
    this -> width = w;
    this -> height = h;
    this -> spacing = spacing;
    this -> pause = false;

    if (type == CIRCLE_EMITTER)
        buildCircleGroup();
    else if (type == SQUARE_EMITTER)
        buildSquareGroup();
}

void Emitter::buildCircleGroup()
{
    Simulation *sim = Simulation::getCurrent();
    Real particleRadius = sim -> getParticleRadius();
    Real dist = 2.0 * particleRadius;
    unsigned int numLevels = floor((width - 0.5 * dist) / dist);

    group.clear();
    group.push_back(Vector3r(0.0, 0.0, 0.0));
    for (unsigned int i = 0; i < numLevels; ++i)
    {
        Real radius = (i + 1) * dist; // Radio en el que se colocaran las particulas
        Real l = 2.0 * M_PI * radius; // Longitud de la circunferencia del radio calculado
        unsigned int nParts = floor (l / dist);

        Real inc = 2 * M_PI / nParts;
        for (unsigned int i = 0; i < nParts; ++i)
        {
            Real angle = i * inc;
            Vector3r position = Vector3r(radius * sin(angle), radius * cos(angle), 0.0);

            group.push_back(position);
        }
    }
}

void Emitter::buildSquareGroup()
{
    Simulation *sim = Simulation::getCurrent();
    Real particleRadius = sim -> getParticleRadius();
    Real dist = particleRadius * 2.0;

    unsigned int numX = floor(width / dist);
    unsigned int numY = floor(height / dist);

    Real realWidth = (numX - 1) * dist;
    Real realHeight = (numY - 1) * dist;

    group.clear();
    for (unsigned int i = 0; i < numX; ++i)
        for (unsigned int j = 0; j < numY; ++j)
        {
            Vector3r pos = Vector3r(realWidth * 0.5 - i * dist, realHeight * 0.5 - j * dist, 0.0);
            group.push_back(pos);
        }
}

// Para la interfaz de opengl
void Emitter::setPause(bool pause) { 
    Simulation *sim = Simulation::getCurrent();
    this -> pause = pause; 
    nextTimeEmit = sim -> getTime();
}

void Emitter::emitParticles()
{
    Simulation *sim = Simulation::getCurrent();
    Real time = sim -> getTime();

    if (time >= nextTimeEmit && !pause)
    {
        if (type == SQUARE_EMITTER)
            emitSquare(time - nextTimeEmit);
        else if (type == CIRCLE_EMITTER)
            emitCircle(time - nextTimeEmit);
    }
}

void Emitter::emitCircle(Real timeOffset)
{
    Simulation *sim = Simulation::getCurrent();
    Real particleRadius = sim -> getParticleRadius();
    Real dist = 2.0 * particleRadius;
    unsigned int numActiveParticles = fm -> getNumActiveParticles();

    unsigned int id = numActiveParticles;
    if (numActiveParticles + group.size() <= numActiveParticles + numParticles) 
    {
        for (unsigned int i = 0; i < group.size(); ++i)
        {
            Vector3r & position = fm -> getPosition(id);
            Vector3r & velocity = fm -> getVelocity(id);

            velocity = Vector3r(0, 0, v);
            position = group[i] + timeOffset * (velocity + sim -> getGravity() * timeOffset);

            position = rot * Vector4r(position, 1.0);
            velocity = rot * Vector4r(velocity, 1.0);

            position += r;

            ++id;
        }

        fm -> setNumActiveParticles(numActiveParticles + group.size());
        numParticles -= group.size();

        nextTimeEmit += spacing * 1.105 * dist / v;
    }
}

void Emitter::emitSquare(Real timeOffset)
{
    Simulation *sim = Simulation::getCurrent();
    Real particleRadius = sim -> getParticleRadius();
    Real dist = particleRadius * 2.0;
    unsigned int numActiveParticles = fm -> getNumActiveParticles();

    unsigned int id = numActiveParticles;
    if (numActiveParticles + group.size() <= numActiveParticles + numParticles) 
    {
        for (unsigned int i = 0; i < group.size(); ++i)
        {
            Vector3r & position = fm -> getPosition(id);
            Vector3r & velocity = fm -> getVelocity(id);

            velocity = Vector3r(0, 0, v);
            position = group[i] + timeOffset * (velocity + sim -> getGravity() * timeOffset);

            position = rot * Vector4r(position, 1.0);
            velocity = rot * Vector4r(velocity, 1.0);

            position += r;

            ++id;
        }

        fm -> setNumActiveParticles(numActiveParticles + group.size());
        numParticles -= group.size();

        nextTimeEmit += spacing * 1.105 * dist / v;
    }
}