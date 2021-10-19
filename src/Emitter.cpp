#include "Emitter.h"
#include "Simulation.h"

#include <iostream>

unsigned int Emitter::SQUARE_EMITTER = 0;
unsigned int Emitter::CIRCLE_EMITTER = 1;

Emitter::Emitter(FluidModel *fm, unsigned int type, unsigned int numParticles, Vector3r r, Real v, Matrix4r rot, Real startTime, Real w, Real h)
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
}

void Emitter::emitParticles()
{
    Simulation *sim = Simulation::getCurrent();
    Real time = sim -> getTime();

    if (time >= nextTimeEmit)
    {
        if (type == SQUARE_EMITTER)
            emitSquare();
        else if (type == CIRCLE_EMITTER)
            emitCircle();
    }
}

void Emitter::emitCircle()
{
    Simulation *sim = Simulation::getCurrent();
    Real particleRadius = sim -> getParticleRadius();
    Real dist = particleRadius * pow(4.0 * M_PI / 3.0, 1.0 / 3.0);
    unsigned int numActiveParticles = fm -> getNumActiveParticles();

    unsigned int numLevels = floor((width - 0.5 * dist) / dist);

    std::vector<Vector3r> pos;
    pos.push_back(Vector3r(0.0, 0.0, 0.0));
    for (unsigned int i = 0; i < numLevels; ++i)
    {
        Real radius = (i + 1) * dist; // Radio en el que se colocaran las particulas
        Real l = 2.0 * M_PI * radius; // Longitud de la circunferencia con el radio anterior
        unsigned int nParts = floor (l / dist);

        Real inc = 2 * M_PI / nParts;
        for (unsigned int i = 0; i < nParts; ++i)
        {
            Real angle = i * inc;
            Vector3r position = Vector3r(radius * sin(angle), radius * cos(angle), 0.0);

            pos.push_back(position);
        }
    }

    unsigned int id = numActiveParticles;

    if (numActiveParticles + pos.size() <= numActiveParticles + numParticles) // Se podría calcular al principio el numero de particulas que se emite cada vez y establecer un numero real del total de particulas que sea multiplo del grupo asi no habria que controlar esto
    {
        for (unsigned int i = 0; i < pos.size(); ++i)
        {
            Vector3r & position = fm -> getPosition(id);
            Vector3r & velocity = fm -> getVelocity(id);

            position = pos[i];
            velocity = Vector3r(0, 0, v);

            position = rot * Vector4r(position, 1.0);
            velocity = rot * Vector4r(velocity, 1.0);

            position += r;

            ++id;
        }

        fm -> setNumActiveParticles(numActiveParticles + pos.size());
        numParticles -= pos.size();

        nextTimeEmit += 1.1 * dist / v;
        // 0.9 abajo
    }
}

void Emitter::emitSquare() // Ver qué se puede dejar calculado en la inicializacion
{
    Simulation *sim = Simulation::getCurrent();
    Real particleRadius = sim -> getParticleRadius();
    Real dist = particleRadius * 2.0 /*pow(4.0 * M_PI / 3.0, 1.0 / 3.0)*/;
    unsigned int numActiveParticles = fm -> getNumActiveParticles();

    unsigned int numX = floor(width / dist);
    unsigned int numY = floor(height / dist);

    Real realWidth = (numX - 1) * dist;
    Real realHeight = (numY - 1) * dist;

    unsigned int id = numActiveParticles;


    if (numActiveParticles + numX * numY <= numActiveParticles + numParticles) // Se podría calcular al principio el numero de particulas que se emite cada vez y establecer un numero real del total de particulas que sea multiplo del grupo asi no habria que controlar esto
    {
        for (unsigned int i = 0; i < numX; ++i)
            for (unsigned int j = 0; j < numY; ++j)
            {
                Vector3r & position = fm -> getPosition(id);
                Vector3r & velocity = fm -> getVelocity(id);

                position = Vector3r(realWidth * 0.5 - i * dist, realHeight * 0.5 - j * dist, 0.0);
                velocity = Vector3r(0, 0, v);

                position = rot * Vector4r(position, 1.0);
                velocity = rot * Vector4r(velocity, 1.0);

                position += r;

                ++id;
            }

        fm -> setNumActiveParticles(numActiveParticles + numX * numY);
        numParticles -= numX * numY;

        nextTimeEmit += 1.1 * dist / v;
    }
}