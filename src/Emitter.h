#ifndef _EMITTER_H_
#define _EMITTER_H_

#include "types.h"
#include <vector>

class FluidModel;

class Emitter
{
    protected:

        unsigned int type;           // Circulo o cuadrado
        unsigned int numParticles;   // Numero de particulas totales que se van a emitir
        std::vector<Vector3r> group; // Grupo de particulas que se emite cada vez

        Vector3r r;
        Real v; 
        Quat4r rot;

        Real startTime;             // Tiempo en el que se empieza a emitir
        Real nextTimeEmit;          // Siguiente tiempo en el que se emite

        Real width;                 // Si el tipo es circular el radio sera width
        Real height;

        Real spacing;

        FluidModel *fm;

        bool pause;

    public:

        static unsigned int SQUARE_EMITTER;
        static unsigned int CIRCLE_EMITTER;

        Emitter(FluidModel *fm, unsigned int type, unsigned int numParticles, Vector3r r, Real v, Quat4r rot, Real startTime, Real w, Real h, Real s);

        void buildCircleGroup();
        void buildSquareGroup();

        void emitParticles(); 
        void emitCircle(Real);
        void emitSquare(Real);

        void setStartTime(Real startTime) { this -> startTime = startTime; }
        void setRotation(Vector3r rot) { this -> rot = rot; }
        void setPosition(Vector3r r) { this -> r = r; }
        void setSize(Real w, Real h) { width = w; height = h; }

        void setPause(bool pause);
};

#endif 
