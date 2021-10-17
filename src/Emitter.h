#ifndef _EMITTER_H_
#define _EMITTER_H_

#include "types.h"
//#include "FluidModel.h"

class FluidModel;

class Emitter
{
    protected:

        unsigned int type; // Circulo o cuadrado
        unsigned int numParticles; // Numero de particulas que se van a emitir

        Vector3r r;
        Real v; 
        Matrix4r rot;

        Real startTime; // Tiempo en el que se empieza a emitir
        Real nextTimeEmit; // Siguiente tiempo en el que se emite

        Real width; // Si el tipo es circular el radio sera width
        Real height;

        FluidModel *fm;

    public:

        static unsigned int SQUARE_EMITTER;
        static unsigned int CIRCLE_EMITTER;

        Emitter(FluidModel *fm, unsigned int type, unsigned int numParticles, Vector3r r, Real v, Matrix4r rot, Real startTime, Real w, Real h);

        void emitParticles(); // En vez de calcular la velocidad a partir de una frecuencia hacer al reves calcular el siguiente tiempo de emision a partir de la velocidad y el tamaño de las particulas
        void emitCircle();
        void emitSquare();

        void setStartTime(Real startTime) { this -> startTime = startTime; }
        void setRotation(Matrix4r rot) { this -> rot = rot; }
        void setPosition(Vector3r r) { this -> r = r; }
        void setSize(Real w, Real h) { width = w; height = h; }


};

#endif 
