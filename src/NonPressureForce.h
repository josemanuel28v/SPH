#ifndef _NON_PRESSURE_FORCES_H_
#define _NON_PRESSURE_FORCES_H_

#include "FluidModel.h"

class NonPressureForce
{
    protected:

        FluidModel *fm;

    public:

        NonPressureForce(FluidModel *fm) { this -> fm = fm; }

        virtual void step() = 0; 
        virtual void resize(const unsigned int) = 0;
};

#endif 

