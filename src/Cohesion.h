#ifndef _COHESION_H_
#define _COHESION_H_

#include "types.h"

#include <iostream>

class Cohesion
{
    private:

        static Real supportRadius;
        static Real supportRadius6;
        static Real supportRadius9;
        static Real valueCoef1;
        static Real valueCoef2;

    public:

        static void setSupportRadius(Real sr)
        {
            supportRadius = sr;
            supportRadius6 = sr * sr * sr * sr * sr * sr;
            supportRadius9 = sr * sr * sr * sr * sr * sr * sr * sr * sr;
            valueCoef1 = 32.0 / (M_PI * supportRadius9);
            valueCoef2 = supportRadius6 / 64.0;
        }

        static Real W(const Vector3r & r)
        {
            Real value = 0;
            Real mod_r = glm::length(r);

            if ((2.0 * mod_r > supportRadius) && (mod_r <= supportRadius))
                value = valueCoef1 * pow(supportRadius - mod_r, 3.0) * pow(mod_r, 3.0);
            else if ((mod_r > 0) && (2.0 * mod_r <= supportRadius))
                value = valueCoef1 * 2.0 * pow(supportRadius - mod_r, 3.0) * pow(mod_r, 3.0) - valueCoef2;

            return value;
        }


        static Real getSupportRadius() { return supportRadius; }
};

#endif 
 
