 
#ifndef _ADHESION_H_
#define _ADHESION_H_

#include "types.h"

#include <iostream>

class Adhesion
{
    private:

        static Real valueCoef;
        static Real supportRadius;

    public:

        static void setSupportRadius(Real sr)
        {
            supportRadius = sr;
            valueCoef = 0.007 / pow(sr, 3.25);
        }

        static Real W(const Vector3r & r)
        {
            Real mod_r = glm::length(r);

            if ((2.0 * mod_r > supportRadius) && (mod_r <= supportRadius))
            {
                Real value = - 4.0 * mod_r * mod_r / supportRadius + 6 * mod_r - 2 * supportRadius;
                return pow(value, 1.0 / 4.0);                
            }
            else
                return 0.0;
        }


        static Real getSupportRadius() { return supportRadius; }
};

#endif 
 
