#ifndef _POLY6_H
#define _POLY6_H_

#include "types.h"

class Poly6
{

    private:

        static Real supportRadius;
        static Real supportRadius2;
        static Real supportRadius9;
        static Real valueCoef;
        static Real gradCoef;

    public:

        Poly6() {}

        static void setSupportRadius(Real sr)
        {
            supportRadius = sr;
            supportRadius2 = sr * sr;
            supportRadius9 = sr * sr * sr * sr * sr * sr * sr * sr * sr;
            valueCoef = 315.0 / (64.0 * M_PI * supportRadius9);
            gradCoef = 945.0 / (32.0 * M_PI * supportRadius9);
        }

        static Real W(const Vector3r & r)
        {
            Real mod_r = length(r);
            Real mod_r2 = mod_r * mod_r;

            if (mod_r2 <= supportRadius2)
            {
                Real diff = supportRadius2 - mod_r2;

                return valueCoef * diff * diff * diff;
            }

            return 0.0;
        }

        static Vector3r gradW(const Vector3r & r)
        {
            Real mod_r = length(r);
            Real mod_r2 = mod_r * mod_r;

            if (mod_r2 <= supportRadius2)
            {
                Real diff = supportRadius2 - mod_r2;
                return - diff * diff * gradCoef * r;
            }

            return Vector3r(0.0);
        }

        static Real getSupportRadius() { return supportRadius; }
};


#endif
