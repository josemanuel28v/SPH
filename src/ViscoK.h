#ifndef _VISCOK_H_
#define _VISCOK_H_

#include "types.h"

class ViscoK
{
    private:

        static Real supportRadius;
        static Real supportRadius2;
        static Real supportRadius3;
        static Real supportRadius6;
        //static Real valueCoef;
        //static Real gradCoef;
        static Real laplCoef;

    public:

        static void setSupportRadius(Real sr)
        {
            supportRadius = sr;
            supportRadius2 = sr * sr;
            supportRadius3 = sr * sr * sr ;
            supportRadius6 = sr * sr * sr * sr * sr * sr;
            //valueCoef = 15 / (2 * M_PI * supportRadius3);
            //gradCoef = 15 / (2 * M_PI * supportRadius3);
            laplCoef = 45 / (M_PI * supportRadius6);
        }

        static Real laplW(const Vector3r & r)
        {
            Real mod_r = length(r);

            if (mod_r <= supportRadius)
                return laplCoef * (supportRadius - mod_r);

            return 0.0;
        }

        static Real getSupportRadius() { return supportRadius; }
};

#endif
 
