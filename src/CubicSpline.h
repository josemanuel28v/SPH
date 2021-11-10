#ifndef _CUBIC_SPLINE_H_
#define _CUBIC_SPLINE_H_

#include "types.h"

class CubicSpline
{
    protected:

        static Real supportRadius;
        static Real k;
        static Real l;

    public:

        static Real getSupportRadius() { return supportRadius; }

        static void setSupportRadius(Real val)
        {
            supportRadius = val;
            const Real pi = static_cast<Real>(M_PI);

            const Real h3 = supportRadius*supportRadius*supportRadius;
            k = static_cast<Real>(8.0) / (pi*h3);
            l = static_cast<Real>(48.0) / (pi*h3);
        }

        static Real W(const Real r)
        {
            Real res = 0.0;
            const Real q = r / supportRadius;
            if (q <= 1.0)
            {
                if (q <= 0.5)
                {
                    const Real q2 = q*q;
                    const Real q3 = q2*q;
                    res = k * (static_cast<Real>(6.0)*q3 - static_cast<Real>(6.0)*q2 + static_cast<Real>(1.0));
                }
                else
                {
                    res = k * (static_cast<Real>(2.0)*pow(static_cast<Real>(1.0) - q, static_cast<Real>(3.0)));
                }
            }
            return res;
        }

        static Real W(const Vector3r &r)
        {
            return W(length(r));
        }

        static Vector3r gradW(const Vector3r &r)
        {
            Vector3r res = Vector3r(0.0, 0.0, 0.0);
            const Real rl = length(r);
            const Real q = rl / supportRadius;
            if ((rl > 1.0e-5) && (q <= 1.0))
            {
                const Vector3r gradq = r * (static_cast<Real>(1.0) / (rl*supportRadius));
                if (q <= 0.5)
                {
                    res = l*q*((Real) 3.0*q - static_cast<Real>(2.0))*gradq;
                }
                else
                {
                    const Real factor = static_cast<Real>(1.0) - q;
                    res = l*(-factor*factor)*gradq;
                }
            }

            return res;
        }
};

#endif
