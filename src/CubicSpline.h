#ifndef _CUBIC_SPLINE_H_
#define _CUBIC_SPLINE_H_

#include "types.h"

class CubicSpline
{
    protected:

        static Real m_radius;
        static Real m_k;
        static Real m_l;

    public:

        static Real getSupportRadius() { return m_radius; }

        static void setSupportRadius(Real val)
        {
            m_radius = val;
            const Real pi = static_cast<Real>(M_PI);

            const Real h3 = m_radius*m_radius*m_radius;
            m_k = static_cast<Real>(8.0) / (pi*h3);
            m_l = static_cast<Real>(48.0) / (pi*h3);
        }

        static Real W(const Real r)
        {
            Real res = 0.0;
            const Real q = r / m_radius;
            if (q <= 1.0)
            {
                if (q <= 0.5)
                {
                    const Real q2 = q*q;
                    const Real q3 = q2*q;
                    res = m_k * (static_cast<Real>(6.0)*q3 - static_cast<Real>(6.0)*q2 + static_cast<Real>(1.0));
                    // 8/(pi*h^3) * 6*x^3 - 6*x^2 + 1
                }
                else
                {
                    res = m_k * (static_cast<Real>(2.0)*pow(static_cast<Real>(1.0) - q, static_cast<Real>(3.0)));
                    // 8/(pi*h^3) * 2*(1 - x) ^ 3
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
            const Real q = rl / m_radius;
            if ((rl > 1.0e-5) && (q <= 1.0))
            {
                const Vector3r gradq = r * (static_cast<Real>(1.0) / (rl*m_radius));
                if (q <= 0.5)
                {
                    res = m_l*q*((Real) 3.0*q - static_cast<Real>(2.0))*gradq;
                }
                else
                {
                    const Real factor = static_cast<Real>(1.0) - q;
                    res = m_l*(-factor*factor)*gradq;
                }
            }

            return res;
        }
};

#endif
