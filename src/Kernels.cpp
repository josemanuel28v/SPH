#include "Poly6.h"
#include "Spiky.h"
#include "ViscoK.h"
#include "Cohesion.h"
#include "CubicSpline.h"
#include "Adhesion.h"
#include "types.h"

Real Poly6::supportRadius;
Real Poly6::supportRadius2;
Real Poly6::supportRadius9;
Real Poly6::valueCoef;
Real Poly6::gradCoef;

Real Spiky::supportRadius;
Real Spiky::supportRadius2;
Real Spiky::supportRadius6;
Real Spiky::valueCoef;
Real Spiky::gradCoef;

Real ViscoK::supportRadius;
Real ViscoK::supportRadius2;
Real ViscoK::supportRadius3;
Real ViscoK::supportRadius6;
Real ViscoK::laplCoef;

Real Cohesion::supportRadius;
Real Cohesion::supportRadius6;
Real Cohesion::supportRadius9;
Real Cohesion::valueCoef1;
Real Cohesion::valueCoef2;

Real CubicSpline::supportRadius;
Real CubicSpline::k;
Real CubicSpline::l;

Real Adhesion::supportRadius;
Real Adhesion::valueCoef;

