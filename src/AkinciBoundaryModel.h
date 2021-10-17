#ifndef _AKINCI_BOUNDARY_MODEL_H_
#define _AKINCI_BOUNDARY_MODEL_H_

#include "BoundaryModel.h"
#include "types.h"

class AkinciBoundaryModel: public BoundaryModel
{
    protected:

        std::vector<Vector3r> v;
        std::vector<Real> volume;

        Real mass;
        Real radius;
        
        Real normalFct;
        Real tangentialFct;

    public:

        void init(std::vector<Vector3r> &);

        void computeVolume();

        void setNormalFct(Real normalFct) { this -> normalFct = normalFct; }
        void setTangentialFct(Real tangentialFct) { this -> tangentialFct = tangentialFct; }
        void setRadius(Real radius) { this -> radius = radius; }
        void setMass(Real mass) { this -> mass = mass; }

        Real getNormalFct() { return normalFct; }
        Real getTangentialFct() { return tangentialFct; }
        Real getRadius() { return radius; }
        Real getMass() { return mass; }

        Real & getVolume(const unsigned int i) { return volume[i]; }

        void sampleCube();
        void sampleSphere();
        void sample();

        void resizeBoundary(const unsigned int);
        unsigned int size() { return r.size(); }
};

#endif  
 
