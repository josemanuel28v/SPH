#ifndef _CUBE_BOUNDARY_MODEL_H_
#define _CUBE_BOUNDARY_MODEL_H_

#include "BoundaryModel.h"
#include "types.h"

class CubeBoundaryModel: public BoundaryModel
{
    protected:

        Vector3r *min;
        Vector3r *max;

        Real normalFct;
        Real tangentialFct;

    public:

        void init(std::vector<Vector3r> &);

        void correctPositionAndVelocity();

        void setNormalFct(Real normalFct) { this -> normalFct = normalFct; }
        void setTangentialFct(Real tangentialFct) { this -> tangentialFct = tangentialFct; }

        Real getNormalFct() { return normalFct; }
        Real getTangentialFct() { return tangentialFct; }


};

#endif 
