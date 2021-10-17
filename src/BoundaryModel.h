#ifndef _BOUNDARY_MODEL_H_
#define _BOUNDARY_MODEL_H_

#include <vector>
#include "types.h"

class BoundaryModel
{
    protected:
        // aqui deberia ir el rigid body que tiene los vertices, de momento sera un array de positiones
        std::vector<Vector3r> r;

    public:

        BoundaryModel() {}
        virtual ~BoundaryModel() {}

        virtual void init(std::vector<Vector3r> &);

        virtual void resizeBoundary(const unsigned int);
        virtual void cleanBoundary();
        unsigned int getNumParticles() { return r.size(); }

        Vector3r & getPosition(const unsigned int i) { return r[i]; }
        std::vector<Vector3r>& getPosition() { return r; }
};

#endif
 
