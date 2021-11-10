#ifndef _BOUNDARY_MODEL_H_
#define _BOUNDARY_MODEL_H_

#include <vector>
#include "types.h"

class BoundaryModel
{
    protected:
        
        std::vector<Vector3r> r;

        // Contendra los rigid body que podran ser estaticos animados o dinamicos

    public:

        BoundaryModel() {}
        virtual ~BoundaryModel() {}

        virtual void init(std::vector<Vector3r> &);

        virtual void resizeBoundary(const unsigned int);
        virtual void cleanBoundary();
        unsigned int getNumParticles() { return r.size(); }

        Vector3r & getPosition(const unsigned int i) { return r[i]; }
        std::vector<Vector3r>& getPosition() { return r; }

        unsigned int size() { return r.size(); }
};

#endif
 
