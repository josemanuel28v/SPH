#ifndef _PCISPH_BOUNDARY_MODEL_H_
#define _PCISPH_BOUNDARY_MODEL_H_

#include "BoundaryModel.h"
#include "types.h"

// De momento esta clase solo sera para geometrias cubicas y esfericas
// Cuando se pueda samplear cualquier geometria se utilizaran otros metodos 
// de bh ya que este requiere un espaciado constante entre particulas para un correcto funcionamiento
class PCISPHBoundaryModel: public BoundaryModel
{
    protected:

        std::vector<Vector3r> v;
        std::vector<Vector3r> n;
        std::vector<Real> density;
        std::vector<Real> pressure;

        Real density0; // Debe ser la misma densidad que el fluido con el que va a colisionar por lo que implica solamente una fase
        Real mass;
        Real radius;
        Real normalFct;
        Real tangentialFct;

    public:

        void init(std::vector<Vector3r> &);

        void setNormalFct(Real normalFct) { this -> normalFct = normalFct; }
        void setTangentialFct(Real tangentialFct) { this -> tangentialFct = tangentialFct; }
        void setRadius(Real radius) { this -> radius = radius; }
        void setMass(Real mass) { this -> mass = mass; }
        void setRefDensity(Real density0) { this -> density0 = density0; }

        Real getNormalFct() { return normalFct; }
        Real getTangentialFct() { return tangentialFct; }
        Real getRadius() { return radius; }
        Real getMass() { return mass; }
        Real getRefDensity() { return density0; }
        
        Real & getDensity(const unsigned int i) { return density[i]; }
        Real & getPressure(const unsigned int i) { return pressure[i]; }

        void sampleCube();
        void sampleSphere();
        void sample();

        void resizeBoundary(const unsigned int);
        unsigned int size() { return r.size(); }

        void correctVelocities();
        void correctPositions();
        void correctPredPositions();
};

#endif  
