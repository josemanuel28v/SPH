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

        Real normalFct;
        Real tangentialFct;

    public:

        void setNormalFct(Real normalFct) { this -> normalFct = normalFct; }
        void setTangentialFct(Real tangentialFct) { this -> tangentialFct = tangentialFct; }
        void setRefDensity(Real density0) { this -> density0 = density0; }
        void setMass(Real mass) { this -> mass = mass; }

        Real getNormalFct() { return normalFct; }
        Real getTangentialFct() { return tangentialFct; }
        Real getRefDensity() { return density0; }
        Real getMass() { return mass; }
        
        Real & getDensity(const unsigned int i) { return density[i]; }
        Real & getPressure(const unsigned int i) { return pressure[i]; }

        void sampleCube(Vector3r min, Vector3r max, std::vector<Vector3r> & points, std::vector<Vector3r> & normals);
        void sampleSphere(Vector3r origen, Real radius, std::vector<Vector3r> & points);

        void addCube(Vector3r min, Vector3r max, bool);
        void addSphere(Vector3r origen, Real radius, bool);

        void resizeBoundary(const unsigned int);

        void correctVelocities();
        void correctPositions();
        void correctPredPositions();
};

#endif  
