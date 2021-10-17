#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <vector>
#include <algorithm>
#include "types.h"
#include "FluidModel.h"
#include "BoundaryModel.h"
#include "NonPressureForce.h"
#include "TimeManager.h"
#include "SPHSolver.h"
#include "HashTable.h"

// Loop over neighborhood

#define forall_fluid_neighbors(code) \
    for (auto & pointInfo: grid -> table[cellId].neighbors) \
    { \
        unsigned int j = pointInfo.id; \
        unsigned int nfmIndex = pointInfo.pointSetId; \
        FluidModel *nfm = sim -> getFluidModel(nfmIndex); \
        code \
    } 

#define forall_fluid_neighbors_in_same_phase(code) \
    for (auto & pointInfo: grid -> table[cellId].neighbors) \
    { \
        unsigned int j = pointInfo.id; \
        code \
    } 

#define forall_boundary_neighbors(code) \
    for (auto & pointInfo: grid -> table[cellId].boundary_neighbors) \
    { \
        unsigned int b = pointInfo.id; \
        unsigned int nbmIndex = pointInfo.pointSetId; \
        code \
    } 

#define forsame_boundary_neighbors(code) \
    for (auto & pointInfo: grid -> table[cellId].boundary_neighbors) \
    { \
        unsigned int b = pointInfo.id; \
        code \
    } 

struct FluidBlockInfo
{
    Vector3r origin;
    unsigned int nx; 
    unsigned int ny; 
    unsigned int nz;
};

/**
 * @brief Clase singleton para controlar la simulacion actual, parametros generales y lectura y escritura en ficheros
 */

class Simulation
{
    private:

        static Simulation *current;

    protected:

        std::vector<FluidModel*> fluidModels;
        std::vector<BoundaryModel*> boundaryModels;

        SPHSolver *solver;
        HashTable *grid; // Ver si es su mejor sitio

        TimeManager tm;

        unsigned int kernelParticles;
        Real particleRadius;
        Real supportRadius;
        Vector3r gravity;

        int current_method;
        int current_bh_method;
        int current_visco_method;
        int current_surften_method;
        int current_adhesion_method;
        int current_kernel_method;
        int current_gradient_method;
        int current_laplacian_method;

    public:

        Real (*value) (const Vector3r &);
        Vector3r (*grad) (const Vector3r &);
        Real (*lapl) (const Vector3r &);

        // Solver methods
        const static int WCSPH_METHOD;
        const static int PCISPH_METHOD;
        const static int DFSPH_METHOD;

        // Boundary methods
        const static int CUBE_BOUNDARY_METHOD;
        const static int PCISPH_BOUNDARY_METHOD;
        const static int AKINCI_BOUNDARY_METHOD;

        // Non pressure forces
        const static int STANDARD_VISCOSITY_METHOD;
        const static int ARTIFICIAL_VISCOSITY_METHOD;
        const static int XSPH_VISCOSITY_METHOD;
        const static int AKINCI_SURFACE_TENSION_METHOD;
        const static int AKINCI_ADHESION_METHOD;

        // Kernel value method
        const static int POLY6_KERNEL_METHOD;
        const static int SPIKY_KERNEL_METHOD;
        
        // Kernel gradient method
        const static int POLY6_GRADIENT_METHOD;
        const static int SPIKY_GRADIENT_METHOD;

        // Kernel laplacian method
        const static int SPLINE_LAPLACIAN_METHOD;

        Simulation();
        ~Simulation();

        static Simulation* getCurrent();
        static void setCurrent(Simulation*);
        static bool hasCurrent();

        void init();
        bool step();
        void run();

        void initKernels();
        void initGrid();

        // Setters
        void setKernelParticles(const int kernelParticles) { this -> kernelParticles = kernelParticles; }
        void setParticleRadius(Real);
        void setSupportRadius(Real);
        void setGravity(Vector3r gravity) { this -> gravity = gravity; }
        void setTime(Real time) { tm.setTime(time); }
        void setTimeStep(Real ts) { tm.setTimeStep(ts); }
        void setFPS(Real fps) { tm.setFPS(fps); }
        void setMinTimeStep(Real minTs) { tm.setMinTimeStep(minTs); }
        void setMaxTimeStep(Real maxTs) { tm.setMaxTimeStep(maxTs); }

        void setSimulationMethod(const int);
        void setBoundaryMethod(const int);
        void setKernelMethod(const int);
        void setGradientMethod(const int);
        void setLaplacianMethod(const int);
        void setViscosityMethod(const int method) { current_visco_method = method; }
        void setSurfaceTensionMethod(const int method) { current_surften_method = method; }
        void setAdhesionMethod(const int method) { current_adhesion_method = method; }

        // Getters
        unsigned int getKernelParticles() { return kernelParticles; }
        Real getParticleRadius() { return particleRadius; }
        Real getSupportRadius() { return supportRadius; }
        Vector3r getGravity() { return gravity; }
        Real getTime() { return tm.getTime(); }
        Real getTimeStep() { return tm.getTimeStep(); }
        void startCounting() { tm.startCounting(); }
        void stopCounting() { tm.stopCounting(); }
        Real getInterval() { return tm.getInterval(); }
        Real getMinTimeStep() { return tm.getMinTimeStep(); }
        Real getMaxTimeStep() { return tm.getMaxTimeStep(); }

        int getSimulationMethod() { return current_method; }
        int getKernelMethod();
        int getGradientMethod();
        int getLaplacianMethod();
        int getViscosityMethod() { return current_visco_method; }
        int getSurfaceTensionMethod() { return current_surften_method; }
        int getAdhesionMethod() { return current_adhesion_method; }

        SPHSolver* getSolver() { return solver; }

        HashTable* getGrid() { return grid; }
        const unsigned int numberFluidModels() { return fluidModels.size(); }
        const unsigned int numberBoundaryModels() { return boundaryModels.size(); }

        void addFluidModel(std::vector<Vector3r>&, std::vector<Vector3r>&);
        FluidModel* getFluidModel(const unsigned int i) { return fluidModels[i]; }

        void addBoundaryModel(BoundaryModel*);
        BoundaryModel* getBoundaryModel(const unsigned int i) { return boundaryModels[i]; }
        int getBoundaryHandlingMethod() { return current_bh_method; }

        __attribute__((always_inline)) Real W(const Vector3r & r) { return value(r); }
        __attribute__((always_inline)) Vector3r gradW(const Vector3r & r) { return grad(r); }
        __attribute__((always_inline)) Real laplW(const Vector3r & r) { return lapl(r); }

        void emitParticles();

        // se necesita para el programa de opengl antiguo
        const unsigned int numberActiveParticles(const unsigned int fmIndex) { return fluidModels[fmIndex] -> getNumActiveParticles(); }
        std::vector<Vector3r> & getPositions(const unsigned int fmIndex) { return fluidModels[fmIndex] -> getPositions(); }
        std::vector<Real> & getPressures(const unsigned int fmIndex) { return fluidModels[fmIndex] -> getPressures(); }
        std::vector<Real> getVelocities(const unsigned int fmIndex) { return fluidModels[fmIndex] -> getVelocities(); }
        std::vector<Vector3r> & getBoundaryPositions(const unsigned int bmIndex) { return boundaryModels[bmIndex] -> getPosition(); }
        std::vector<Real> & getDivergenceError(const unsigned int fmIndex);

        // Construir bloque de fluido
        FluidModel* buildFluidBlock(Real, std::vector<FluidBlockInfo>);

        // Añadir boundary model a partir de un modelo 3d obj
        void addBoundaryModelFromOBJ(std::string, Real, Vector3r, Vector3r, Vector3r);
};

#endif