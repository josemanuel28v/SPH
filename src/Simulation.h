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

// Structs para exportar e importar escenas en formato json
struct SceneInfo
{
    Real startTime;
    Real endTime;
    Real timeStep;
    Real fps;
    Real minTimeStep;
    Real maxTimeStep;
    Vector3r gravity;
    Real particleRadius;
    int simulationMethod;
    int boundaryMethod;
    Real stiffness;
    Real gamma;
    Real eta;
    Real etaV;
    Real cflFactor;
};

struct BlockInfo
{
    Vector3r min; // fluid block
    Vector3r max; // fluid block
};

struct EmitterInfo
{
    unsigned int type; 
    unsigned int numParticles; 
    Vector3r r;
    Real v; 
    Quat4r rot; // Cuando se aclare la forma de dar la orientacion desde blender
    Real startTime; 
    Real width; // Si el tipo es circular el radio sera width
    Real height;
    Real spacing;
};

struct Fluid // Representa un fluidModel en la misma fase que puede inicializarse con fluidBlocks geometrias de volumen y emisores
{
    std::vector<BlockInfo> fluidBlocks;
    std::vector<EmitterInfo> emitters;
    //std::vector<> geometries; // Conjunto de geometrias que seran sampleadas (volumen)

    Real viscosity;
    Real boundaryViscosity;
    Real surfaceTension;
    Real adhesion;
    Real density0;
};

struct FluidInfo
{
    std::vector<Fluid> fluids;

    int viscosityMethod;
    int surfaceTensionMethod;
    int adhesionMethod;
};

struct Sphere
{
    Vector3r pos;
    Real radius;
};

struct Geometry
{
    std::string path;
    Real spacing;
};

struct Boundary // Representa un boundary model
{
    std::vector<std::pair<BlockInfo, bool>> box;
    std::vector<std::pair<Sphere, bool>> sphere;
    std::vector<Geometry> geometry;

    Real normalFct;
    Real tangFct;
};

struct BoundaryInfo
{
    std::vector<Boundary> boundaries;
};

struct SimulationInfo
{
    SceneInfo sceneData;
    FluidInfo fluidData;
    BoundaryInfo boundaryData;
};

/**
 * @brief Clase singleton para controlar la simulacion actual, parametros generales y lectura y escritura en ficheros
 */

class Simulation
{
    private:

        static Simulation *current;
        std::string name;

    protected:

        std::vector<FluidModel*> fluidModels;
        std::vector<BoundaryModel*> boundaryModels;

        SPHSolver *solver;
        HashTable *grid;

        TimeManager tm;

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

        bool activeSave;

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
        const static int CUBIC_KERNEL_METHOD;
        
        // Kernel gradient method
        const static int POLY6_GRADIENT_METHOD;
        const static int SPIKY_GRADIENT_METHOD;
        const static int CUBIC_GRADIENT_METHOD;

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
        void initMasses();

        // Setters
        void setParticleRadius(Real);
        void setSupportRadius(Real);
        void setGravity(Vector3r gravity) { this -> gravity = gravity; }
        void setTime(Real time) { tm.setTime(time); }
        void setStartTime(Real startTime) { tm.setStartTime(startTime); }
        void setEndTime(Real endTime) { tm.setEndTime(endTime); }
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
        Real getParticleRadius() { return particleRadius; }
        Real getSupportRadius() { return supportRadius; }
        Vector3r getGravity() { return gravity; }
        Real getTime() { return tm.getTime(); }
        unsigned int getFrame() { return tm.getFrame(); }
        Real getStartTime() { return tm.getStartTime(); }
        Real getEndTime() { return tm.getEndTime(); }
        Real getTimeStep() { return tm.getTimeStep(); }
        Real getFPS() { return tm.getFPS(); }
        void startCounting(std::string name) { tm.startCounting(name); }
        void stopCounting(std::string name) { tm.stopCounting(name); }
        Real getInterval(std::string name) { return tm.getInterval(name); }
        Real getMinTimeStep() { return tm.getMinTimeStep(); }
        Real getMaxTimeStep() { return tm.getMaxTimeStep(); }
        std::string getName() { return name; }

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

        FluidModel* addFluidModel(std::vector<Vector3r>&, std::vector<Vector3r>&);
        FluidModel* addFluidModel();
        FluidModel* getFluidModel(const unsigned int i) { return fluidModels[i]; }

        void addBoundaryModel(BoundaryModel*);
        BoundaryModel* getBoundaryModel(const unsigned int i) { return boundaryModels[i]; }
        int getBoundaryHandlingMethod() { return current_bh_method; }

        __attribute__((always_inline)) Real W(const Vector3r & r) { return value(r); }
        __attribute__((always_inline)) Vector3r gradW(const Vector3r & r) { return grad(r); }
        __attribute__((always_inline)) Real laplW(const Vector3r & r) { return lapl(r); }

        void computeNonPressureForces();

        void printInfo();

        void emitParticles();

        bool importScene(std::string);

        void activateSave(bool active) { activeSave = active; }

        // Para pintar particulas con OpenGL
        const unsigned int numberActiveParticles(const unsigned int fmIndex) { return fluidModels[fmIndex] -> getNumActiveParticles(); }
        std::vector<Vector3r> & getPositions(const unsigned int fmIndex) { return fluidModels[fmIndex] -> getPositions(); }
        std::vector<Real> & getPressures(const unsigned int fmIndex) { return fluidModels[fmIndex] -> getPressures(); }
        std::vector<Real> getVelocities(const unsigned int fmIndex) { return fluidModels[fmIndex] -> getVelocities(); }
        std::vector<Vector3r> & getBoundaryPositions(const unsigned int bmIndex) { return boundaryModels[bmIndex] -> getPosition(); }
        std::vector<Real> & getDivergenceError(const unsigned int fmIndex);

        // Construir bloque de fluido
        FluidModel* buildFluidBlock(const std::vector<BlockInfo> &);

        void setEmittersPause(bool pause) { for (auto & fm: fluidModels) { fm -> setEmittersPause(pause); } }
};

#endif