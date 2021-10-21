#ifndef _SCENE_LOADER_H_
#define _SCENE_LOADER_H_

#include <string>
#include <vector>
#include "types.h"

class SceneLoader
{
    public:

        struct FluidModelInfo
        {
            std::vector<Vector3r> position;
            std::vector<Vector3r> velocity;
            Real density0;
            Real eta;
            Real viscosity;
            Real stiffness;
            Real gamma;
            Real surften;
            Vector3r gravity;
        };

        struct SceneInfo
        {
            Real particleVolume;
            Real startTime;
            Real endTime;
            Real timeStep;
            Real fps;
        };

        SceneLoader(std::string path) { this -> path = path; }

        // Problema con las funciones que leen con un if else if
        // parametros que no se lean pero empiecen igual que otro 
        // que si se lea sobreescriben el parametro que se queria
        // leer en un principio
        bool readDomain();
        bool readScene();
        bool readState(std::string);
        bool readParameters();
        void read();

        std::vector<Vector3r> & getPositions() { return fluidData.position; }
        std::vector<Vector3r> & getVelocities() { return fluidData.velocity; }
        Real getRefDensity() { return fluidData.density0; }
        Real getMaxError() { return fluidData.eta; } // Mover a sceneData si se utiliza el mismo eta para todos los fluidModels
        Real getViscosity() { return fluidData.viscosity; }
        Real getStiffness() { return fluidData.stiffness; }
        Real getGamma() { return fluidData.gamma; }
        Real getSurfaceTension() { return fluidData.surften; }
        Vector3r getGravity() { return fluidData.gravity; }
        Real getParticleVolume() { return sceneData.particleVolume; }
        Real getStartTime() { return sceneData.startTime; }
        Real getEndTime() { return sceneData.endTime; }
        Real getTimeStep() { return sceneData.timeStep; }
        Real getFPS() { return sceneData.fps; }
        std::vector<Vector3r> & getDomain() { return domain; }

    private:

        std::string path;

        FluidModelInfo fluidData;
        SceneInfo sceneData;
        std::vector<Vector3r> domain;
};

#endif
