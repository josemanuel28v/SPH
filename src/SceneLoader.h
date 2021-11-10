#ifndef _SCENE_LOADER_H_
#define _SCENE_LOADER_H_

#include <string>
#include <vector>
#include "types.h"
#include "Simulation.h"
#include "../extern/json.hpp"

using nlohmann::json;

class SceneLoader
{
    private:

        static unsigned int maxDigits;

    public:

        static bool readConfiguration(SimulationInfo &, std::string path);
        static void readSceneInfo(SimulationInfo &, json file);
        static void readFluidInfo(SimulationInfo &, json file);
        static void readBoundaryInfo(SimulationInfo &, json file);
        static void writeFluid();
};

#endif
