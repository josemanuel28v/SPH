#include "OGL.h"
#include "SceneLoader.h"
#include "Simulation.h"
#include "AkinciBoundaryModel.h"
#include "PCISPHBoundaryModel.h"
#include "DFSPHSolver.h"
#include "WCSPHSolver.h"
#include "PCISPHSolver.h"
#include "InputParser.h"
#include <glm/gtc/matrix_transform.hpp>

void initSPHSimulation(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    initSPHSimulation(argc, argv);

	return 0;
}

void initSPHSimulation(int argc, char *argv[])
{
	InputParser(argc, argv);

	if (InputParser::cmdOptionExists("-f"))
	{
		Simulation *sim = Simulation::getCurrent();
		OGL ogl;
		std::string path = InputParser::getCmdOption("-f");

		if (InputParser::cmdOptionExists("-o"))
		{
			std::string activeOutput = InputParser::getCmdOption("-o");

			if (activeOutput == "true")
				sim -> activateSave(true);
			else if (activeOutput == "false")
				sim -> activateSave(false);
		}

		if (path.find("json") != string::npos)
		{
			if (sim -> importScene(path))
				ogl.mainLoop(argc, argv);
			else
				ERROR("Cannot import scene!");
		}
		else
			ERROR("File must be .json");

	}
}

