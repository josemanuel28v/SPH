#include "OGL.h"
#include "SceneLoader.h"
#include "Simulation.h"
#include "AkinciBoundaryModel.h"
#include "PCISPHBoundaryModel.h"
#include "DFSPHSolver.h"
#include "WCSPHSolver.h"
#include "PCISPHSolver.h"
#include "InputParser.h"

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
		SceneLoader sceneLoader(InputParser::getCmdOption("-f"));
		sceneLoader.read();

		Simulation *sim = Simulation::getCurrent();

		Real maxTimeStep = 0.005;
		Real minTimeStep = 0.000;
		Real eta = 0.0005;  // 0.05%
		Real etaV = 0.001; // 0.1%
		Real bViscosity = 0.0;
		Real beta = 1.0; // 10000
		Real side = pow(sceneLoader.getParticleVolume(), 1.0 / 3.0); // Lado del cubo con el volumen recogido de blender
		Real radius = side / 2.0;


		////////// NOTA INICIALIZACION FLUIDO - RADIO Y SUPPORT RADIUS //////////////////////
		//
		//	El espaciado entre particulas sera de 2 * radius
		//	El supportRadius sera 4 * radius
		//	El volumen sera un cubo con lado 2 * radius
		//
		// 	El particleVolume que se coge de los ficheros de blender no es el volumen de la particula sino del cubo
		////////////////////// General parameters //////////////////////////
		sim -> setParticleRadius(radius);
		sim -> setGravity(sceneLoader.getGravity());
		sim -> setTimeStep(sceneLoader.getTimeStep());
		sim -> setStartTime(sceneLoader.getStartTime());
		sim -> setEndTime(sceneLoader.getEndTime());
		sim -> setFPS(sceneLoader.getFPS());
		sim -> setMinTimeStep(minTimeStep);
		sim -> setMaxTimeStep(maxTimeStep);
		sim -> setSurfaceTensionMethod(Simulation::AKINCI_SURFACE_TENSION_METHOD);
		sim -> setAdhesionMethod(Simulation::AKINCI_ADHESION_METHOD);
		sim -> setViscosityMethod(Simulation::ARTIFICIAL_VISCOSITY_METHOD); // 0.01 va bien
		//sim -> setViscosityMethod(Simulation::STANDARD_VISCOSITY_METHOD);
		//sim -> setViscosityMethod(Simulation::XSPH_VISCOSITY_METHOD); // entre 0 y 1
		/*sim -> setKernelMethod(Simulation::POLY6_KERNEL_METHOD);
		sim -> setGradientMethod(Simulation::SPIKY_GRADIENT_METHOD);
		sim -> setLaplacianMethod(Simulation::SPLINE_LAPLACIAN_METHOD);*/
		////////////////////////////////////////////////////////////////////

		////////////////////// Simulation method ///////////////////////////
		sim -> setSimulationMethod(Simulation::WCSPH_METHOD);
		WCSPHSolver *solver = static_cast<WCSPHSolver*>(sim -> getSolver());
		solver -> setStiffness(200);
		solver -> setGamma(1);

		/*sim -> setSimulationMethod(Simulation::PCISPH_METHOD);
		PCISPHSolver *solver = static_cast<PCISPHSolver*>(sim -> getSolver());
		solver -> setMaxError(sceneLoader.getMaxError());*/

		/*sim -> setSimulationMethod(Simulation::DFSPH_METHOD);
		DFSPHSolver *solver = static_cast<DFSPHSolver*>(sim -> getSolver());
		solver -> setMaxError(eta);
		solver -> setMaxErrorV(etaV);*/
		////////////////////////////////////////////////////////////////////

		////////////////////// Fluid Model /////////////////////////////////
		//fluidPoints.resize(0);
		//fluidVelocities.resize(0);
		FluidModel *fm = sim -> addFluidModel(sceneLoader.getPositions(), sceneLoader.getVelocities());

		/*std::vector<FluidBlockInfo> fbi_vec = {{Vector3r(0.0), 26, 26, 26}};
		FluidModel *fm = sim -> buildFluidBlock(radius, fbi_vec);*/

		fm -> setRefDensity(sceneLoader.getRefDensity());
		fm -> setMasses(sceneLoader.getRefDensity() * sceneLoader.getParticleVolume());
		fm -> setViscosityForce(sceneLoader.getViscosity(), bViscosity);
		fm -> setSurfaceTensionForce(sceneLoader.getSurfaceTension());
		//fm -> setAdhesionForce(beta);
		////////////////////////////////////////////////////////////////////

		////////////////////// Emitter /////////////////////////////////
		/*Matrix4r rot(1.0);
		//rot = glm::rotate(rot, 0.5 * M_PI, Vector3r(1.0, 0.0, 0.0));
		rot = glm::rotate(rot, M_PI * 0.5, Vector3r(0.0, 1.0, 0.0));

		fm -> addEmitter(Emitter::SQUARE_EMITTER, 
						40000, 
						Vector3r(0.01, 0.0, 0.3), 
						3.0, 
						rot, 
						0.0, 
						0.03, 
						0.03);*/
		////////////////////////////////////////////////////////////////////

		//////////////////// Boundary Model ////////////////////////////////
		/*CubeBoundaryModel *cbm = new CubeBoundaryModel();
		cbm -> init(boundaryPoints);
		sim -> addBoundaryModel(cbm);*/

		/*PCISPHBoundaryModel *pcibm = new PCISPHBoundaryModel();
		pcibm -> setMass(sceneLoader.getRefDensity() * sceneLoader.getParticleVolume());
		pcibm -> setRefDensity(sceneLoader.getRefDensity());
		pcibm -> init(sceneLoader.getDomain());
		sim -> setBoundaryMethod(Simulation::PCISPH_BOUNDARY_METHOD);
		sim -> addBoundaryModel(pcibm);*/

		// domain
		sim -> setBoundaryMethod(Simulation::AKINCI_BOUNDARY_METHOD);
		AkinciBoundaryModel *abm = new AkinciBoundaryModel();
		abm -> init(sceneLoader.getDomain());
		sim -> addBoundaryModel(abm);

		////////////////////////////////////////////////////////////////////
		
		///////////// Init simulation //////////////////////////////////////
		sim -> init();
		////////////////////////////////////////////////////////////////////

		//////////// Launch simulation in openGL ///////////////////////////
		OGL ogl;
    	ogl.mainLoop(argc, argv);
		////////////////////////////////////////////////////////////////////
	}
}

