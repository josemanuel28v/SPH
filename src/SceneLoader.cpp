#include "SceneLoader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

void SceneLoader::read()
{
    std::cout << "Reading params..." << std::endl;
    if (readParameters())
        std::cout << "Done!" << std::endl << std::endl;
    else
        std::cerr << "Cannot open file " << path + ".prm" << std::endl << std::endl;

    std::cout << "Reading scene..." << std::endl;
    if (readScene())
        std::cout << "Done!" << std::endl << std::endl;
    else
        std::cerr << "Cannot open file " << path + ".scn" << std::endl << std::endl;

    std::cout << "Reading domain..." << std::endl;
    if (readDomain())
        std::cout << "Done!" << std::endl << std::endl;
    else
        std::cerr << "Cannot open file " << path + ".dom" << std::endl << std::endl;

    std::cout << "Reading positions and velocities..." << std::endl;
    if (readState(path))
        std::cout << "Done!" << std::endl << std::endl;
    else
        std::cerr << "Cannot open file " << path + ".dat" << std::endl << std::endl;
}

bool SceneLoader::readState(std::string path)
{
    std::ifstream in_file(path + ".dat");
    unsigned n;
    std::string tmp;
    std::string x, y, z;

    if (!in_file.is_open())
        return false;
    else
    {
        std::stringstream buffer;
        buffer << in_file.rdbuf();
        in_file.close();
        
        // Leer y guardar el ultimo paso de simulacion calculado para ese frame
        buffer >> tmp;     //time = stod(tmp); // Time
        buffer >> tmp;     // Num parts

        n = (unsigned) stoul(tmp);

        fluidData.position.resize(n);
        fluidData.velocity.resize(n);

        unsigned i = 0;
        while (buffer >> x)
        {
            buffer >> y;
            buffer >> z;

            if (i < n)
            {
                fluidData.position[i].x = (Real) stod(x);
                fluidData.position[i].y = (Real) stod(y);
                fluidData.position[i].z = (Real) stod(z);
            }
            else
            {
                fluidData.velocity[i - n].x = (Real) stod(x);
                fluidData.velocity[i - n].y = (Real) stod(y);
                fluidData.velocity[i - n].z = (Real) stod(z);
            }

            i++;
        }

        return true;
    }
}

bool SceneLoader::readScene()
{
    std::ifstream in_file(path + ".scn");
    std::string line;

    if (!in_file.is_open())
    {
        return false;  
    }
    else
    {
        while (getline(in_file, line))
        {
            size_t i = line.find("=") + 1;
            std::string param = line.substr(i, line.size());

            /*if (line.find("nparts") != line.npos)
                sys.resize(stoul(param));*/
            if (line.find("particle_volume") != line.npos)
                sceneData.particleVolume = static_cast<Real>(stod(param));
            /*else if (line.find("fluid_volume") != line.npos)
                sys.fluid_volume = stod(param);
            else if (line.find("kernel_parts") != line.npos)
                sys.kernel_parts = stoul(param);*/
            else if (line.find("tini") != line.npos)
                sceneData.startTime = stod(param);
            else if (line.find("tfin") != line.npos)
                sceneData.endTime = stod(param);
            else if (line.find("tstep") != line.npos)
                sceneData.timeStep = stod(param);
            else if (line.find("save_freq") != line.npos)
                sceneData.fps = stoul(param);
            else
                std::cout << "Cannot read line: " << line << std::endl;
        }

        in_file.close();

        return true;
    }
}

bool SceneLoader::readDomain()
{
    std::ifstream in_file(path + ".dom");
    std::string word;
    std::vector<std::string> dom_pts;

    domain.resize(2);
    
    if (!in_file.is_open())
        return false;
    else
    {
        while (in_file >> word)
            dom_pts.push_back(word);

        in_file.close();

        domain[0].x = (Real) stod(dom_pts[0]);
        domain[1].x = (Real) stod(dom_pts[1]);
        domain[0].y = (Real) stod(dom_pts[2]);
        domain[1].y = (Real) stod(dom_pts[3]);
        domain[0].z = (Real) stod(dom_pts[4]);
        domain[1].z = (Real) stod(dom_pts[5]);

        return true;
    }
} 

bool SceneLoader::readParameters()
{
    std::ifstream in_file(path + ".prm");
    std::string line;

    if (!in_file.is_open())
        return false; 
    else
    {
        while (getline(in_file, line))
        {
            size_t i = line.find("=") + 1;
            std::string param = line.substr(i, line.size()); // probar a quitar line.size de la llamada

            if (line.find("density") != line.npos)
                fluidData.density0 = stod(param);
            else if (line.find("visco") != line.npos)
                fluidData.viscosity = stod(param);
            /*else if (line.find("stiff") != line.npos)
                sys.stiffness = stod(param);*/
            else if (line.find("wsph_b") != line.npos)
                fluidData.stiffness = stod(param);
            else if (line.find("wsph_g") != line.npos)
                fluidData.gamma = stod(param);
            /*else if (line.find("wsph") != line.npos)
            {
                if (param == "True" || param == "true")
                    sys.model = SPHSystem::WCSPH;
            }*/
            else if (line.find("pcisph_fluct") != line.npos)
                fluidData.eta = stod(param);
            /*else if (line.find("pcisph") != line.npos)
            {
                if (param == "True" || param == "true")
                    sys.model = SPHSystem::PCISPH;
            }*/
            else if (line.find("surften_threshold") != line.npos)
                stod(param);
            else if (line.find("surften") != line.npos)
            {
                fluidData.surften = stod(param);
            }
            /*else if (line.find("collision_restitution") != line.npos)
                sys.restitution = stod(param);
            else if (line.find("collision_tang_conservation") != line.npos)
                sys.tang_conservation = stod(param);*/
            else if (line.find("gravity") != line.npos)
            {
                param.erase(remove(param.begin(), param.end(), ' '), param.end()); // Eliminar los espacios

                std::string x = param.substr(1, param.find_first_of(",") - 1);
                param = param.substr(param.find_first_of(",") + 1);

                std::string y = param.substr(0, param.find_first_of(","));
                param = param.substr(param.find_first_of(",") + 1);

                param.erase(param.end() - 1, param.end());
                std::string z = param;

                fluidData.gravity.x = (Real) stod(x);
                fluidData.gravity.y = (Real) stod(y);
                fluidData.gravity.z = (Real) stod(z);
            }
            else
                std::cout << "Cannot read line: " << line << std::endl;
        }

        in_file.close();

        return true;
    }
}
