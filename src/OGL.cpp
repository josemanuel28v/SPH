#include "OGL.h"
#include "lodepng.h"
#include "AkinciBoundaryModel.h"
#include "DFSPHSolver.h"
#include "SceneLoader.h"
#include <iostream>
#include <chrono>
#include <sys/stat.h>
#include <glm/gtc/matrix_transform.hpp>

OGL::MaterialInfo gold = {glm::vec3(0.24725f, 0.1995f, 0.0745f), glm::vec3(0.75164f, 0.60648f, 0.22648f), glm::vec3(0.628281f, 0.555802f, 0.366065f), 52.0f};
OGL::MaterialInfo perl = {glm::vec3(0.25f, 0.20725f, 0.20725f), glm::vec3(1.0f, 0.829f, 0.829f), glm::vec3(0.296648f, 0.296648f, 0.296648f), 12.0f};
OGL::MaterialInfo bronze = {glm::vec3(0.2125f, 0.1275f, 0.054f), glm::vec3(0.714f, 0.4284f, 0.18144f), glm::vec3(0.393548f, 0.271906f, 0.166721f), 25.0f};
OGL::MaterialInfo brass = {glm::vec3(0.329412f, 0.223529f, 0.027451f), glm::vec3(0.780392f, 0.568627f, 0.113725f), glm::vec3(0.992157f, 0.941176f, 0.807843f), 28.0f};
OGL::MaterialInfo emerald = {glm::vec3(0.0215f, 0.1745f, 0.0215f), glm::vec3(0.07568f, 0.61424f, 0.07568f), glm::vec3(0.633f, 0.727811f, 0.633f), 28.0f};
OGL::MaterialInfo tin = {glm::vec3(0.105882f, 0.058824f, 0.113725f), glm::vec3(0.427451f, 0.470588f, 0.541176f), glm::vec3(0.333333f, 0.333333f, 0.521569f), 9.84615f};

GLint OGL::g_Width = 800;
GLint OGL::g_Height = 800;
GLboolean OGL::pause = false;
GLboolean OGL::cacheMode = false;
GLboolean OGL::fullscreen = false;
GLboolean OGL::mouseDown = false;
GLboolean OGL::drawBoundary = false;
GLboolean OGL::saveMode = false;
GLfloat OGL::xdiff;
GLfloat OGL::ydiff;
GLfloat OGL::xrot;
GLfloat OGL::yrot;
GLfloat OGL::desiredFps = 60;
GLfloat OGL::scaleFac = 6.0;
GLint OGL::colorMode = 1;
GLint OGL::m_index = 0;
std::string OGL::png_path;
double OGL::lastRender = 0;

glm::vec3 OGL::cameraPos = glm::vec3(0.0f, 0.5f, 5.0f);
glm::vec3 OGL::look = glm::vec3(0.0, 0.0, 0.0);

OGL::LightInfo OGL::pointLight  = { glm::vec4(0.0, 1.0, 5.0, 1.0), 1.0f * glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(1.0f, 0.00f, 0.000f), };
OGL::Cache OGL::cache;
OGL::MaterialInfo OGL::materials[] = {gold, perl, bronze, brass, emerald, tin};
OGL::MaterialInfo OGL::currentMaterial = brass;

ShaderLoader OGL::spriteShader;

void OGL::init()
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glClearDepth(1.0f);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glShadeModel(GL_SMOOTH);

	initSpriteProgram();
}

void OGL::initSpriteProgram()
{
	GLint programID = glCreateProgram();

	spriteShader.setProgramID(programID);
    spriteShader.addUniform("mvpm", 0);
    spriteShader.addUniform("mvm", 1);
    spriteShader.addUniform("matAmb", 7);
    spriteShader.addUniform("matDiff", 8);
    spriteShader.addUniform("matSpec", 9);
    spriteShader.addUniform("matShine", 10);
    spriteShader.addUniform("lPos", 11);
    spriteShader.addUniform("lIntensity", 12);
    spriteShader.addUniform("lK", 13);
    spriteShader.addUniform("color", 14);
    spriteShader.addUniform("colorMode", 15);
    spriteShader.addUniform("proyM", 3);
    spriteShader.addUniform("spriteTex", 4);
    spriteShader.addUniform("scale", 5);
    spriteShader.addUniform("pos", 6);
	spriteShader.addShader("shaders/sprite.vert");
	spriteShader.addShader("shaders/sprite.geom");
	spriteShader.addShader("shaders/sprite.frag");
	spriteShader.loadShaders();

	// Sprite texture
	std::vector<GLubyte> img_data;
	GLuint img_width, img_height;
	GLuint error;
	const GLchar img_filename[] = "textures/white_sphere.png";

	GLuint textId;
	glGenTextures(1, &textId);

	error = lodepng::decode(img_data, img_width, img_height, img_filename);
	if (!error)
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, textId);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, &img_data[0]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

		std::cout << "Textura cargada " << img_filename << std::endl;
	}
	else
		std::cout << "Error al cargar la textura " << img_filename << std::endl;
	img_data.clear();
}

int OGL::mainLoop(int argc, char *argv[])
{
	glutInit(&argc, argv); 
	glutInitWindowPosition(50, 50);
	glutInitWindowSize(g_Width, g_Height);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Particles");
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
	  // Problem: glewInit failed, something is seriously wrong. 
	  std::cerr << "Error: " << glewGetErrorString(err) << std::endl;
	  //system("pause");
	  exit(-1);
	}
	init();

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(specialKeyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(mouseMotion);
	glutReshapeFunc(resize);
	glutIdleFunc(idle);
 
	glutMainLoop();
 
	return EXIT_SUCCESS;
}

void OGL::drawSprites()
{
	glDrawArrays(GL_POINTS, 0, 1);
}

void OGL::displaySpriteScene(const glm::mat4 & Projection, const glm::mat4 & View, const glm::mat4 & Model)
{
	Simulation *sim = Simulation::getCurrent();

	glm::mat4 mv = View * Model;  
    glm::vec4 lpos = View * Model * pointLight.lightPos;
	glm::vec4 pos;

	vector<Vector3r> * positions;
	vector<Vector3r> * bpositions;
	vector<Real> * parameter;
	vector<Real> vel; 
	uint n, n_boundary = 0;
	float max_param;

    spriteShader.begin();

    glUniform1i(spriteShader.getUniformLocation("spriteTex"), 0);
	glUniform1f(spriteShader.getUniformLocation("scale"), sim -> getParticleRadius() * scaleFac * 1.3); // * 1.3 tamaño real
	glUniformMatrix4fv(spriteShader.getUniformLocation("mvm"), 1, GL_FALSE, &mv[0][0]);
	glUniformMatrix4fv(spriteShader.getUniformLocation("proyM"), 1, GL_FALSE, &Projection[0][0]);
	glUniform3fv(spriteShader.getUniformLocation("matAmb"), 1, &(currentMaterial.ambient.x));
	glUniform3fv(spriteShader.getUniformLocation("matDiff"), 1, &(currentMaterial.diffuse.x));
	glUniform3fv(spriteShader.getUniformLocation("matSpec"), 1, &(currentMaterial.specular.x));
	glUniform1f(spriteShader.getUniformLocation("matShine"), currentMaterial.shininess);
	glUniform4fv(spriteShader.getUniformLocation("lPos"), 1, &(lpos.x));
	glUniform3fv(spriteShader.getUniformLocation("lIntensity"), 1, &(pointLight.intensity.r));
	glUniform3fv(spriteShader.getUniformLocation("lK"), 1, &(pointLight.k[0]));
	glUniform1i(spriteShader.getUniformLocation("colorMode"), colorMode);

	// Posiciones guardadas
	if (cacheMode && cache.position.size() > 0)
	{
		n = cache.numParts[cache.frameCount];
		positions = &cache.position[cache.frameCount];
		
		if (colorMode == 2 && cache.pressure[cache.frameCount].size() != 0)
		{
			parameter = &cache.pressure[cache.frameCount];
			max_param = 3000;
		}
		else
		{
			parameter = &cache.velocity[cache.frameCount];
			max_param = 2;
		}
	}
	// Posiciones actuales del sistema de particulas
	else
	{
		n = sim -> numberActiveParticles(0);
		positions = &sim -> getPositions(0);

		if (colorMode == 2)
		{
			parameter = &sim -> getDivergenceError(0);
			max_param = 3000;
			//max_param = sim -> getFluidModel(0) -> getRefDensity();
		}
		else 
		{
			vel = sim -> getVelocities(0);
			parameter = &vel;
			max_param = 2;			
		}
	}

	if (drawBoundary)
	{
		// Boundary
		for (unsigned int bmi = 0; bmi < sim -> numberBoundaryModels(); ++bmi)
		{
			bpositions = &sim -> getBoundaryPositions(bmi);
			n_boundary = bpositions -> size();
			for (uint i = 0; i < n_boundary; ++i)
			//for (uint i: sim.sys.movingWallIndexPX)
			{		
				pos = glm::vec4((bpositions->at(i)), 1.0);
				glUniform4fv(spriteShader.getUniformLocation("pos"), 1, &pos.x);
				glUniform1f(spriteShader.getUniformLocation("color"), 0.5);

				drawSprites();
			}
		}
	}

	// Fluid
	for (uint i = 0; i < n; ++i)
	{	
		pos = glm::vec4((positions->at(i)), 1.0);
		glUniform4fv(spriteShader.getUniformLocation("pos"), 1, &pos.x);
		glUniform1f(spriteShader.getUniformLocation("color"), glm::clamp(sqrt((float) parameter->at(i) / max_param), 0.0, 1.0));

		drawSprites();
	}

	spriteShader.end();
}

void OGL::display()
{
	double frame_delay = 1.0 / desiredFps; 
	double currentTime = getCurrentTime();

	if ((currentTime - lastRender) > frame_delay)
	{
		updateTitle(currentTime - lastRender);

		lastRender = currentTime;

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		

        glm::mat4 Projection = glm::perspective(45.0f, 1.0f * g_Width / g_Height, 1.0f, 100.0f);
        glm::mat4 View = glm::lookAt(cameraPos, look, glm::vec3(0.0f, 1.0f, 0.0f));
        glm::mat4 Model = glm::scale(glm::mat4(1.0f), glm::vec3(scaleFac, scaleFac, scaleFac));

        // Rotacion del modelo mediante el raton
        Model = glm::rotate(Model, xrot * 0.0025f, glm::vec3(1.0f, 0.0f, 0.0f));
        Model = glm::rotate(Model, yrot * 0.0025f, glm::vec3(0.0f, 1.0f, 0.0f));

        // Rotacion para pasar del s.c. de Blender (z -> altura y -> profundidad) al de OpenGL (y -> altura z -> profundidad)
        Model = glm::rotate(Model, - (GLfloat) (M_PI * /*0.35*/ 0.5), glm::vec3(1, 0, 0));
        //Model = rotate(Model, (GLfloat) (M_PI * 0.25), vec3(0, 0, 1));

        // Dibujar la escena
        displaySpriteScene(Projection, View, Model);

        glutSwapBuffers();

        // Guardar frames
        /*if (saveMode)
        {
            savePNGFrame();
            if (cache.frameCount == cache.position.size() - 1)
                saveMode = false;
        }*/

        // Actualizar contador de la cache
        if (!pause && cacheMode && cache.frameCount < cache.position.size() - 1)
            ++cache.frameCount;
	}
}

void OGL::keyboard(GLubyte key, GLint x, GLint y)
{
	switch(key)
	{
		case 27 : 
		case 'q': 
		case 'Q':
			exit(1); 
			break;

		case 'i':
		case 'I':
			pointLight.lightPos.z += 0.1;
			break;

		case 'k':
		case 'K':
			pointLight.lightPos.z -= 0.1;
			break;

		case 'c':
		case 'C':
			cache.frameCount = 0;
			cacheMode = !cacheMode;
			break;

		case 'p':
		case 'P':
			pause = !pause;
			break;

		case 'n':
		case 'N':
			colorMode = (colorMode + 1) % 5; // Material -> Normal -> Velocities -> Pressure -> Sin iluminacion
			break;

		case '0':
			cache.frameCount = 0;
			break;

		case '+':
			scaleFac += 0.1;
			break;

		case '-':
			scaleFac -= 0.1;
			break;

		case 'b':
		case 'B':
			drawBoundary = !drawBoundary;
			break;	

		case '8':	// UP
			look.z += 0.1;
			cameraPos.z += 0.1;
			break;

		case '2':	// DOWN
			look.z -= 0.1;
			cameraPos.z -= 0.1;
			break;

		case '4':	// LEFT
			look.x -= 0.1;
			cameraPos.x -= 0.1;
			break;

		case '6':	// RIGHT
			look.x += 0.1;
			cameraPos.x += 0.1;
			break;

		case '9':	// BEHIND
			look.y -= 0.1;
			cameraPos.y -= 0.1;
			break;

		case '1':	// FRONT
			look.y += 0.1;
			cameraPos.y += 0.1;
			break;

		case 'm': // Change material
			m_index = (m_index + 1) % 6;
			currentMaterial = materials[m_index];
			break;

		case 'w': // Save frames
		{
			// Comprobar si existe el directorio frames, si no existe, se crea (Linux)
			std::string png_dir_path = png_path.substr(0, png_path.find_last_of("/"));
			struct stat st;
			if(stat(png_dir_path.c_str(), &st) != 0)
			{
				std::cout << std::endl << "Se ha creado el directorio " + png_dir_path << std::endl;
				std::cout << mkdir(png_dir_path.c_str(), 0777) << std::endl;	
			}

			saveMode = !saveMode;
			cacheMode = saveMode;
			cache.frameCount = 0;
		}
			break;

	}
}

void OGL::specialKeyboard(GLint key, GLint x, GLint y)
{
	if (key == GLUT_KEY_F1)
	{
		fullscreen = !fullscreen;
 
		if (fullscreen)
			glutFullScreen();
		else
		{
			glutReshapeWindow(g_Width, g_Height);
			glutPositionWindow(50, 50);
		}
	}
	else if (key == GLUT_KEY_LEFT)
	{
		desiredFps -= 1;
		if (desiredFps < 1)
			desiredFps = 1;
	}
	else if (key == GLUT_KEY_RIGHT)
		desiredFps += 1;
}

void OGL::mouse(GLint button, GLint state, GLint x, GLint y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		mouseDown = true;
 
		xdiff = x - yrot;
		ydiff = -y + xrot;
	}
	else
		mouseDown = false;
}

void OGL::mouseMotion(GLint x, GLint y)
{
	if (mouseDown)
	{
		yrot = x - xdiff;
		xrot = y + ydiff;
 
		glutPostRedisplay();
	}
}

void OGL::resize(GLint w, GLint h)
{
	g_Width = w;
	g_Height = h;
	glViewport(0, 0, g_Width, g_Height);
}

void OGL::idle()
{	
	Simulation *sim = Simulation::getCurrent();

    if (!pause && !cacheMode)
			if (sim -> step())
			{
				cache.position.push_back(sim -> getPositions(0));
				cache.numParts.push_back(sim -> numberActiveParticles(0));
				cache.time.push_back(sim -> getTime());

				//cacheParam3D.push_back(sim.getVelocities3D());
				cache.velocity.push_back(sim -> getVelocities(0));
				cache.pressure.push_back(sim -> getPressures(0));

				cache.size ++;
			}

	glutPostRedisplay();
}

double OGL::getCurrentTime()
{
	using Duration = std::chrono::duration<double>;
	return std::chrono::duration_cast<Duration>( 
		std::chrono::high_resolution_clock::now().time_since_epoch() 
	).count();
}

void OGL::updateTitle(double delay)
{
	static unsigned int count = 0;

	if (count % 10 == 0)
	{
		double current_fps = round(1.0 / delay);
		std::string current_fps_s = std::to_string(current_fps);
		std::string desiredFps_s =  std::to_string(desiredFps);
		std::string current_frame = std::to_string(cache.frameCount + 1);
		std::string total_frames  = std::to_string(cache.position.size());
		std::string title = current_fps_s.substr(0, current_fps_s.find(".") + 3) + "/" + 
							desiredFps_s.substr(0, desiredFps_s.find(".")) + "		" + 
							current_frame + "/" + total_frames;
		glutSetWindowTitle(title.c_str());
	}

	count ++;
}
