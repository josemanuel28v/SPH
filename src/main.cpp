#include <GL/glew.h>
#include <GL/glut.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtx/string_cast.hpp>

#include "lodepng.h"
#include "Simulation.h"
//#include "Domain.h"
//#include "InputParser.h"
#include "ShaderLoader.h"
//#include "FluidBuilder.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <vector>
#include <sys/stat.h>

using std::cout;
using std::endl;
using std::string;

using glm::vec3;
using glm::vec4;
using glm::mat3;
using glm::mat4;
using glm::clamp;
using glm::rotate;
using glm::translate;
using glm::scale;
using glm::lookAt;
using glm::perspective;

/********NEW********/
using namespace std;

bool readState(string, vector<Vector3r>&, vector<Vector3r>&);
bool readScene(string);

Real particleVolume;
Simulation *sim = Simulation::getCurrent();

inline bool solve();
/*******************/

// Inicializa localizadores uniform
void initUniform();

// Sprites
void initSpriteProgram();

// Funciones para pintar las particulas como esferas
void initSphereProgram();
GLint buildSphere(GLfloat, GLuint, GLuint);
GLint buildPlane(GLfloat, GLfloat, GLint, GLint);
void drawSphere();
void drawSPlane();
void displayGeometryScene(const mat4 &, const mat4 &, const mat4 &);

// Funciones de callback que se van a registrar
GLboolean init();
void display();
void resize(GLint, GLint);
void idle();
void keyboard(GLubyte, GLint, GLint);
void specialKeyboard(GLint, GLint, GLint);
void mouse(GLint, GLint, GLint, GLint);
void mouseMotion(GLint, GLint);

// Funcion principal que inicializa la simulacion
void simulationMain(int, char* []);

// Funciones de tiempo
bool printFPS(float &);
double getCurrentTime();

// Estructura para la informacion de la luz
struct LightInfo 
{
	vec4 lightPos; 	
	vec3 intensity;
	vec3 k;
};

// Estructura para la informacion de un material
struct MaterialInfo 
{
	vec3 ambient;
	vec3 diffuse;
	vec3 specular;
	GLfloat shininess;
};

// Estructura para almacenar las posiciones, velocidades, presiones y numero de particulas de cada frame de la simulacion
struct Cache 
{
	vector<vector<Vector3r>> position;
	vector<vector<Real>> velocity;
	vector<vector<Real>> pressure;

	vector<GLuint> numParts;
	vector<Real> time;

	GLuint frameCount = 0;
	GLuint size = 0;
};

enum PROGRAM
{
	SPRITE = 0,
	SPHERE = 1
};

LightInfo pointLight = 
{
	vec4(0.0, 1.0, 5.0, 1.0),
	1.0f * vec3(1.0f, 1.0f, 1.0f), 
	vec3(1.0f, 0.00f, 0.000f), 
};

LightInfo dirLight = 
{
	vec4(0.0, -1.0, 0.0, 0.0),
	0.5f * vec3(1.0f, 1.0f, 1.0f), 
	vec3(0.0f, 0.0f, 0.0f), 
};

// Materiales
MaterialInfo gold = {vec3(0.24725f, 0.1995f, 0.0745f), vec3(0.75164f, 0.60648f, 0.22648f), vec3(0.628281f, 0.555802f, 0.366065f), 52.0f};
MaterialInfo perl = {vec3(0.25f, 0.20725f, 0.20725f), vec3(1.0f, 0.829f, 0.829f), vec3(0.296648f, 0.296648f, 0.296648f), 12.0f};
MaterialInfo bronze = {vec3(0.2125f, 0.1275f, 0.054f), vec3(0.714f, 0.4284f, 0.18144f), vec3(0.393548f, 0.271906f, 0.166721f), 25.0f};
MaterialInfo brass = {vec3(0.329412f, 0.223529f, 0.027451f), vec3(0.780392f, 0.568627f, 0.113725f), vec3(0.992157f, 0.941176f, 0.807843f), 28.0f};
MaterialInfo emerald = {vec3(0.0215f, 0.1745f, 0.0215f), vec3(0.07568f, 0.61424f, 0.07568f), vec3(0.633f, 0.727811f, 0.633f), 28.0f};
MaterialInfo tin = {vec3(0.105882f, 0.058824f, 0.113725f), vec3(0.427451f, 0.470588f, 0.541176f), vec3(0.333333f, 0.333333f, 0.521569f), 9.84615f};

MaterialInfo materials[] = {gold, perl, bronze, brass, emerald, tin};
MaterialInfo * currentMaterial = &gold;
const GLuint NUM_MATERIALS = 6;
GLuint m_index = 0;
 
// Rotacion del modelo
GLfloat xrot = 0.0f;
GLfloat yrot = 0.0f;
GLfloat xdiff = 0.0f;
GLfloat ydiff = 0.0f;

GLint g_Width = 800;                          // Ancho inicial de la ventana
GLint g_Height = 800;                         // Altura incial de la ventana

// Parametros de camara
vec3 look = vec3(0.0, 0.0, 0.0);
vec3 cameraPos = vec3(0.0f, 0.5f, 5.0f);

// Variables de control la interfaz
GLboolean cache_mode = false;
GLboolean pause = false;
GLboolean drawBoundary = false;
GLboolean fullscreen = false;
GLboolean mouseDown = false;
GLboolean animation = false;
GLboolean sprite = true;
GLboolean graphic_mode = true;
GLboolean save_mode = false;
GLuint colorMode = 3;

// Frecuencia de dibujado deseada
GLfloat desired_fps = 60;
GLfloat current_fps = 60;

// Radio y escalado de las particulas
GLfloat particleRad;
GLfloat particleScale = 1;

// Sphere and plane
GLuint sphereVAO, planeVAO, spriteVAO;
GLuint numVertSphere, numVertPlane;
GLfloat scaleFac = 6;
GLuint programID[2];

// Localizadores globales
GLuint locUniformMVPM, locUniformMVM, locUniformNM, locUniformPM;

// Localizadorese las variables uniform de la luz y los materiales (Geometry Program)
GLuint locUniformLightPos, locUniformLightIntensity, locUniformLightK;
GLuint locUniformDirLightPos, locUniformDirLightIntensity, locUniformDirLightK;
GLuint locUniformMaterialAmbient, locUniformMaterialDiffuse, locUniformMaterialSpecular, locUniformMaterialShininess;
GLuint locUniformObject, locUniformColorMode, locUniformParamColor, locUniformParamColor3D;

// Loclizador sprite texture
GLuint locUniformSpriteTex, locUniformScale, locUniformPos, locUniformColor;
GLuint spriteHandle[2];

// Ruta para los ficheros .png
string png_path;

Cache cache;

/**
 * @brief Lee los pixels del frame buffer y los exporta a un fichero .png utilizando lodepng
 * 
 */
void savePNGFrame()
{
	cout << "Frame " << cache.frameCount + 1;

	vector<GLubyte> imageData(g_Height * g_Width * 4);

	glFinish();
	glReadPixels(0, 0, g_Width, g_Height, GL_RGBA, GL_UNSIGNED_BYTE, imageData.data());

	std::vector<GLubyte> pngBuffer(imageData.size());

	for(GLint i = 0; i < g_Height; ++i)
	{
		for(GLint j = 0; j < g_Width; ++j)
		{
			size_t OldPos = (g_Height - i - 1) * (g_Width * 4) + 4 * j;
			size_t NewPos = i * (g_Width * 4) + 4 * j;
			pngBuffer[NewPos + 0] = imageData[OldPos + 0]; 
			pngBuffer[NewPos + 1] = imageData[OldPos + 1]; 
			pngBuffer[NewPos + 2] = imageData[OldPos + 2]; 
			pngBuffer[NewPos + 3] = imageData[OldPos + 3]; 
		}
	}

	vector<GLubyte> outBuffer;
	lodepng::encode(outBuffer, pngBuffer, g_Width, g_Height);
	lodepng::save_file(outBuffer, png_path + to_string(cache.frameCount + 1) + ".png");

	cout << " saved" << endl;
}

/**
 * @brief Construye la geometria de una esfera preparando los VBO y su VAO
 * 
 * @param radius Radio de la esfera
 * @param rings Anillos de la esfera
 * @param sectors Sectores de la esfera
 * 
 * @return Numero de vertices de la esfera
 */
GLint buildSphere(GLfloat radius, GLuint rings, GLuint sectors)
{
    const GLfloat R = 1.0f/(GLfloat)(rings-1);
    const GLfloat S = 1.0f/(GLfloat)(sectors-1);
	const GLfloat PI = 3.14159265358979323846;

    GLfloat *sphere_vertices = new GLfloat[rings * sectors * 3];
    GLfloat *sphere_normals = new GLfloat[rings * sectors * 3];
    GLfloat *sphere_texcoords = new GLfloat[rings * sectors * 2];
    GLfloat *v = sphere_vertices;
    GLfloat *n = sphere_normals;
    GLfloat *t = sphere_texcoords;
    for(GLuint r = 0; r < rings; r++) for(GLuint s = 0; s < sectors; s++) {
            GLfloat const y = GLfloat( sin( -PI/2 + PI * r * R ) );
            GLfloat const x = GLfloat( cos(2*PI * s * S) * sin( PI * r * R ) );
            GLfloat const z = GLfloat( sin(2*PI * s * S) * sin( PI * r * R ) );

            *t++ = s*S;
            *t++ = r*R;

            *v++ = x * radius;
            *v++ = y * radius;
            *v++ = z * radius;

            *n++ = x;
            *n++ = y;
            *n++ = z;
    }

    GLushort *sphere_indices = new GLushort[rings * sectors * 4];
    GLushort *i = sphere_indices;
    for(GLuint r = 0; r < rings; r++) for(GLuint s = 0; s < sectors; s++) {
            *i++ = r * sectors + s;
            *i++ = r * sectors + (s+1);
            *i++ = (r+1) * sectors + (s+1);
            *i++ = (r+1) * sectors + s;
    }

    glGenVertexArrays( 1, &sphereVAO );
    glBindVertexArray(sphereVAO);

    GLuint handle[4];
    glGenBuffers(4, handle);

    glBindBuffer(GL_ARRAY_BUFFER, handle[0]);
    glBufferData(GL_ARRAY_BUFFER, (rings * sectors * 3) * sizeof(GLfloat), sphere_vertices, GL_STATIC_DRAW); // Datos de la posicion de los vertices
	GLuint loc1 = glGetAttribLocation(programID[SPHERE], "aPosition");   
	glEnableVertexAttribArray(loc1); // Vertex position
	glVertexAttribPointer( loc1, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL + 0 ); 

    glBindBuffer(GL_ARRAY_BUFFER, handle[1]);
    glBufferData(GL_ARRAY_BUFFER, (rings * sectors * 3) * sizeof(GLfloat), sphere_normals, GL_STATIC_DRAW); // Datos de las normales de los vertices
	GLuint loc2 = glGetAttribLocation(programID[SPHERE], "aNormal");   
	glEnableVertexAttribArray(loc2); // Vertex normal
	glVertexAttribPointer( loc2, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL + 0 ); 

    glBindBuffer(GL_ARRAY_BUFFER, handle[2]);
    glBufferData(GL_ARRAY_BUFFER, (rings * sectors * 2) * sizeof(GLfloat), sphere_texcoords, GL_STATIC_DRAW); // Datos de las coordenadas de textura
	GLuint loc3 = glGetAttribLocation(programID[SPHERE], "aTexCoord");   
	glEnableVertexAttribArray(loc3); // Texture coords
	glVertexAttribPointer( loc3, 2, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL + 0 ); 

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, handle[3]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, rings * sectors * 4 * sizeof(GLushort), sphere_indices, GL_STATIC_DRAW); // Array de indices

    delete [] sphere_vertices;
    delete [] sphere_normals;
    delete [] sphere_texcoords;
    delete [] sphere_indices;

    glBindVertexArray(0);

	return rings * sectors * 4;
}

/**
 * @brief Construye la geometria de una plano preparando los VBO y su VAO
 * 
 * @param xsize Longitud en el eje x
 * @param zsize Longitud en el eje z
 * @param xdivs Divisiones en el eje x
 * @param zdivs Divisiones en el eje z
 * 
 * @return Numero de vertices del plano
 */
GLint buildPlane(GLfloat xsize, GLfloat zsize, GLint xdivs, GLint zdivs)
{
    
    GLfloat * v = new GLfloat[3 * (xdivs + 1) * (zdivs + 1)];
	GLfloat * n = new GLfloat[3 * (xdivs + 1) * (zdivs + 1)];
    GLfloat * tex = new GLfloat[2 * (xdivs + 1) * (zdivs + 1)];
    GLuint * el = new GLuint[6 * xdivs * zdivs];

    GLfloat x2 = xsize / 2.0f;
    GLfloat z2 = zsize / 2.0f;
    GLfloat iFactor = (GLfloat)zsize / zdivs;
    GLfloat jFactor = (GLfloat)xsize / xdivs;
    GLfloat texi = 1.0f / zdivs;
    GLfloat texj = 1.0f / xdivs;
    GLfloat x, z;
    GLint vidx = 0, tidx = 0;
    for( GLint i = 0; i <= zdivs; i++ ) {
        z = iFactor * i - z2;
        for( GLint j = 0; j <= xdivs; j++ ) {
            x = jFactor * j - x2;
            v[vidx] = x;
            v[vidx+1] = 0.0f;
            v[vidx+2] = z;
			n[vidx] = 0.0f;
			n[vidx+1] = 1.0f;
			n[vidx+2] = 0.0f;
            vidx += 3;
            tex[tidx] = j * texi;
            tex[tidx+1] = i * texj;
            tidx += 2;
        }
    }

    GLuint rowStart, nextRowStart;
    GLint idx = 0;
    for( GLint i = 0; i < zdivs; i++ ) {
        rowStart = i * (xdivs+1);
        nextRowStart = (i+1) * (xdivs+1);
        for( GLint j = 0; j < xdivs; j++ ) {
            el[idx] = rowStart + j;
            el[idx+1] = nextRowStart + j;
            el[idx+2] = nextRowStart + j + 1;
            el[idx+3] = rowStart + j;
            el[idx+4] = nextRowStart + j + 1;
            el[idx+5] = rowStart + j + 1;
            idx += 6;
        }
    }

    GLuint handle[4];
    glGenBuffers(4, handle);

	glGenVertexArrays( 1, &planeVAO );
    glBindVertexArray(planeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, handle[0]);
    glBufferData(GL_ARRAY_BUFFER, 3 * (xdivs+1) * (zdivs+1) * sizeof(GLfloat), v, GL_STATIC_DRAW);
	GLuint loc1 = glGetAttribLocation(programID[SPHERE], "aPosition");   
    glVertexAttribPointer( loc1, 3, GL_FLOAT, GL_FALSE, 0, ((GLubyte *)NULL + (0)) );
    glEnableVertexAttribArray(loc1);  // Vertex position

	glBindBuffer(GL_ARRAY_BUFFER, handle[1]);
    glBufferData(GL_ARRAY_BUFFER, 3 * (xdivs+1) * (zdivs+1) * sizeof(GLfloat), n, GL_STATIC_DRAW);
	GLuint loc2 = glGetAttribLocation(programID[SPHERE], "aNormal");   
    glVertexAttribPointer( loc2, 3, GL_FLOAT, GL_FALSE, 0, ((GLubyte *)NULL + (0)) );
    glEnableVertexAttribArray(loc2);  // Vertex normal

    glBindBuffer(GL_ARRAY_BUFFER, handle[2]);
    glBufferData(GL_ARRAY_BUFFER, 2 * (xdivs+1) * (zdivs+1) * sizeof(GLfloat), tex, GL_STATIC_DRAW);
	GLuint loc3 = glGetAttribLocation(programID[SPHERE], "aTexCoord");   
    glVertexAttribPointer( loc3, 2, GL_FLOAT, GL_FALSE, 0, ((GLubyte *)NULL + (0)) );
    glEnableVertexAttribArray(loc3);  // texture coords

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, handle[3]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 6 * xdivs * zdivs * sizeof(GLuint), el, GL_STATIC_DRAW);

    glBindVertexArray(0);
    
    delete [] v;
	delete [] n;
    delete [] tex;
    delete [] el;

	return 6 * xdivs * zdivs;
}

int main(int argc, char *argv[])
{
	simulationMain(argc, argv);

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

/**
 * @brief Inicializa los estados de la OpenGL y los shaders de los programas graficos
 * 
 * @retval True
 */
GLboolean init()
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glClearDepth(1.0f);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glShadeModel(GL_SMOOTH);

	initUniform();
	initSphereProgram();
	initSpriteProgram();

	return true;
}

/**
 * @brief Inicializa los localizadores de las variables uniform de los shaders
 * 
 */
void initUniform()
{
	// Comunes
	locUniformMVPM = 0;
	locUniformMVM = 1;

	locUniformMaterialAmbient = 7;
	locUniformMaterialDiffuse = 8;
	locUniformMaterialSpecular = 9;
	locUniformMaterialShininess = 10;

	locUniformLightPos = 11;
	locUniformLightIntensity = 12;
	locUniformLightK = 13;

	locUniformColor = 14;
	locUniformColorMode = 15;

	// Sphere
	locUniformNM = 2;

	// Sprites
	locUniformPM = 3;
	locUniformSpriteTex = 4;
	locUniformScale = 5;
	locUniformPos = 6;	
}

/**
 * @brief Carga los shaders y las gemoterias para dibujar las particulas como esferas sobre un plano
 *
 */
void initSphereProgram()
{
	programID[SPHERE] = glCreateProgram();
	ShaderLoader shaderLoader(programID[SPHERE]);
	shaderLoader.addShader("shaders/sphere.vert");
	shaderLoader.addShader("shaders/sphere.frag");
	shaderLoader.loadShaders();

	numVertSphere = buildSphere(particleRad, 20 / 2, 30 / 3);
	numVertPlane = buildPlane(10.0f, 10.0f, 20, 20);
	
	// Variables de control
	locUniformObject = glGetUniformLocation(programID[SPHERE], "uObject");
}

/**
 * @brief Carga los shaders y la textura para dibujar las particulas como sprites orientados
 *
 */
void initSpriteProgram()
{
	programID[SPRITE] = glCreateProgram();

	ShaderLoader shaderLoader(programID[SPRITE]);
	shaderLoader.addShader("shaders/sprite.vert");
	shaderLoader.addShader("shaders/sprite.geom");
	shaderLoader.addShader("shaders/sprite.frag");
	shaderLoader.loadShaders();

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

		cout << "Textura cargada " << img_filename << std::endl;
	}
	else
		cout << "Error al cargar la textura " << img_filename << std::endl;
	img_data.clear();
}

/**
 * @brief Dibuja una particula como sprite
 *
 */
void drawSprites()
{
	glDrawArrays(GL_POINTS, 0, 1);
}

/**
 * @brief Dibuja una particula como esfera
 *
 */
void drawSphere()  
{
    glBindVertexArray(sphereVAO);
    glDrawElements(GL_QUADS, numVertSphere, GL_UNSIGNED_SHORT, ((GLubyte *)NULL + (0)));
	glBindVertexArray(0);
}

/**
 * @brief Dibuja un plano
 *
 */
void drawPlane() 
{
    glBindVertexArray(planeVAO);
    glDrawElements(GL_TRIANGLES, numVertPlane, GL_UNSIGNED_INT, ((GLubyte *)NULL + (0)));
	glBindVertexArray(0);
}

/**
 * @brief Funcion de dibujado de toda la escena representando las particulas como sprites
 *
 * Envia los valores de las variables uniform y ejecuta el correspondiente programa grafico
 * 
 * Si esta el modo cache activado se dibujan los frames ya simulados y almacenados en la aplicacion, 
 * si no esta activado, se dibuja el estado actual de la simulacion 
 */
void displaySpriteScene(const mat4 & Projection, const mat4 & View, const mat4 & Model)
{
	mat4 mv = View * Model;  
	vec4 pos;

	vector<Vector3r> * positions;
	vector<Vector3r> * bpositions;
	vector<Real> * parameter;
	vector<Real> vel; 
	uint n, n_boundary = 0;
	float max_param;

	glUseProgram(programID[SPRITE]);

	glUniform1i(locUniformSpriteTex, 0);
	glUniform1f(locUniformScale, sim -> getParticleRadius() * scaleFac * 1.5); // * 1.3 tamaño real
	glUniformMatrix4fv(locUniformMVM, 1, GL_FALSE, &mv[0][0]);
	glUniformMatrix4fv(locUniformPM, 1, GL_FALSE, &Projection[0][0]);

	// Material
	glUniform3fv(locUniformMaterialAmbient, 1, &(currentMaterial -> ambient.x));
	glUniform3fv(locUniformMaterialDiffuse, 1, &(currentMaterial -> diffuse.x));
	glUniform3fv(locUniformMaterialSpecular, 1, &(currentMaterial -> specular.x));
	glUniform1f(locUniformMaterialShininess, currentMaterial -> shininess);

	// Lights
	vec4 lpos = View * Model * pointLight.lightPos;
	glUniform4fv(locUniformLightPos, 1, &(lpos.x));
	glUniform3fv(locUniformLightIntensity, 1, &(pointLight.intensity.r));
	glUniform3fv(locUniformLightK, 1, &(pointLight.k[0]));

	// Color Mode
	glUniform1i(locUniformColorMode, colorMode);

	// Posiciones guardadas
	if (cache_mode && cache.position.size() > 0)
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
				pos = vec4((bpositions->at(i)), 1.0);
				glUniform4fv(locUniformPos, 1, &pos.x);
				glUniform1f(locUniformColor, 0.5);

				drawSprites();
			}
		}
	}

	// Fluid
	for (uint i = 0; i < n; ++i)
	{	
		pos = vec4((positions->at(i)), 1.0);
		glUniform4fv(locUniformPos, 1, &pos.x);
		glUniform1f(locUniformColor, clamp(sqrt((float) parameter->at(i) / max_param), 0.0f, 1.0f));

		drawSprites();
	}

	glUseProgram(0);
}

/**
 * @brief Funcion de dibujado de toda la escena representando las particulas como esferas
 *
 * Envia los valores de las variables uniform y ejecuta el correspondiente programa grafico
 * 
 * Si esta el modo cache activado se dibujan los frames ya simulados y almacenados en la aplicacion, 
 * si no esta activado, se dibuja el estado actual de la simulacion 
 */
void displayGeometryScene(const mat4 & Projection, const mat4 & View, const mat4 & Model)
{
	mat4 mv;
	mat4 mvp;
	mat3 nm;  
	vec4 pos;

	vector<Vector3r> * positions;
	vector<Real> * parameter;
	vector<Real> vel; 
	uint n/*, n_boundary */= 0;
	float max_param = 2.5;

	glUseProgram(programID[SPHERE]);
	glUniform1i(locUniformObject, 1);

	// Material
	glUniform3fv(locUniformMaterialAmbient, 1, &(currentMaterial -> ambient.x));
	glUniform3fv(locUniformMaterialDiffuse, 1, &(currentMaterial -> diffuse.x));
	glUniform3fv(locUniformMaterialSpecular, 1, &(currentMaterial -> specular.x));
	glUniform1f(locUniformMaterialShininess, currentMaterial -> shininess);

	// Point light
	vec4 lpos = View * Model * pointLight.lightPos;
	glUniform4fv(locUniformLightPos, 1, &(lpos.x));
	glUniform3fv(locUniformLightIntensity, 1, &(pointLight.intensity.r));
	glUniform3fv(locUniformLightK, 1, &(pointLight.k[0]));

	// Directional Light
	mv = View * Model;
	nm = mat3(transpose(inverse(mv)));
	lpos = mv * dirLight.lightPos;	// S.R. vista
	glUniform4fv(locUniformDirLightPos, 1, &(lpos.x));
	glUniform3fv(locUniformDirLightIntensity, 1, &(dirLight.intensity.r));
	glUniform3fv(locUniformDirLightK, 1, &(dirLight.k[0]));

	// Color Mode
	glUniform1i(locUniformColorMode, colorMode);

	// Posiciones guardadas
	if (cache_mode && cache.position.size() > 0)
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
			max_param = 5;
		}
	}
	// Posiciones actuales del sistema de particulas
	else
	{
		n = sim -> numberActiveParticles(0);
		positions = &sim -> getPositions(0);

		if (colorMode == 2)
		{
			parameter = &sim -> getPressures(0);
			max_param = 3000;
		}
		else 
		{
			vel = sim -> getVelocities(0);
			parameter = &vel;
			max_param = 5;			
		}
	}

	/*if (drawBoundary)
		n_boundary = sim.sys.n_boundary;*/

	// Fluid
	for (uint i = 0; i < n; ++i)
	{
		pos = vec4((positions->at(i)), 1.0);
		mv = View * scale(translate(Model, vec3(pos)), vec3(particleScale));
		mvp = Projection * mv;
		nm = mat3(transpose(inverse(mv)));

		glUniformMatrix4fv(locUniformMVPM, 1, GL_FALSE, &mvp[0][0]);
		glUniformMatrix4fv(locUniformMVM, 1, GL_FALSE, &mv[0][0]);
		glUniformMatrix3fv(locUniformNM, 1, GL_FALSE, &nm[0][0]);
		glUniform1f(locUniformColor, clamp(sqrt((float) parameter->at(i)/ max_param), 0.0f, 1.0f));

		drawSphere();
	}

	// Boundary
	/*for (uint i = sim -> numberActiveParticles(0); i < sim -> numberActiveParticles + n_boundary; ++i)
	//for (uint i: sim.sys.movingWallIndexPX)
	{		
		pos = vec4((positions->at(i)), 1.0);
		mv = View * scale(translate(Model, vec3(pos)), vec3(particleScale));
		mvp = Projection * mv;
		nm = mat3(transpose(inverse(mv)));

		glUniformMatrix4fv(locUniformMVPM, 1, GL_FALSE, &mvp[0][0]);
		glUniformMatrix4fv(locUniformMVM, 1, GL_FALSE, &mv[0][0]);
		glUniformMatrix3fv(locUniformNM, 1, GL_FALSE, &nm[0][0]);
		glUniform1f(locUniformColor, clamp(sqrt((float) parameter->at(i)/ max_param), 0.0f, 1.0f));

		drawSphere();
	}*/

	glUseProgram(0);
}

/**
 * @brief Funcion de callback de dibujado de la escena
 * 
 * Crea las matrices Model, View y Projection para los shaders y llama al dibujado de la escena correspondiente (Sprites o Spheres)
 */
void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		

	mat4 Projection = perspective(45.0f, 1.0f * g_Width / g_Height, 1.0f, 100.0f);
	mat4 View = lookAt(cameraPos, look, vec3(0.0f, 1.0f, 0.0f));
	mat4 Model = scale(mat4(1.0f), vec3(scaleFac, scaleFac, scaleFac));

	// Rotacion del modelo mediante el raton
	Model = rotate(Model, xrot * 0.0025f, vec3(1.0f, 0.0f, 0.0f));
	Model = rotate(Model, yrot * 0.0025f, vec3(0.0f, 1.0f, 0.0f));

	// Rotacion para pasar del s.c. de Blender (z -> altura y -> profundidad) al de OpenGL (y -> altura z -> profundidad)
	Model = rotate(Model, - (GLfloat) (M_PI * /*0.35*/ 0.5), vec3(1, 0, 0));
	//Model = rotate(Model, (GLfloat) (M_PI * 0.25), vec3(0, 0, 1));

	// Dibujar la escena
	if (sprite)
		displaySpriteScene(Projection, View, Model);
	else
		displayGeometryScene(Projection, View, Model);

	glutSwapBuffers();

	// Mostrar frecuencia de dibujado
	printFPS(current_fps);
	string current_fps_s = to_string(current_fps);
	string desired_fps_s = to_string(desired_fps);
	string current_frame = to_string(cache.frameCount + 1);
	string total_frames  = to_string(cache.position.size());
	string title = current_fps_s.substr(0, current_fps_s.find(".") + 3) + "/" + 
					desired_fps_s.substr(0, desired_fps_s.find(".")) + "		" + 
					current_frame + "/" + total_frames;
	glutSetWindowTitle(title.c_str());

	// Guardar frames
	if (save_mode)
	{
		savePNGFrame();
		if (cache.frameCount == cache.position.size() - 1)
			save_mode = false;
	}

	// Actualizar contador de la cache
	if (!pause && cache_mode && cache.frameCount < cache.position.size() - 1)
		++cache.frameCount;
}

/**
 * @brief Funcion de callback de redimensionado de la ventana
 */
void resize(GLint w, GLint h)
{
	g_Width = w;
	g_Height = h;
	glViewport(0, 0, g_Width, g_Height);
}
 
/**
 * @brief Funcion de callback de idle
 * 
 * Ejecuta un paso de simulacion a la frecuencia indicada en desired_fps y llama a la funcion de redibujado
 */
void idle()
{
	static double last_render = 0;

	double frame_delay = 1.0 / desired_fps; 
	double current_time = getCurrentTime();

	if ((current_time - last_render) > frame_delay)
	{
		last_render = current_time;

		if (!pause && !cache_mode)
			if (solve())
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
}

/**
 * @brief Funcion de callback de teclado
 * 
 */
void keyboard(GLubyte key, GLint x, GLint y)
{
	switch(key)
	{
		case 27 : 
		case 'q': 
		case 'Q':
			exit(1); 
			break;

		case 'a': 
		case 'A':
			animation = !animation;
			break;

		case 'z':
		case 'Z':
			particleScale += 0.1;
			break;
		
		case 'x':
		case 'X':
			particleScale -= 0.1;
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
			cache_mode = !cache_mode;
			break;

		case 'p':
		case 'P':
			pause = !pause;
			break;

		case 'n':
		case 'N':
			colorMode = (colorMode + 1) % 6; // Material -> Normal -> Velocities -> Pressure -> Sin iluminacion
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

		case 's':
		case 'S':
			sprite = !sprite;
			break;	

		case 'm': // Change material
			m_index = (m_index + 1) % NUM_MATERIALS;
			currentMaterial = &materials[m_index];
			break;

		case 'w': // Save frames
		{
			// Comprobar si existe el directorio frames, si no existe, se crea (Linux)
			string png_dir_path = png_path.substr(0, png_path.find_last_of("/"));
			struct stat st;
			if(stat(png_dir_path.c_str(), &st) != 0)
			{
				cout << endl << "Se ha creado el directorio " + png_dir_path << endl;
				cout << mkdir(png_dir_path.c_str(), 0777) << endl;	
			}

			save_mode = !save_mode;
			cache_mode = save_mode;
			cache.frameCount = 0;
		}
			break;

		case 'o':
		case 'O':
			//sim.sys.saveStepData();
			break;
	}
}

/**
 * @brief Funcion de callback de teclado (especial)
 * 
 */
void specialKeyboard(GLint key, GLint x, GLint y)
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
		desired_fps -= 1;
			if (desired_fps < 1)
				desired_fps = 1;
	}
	else if (key == GLUT_KEY_RIGHT)
		desired_fps += 1;
}
 
/**
 * @brief Funcion de callback del raton
 * 
 */
void mouse(GLint button, GLint state, GLint x, GLint y)
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
 
/**
 * @brief Funcion de callback de movimiento de raton
 * 
 */
void mouseMotion(GLint x, GLint y)
{
	if (mouseDown)
	{
		yrot = x - xdiff;
		xrot = y + ydiff;
 
		glutPostRedisplay();
	}
}

/**
 * @brief Devuelve el tiempo actual usando el reloj de alta resolucion
 * 
 */
double getCurrentTime()
{
	using Duration = std::chrono::duration<double>;
	return std::chrono::duration_cast<Duration>( 
		std::chrono::high_resolution_clock::now().time_since_epoch() 
	).count();
}

/**
 * @brief Muestra por consola cada segundo la frecuencia a la que se esta dibujando 
 * 
 * @retval Si se ha calculado la frecuencia
 */
bool printFPS(float & fps)
{
	static GLint count = 0;			//  Numero de frames
	static GLint currentTime = 0, previousTime = 0;

	count++;

	currentTime = glutGet(GLUT_ELAPSED_TIME);

	GLuint timeInterval = currentTime - previousTime;

	if (timeInterval > 1000)
	{
		fps = count / (timeInterval / 1000.0f);
		previousTime = currentTime;
		count = 0;

		return true;
	}

	return false;
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "CubeBoundaryModel.h"
#include "PCISPHBoundaryModel.h"
#include "AkinciBoundaryModel.h"
#include "WCSPHSolver.h"
#include "PCISPHSolver.h"
#include "DFSPHSolver.h"
#include "Emitter.h"

bool readState(string filename, vector<Vector3r> & points, vector<Vector3r> & velocities)
{
    ifstream in_file(filename);
    unsigned n;
    string tmp;
    string x, y, z;

    if (!in_file.is_open())
        return false;
    else
    {
        std::stringstream buffer;
        buffer << in_file.rdbuf();
        in_file.close();
        
        // Leer y guardar el ultimo paso de simulacion calculado para ese frame
        buffer >> tmp;     //time = stod(tmp); // Time
        buffer >> tmp; // Num parts

        n = (unsigned) stoul(tmp);

        points.resize(n);
        velocities.resize(n);

        unsigned i = 0;
        while (buffer >> x)
        {
            buffer >> y;
            buffer >> z;

            if (i < n)
            {
                points[i].x = (Real) stod(x);
                points[i].y = (Real) stod(y);
                points[i].z = (Real) stod(z);
            }
            else
            {
                velocities[i - n].x = (Real) stod(x);
                velocities[i - n].y = (Real) stod(y);
                velocities[i - n].z = (Real) stod(z);
            }

            i++;
        }

        return true;
    }
}

bool readScene(string filename)
{
    ifstream in_file(filename);
    string line;

    if (!in_file.is_open())
    {
        return false;  
    }
    else
    {
        while (getline(in_file, line))
        {
            size_t i = line.find("=") + 1;
            string param = line.substr(i, line.size());

            /*if (line.find("nparts") != line.npos)
                sys.resize(stoul(param));*/
            if (line.find("particle_volume") != line.npos)
                particleVolume = static_cast<Real>(stod(param));
            /*else if (line.find("fluid_volume") != line.npos)
                sys.fluid_volume = stod(param);
            else if (line.find("kernel_parts") != line.npos)
                sys.kernel_parts = stoul(param);
            else if (line.find("tini") != line.npos)
                t_ini = stod(param);
            else if (line.find("tfin") != line.npos)
                t_fin = stod(param);
            else if (line.find("tstep") != line.npos)
                step = stod(param);
            else if (line.find("save_freq") != line.npos)
                fps = stoul(param);
            else
                cout << "Cannot read line: " << line << endl;*/
        }

        in_file.close();

        return true;
    }
}

bool readDomain(string filename, vector<Vector3r> & points)
{
    ifstream in_file(filename);
    string word;
    vector<string> dom_pts;
    
    points.resize(2);

    if (!in_file.is_open())
        return false;
    else
    {
        while (in_file >> word)
            dom_pts.push_back(word);

        in_file.close();

        points[0].x = (Real) stod(dom_pts[0]);
        points[1].x = (Real) stod(dom_pts[1]);
        points[0].y = (Real) stod(dom_pts[2]);
        points[1].y = (Real) stod(dom_pts[3]);
        points[0].z = (Real) stod(dom_pts[4]);
        points[1].z = (Real) stod(dom_pts[5]);

        return true;
    }
}

// Crea una escena con dos cubos y un dragon 
void createBoundaryScene(Real radius)
{
	// Cambiar y por z en la escala al coger escena de blender

	sim -> setBoundaryMethod(Simulation::AKINCI_BOUNDARY_METHOD);

	// Dragon obstaculo
	/*sim -> addBoundaryModelFromOBJ("data/Models/Dragon_50k.obj", radius * 2.0, Vector3r(0.642, 0.642, 0.642), Vector3r(0.588891, 0, -0.5), Vector3r(M_PI_2, 0.0, 0.0));
	// Caja obstaculo
	sim -> addBoundaryModelFromOBJ("data/Models/UnitBox.obj", radius * 1.0, Vector3r(0.277), Vector3r(-0.41778, -0.22155, -0.35765), Vector3r(0.0));
	// Caja obstaculo
	sim -> addBoundaryModelFromOBJ("data/Models/UnitBox.obj", radius * 1.0, Vector3r(0.277), Vector3r(-0.034815, 0.221551, -0.357647), Vector3r(0.0));
	// Caja dominio
	sim -> addBoundaryModelFromOBJ("data/Models/UnitBox.obj", radius * 1.5, Vector3r(2.1866, 0.831253, 1.0), Vector3r(0.0, 0.0, 0.0), Vector3r(0.0));*/

	sim -> addBoundaryModelFromOBJ("data/Models/test5.obj", radius * 1.5, Vector3r(1.0), Vector3r(0.0), Vector3r(M_PI_2, 0.0, 0.0));
}

FluidModel* createCrownScene(Real radius)
{
	FluidBlockInfo fbi1 = {Vector3r(0.0, 0.0, 0.024537), 75, 75, 20};
	FluidBlockInfo fbi2 = {Vector3r(0.0, 0.0, 0.972614), 15, 15, 15};
	std::vector<FluidBlockInfo> fbi_vec;

	fbi_vec.push_back(fbi1);
	fbi_vec.push_back(fbi2);

	//radius = 0.00099;

	FluidModel *fm = sim -> buildFluidBlock(radius, fbi_vec);

	return fm;
}

void simulationMain(int argc, char* argv[])
{
    vector<Vector3r> fluidPoints;
    vector<Vector3r> fluidVelocities;
    vector<Vector3r> boundaryPoints;

	string filename = "data/Scenes/test4";

    cout << "Reading..." << endl;
    if (readState(filename + ".dat", fluidPoints, fluidVelocities))
        cout << "Done!" << endl;
    else
        cout << "Error reading" << endl;

    cout << "Reading..." << endl;
    if (readScene(filename + ".scn"))
        cout << "Done!" << endl;
    else
        cout << "Error reading" << endl;

    cout << "Reading..." << endl;
    if (readDomain(filename + ".dom", boundaryPoints))
        cout << "Done!" << endl;
    else
        cout << "Error reading" << endl;

    Simulation *sim = Simulation::getCurrent();

	unsigned int kernelParticles = 27;
	Real fps = 120;
	Real timeStep = 0.0001;
	Real maxTimeStep = 0.002;
	Real minTimeStep = 0.000;
	Real density0 = 1000;
	Real stiffness = 100;
	Real gamma = 1;
	Real eta = 0.0005;  // 0.05%
	Real etaV = 0.001; // 0.1%
	Real etaPCI = 0.01; // 1 %
	Real viscosity = 0.01;
	Real bViscosity = 0.0;
	Real surften = 0.175;
	Real beta = 1.0; // 10000
	Real side = pow(particleVolume, 1.0 / 3.0); // Lado del cubo con el volumen recogido de blender
	Real radius = side / 2.0;
	particleVolume = pow(2.0 * radius, 3.0);
	Vector3r gravity(0.0, 0.0, -9.8);


	////////// NOTA INICIALIZACION FLUIDO RADIO Y SUPPORT RADIUS //////////////////////
	//
	//	El espaciado entre particulas sera de 2 * radius
	//	El supportRadius sera 4 * radius
	//	El volumen sera un cubo con lado 2 * radius
	//
	// 	El particleVolume que se coge de los ficheros de blender no es el volumen de la particula sino del cubo

	////////////////////// General parameters //////////////////////////
	sim -> setKernelParticles(kernelParticles);
    sim -> setParticleRadius(radius);
	sim -> setGravity(gravity);
    sim -> setTimeStep(timeStep);
	sim -> setFPS(fps);
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
    /*sim -> setSimulationMethod(Simulation::WCSPH_METHOD);
	WCSPHSolver *solver = static_cast<WCSPHSolver*>(sim -> getSolver());
	solver -> setStiffness(stiffness);
	solver -> setGamma(gamma);*/

	/*sim -> setSimulationMethod(Simulation::PCISPH_METHOD);
	PCISPHSolver *solver = static_cast<PCISPHSolver*>(sim -> getSolver());
	solver -> setMaxError(etaPCI);*/

	sim -> setSimulationMethod(Simulation::DFSPH_METHOD);
	DFSPHSolver *solver = static_cast<DFSPHSolver*>(sim -> getSolver());
	solver -> setMaxError(eta);
	solver -> setMaxErrorV(etaV);
	////////////////////////////////////////////////////////////////////

	////////////////////// Fluid Model /////////////////////////////////
	//fluidPoints.resize(0);
	//fluidVelocities.resize(0);
    sim -> addFluidModel(fluidPoints, fluidVelocities);
	FluidModel *fm = sim -> getFluidModel(0);

	/*std::vector<FluidBlockInfo> fbi_vec = {{Vector3r(0.0), 26, 26, 26}};
	FluidModel *fm = sim -> buildFluidBlock(radius, fbi_vec);*/

	//FluidModel *fm = createCrownScene(radius);

    fm -> setRefDensity(density0);
    fm -> setMasses(density0 * particleVolume);
	fm -> setViscosityForce(viscosity, bViscosity);
	fm -> setSurfaceTensionForce(surften);
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
					 0.03);

	rot = glm::rotate(rot, M_PI, Vector3r(0.0, 1.0, 0.0));

	fm -> addEmitter(Emitter::SQUARE_EMITTER, 
					 40000, 
					 Vector3r(-0.01, 0.0, 0.3), 
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
	pcibm -> setMass(density0 * particleVolume);
	pcibm -> setRefDensity(density0);
	pcibm -> init(boundaryPoints);
	sim -> setBoundaryMethod(Simulation::PCISPH_BOUNDARY_METHOD);
	sim -> addBoundaryModel(pcibm);*/

	// Escena personalizada
	//createBoundaryScene(radius);

	

	// Sphere
	/*std::vector<Vector3r> sphereBoundaryPoints = {Vector3r(0.0, 0.0, 0.18838)};
	sim -> setBoundaryMethod(Simulation::AKINCI_BOUNDARY_METHOD);
	AkinciBoundaryModel *abm2 = new AkinciBoundaryModel();
	abm2 -> setRadius(0.53/2);
	abm2 -> init(sphereBoundaryPoints);
	sim -> addBoundaryModel(abm2);*/

	// domain
	sim -> setBoundaryMethod(Simulation::AKINCI_BOUNDARY_METHOD);
	AkinciBoundaryModel *abm = new AkinciBoundaryModel();
	abm -> init(boundaryPoints);
	sim -> addBoundaryModel(abm);

	////////////////////////////////////////////////////////////////////
    
	///////////// Init simulation //////////////////////////////////////
    sim -> init();
	////////////////////////////////////////////////////////////////////

	std::cout << "Num fluid particles " << sim -> numberActiveParticles(0) << std::endl;
}

inline bool solve()
{
    return sim -> step();
}

