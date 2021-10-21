#ifndef _OGL_H_
#define _OGL_H_

#include <glm/glm.hpp>
#include <GL/glew.h>
#include <GL/glut.h>
#include <vector>
#include <string>
#include "types.h"
#include "Simulation.h"
#include "ShaderLoader.h"

class OGL
{
    public:

        struct LightInfo 
        {
            glm::vec4 lightPos; 	
            glm::vec3 intensity;
            glm::vec3 k;
        };

        struct Cache 
        {
            std::vector<std::vector<Vector3r>> position;
            std::vector<std::vector<Real>> velocity;
            std::vector<std::vector<Real>> pressure;

            std::vector<GLuint> numParts;
            std::vector<Real> time;

            GLuint frameCount = 0;
            GLuint size = 0;
        };

        struct MaterialInfo 
        {
            glm::vec3 ambient;
            glm::vec3 diffuse;
            glm::vec3 specular;
            GLfloat shininess;
        };

        
        int mainLoop(int argc, char *argv[]);

    private:

        static GLint g_Width;
        static GLint g_Height;
        static GLboolean pause;
        static GLboolean cacheMode;
        static GLboolean fullscreen;
        static GLboolean mouseDown;
        static GLboolean drawBoundary;
        static GLboolean saveMode;
        static GLfloat xdiff;
        static GLfloat ydiff;
        static GLfloat xrot;
        static GLfloat yrot;
        static GLfloat desiredFps;
        static GLfloat scaleFac;
        static GLint colorMode;
        static GLint m_index;
        static std::string png_path;
        static double lastRender;

        static glm::vec3 cameraPos;
        static glm::vec3 look;

        static ShaderLoader spriteShader;

        static LightInfo pointLight;
        static Cache cache;
        static MaterialInfo materials[];
        static MaterialInfo currentMaterial;

        // Auxiliar fucntions
        static double getCurrentTime();
	    static void initSpriteProgram();
        static void updateTitle(double);
        void init();
        static void drawSprites();
        static void displaySpriteScene(const glm::mat4 & Projection, const glm::mat4 & View, const glm::mat4 & Model);

        // Callback functions
        static void display();
        static void keyboard(GLubyte key, GLint x, GLint y);
        static void specialKeyboard(GLint key, GLint x, GLint y);
        static void mouse(GLint button, GLint state, GLint x, GLint y);
        static void mouseMotion(GLint, GLint);
        static void resize(GLint, GLint);
        static void idle();
};

#endif
