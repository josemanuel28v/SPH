#ifndef _SHADER_LOADER_
#define _SHADER_LOADER_

#include <iostream>
#include <vector>
#include <fstream>
#include <GL/glew.h>

class ShaderLoader
{
    GLuint programID;

    std::vector<std::string> filenames;
    std::vector<GLuint> shaderIDs;

    public:

        ShaderLoader(GLuint);

        void createShader();
        void loadSource(std::string, GLuint);
        void printCompileInfoLog(GLuint);
        void printLinkInfoLog();
        void validateProgram();

        void loadShaders();

        void setProgramID(GLuint);
        void addShader(std::string);
};

#endif
