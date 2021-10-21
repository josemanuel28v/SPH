#ifndef _SHADER_LOADER_
#define _SHADER_LOADER_

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <GL/glew.h>

class ShaderLoader
{
    GLuint programID;

    std::map<std::string, GLint> attributes;
    std::map<std::string, GLint> uniforms;

    std::vector<std::string> filenames;
    std::vector<GLuint> shaderIDs;

    void createShader();
    void loadSource(std::string, GLuint);
    void printCompileInfoLog(GLuint);
    void printLinkInfoLog();
    void validateProgram();

    public:

        ShaderLoader() {}
        ShaderLoader(GLuint);

        void loadShaders();

        void setProgramID(GLuint);
        void addShader(std::string);
        
        void addAttribute(std::string);
        void addUniform(std::string, GLint);

        GLint getUniformLocation(std::string name);
        GLint getAttributeLocation(std::string name);

        void begin();
        void end();
};

#endif
