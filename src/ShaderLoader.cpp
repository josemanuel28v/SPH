#include "ShaderLoader.h"
#include "Logger.h"

ShaderLoader::ShaderLoader(GLuint programID)
{
    this -> programID = programID;
}

void ShaderLoader::createShader()
{
    for (std::string filename: filenames)
    {
        if (filename.find(".vert") != std::string::npos)
            shaderIDs.push_back(glCreateShader(GL_VERTEX_SHADER));
        else if (filename.find(".geom") != std::string::npos)
            shaderIDs.push_back(glCreateShader(GL_GEOMETRY_SHADER));
        else if (filename.find(".frag") != std::string::npos)
            shaderIDs.push_back(glCreateShader(GL_FRAGMENT_SHADER));
    }
}

void ShaderLoader::loadSource(std::string filename, GLuint id) 
{
    std::ifstream f(filename.c_str());
    if (!f.is_open()) 
    {
        ERROR("File not found ", filename.c_str());
        exit(EXIT_FAILURE);
    }

    // now read in the data
    std::string *source;
    source = new std::string( std::istreambuf_iterator<char>(f),   
                        std::istreambuf_iterator<char>() );
    f.close();

    // add a null to the string
    *source += "\0";
    const GLchar * data = source->c_str();
    glShaderSource(id, 1, &data, NULL);
    delete source;
}

void ShaderLoader::printCompileInfoLog(GLuint id) 
{
    GLint compiled;
    glGetShaderiv( id, GL_COMPILE_STATUS, &compiled );
    if (compiled == GL_FALSE)
    {
        GLint infoLength = 0;
        glGetShaderiv( id, GL_INFO_LOG_LENGTH, &infoLength );

        GLchar *infoLog = new GLchar[infoLength];
        GLint chsWritten = 0;
        glGetShaderInfoLog( id, infoLength, &chsWritten, infoLog );

        ERROR("Shader compiling failed: ", infoLog);

        //system("pause");
        delete [] infoLog;

        exit(EXIT_FAILURE);
    }
}

void ShaderLoader::printLinkInfoLog()
{
    GLint linked;
    glGetProgramiv( programID, GL_LINK_STATUS, &linked );
    if(! linked)
    {
        GLint infoLength = 0;
        glGetProgramiv( programID, GL_INFO_LOG_LENGTH, &infoLength );

        GLchar *infoLog = new GLchar[infoLength];
        GLint chsWritten = 0;
        glGetProgramInfoLog( programID, infoLength, &chsWritten, infoLog );

        ERROR("Shader linking failed: ", infoLog);
        //system("pause");
        delete [] infoLog;

        exit(EXIT_FAILURE);
    }
}

void ShaderLoader::validateProgram()
{
    GLint status;
    glValidateProgram( programID );
    glGetProgramiv( programID, GL_VALIDATE_STATUS, &status );

    if( status == GL_FALSE ) 
    {
        GLint infoLength = 0;
        glGetProgramiv( programID, GL_INFO_LOG_LENGTH, &infoLength );

        if( infoLength > 0 ) 
        {
            GLchar *infoLog = new GLchar[infoLength];
            GLint chsWritten = 0;
            glGetProgramInfoLog( programID, infoLength, &chsWritten, infoLog );
            ERROR("Program validating failed: ",  infoLog);
            //system("pause");
            delete [] infoLog;

            exit(EXIT_FAILURE);
        }
    }
}

void ShaderLoader::setProgramID(GLuint programID)
{
    this -> programID = programID;
}

void ShaderLoader::addShader(std::string filename)
{
    filenames.push_back(filename);
}

void ShaderLoader::loadShaders()
{
    createShader();                                         // Crear ID de cada shader
    for (unsigned i = 0; i < filenames.size(); ++i)
    {
        loadSource(filenames[i], shaderIDs[i]);             // Leer cada shader
        LOG("Compilando shader ", filenames[i]);
        glCompileShader(shaderIDs[i]);                      // Compilar shader
        printCompileInfoLog(shaderIDs[i]);                  // Imprimir log de compilacion
        glAttachShader(programID, shaderIDs[i]);            // Añadir el shader al programa
    }
    glLinkProgram(programID);                               // Enlazar programa
    LOG("Linkando objeto programa ...");
    printLinkInfoLog();
    validateProgram();                                  
}

void ShaderLoader::addAttribute(std::string name)
{
    attributes[name] = glGetAttribLocation(programID, name.c_str());
}

void ShaderLoader::addUniform(std::string name, GLint id)
{
    uniforms[name] = id;
}

GLint ShaderLoader::getUniformLocation(std::string name)
{
    return uniforms[name];
}

GLint ShaderLoader::getAttributeLocation(std::string name)
{
    return attributes[name];
}

void ShaderLoader::begin()
{
    glUseProgram(programID);
}

void ShaderLoader::end()
{
    glUseProgram(0);
}