#include "ShaderLoader.h"

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
        std::cerr << "File not found " << filename.c_str() << std::endl;
        //system("pause");
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

        std::cerr << "Shader compiling failed:" << infoLog << std::endl;
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

        std::cerr << "Shader linking failed:" << infoLog << std::endl;
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
            std::cerr << "Program validating failed:" << infoLog << std::endl;
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
        std::cout << "Compilando shader " << filenames[i] << std::endl;
        glCompileShader(shaderIDs[i]);                      // Compilar shader
        printCompileInfoLog(shaderIDs[i]);                  // Imprimir log de compilacion
        glAttachShader(programID, shaderIDs[i]);            // Añadir el shader al programa
    }
    glLinkProgram(programID);                               // Enlazar programa
    std::cout << "Linkando objeto programa ..." << std::endl;
    printLinkInfoLog();
    validateProgram();                                  
}