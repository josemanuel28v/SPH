#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <sstream>
#include <iostream>

class Logger
{
    private:

        // Añadir niveles de log para seleccionar a partir de que nivel se muestran cosas
        // ERROR, WARN, INFO, TIME (importantes info y time para poder ver o no los tiempos)

        static std::ostringstream stream;

    public:

        template<typename T, typename... Types>
        static void write(T s, Types... args)
        {    
            std::cout << s;      

            write(args...);   
        }

        static void write()
        {
            std::cout << std::endl;
        }

        template<typename T, typename... Types>
        static void writeError(T s, Types... args)
        {  
            std::cerr << s;      

            write(args...);   
        }

        static void writeError()
        {
            std::cerr << std::endl;
        }

        static std::string str()
        {
            return stream.str();
        }
};

#define LOG Logger::write
#define ERROR Logger::writeError
//#define LOG Logger::write(ignored) para que el compilador ignore los LOG al compilar

#endif
